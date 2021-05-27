module Model

using LightGraphs
using GraphIO

using DataStructures
using SparseArrays
using Optim
using LineSearches
using StatsBase
using Statistics
using JSON
using Dates
using Serialization
import Optim: converged


include("./structures.jl")
include("./nodemap.jl")

setprecision(BigFloat, 64)

mutable struct GraphModel
    dataset::String # constant
    structure_alphabet_size::Int64 # constant
    G::LightGraphs.SimpleGraph{Int64} # constant
    A::SparseArrays.SparseMatrixCSC{Int64,Int64} # constant
    n::Int64 # constant
    m::Int64 # constant
    maximum_m_simple::Int64 # constant
    maximum_m_loopy::Int64 # constant
    lambdas::Array{BigFloat,1} # growing at end (and values changing by optimization)
    constraint_values::Array{BigFloat,1} # growing at end
    micro_structures::Array{Union{Bool, Set{Tuple{Int64,Int64}}},1} # growing at end
    edge_lambda_indices::SortedDict{Tuple{Int64,Int64},Array{Int64,1}} # new edges inserted, growing at selected value ends
    edge_equivalence_classes::Dict{Set{Int64},Array{Tuple{Int64,Int64},1}} # currently reconstructed but could be smarter
    backup_optimization_result::Optim.MultivariateOptimizationResults # reconstructed
    last_optimization_result::Optim.MultivariateOptimizationResults # reconstructed
    entropy::BigFloat # changing
    equivalence_class_edge_probabilities::Dict{Set{Int64},BigFloat} # reconstructed
    surprise_matrix::SparseArrays.SparseMatrixCSC{BigFloat,Int64} # reconstructed
    surprise_matrix_transpose::SparseArrays.SparseMatrixCSC{BigFloat,Int64} # reconstructed
    macro_structures::Array{Structures.Structure} # growing at end
    macro_structure_description_lengths::Array{BigFloat} # growing at end
    edge_probability_description_length_maxent::BigFloat # changing
    structure_type_description_length::BigFloat # changing
    node_surprises::Array{Tuple{Int64,BigFloat},1} # reconstructed
    surprising_edges::Array{Tuple{Int64,Int64},1} # reconstructed
    start_description_length_maxent::BigFloat # constant
    total_description_length_maxent::BigFloat # changing
    description_lengths_over_time::Array{BigFloat,1} # growing at end
    parameters::Dict{String,Float64} # constant resp. set only upon saving
    io::IO # constant
    hub_spoke_edges::Set{Tuple{Int64,Int64}} # changing
    star_spoke_edges::Set{Tuple{Int64,Int64}} # changing
    zero_edges::Set{Tuple{Int64,Int64}} # changing
    last_hub_spoke_edges::Set{Tuple{Int64,Int64}} # changing
    last_star_spoke_edges::Set{Tuple{Int64,Int64}} # changing
    last_zero_edges::Set{Tuple{Int64,Int64}} # changing - note that we'll be a bit imprecise here: if a zero edge is discovered after already having been covered, it will remain also in its previous equivalence class, if any
    covered_nodes::Set{Int64} # growing
    function GraphModel(dataset::String, io=stdout)
        setup_graph!(new(), dataset, io)
    end
    function GraphModel(kwargs::Dict)
        M = new()
        for (k,v) in kwargs
            setproperty!(M, k, v)
        end
        return M
    end
end

Base.show(io::IO, M::GraphModel) = print(io, "Model for $(M.dataset):\n$(typeof(M.G)) with n = $(M.n), m=$(M.m)")
JSON.lower(M::GraphModel) = Dict(key=> getfield(M, key) for key in fieldnames(GraphModel) if key in
            (:dataset, :n, :m, :macro_structures, :macro_structure_description_lengths,
                :edge_probability_description_length_maxent,
                :start_description_length_maxent,
                :total_description_length_maxent,
                :description_lengths_over_time,
                :lambdas,
                :parameters,
                :structure_alphabet_size,
                :structure_type_description_length
            )
)

"""
helper providing basic logging
"""
function log_progress!(M::GraphModel, msg::String)
    println(M.io, now(), " - ", msg)
    flush(M.io)
end

"""
save model as jldb (serialization format) or json
"""
function save_model(M::GraphModel, directory::String; final=true, overwrite=true)
    path, file = splitdir(M.dataset)
    last_folder = splitdir(path)[end]
    file_base, ext = splitext(file)
    size_threshold = M.parameters["size_threshold"]
    stop_after_n_structures = M.parameters["stop_after_n_structures"]
    filepath = "$(directory)/$(last_folder)/$(file_base)_size-$(size_threshold)_max-$(stop_after_n_structures)"
    mkpath(splitdir(filepath)[1])
    if !final
        open(filepath*".jldb", "w") do io
            serialize(io, M)
        end
    else
        open(filepath*".json", "w") do io
            write(io, json(M, 4))
        end
        log_progress!(M, "Saved model to $(filepath).json")
    end
    mv("$(path)/$(file_base)-nodemap.csv", "$(directory)/$(last_folder)/$(file_base)-nodemap.csv", force=overwrite)
end

"""
restore model from json
"""
function restore_model(path_to_file::String)
    model_dict = JSON.parsefile(path_to_file)
    M = Model.GraphModel(model_dict["dataset"])
    structures = Array{Structures.Structure,1}()
    for sd in model_dict["macro_structures"]
        if sd["structure_type"] == "clique"
            structure = Structures.Clique(sd["n_edges_total"], Array{Int64,1}(sd["nodes"]))
        elseif sd["structure_type"] == "star"
            structure = Structures.Star(sd["n_edges_total"], sd["hub"], Array{Int64,1}(sd["spokes"]))
        elseif sd["structure_type"] == "biclique"
            structure = Structures.Biclique(sd["n_edges_total"], Array{Int64,1}(sd["left_nodes"]), Array{Int64,1}(sd["right_nodes"]), sd["n_edges_in_left"], sd["n_edges_in_right"])
        elseif sd["structure_type"] == "starclique"
            structure = Structures.Starclique(sd["n_edges_total"], Array{Int64,1}(sd["left_nodes"]), Array{Int64,1}(sd["right_nodes"]), sd["n_edges_in_left"], sd["n_edges_in_right"])
        # generic structures are already present
        else
            continue
        end
        push!(structures, structure)
    end
    add_macro_structures!(M, structures; lambdas=get!(model_dict, "lambdas", []))
    if model_dict["lambdas"] == []
        println("WARNING: Recomputing lambdas when restoring model - true edge probability description length for computed lambdas might differ from restored edge probability description length due to different optimization process...")
    end
    M.parameters["surprise_threshold"] = model_dict["parameters"]["surprise_threshold"]  === nothing ? Inf : model_dict["parameters"]["surprise_threshold"]
    M.parameters["size_threshold"] = model_dict["parameters"]["size_threshold"]  === nothing ? Inf : model_dict["parameters"]["size_threshold"]
    M.description_lengths_over_time = model_dict["description_lengths_over_time"]
    M.macro_structure_description_lengths = model_dict["macro_structure_description_lengths"]
    M.edge_probability_description_length_maxent = model_dict["edge_probability_description_length_maxent"]
    M.total_description_length_maxent = M.edge_probability_description_length_maxent + M.structure_type_description_length + sum(M.macro_structure_description_lengths)
    return M
end

"""
initial graph setup with first optimization
"""
function setup_graph!(M::GraphModel, dataset::String, io::IO)
    set_constant_terms!(M, dataset, io)
    M.lambdas = zeros(2)
    M.constraint_values = [M.m, 0]
    M.micro_structures = [true, Set([(i,i) for i=1:M.n])]
    M.edge_lambda_indices = SortedDict((i,i)=>[1,2] for i=1:M.n)
    M.edge_equivalence_classes = Dict(Set([1,2])=>[(i,i) for i=1:M.n])
    M.macro_structures = [Structures.GenericStructure("Overall number of edges", M.m), Structures.GenericStructure("Loop-freeness", 0)]
    # below, we add the cost of transmitting n to the cost of transmitting m (since it is not a proper structure but necessary to specify m)
    M.macro_structure_description_lengths = [Structures.MDL.universal_integer(M.n) + Structures.MDL.log2(M.maximum_m_loopy), 1.]
    M.structure_type_description_length = 0
    M.description_lengths_over_time = []
    M.star_spoke_edges = Set()
    M.hub_spoke_edges = Set()
    M.zero_edges = Set()
    M.last_star_spoke_edges = Set()
    M.last_hub_spoke_edges = Set()
    M.last_zero_edges = Set()
    M.covered_nodes = Set()
    update_variable_terms!(M)
    M.backup_optimization_result = deepcopy(M.last_optimization_result)
    M.start_description_length_maxent = M.total_description_length_maxent
    # the below is to account for the fact that we have two initial constraints but perform only one optimization
    append!(M.description_lengths_over_time, [M.start_description_length_maxent])
    M.parameters = Dict("surprise_threshold" => Inf, "size_threshold" => Inf)
    return M
end

"""
initial graph setup: set constant terms - NB: dataset now needs to be full path with extension
"""
function set_constant_terms!(M::GraphModel, dataset::String, io::IO)
    M.dataset = dataset
    M.structure_alphabet_size = 4 # cliques, stars, bicliques, starcliques
    M.io = io
    M.G = SimpleGraph(NodeMap.loadgraph(dataset, GraphIO.EdgeList.EdgeListFormat())[1])
    NodeMap.create_node_map(dataset)
    remove_loops!(M.G)
    M.A = LightGraphs.LinAlg.adjacency_matrix(M.G)
    M.n = nv(M.G)
    M.m = ne(M.G)
    M.maximum_m_simple = Structures.MDL.choose(M.n, 2)
    M.maximum_m_loopy = M.maximum_m_simple + M.n
end

"""
make graph loop-free
"""
function remove_loops!(G::SimpleGraph)
    for e in [e for e in edges(G) if src(e) == dst(e)]
        rem_edge!(G, e)
    end
end

"""
each step: update model (re-optimize, update description length, update surprises)
"""
function update_variable_terms!(M::GraphModel)
    log_progress!(M, "starting optimization")
    optimize_maxent!(M)
    log_progress!(M, "finished optimization")
    log_progress!(M, "updating description length")
    update_description_length!(M)
end

"""
optimization wrapper function
"""
function optimize_maxent!(M::Model.GraphModel)
    n_iter = 50
    td = TwiceDifferentiable(x -> lagrange_dual_equiv(x, M), M.lambdas; autodiff = :forward)
        result = optimize(
        td, M.lambdas,
        Newton(linesearch=LineSearches.BackTracking()), # the solver - using BackTracking instead of the default HagerZhang because that sometimes throws a weird assertion error
        Optim.Options(allow_f_increases=true,
            show_trace=false, show_every=100, iterations=n_iter, x_tol=1e-8, f_tol=1e-8, g_tol=1e-6)
    )
    M.last_optimization_result = result
    M.lambdas = M.last_optimization_result.minimizer
    M.entropy = M.last_optimization_result.minimum
    if !converged(result)
        log_progress!(M, "NO CONVERGENCE AFTER $(n_iter) ITERATIONS")
    end
end

"""
optimization computation: wrapper
"""
function lagrange_dual_equiv(lambdas, M::Model.GraphModel)
    return (get_first_term(lambdas, M) - get_second_term(lambdas, M))
end

"""
optimization computation: first term
"""
function get_first_term(lambdas, M::Model.GraphModel)
    #log_progress!(M, "computing first term in optimization")
    rest_lambda = lambdas[1]
    if isempty(M.edge_equivalence_classes)
        return M.maximum_m_loopy * log(get_normalizer(rest_lambda))
    else
        lambda_total_contribution = BigFloat(0)
        lambda_total_sizes = 0
        for equivalence_class_key in keys(M.edge_equivalence_classes)
            size = get_equiv_size(equivalence_class_key, M.edge_equivalence_classes)
            lambda_sum = get_lambda_sum_from_equiv(equivalence_class_key, lambdas)
            normalized_sum = get_normalizer(lambda_sum)
            contribution = log(normalized_sum)
            lambda_total_contribution += size * contribution
            lambda_total_sizes += size
        end
        rest_class_n = M.maximum_m_loopy - sum(lambda_total_sizes)
        rest_contribution = rest_class_n * log(get_normalizer(rest_lambda))
        return lambda_total_contribution + rest_contribution
    end
end

"""
optimization helper (normalization)
"""
function get_normalizer(lambda_sum)
    return (1 + exp(lambda_sum))
end

"""
optimization helper (size of equivalence class)
"""
function get_equiv_size(equivalence_class_key::Set{Int64}, equivalence_classes::Dict{Set{Int64},Array{Tuple{Int64,Int64},1}})
    return length(equivalence_classes[equivalence_class_key])
end

"""
optimization helper (lambda sum for equivalence class)
"""
function get_lambda_sum_from_equiv(equivalence_class_key::Set{Int64}, lambdas)
    lambda_sum = BigFloat(0)
    for idx in equivalence_class_key
        lambda_sum += lambdas[idx]
    end
    return lambda_sum 
end

"""
optimization computation: second term
"""
function get_second_term(lambdas, M::Model.GraphModel)
    #log_progress!(M, "computing second term in optimization")
    second_term = BigFloat(0)
    for (lambda, constraint_value) in zip(lambdas, M.constraint_values)
        second_term += lambda * constraint_value
    end
    return second_term
end

"""
description length computation wrapper
"""
function update_description_length!(M::Model.GraphModel)
    update_equivalence_class_edge_probabilities!(M)
    update_edge_probability_description_length!(M)
    update_structure_type_description_length!(M)
    update_total_description_length!(M)
    # NB: the below gives misleading results when multiple structures are added at a time, e.g., during replay -> we have taken care of that in restore_model
    append!(M.description_lengths_over_time, M.total_description_length_maxent)
end

"""
description length helper: structure type encoding costs
"""
function update_structure_type_description_length!(M::Model.GraphModel)
    n_different_structures = length(M.macro_structures) - 2 # we don't count the generic macro structures as their type is equivalence_class_edge_probabilities
    M.structure_type_description_length = (Structures.MDL.universal_integer(n_different_structures + 1) + Structures.MDL.log2_choose(M.structure_alphabet_size + n_different_structures - 1, n_different_structures - 1)
    )
end

"""
description length helper: edge probabilities for equivalence classes
"""
function update_equivalence_class_edge_probabilities!(M::Model.GraphModel)
    probabilities = Dict{Array{Int64,1},BigFloat}(
        equivalence_class=>edge_probability_for_equiv(equivalence_class, M.lambdas) for equivalence_class in keys(M.edge_equivalence_classes))
    probabilities[[1,]] = edge_probability_for_equiv(Set([1,]), M.lambdas)
    M.equivalence_class_edge_probabilities = probabilities
end

"""
description length helper: edge probability for one equivalence class
"""
function edge_probability_for_equiv(equivalence_class, lambdas)
    if equivalence_class == Set([1,])
        lambda_sum = lambdas[1]
    else
        lambda_sum = get_lambda_sum_from_equiv(equivalence_class, lambdas)
    end
    edge_prop = exp(lambda_sum) / get_normalizer(lambda_sum)
    return edge_prop
end

"""
description length helper: surprise matrices update
"""
function update_surprise_matrices!(M::Model.GraphModel)
    # NB: we currently keep track only of surprising EDGES (else we get a non-sparse matrix :( )
    surprise_matrix = SparseMatrixCSC{BigFloat,Int64}(M.A)
    for idx in filter!(x->(x[1] < x[2]), findall(isone, M.A))
        edge = (idx[1],idx[2])
        # due to symmetry
        eq_class = Set(get(M.edge_lambda_indices, edge, [1,]))
        rest_class_probability = M.equivalence_class_edge_probabilities[Set([1,])]
        if eq_class == Set([1,])
            surprise_matrix[edge...] = 1 - rest_class_probability
        else
            surprise_matrix[edge...] = 0
        end
        surprise_matrix[reverse(edge)...] = 0
    end
    M.surprise_matrix = surprise_matrix
    M.surprise_matrix_transpose = sparse(transpose(M.surprise_matrix))
end

"""
description length helper: total description length of edge probabilities (error in MDL, L(D|M))
"""
function update_edge_probability_description_length!(M::Model.GraphModel)
    dl = 0.
    rest_size = M.maximum_m_loopy # because we have loop-freeness as a separate constraint
    n_guaranteed_edges = length(M.hub_spoke_edges) # guaranteed to be there and not appearing in any equivalence class
    n_zero_edges = length(M.zero_edges) # guaranteed to be not there and not appearing in any equivalence class
    rest_size -=  (n_guaranteed_edges + n_zero_edges)
    rest_edges = M.m - n_guaranteed_edges
    for (eq_class, edges) in M.edge_equivalence_classes
        size = length(edges)
        size_ones = length(filter(edge -> has_edge(M.G, edge...), edges))
        size_zeros = size - size_ones
        eq_class_probability = M.equivalence_class_edge_probabilities[eq_class]
        contribution = size_ones * -Structures.MDL.log2_zero(eq_class_probability) + size_zeros * -Structures.MDL.log2_zero(1 - eq_class_probability)
        dl += contribution
        rest_size -= size
        rest_edges -= size_ones
    end
    rest_probability = M.equivalence_class_edge_probabilities[Set([1,])]
    rest_contribution = rest_edges * -Structures.MDL.log2_zero(rest_probability) + (rest_size - rest_edges) * -Structures.MDL.log2_zero(1 - rest_probability)
    dl += rest_contribution
    M.edge_probability_description_length_maxent = dl
end

"""
description length helper: total description length (L(M) + L(D|M))
"""
function update_total_description_length!(M::Model.GraphModel)
    M.total_description_length_maxent = sum(M.macro_structure_description_lengths) + M.edge_probability_description_length_maxent + M.structure_type_description_length
end

"""
sanity check whether probability mass approximately sums to 1
"""
function sanity_check(M::Model.GraphModel)
    if !isempty(M.edge_equivalence_classes)
        equiv = sum([edge_probability_for_equiv(equivalence_class, M.lambdas) * get_equiv_size(equivalence_class, M.edge_equivalence_classes)
            for equivalence_class in keys(M.edge_equivalence_classes)])
        rest_size = M.maximum_m_loopy - sum([get_equiv_size(equivalence_class, M.edge_equivalence_classes)
            for equivalence_class in keys(M.edge_equivalence_classes)])
    else
        equiv = 0
        rest_size = M.maximum_m_loopy
    end
    rest = edge_probability_for_equiv([1,], M.lambdas)
    probability_mass_sum = equiv + rest * rest_size
    print("Sanity Check: probability mass sums to $(probability_mass_sum), expected $(M.m)\n")
    return abs(probability_mass_sum-M.m) <= 1.0e-2
end

"""
generate micro structures and constraints from clique candidate
"""
function generate_micro_structures(structure::Structures.Clique, M::Model.GraphModel)
    micro_structures = [Set([(n1,n2) for n1 in structure.nodes for n2 in structure.nodes if n1 < n2])]
    clique_edges_in_hs_edges = intersect(micro_structures[1], M.hub_spoke_edges)
    setdiff!(micro_structures[1], clique_edges_in_hs_edges)
    constraint_values = [structure.n_edges_total - length(clique_edges_in_hs_edges)]
    return micro_structures, constraint_values
end

"""
generate micro structures and constraints from star candidate
"""
function generate_micro_structures(structure::Structures.Star, M::Model.GraphModel)
    hub = structure.hub 
    micro_structures = [Set([(hub < n ? (hub,n) : (n, hub)) for n in structure.spokes]),
        Set([(n1,n2) for n1 in structure.spokes for n2 in structure.spokes if (n1 < n2) && (n1,n2) ∉ M.star_spoke_edges])]
    hs_edges_in_hs_edges = intersect(micro_structures[1], M.hub_spoke_edges)
    setdiff!(micro_structures[1],hs_edges_in_hs_edges)
    ss_edges_in_hs_edges = intersect(micro_structures[end], M.hub_spoke_edges)
    setdiff!(micro_structures[end],ss_edges_in_hs_edges)
    constraint_values = [structure.n_nodes_in_spokes-length(hs_edges_in_hs_edges), structure.n_edges_in_spokes-length(ss_edges_in_hs_edges)] # note that this is theoretically wrong - we do it for performance reasons
    return micro_structures, constraint_values
end

"""
generate micro structures and constraints from biclique candidate
"""
function generate_micro_structures(structure::Union{Structures.Biclique,Structures.Starclique}, M::Model.GraphModel)
    if structure isa Structures.Biclique
        left = Set([(n1, n2) for n1 in structure.left_nodes for n2 in structure.left_nodes if n1 < n2 && (n1,n2) ∉ M.star_spoke_edges])
    else # assuming Starclique
        left = Set([(n1, n2) for n1 in structure.left_nodes for n2 in structure.left_nodes if n1 < n2])
    end
    left_edges_in_hs_edges = intersect(left, M.hub_spoke_edges)
    setdiff!(left,left_edges_in_hs_edges)
    right = Set([(n1, n2) for n1 in structure.right_nodes for n2 in structure.right_nodes if n1 < n2 && (n1,n2) ∉ M.star_spoke_edges])
    right_edges_in_hs_edges = intersect(right, M.hub_spoke_edges)
    setdiff!(right,right_edges_in_hs_edges)
    across = Set([(n1 < n2) ? (n1, n2) : (n2, n1) for n1 in structure.left_nodes for n2 in structure.right_nodes])
    across_edges_in_hs_edges = intersect(across, M.hub_spoke_edges)
    setdiff!(across, across_edges_in_hs_edges)
    micro_structures = [left, right, across]
    constraint_values = [   structure.n_edges_in_left - length(left_edges_in_hs_edges), 
                            structure.n_edges_in_right - length(right_edges_in_hs_edges), 
                            structure.n_edges_across - length(across_edges_in_hs_edges)
                        ]
    return micro_structures, constraint_values
end

"""
try adding structure to model, keeping the result only if it reduces the description length
"""
function test_adding_structure!(M::Model.GraphModel, structure::Structures.Structure, size_threshold::Float64=1.)
    log_progress!(M, "saving previous state...")
    M.backup_optimization_result = deepcopy(M.last_optimization_result)
    log_progress!(M, "testing structure $(structure)...")
    if structure.n_nodes_total < size_threshold
        log_progress!(M, "failure: structure too small")
        return false
    end
    add_macro_structures!(M, Array{Structures.Structure,1}([structure]))
    if converged(M.last_optimization_result) && M.description_lengths_over_time[end] <= M.description_lengths_over_time[end-1]
        log_progress!(M, "success: new description length $(M.total_description_length_maxent)")
        return true
    else
        roll_back_model!(M)
        log_progress!(M, "failure: no decrease in description length")
        return false
    end
end

"""
add macro structures to model (assuming no connector knowledge)
we allow for adding multiple structures at once to enable faster replay of macro_structure_description_lengths
(in the model building phase, we test adding one structure at a time only)
"""
function add_macro_structures!(M::Model.GraphModel, structures::Array{Structures.Structure,1}; lambdas=[])
    log_progress!(M, "updating micro structures et al")
    for structure in structures
        micro_structures, constraint_values = generate_micro_structures(structure, M)
        if structure isa Structures.Star
            M.last_hub_spoke_edges = deepcopy(M.hub_spoke_edges)
            M.last_star_spoke_edges = deepcopy(M.star_spoke_edges)
            n_edges_in_spokes = structure.n_edges_in_spokes
            union!(M.star_spoke_edges, micro_structures[end])
            union!(M.hub_spoke_edges, micro_structures[1])
            if n_edges_in_spokes > 0 # if the star is perfect, we know exactly what its associated probabilities are
                append!(M.micro_structures, [micro_structures[end]])
                append!(M.constraint_values, [constraint_values[end]])
                append!(M.lambdas, zeros(1))
                update_edge_lambda_indices!(M, [micro_structures[end]], [constraint_values[end]])
                update_equivalence_classes_from_edge_lambda_indices!(M)
            else
                M.last_zero_edges = deepcopy(M.last_zero_edges)
                union!(M.zero_edges, micro_structures[1])
            end
        elseif structure isa Structures.Starclique || structure isa Structures.Biclique
            M.last_star_spoke_edges = deepcopy(M.star_spoke_edges)
            union!(M.star_spoke_edges, micro_structures[2])
            if structure isa Structures.Biclique
                union!(M.star_spoke_edges, micro_structures[1])
            end
            append!(M.micro_structures, micro_structures)
            append!(M.constraint_values, constraint_values)
            append!(M.lambdas, zeros(length(micro_structures)))
            update_edge_lambda_indices!(M, micro_structures, constraint_values)
            update_equivalence_classes_from_edge_lambda_indices!(M)
        else
            append!(M.micro_structures, micro_structures)
            append!(M.constraint_values, constraint_values)
            append!(M.lambdas, zeros(length(micro_structures)))
            update_edge_lambda_indices!(M, micro_structures, constraint_values)
            update_equivalence_classes_from_edge_lambda_indices!(M)
        end  
    end
    if !isempty(lambdas)
        M.lambdas = lambdas
    end
    update_variable_terms!(M)
    for structure in structures
        append!(M.macro_structures, [structure])
        append!(M.macro_structure_description_lengths, [Structures.compute_description_length(structure, M.n, false)])
    end
end

function update_edge_lambda_indices!(M::Model.GraphModel, new_micro_structures, new_constraint_values)
    log_progress!(M, "updating edge lambda indices")
    start_lambda_index = length(M.lambdas) - length(new_constraint_values)
    for (idx, structure) in enumerate(new_micro_structures)
        new_idx = start_lambda_index+idx
        present_edges = intersect(keys(M.edge_lambda_indices), structure)
        for edge in present_edges
            append!(M.edge_lambda_indices[edge], new_idx)
        end
        absent_edge_dict = Dict(e=>[1,new_idx] for e in setdiff(structure,present_edges))
        merge!(M.edge_lambda_indices,absent_edge_dict)
    end
end

"""
roll back the model if the added structure is no good
eliminates last macro structure added and all its consequences
"""
function roll_back_model!(M::Model.GraphModel)
    structure_to_remove = M.macro_structures[end]
    micros_to_remove = determine_micros_to_remove(structure_to_remove)
    if structure_to_remove isa Structures.Star
        M.hub_spoke_edges = deepcopy(M.last_hub_spoke_edges)
        M.star_spoke_edges = deepcopy(M.last_star_spoke_edges)
        if structure_to_remove.n_edges_in_spokes == 0
            micros_to_remove = 0
            M.last_zero_edges = deepcopy(M.last_zero_edges)
        end
    elseif structure_to_remove isa Structures.Starclique || structure_to_remove isa Structures.Biclique
        M.star_spoke_edges = deepcopy(M.last_star_spoke_edges)
    end
    if micros_to_remove > 0
        roll_back_edge_lambda_indices!(M, micros_to_remove)
        roll_back_lambda_related!(M, micros_to_remove)
        update_equivalence_classes_from_edge_lambda_indices!(M)
    end
    update_equivalence_class_edge_probabilities!(M)
    update_edge_probability_description_length!(M)
    pop!(M.macro_structures)
    pop!(M.macro_structure_description_lengths)
    pop!(M.description_lengths_over_time)
    roll_back_total_description_length!(M)
    roll_back_optimization_result!(M)
end

function determine_micros_to_remove(structure_to_remove::Structures.Structure)
    if structure_to_remove isa Structures.Clique || structure_to_remove isa Structures.Star
        micros_to_remove = 1
    elseif structure_to_remove isa Structures.Biclique || structure_to_remove isa Structures.Starclique
        micros_to_remove = 3
    else
        throw(ErrorException("Cannot roll back for structure of type $(typeof(structure_to_remove))!"))
    end
    return micros_to_remove
end

function roll_back_edge_lambda_indices!(M::Model.GraphModel, micros_to_remove::Int64)
    for micro_structure in M.micro_structures[end-micros_to_remove+1:end]
        for edge in micro_structure
            pop!(M.edge_lambda_indices[edge])
            if M.edge_lambda_indices[edge] == [1]
                delete!(M.edge_lambda_indices, edge)
            end
        end
    end
end

function roll_back_lambda_related!(M::Model.GraphModel, micros_to_remove::Int64)
    first_to_delete = length(M.micro_structures)-micros_to_remove+1
    last_to_delete = length(M.micro_structures)
    deleteat!(M.lambdas, first_to_delete:last_to_delete)
    deleteat!(M.constraint_values, first_to_delete:last_to_delete)
    deleteat!(M.micro_structures, first_to_delete:last_to_delete)
end

function update_equivalence_classes_from_edge_lambda_indices!(M::Model.GraphModel)
    equivalence_classes = Dict{Set{Int64},Array{Tuple{Int64,Int64},1}}()
    for (k,v) in M.edge_lambda_indices
        append!(get!(equivalence_classes, Set(v), Array{Int64,1}()), [k])
    end
    M.edge_equivalence_classes = equivalence_classes
end

function update_equivalence_class_edge_probabilities!(M::Model.GraphModel)
    probabilities = Dict{Set{Int64},BigFloat}(
        equivalence_class=>edge_probability_for_equiv(equivalence_class, M.lambdas) for equivalence_class in keys(M.edge_equivalence_classes))
    probabilities[Set([1,])] = edge_probability_for_equiv(Set([1,]), M.lambdas)
    M.equivalence_class_edge_probabilities = probabilities
end

function roll_back_total_description_length!(M::Model.GraphModel)
    M.total_description_length_maxent = M.description_lengths_over_time[end]
end

function roll_back_optimization_result!(M::Model.GraphModel)
    M.last_optimization_result = deepcopy(M.backup_optimization_result)
    M.entropy = M.last_optimization_result.minimum
end

end