module Structures

include("mdl.jl")

mutable struct GenericStructure
    structure_type::String
    structure_role::String
    n_edges_total::Int64
    function GenericStructure(structure_role::String, n_edges_total::Int64)
        return new("generic", structure_role, n_edges_total)
    end
end

mutable struct Clique
    structure_type::String
    n_nodes_total::Int64
    n_edges_total::Int64
    nodes::Array{Int64,1}
    function Clique(n_edges_total::Int64, nodes::Array{Int64,1})
        length_nodes = length(nodes)
        @assert n_edges_total <= (length_nodes * (length_nodes - 1)) / 2 "Impossible number of edges! nv=$(length_nodes), ne=$(n_edges_total)"
        @assert length_nodes == length(unique(nodes)) "Input nodes must be unique!"
        return new("clique", length_nodes, n_edges_total, nodes)
    end
end

mutable struct Star
    structure_type::String
    n_nodes_total::Int64
    n_edges_total::Int64
    n_nodes_in_spokes::Int64
    n_edges_in_spokes::Int64
    hub::Int64
    spokes::Array{Int64,1}
    function Star(n_edges_total::Int64, hub::Int64, spokes::Array{Int64,1})
        length_spokes = length(spokes)
        @assert n_edges_total >= length_spokes "Too few edges!"
        @assert n_edges_total <= length_spokes + length_spokes * (length_spokes - 1) / 2 "Too many edges!"
        @assert length_spokes + 1 == length(unique([spokes...,hub])) "Input nodes must be unique!"
        return new("star", length_spokes + 1, n_edges_total, length_spokes, n_edges_total - length_spokes, hub, spokes)
    end
end

mutable struct Biclique
    structure_type::String
    n_nodes_total::Int64
    n_edges_total::Int64
    n_nodes_in_left::Int64
    n_nodes_in_right::Int64
    n_edges_in_left::Int64
    n_edges_in_right::Int64
    n_edges_across::Int64
    left_nodes::Array{Int64,1}
    right_nodes::Array{Int64,1}
    function Biclique(n_edges_total::Int64, left_nodes::Array{Int64,1}, right_nodes::Array{Int64,1}, n_edges_in_left, n_edges_in_right)
        left_length = length(left_nodes)
        right_length = length(right_nodes)
        n_edges_across = n_edges_total - n_edges_in_left - n_edges_in_right
        max_edges_in_left = left_length * (left_length - 1) / 2
        max_edges_in_right = right_length * (right_length - 1) / 2
        max_edges_across = left_length * right_length
        @assert n_edges_in_left <= max_edges_in_left "Too many edges in left!"
        @assert n_edges_in_right <= max_edges_in_right "Too many edges in right!"
        @assert n_edges_across <= max_edges_across "Too many edges in across!"
        @assert left_length + right_length == length(unique([left_nodes..., right_nodes...])) "Input nodes must be unique!"
        return new("biclique", left_length + right_length, n_edges_total, left_length, right_length,
            n_edges_in_left, n_edges_in_right, n_edges_across, left_nodes, right_nodes)
    end
end

mutable struct Starclique
    structure_type::String
    n_nodes_total::Int64
    n_edges_total::Int64
    n_nodes_in_left::Int64
    n_nodes_in_right::Int64
    n_edges_in_left::Int64
    n_edges_in_right::Int64
    n_edges_across::Int64
    left_nodes::Array{Int64,1}
    right_nodes::Array{Int64,1}
    function Starclique(n_edges_total::Int64, left_nodes::Array{Int64,1}, right_nodes::Array{Int64,1}, n_edges_in_left, n_edges_in_right)
        left_length = length(left_nodes)
        right_length = length(right_nodes)
        n_edges_across = n_edges_total - n_edges_in_left - n_edges_in_right
        max_edges_in_left = left_length * (left_length - 1) / 2
        max_edges_in_right = right_length * (right_length - 1) / 2
        max_edges_across = left_length * right_length
        @assert n_edges_in_left <= max_edges_in_left "Too many edges in left!"
        @assert n_edges_in_right <= max_edges_in_right "Too many edges in right!"
        @assert n_edges_across <= max_edges_across "Too many edges in across!"
        @assert left_length + right_length == length(unique([left_nodes..., right_nodes...])) "Input nodes must be unique!"
        return new("starclique", left_length + right_length, n_edges_total, left_length, right_length,
            n_edges_in_left, n_edges_in_right, n_edges_across, left_nodes, right_nodes)
    end
end

Structure = Union{GenericStructure, Clique, Star, Biclique, Starclique}
GenericBipartite = Union{Biclique, Starclique}

Base.show(io::IO, GS::GenericStructure) = print(io, "Generic Structure: $(GS.structure_role) with $(GS.n_edges_total) edges")
Base.show(io::IO, C::Clique) = print(io, "Clique with $(C.n_nodes_total) nodes and $(C.n_edges_total) edges")
Base.show(io::IO, S::Star) = print(io, "Star with $(S.n_nodes_in_spokes) spokes and $(S.n_edges_in_spokes) edges between spokes")
Base.show(io::IO, B::Biclique) = print(io, "Biclique with $(B.n_nodes_in_left) left nodes, $(B.n_nodes_in_right) right nodes, $(B.n_edges_in_left) edges in left, $(B.n_edges_in_right) edges in right, and $(B.n_edges_across) edges across")
Base.show(io::IO, B::Starclique) = print(io, "Starclique with $(B.n_nodes_in_left) left nodes, $(B.n_nodes_in_right) right nodes, $(B.n_edges_in_left) edges in left, $(B.n_edges_in_right) edges in right, and $(B.n_edges_across) edges across")

"""
generic structure description length
"""
function compute_description_length(structure::GenericStructure, maximum_n_loopy::Int64)
    @assert structure.structure_role in ["Overall number of edges", "Loop-freeness"] "Unexpected generic structure with no proper description length computation!"
    if structure.structure_role == "Loop-freeness"
        return 1.
    else
        return MDL.log2_zero(maximum_n_loopy)
    end
end

function compute_rest(clique::Clique)
    max_n_edges = convert(Int64, MDL.choose(clique.n_nodes_total, 2))
    cost_sparse_or_dense = 1
    cost_n_edges = MDL.log2_zero(MDL.log2_zero(floor(max_n_edges/2))) + MDL.log2_zero(min(clique.n_edges_total, max_n_edges - clique.n_edges_total))
    #cost_n_edges = MDL.universal_integer(min(clique.n_edges_total, max_n_edges - clique.n_edges_total) + 1)
    return cost_sparse_or_dense + cost_n_edges
end

### CLIQUES #######################################################################################################################################

"""
clique description length without connector or with path connector
(note that total_n_remaining is computed differently depending on the connector)
"""
function compute_description_length(clique::Clique, total_n_remaining::Int64, path::Bool; with_node_ids_cost::Bool=true)
    cost_n_nodes = path ? MDL.universal_integer(clique.n_nodes_total - 1) : MDL.universal_integer(clique.n_nodes_total)
    if with_node_ids_cost
        cost_node_ids = path ? MDL.log2_choose(total_n_remaining, clique.n_nodes_total - 1) : MDL.log2_choose(total_n_remaining, clique.n_nodes_total)
    else
        cost_node_ids = 0
    end
    rest = compute_rest(clique)
    return cost_n_nodes + cost_node_ids + rest
end

"""
clique description length with overlap connector
(note that total_n_remaining is computed differently depending on the connector)
"""
function compute_description_length(clique::Clique, parent_size::Int64, overlap_size::Int64, total_n_remaining::Int64)
    cost_n_nodes = MDL.universal_integer(clique.n_nodes_total - overlap_size)
    cost_node_ids = MDL.log2_choose(parent_size, overlap_size) + MDL.log2_choose(total_n_remaining, clique.n_nodes_total - overlap_size)
    rest = compute_rest(clique)
    return cost_n_nodes + cost_node_ids + rest
end

### STARS #########################################################################################################################################

"""
helper to remove redundancy in star description length computation
"""
function compute_rest(star::Star)
    cost_of_hub = MDL.log2_zero(star.n_nodes_total)
    cost_n_edges = MDL.log2_zero(MDL.log2_zero((star.n_nodes_in_spokes*(star.n_nodes_in_spokes-1)/2))) + MDL.log2_zero(star.n_edges_in_spokes)
    #cost_n_edges = MDL.universal_integer(star.n_edges_in_spokes + 1)
    return cost_of_hub + cost_n_edges
end

"""
star description length without connector or with path connector
(note that total_n_remaining is computed differently depending on the connector)
"""
function compute_description_length(star::Star, total_n_remaining::Int64, path::Bool; with_node_ids_cost::Bool=true)
    cost_n_spokes = path ? MDL.universal_integer(star.n_nodes_in_spokes - 1) : MDL.universal_integer(star.n_nodes_in_spokes)
    if with_node_ids_cost
        cost_node_ids = path ? MDL.log2_choose(total_n_remaining, star.n_nodes_total - 1) : MDL.log2_choose(total_n_remaining, star.n_nodes_total)
    else
        cost_node_ids = 0
    end
    rest = compute_rest(star)
    return cost_n_spokes + cost_node_ids + rest
end

"""
star description length with overlap connector
(note that total_n_remaining is computed differently depending on the connector)
"""
function compute_description_length(star::Star, parent_size::Int64, overlap_size::Int64, total_n_remaining::Int64)
    cost_n_spokes = MDL.universal_integer(star.n_nodes_in_spokes)
    cost_node_ids = MDL.log2_choose(parent_size, overlap_size) + MDL.log2_choose(total_n_remaining, star.n_nodes_total - overlap_size)
    rest = compute_rest(star)
    return cost_n_spokes + cost_node_ids + rest
end

### BICLIQUES AND STARCLIQUES #####################################################################################################################


"""
helper to remove redundancy in biclique description length computation
"""
function compute_rest(biclique::GenericBipartite)
    cost_n_left = MDL.log2_zero(biclique.n_nodes_total)
    cost_left_ids = MDL.log2_choose(biclique.n_nodes_total, biclique.n_nodes_in_left)
    # across, we definitely want density
    max_across_edges = biclique.n_nodes_in_left * biclique.n_nodes_in_right
    cost_n_across = MDL.log2_zero(MDL.log2_zero(max_across_edges)) + MDL.log2_zero(biclique.n_edges_across)
    max_left_edges = convert(Int64, MDL.choose(biclique.n_nodes_in_left, 2))
    max_right_edges = convert(Int64, MDL.choose(biclique.n_nodes_in_right, 2))
    # left we can have sparsity or density
    cost_n_left_edges = MDL.log2_zero(MDL.log2_zero(floor(max_left_edges/2))) + MDL.log2_zero(min(biclique.n_edges_in_left, max_left_edges - biclique.n_edges_in_left))
    # right we want sparsity
    cost_n_right_edges = MDL.log2_zero(MDL.log2_zero(floor(max_right_edges/2))) + MDL.log2_zero(biclique.n_edges_in_right)
    return (cost_n_left + cost_left_ids
            + cost_n_across + cost_n_left_edges + cost_n_right_edges)
end

"""
biclique description length without connector or with path connector
(note that total_n_remaining is computed differently depending on the connector)
"""
function compute_description_length(biclique::GenericBipartite, total_n_remaining::Int64, path::Bool; with_node_ids_cost::Bool=true)
    cost_n_nodes = path ? MDL.universal_integer(biclique.n_nodes_total - 1) : MDL.universal_integer(biclique.n_nodes_total)
    if with_node_ids_cost
        cost_node_ids = path ? MDL.log2_choose(total_n_remaining, biclique.n_nodes_total - 1) : MDL.log2_choose(total_n_remaining, biclique.n_nodes_total)
    else
        cost_node_ids = 0
    end
    rest = compute_rest(biclique)
    return (cost_n_nodes + cost_node_ids + rest)
end

"""
biclique description length with overlap connector
(note that total_n_remaining is computed differently depending on the connector)
"""
function compute_description_length(biclique::GenericBipartite, parent_size::Int64, overlap_size::Int64, total_n_remaining::Int64)
    cost_n_nodes = MDL.universal_integer(biclique.n_nodes_total - overlap_size)
    cost_node_ids = MDL.log2_choose(parent_size, overlap_size) + MDL.log2_choose(total_n_remaining, biclique.n_nodes_total - overlap_size)
    rest = compute_rest(biclique)
    return (cost_n_nodes + cost_node_ids + rest)
end

"""
get all nodes from a clique
"""
function get_node_set(clique::Clique)
    return clique.nodes
end

"""
get all nodes from a star
"""
function get_node_set(star::Star)
    return [star.hub, star.spokes...] 
end

"""
get all nodes from a biclique
"""
function get_node_set(biclique::GenericBipartite)
    return [biclique.left_nodes..., biclique.right_nodes...]
end

end