module Mining

ccall(:jl_exit_on_sigint, Cvoid, (Cint,), 0)

include("model.jl")

using LightGraphs
using StatsBase


"""
Decompose graph into components to serve as seeds for structure candidates
"""
function decompose_graph(G, minimum_component_size::Int64=10, neighborhood_size::Int64=1, return_slashed::Bool=false)
    components = []
    slashed = []

    ccs = sort!(connected_components(G), by=length, rev=true)
    lcc = ccs[1]
    lccG, vmap = induced_subgraph(G, lcc)
    largest_degree = Δ(lccG)

    while length(lcc) >= minimum_component_size && largest_degree >= minimum_component_size
        new_hub = sort!(filter(x->x[end] == largest_degree, collect(zip(vertices(lccG),degree(lccG)))), by=x->x[1])[1][1]
        new_hub_original = vmap[new_hub]
        push!(slashed, new_hub_original)
        hub_neighbors = neighborhood(G, new_hub_original, neighborhood_size)
        push!(components, hub_neighbors)
        for edge in filter(x -> src(x) == new_hub_original || dst(x) == new_hub_original, collect(edges(G)))
            rem_edge!(G, edge)
        end
        ccs = sort!(connected_components(G), by=length, rev=true)
        lcc = ccs[1]
        lccG, vmap = induced_subgraph(G, lcc)
        largest_degree = Δ(lccG)
    end

    filter!(x->length(x) >= minimum_component_size, ccs)
    append!(components, ccs)
    append!(slashed, repeat([nothing],length(ccs)))
    if !return_slashed
        return components
    else
        return components, slashed
    end
end


"""
Create a star candidate from a component
"""
function find_star_from_component(component, G)
    cG, vmap = induced_subgraph(G, component)
    largest_degree = Δ(cG)
    top_nodes = sort!(filter(x->x[end] == largest_degree, collect(zip(vertices(cG),degree(cG)))), by=x->x[1])
    if length(top_nodes) > 1
        return nothing, [], nothing
    else
        new_hub = top_nodes[1][1]
        new_hub_original = vmap[new_hub]
        sG, vmap = induced_subgraph(G, neighborhood(G, new_hub_original, 1))
        nodes_original = [vmap[x] for x in vertices(sG)]
        n_spokes = nv(sG) - 1
        fraction = 0.1
        candidates = sort!(filter(x -> n_spokes - 1 > x[end] - 1 > 0.05 * n_spokes, collect(zip(vertices(sG),degree(sG)))), by=x->x[end], rev=true)
        cut_point = ceil(fraction * length(candidates))
        nodes_to_cut = candidates[1:convert(Int, cut_point)]
        while !isempty(nodes_to_cut) && nv(sG) >= 10
            nodes_to_cut_original = [vmap[x[1]] for x in nodes_to_cut]
            nodes_to_keep = filter(x -> x ∉ nodes_to_cut_original, nodes_original)
            sG, vmap = induced_subgraph(G, nodes_to_keep)
            nodes_original = [vmap[x] for x in vertices(sG)]
            n_spokes = nv(sG) - 1
            fraction = min(fraction + 0.01, 1.0)
            candidates = sort!(filter(x -> n_spokes - 1 > x[end] - 1 > 0.05 * n_spokes, collect(zip(vertices(sG),degree(sG)))), by=x->x[end], rev=true)
            cut_point = ceil(fraction * length(candidates))
            nodes_to_cut = candidates[1:convert(Int, cut_point)]
        end
        if nv(sG) < 10
            return nothing, [], nothing
        else
            return sG, [vmap[x] for x in vertices(sG) if vmap[x] != new_hub_original], new_hub_original
        end
    end
end


"""
Create a clique candidate from a component
"""
function find_clique_from_component(component, G, min_connectivity_fraction=0.5)
    cG, vmap = induced_subgraph(G, component)

    maxcliques = sort!(maximal_cliques(cG), by=x->length(x), rev=true)
    nodes_original = [vmap[x] for x in maxcliques[1]]
    cG, vmap = induced_subgraph(G, nodes_original)

    while true
        cutoff = min_connectivity_fraction * length(nodes_original)
        nodes_to_add = sort!([x[1] for x in filter(x -> x[end] >= cutoff && x[1] ∉ nodes_original, 
                                collect(countmap(vcat([neighbors(G, x) for x in nodes_original]...))))], by=x->degree(G,x), rev=true)
        if !isempty(nodes_to_add)
            node_to_add = nodes_to_add[1]
        else
            break
        end
        union!(nodes_original, [node_to_add])
        cG, vmap = induced_subgraph(G, nodes_original)
    end

    if nv(cG) < 10 || !is_connected(cG)
        return nothing, []
    else
        return cG, nodes_original
    end
end


"""
Create a biclique candidate from a component
"""
function find_biclique_from_component(component, G)
    max_intra_connectivity = 0.05
    min_inter_connectivity = 0.5
    cG, vmap = induced_subgraph(G, component)
    
    right = independent_set(cG, MaximalIndependentSet(), seed=1234)
    original_right = [vmap[x] for x in right]
    sort!(original_right, by=x->degree(G, x), rev=true)
    original_right = original_right[1:min(length(original_right),5)] 
    cG, vmap = induced_subgraph(G, original_right)

    left_candidates = [x[1] for x in filter(x -> x[end] >= min_inter_connectivity * length(original_right) && x[1] ∉ original_right, collect(countmap(vcat([neighbors(G, x) for x in original_right]...))))]
    lG, vmap2 = induced_subgraph(G, left_candidates)
    left = independent_set(lG, MaximalIndependentSet(), seed=1234)
    original_left = [vmap2[x] for x in left]
    sort!(original_left, by=x->degree(G, x), rev=true)
    original_left = original_left[1:min(length(original_left),5)]
    left_to_add = [x[1] for x in filter(x -> x[end] <= max_intra_connectivity * length(original_left) && x in left_candidates, collect(countmap(vcat([neighbors(G, x) for x in original_left]...))))]
    if !isempty(left_to_add)
        union!(original_left,left_to_add)
    end
    
    if length(original_left) < 3 || length(original_right) < 5
        return nothing, [], []
    end
    
    while true
        added_left = filter(
                    x -> 
                    x[1] ∉ original_left 
                    && x[1] ∉ original_right 
                    && x[end] >= min_inter_connectivity * length(original_right) 
                    && length(intersect(neighbors(G, x[1]), original_left)) <= max_intra_connectivity * length(original_left), 
                        collect(countmap(vcat([neighbors(G, x) for x in original_right]...)))
                    )
        if !isempty(added_left)
            added_left = map(x->x[1],sort!(added_left, by=x->x[end], rev=true))
            append!(original_left, [added_left[1]])
        end

        added_right = filter(
                    x -> 
                    x[1] ∉ original_left
                    && x[1] ∉ original_right
                    && x[end] >= min_inter_connectivity * length(original_left) 
                    && length(intersect(neighbors(G, x[1]), original_right)) <= max_intra_connectivity * length(original_right), 
                        collect(countmap(vcat([neighbors(G, x) for x in original_left]...)))
                    )
        if !isempty(added_right)
            added_right = map(x->x[1],sort!(added_right, by=x->x[end], rev=true))
            append!(original_right, [added_right[1]])
        end
        
        if isempty(added_left) && isempty(added_right)
            break
        end
    end

    cG, vmap = induced_subgraph(G, union(original_left, original_right))
    isolates = [vmap[n] for n in vertices(cG) if degree(cG, n) == 0]
    setdiff!(original_left, isolates)
    setdiff!(original_right, isolates)
    cG, vmap = induced_subgraph(G, union(original_left, original_right))

    if nv(cG) < 10 || !is_connected(cG) #|| length(original_left) < 3 || length(original_right) < 3
        return nothing, [], []
    else
        return cG, original_left, original_right
    end
end


"""
Create a starclique candidate from a component
"""
function find_starclique_from_component(component, G)
    max_intra_connectivity = 0.05
    min_intra_connectivity = 0.5
    min_inter_connectivity = 0.5
    cG, vmap = induced_subgraph(G, component)
    
    maxcliques = sort!(maximal_cliques(cG), by=x->length(x), rev=true)
    original_left = [vmap[x] for x in maxcliques[1]]
    cG, vmap = induced_subgraph(G, original_left)

    right_candidates = [x[1] for x in filter(x -> x[end] >= min_inter_connectivity * length(original_left) && x[1] ∉ original_left, collect(countmap(vcat([neighbors(G, x) for x in original_left]...))))]
    rG, rvmap = induced_subgraph(G, right_candidates)
    right = independent_set(rG, MaximalIndependentSet(), seed=1234)
    original_right = [rvmap[x] for x in right]

    while true
        added_left = filter(
                    x -> 
                    x[1] ∉ original_left 
                    && x[1] ∉ original_right 
                    && x[end] >= min_inter_connectivity * length(original_right) 
                    && length(intersect(neighbors(G, x[1]), original_left)) >= min_intra_connectivity * length(original_left), 
                        collect(countmap(vcat([neighbors(G, x) for x in original_right]...)))
                    )
        if !isempty(added_left)
            added_left = map(x->x[1],sort!(added_left, by=x->x[end], rev=true))
            append!(original_left, [added_left[1]])
        end

        added_right = filter(
                    x -> 
                    x[1] ∉ original_left
                    && x[1] ∉ original_right
                    && x[end] >= min_inter_connectivity * length(original_left) 
                    && length(intersect(neighbors(G, x[1]), original_right)) <= max_intra_connectivity * length(original_right), 
                        collect(countmap(vcat([neighbors(G, x) for x in original_left]...)))
                    )
        if !isempty(added_right)
            added_right = map(x->x[1],sort!(added_right, by=x->x[end], rev=true))
            append!(original_right, [added_right[1]])
        end
        
        if isempty(added_left) && isempty(added_right)
            break
        end
    end

    cG, vmap = induced_subgraph(G, union(original_left, original_right))
    isolates = [vmap[n] for n in vertices(cG) if degree(cG, n) == 0]
    setdiff!(original_left, isolates)
    setdiff!(original_right, isolates)
    cG, vmap = induced_subgraph(G, union(original_left, original_right))
    if nv(cG) < 10 || !is_connected(cG) || length(original_left) < 3 || length(original_right) < 3
        return nothing, [], []
    else
        return cG, original_left, original_right
    end
end


"""
Helper for sorting candidates
"""
function get_structure_type_priority(structure_type)
    if structure_type == "biclique"
        return 4
    elseif structure_type == "starclique"
        return 3
    elseif structure_type == "star"
        return 2
    else # "clique"
        return 1
    end
end


"""
Model building function which does the heavy lifting
"""
function build_decomposition_based_model!(M::Model.GraphModel, size_threshold::Float64=10., stop_threshold::Float64=10., stop_after_n_structures::Float64=Inf, sort_heuristic::String="plain")
    minimum_component_size = convert(Int64, ceil(size_threshold))
    
    G = deepcopy(M.G)
    Model.log_progress!(M, "Created graph with $(nv(G)) nodes and $(ne(G)) edges...")
    components = decompose_graph(G, minimum_component_size, 1, false)
    Model.log_progress!(M, "Decomposed graph")

    stars, star_spokes_original, hubs = generate_star_component_candidates(M, minimum_component_size, components)
    Model.log_progress!(M, "Generated star candidates")
    cliques, clique_nodes_original = generate_clique_component_candidates(M, minimum_component_size, components)
    Model.log_progress!(M, "Generated clique candidates")
    bicliques, left_nodes_original, right_nodes_original = generate_biclique_component_candidates(M, minimum_component_size, components)
    Model.log_progress!(M, "Generated biclique candidates")
    starcliques, s_left_nodes_original, s_right_nodes_original = generate_starclique_component_candidates(M, minimum_component_size, components)
    Model.log_progress!(M, "Generated starclique candidates")

    star_structures = generate_all_star_structures(hubs, star_spokes_original, stars)
    clique_structures = generate_all_clique_structures(clique_nodes_original, cliques)
    biclique_structures = generate_all_biclique_structures(left_nodes_original, right_nodes_original, bicliques, M.G, "biclique")
    starclique_structures = generate_all_biclique_structures(s_left_nodes_original, s_right_nodes_original, starcliques, M.G, "starclique")
    all_structures = vcat(clique_structures, star_structures, biclique_structures, starclique_structures)
    sort!(all_structures, by=x -> (x.n_nodes_total, x.n_edges_total, get_structure_type_priority(x.structure_type)), rev=true)
    Model.log_progress!(M, "List of candidate structures, ordered by size:")
    for structure in all_structures
        Model.log_progress!(M, "$(structure)")
    end
    Model.log_progress!(M, "")

    successive_failures = 0
    total_failures = 0
    while !isempty(all_structures)
        Model.log_progress!(M, "$(length(all_structures)) structures remaining")
        structure = popfirst!(all_structures)
        if Model.test_adding_structure!(M, structure, size_threshold)
            Model.log_progress!(M, "added: $(structure)")
            Model.log_progress!(M, "number of equivalence classes: $(length(M.equivalence_class_edge_probabilities))")
            Model.log_progress!(M, "number of lambdas: $(length(M.lambdas))")
            successive_failures = 0
            union!(M.covered_nodes,Mining.Model.Structures.get_node_set(structure))
            if sort_heuristic == "uncovered_nodes"
                sort!(all_structures, by=x ->  (get_uncovered_score(x, M), x.n_nodes_total, x.n_edges_total), rev=true)
            end
        else
            successive_failures += 1
            total_failures += 1
            if total_failures >= stop_threshold
                Model.log_progress!(M, "Stopping model computation - $(total_failures) successive failures...")
                break
            end
        end
        if length(M.macro_structures) - 2 >= stop_after_n_structures
            Model.log_progress!(M, "Stopping model computation - $(stop_after_n_structures) structures found...")
            break
        end
    end
end


"""
Helper for sort heuristic 'uncovered_nodes'
"""
function get_uncovered_score(structure, M)
    nodes = Mining.Model.Structures.get_node_set(structure)
    n_uncovered_nodes = length(setdiff(nodes, M.covered_nodes))
    n_covered_nodes = length(nodes) - n_uncovered_nodes
    return n_uncovered_nodes * (n_uncovered_nodes - 1) / 2 + n_uncovered_nodes * n_covered_nodes # this can be optimistic or pessimistic depending on the structure
end


"""
Create *all* star candidates (raw)
"""
function generate_star_component_candidates(M::Model.GraphModel, minimum_component_size::Int64=10, components=[])
    stars = []
    star_spokes_original = []
    hubs = []
    for c in components
        star, original_nodes, hub = find_star_from_component(c, M.G)
        if star !== nothing
            push!(stars,star)
            push!(star_spokes_original, original_nodes)
            push!(hubs, hub)
        end
    end
    return stars, star_spokes_original, hubs
end


"""
Create *all* clique candidates (raw)
"""
function generate_clique_component_candidates(M::Model.GraphModel, minimum_component_size::Int64=10, components=[])
    cliques = []
    clique_nodes_original = []
    for c in components
        clique, original_nodes = find_clique_from_component(c, M.G)
        clique_length = length(original_nodes)
        found = false
        for (idx,other_nodes) in enumerate(clique_nodes_original)
            if length(intersect(original_nodes, other_nodes))/min(clique_length,length(other_nodes)) >= 0.9 # was: 0.95
                clique_nodes_original[idx] = union(original_nodes, other_nodes)
                cliques[idx], _ = induced_subgraph(M.G, clique_nodes_original[idx])
                found = true
                break
            end
        end
        if !found && clique !== nothing
            push!(cliques, clique)
            push!(clique_nodes_original, original_nodes)
        end
    end
    return cliques, clique_nodes_original
end


"""
Create *all* biclique candidates (raw)
"""
function generate_biclique_component_candidates(M::Model.GraphModel, minimum_component_size::Int64=10, components=[], bipartite_type="biclique")
    bicliques = []
    biclique_left_nodes = []
    biclique_right_nodes = []
    for c in components
        biclique, original_nodes_left, original_nodes_right = find_biclique_from_component(c, M.G)
        left_length = length(original_nodes_left)
        right_length = length(original_nodes_right)
        found = false
        for (idx,(left,right)) in enumerate(zip(biclique_left_nodes, biclique_right_nodes))
            if length(intersect(original_nodes_left, left))/min(left_length,length(left)) >= 0.90 && length(intersect(original_nodes_right, right))/min(right_length,length(right)) >= 0.90 # was: 0.9
                biclique_left_nodes[idx] = union(original_nodes_left, left)
                biclique_right_nodes[idx] = setdiff(union(original_nodes_right, right), biclique_left_nodes[idx])
                bicliques[idx], _ = induced_subgraph(M.G, union(biclique_left_nodes[idx], biclique_right_nodes[idx]))
                found = true
                break
            end
        end
        if !found && biclique !== nothing
            push!(bicliques, biclique)
            push!(biclique_left_nodes, original_nodes_left)
            push!(biclique_right_nodes, original_nodes_right)
        end
    end
    return bicliques, biclique_left_nodes, biclique_right_nodes
end


"""
Create *all* starclique candidates (raw)
"""
function generate_starclique_component_candidates(M::Model.GraphModel, minimum_component_size::Int64=10, components=[])
    starcliques = []
    starclique_left_nodes = []
    starclique_right_nodes = []
    for c in components
        starclique, original_nodes_left, original_nodes_right = find_starclique_from_component(c, M.G)
        left_length = length(original_nodes_left)
        right_length = length(original_nodes_right)
        found = false
        for (idx,(left,right)) in enumerate(zip(starclique_left_nodes, starclique_right_nodes))
            if length(intersect(original_nodes_left, left))/min(left_length,length(left)) >= 0.9 && length(intersect(original_nodes_right, right))/min(right_length,length(right)) >= 0.9
                starclique_left_nodes[idx] = union(original_nodes_left, left)
                starclique_right_nodes[idx] = setdiff(union(original_nodes_right, right), starclique_left_nodes[idx])
                starcliques[idx], _ = induced_subgraph(M.G, union(starclique_left_nodes[idx], starclique_right_nodes[idx]))
                found = true
                break
            end
        end
        if !found && starclique !== nothing
            push!(starcliques, starclique)
            push!(starclique_left_nodes, original_nodes_left)
            push!(starclique_right_nodes, original_nodes_right)
        end
    end
    return starcliques, starclique_left_nodes, starclique_right_nodes
end



"""
Create *all* star candidates (as objects)
"""
function generate_all_star_structures(hubs, star_spokes_original, stars)
    star_structures = []
    for (hub, spokes, sG) in zip(hubs, star_spokes_original, stars)
        if !isempty(spokes)
            s = Mining.Model.Structures.Star(ne(sG), hub, spokes)
            push!(star_structures, s)
        end
    end
    return star_structures
end


"""
Create *all* clique candidates (as objects)
"""
function generate_all_clique_structures(clique_nodes_original, cliques)
    clique_structures = []
    clique_nodes = []
    for (nodes, cG) in zip(clique_nodes_original, cliques)
        for (idx,elem) in enumerate(map(elem -> elem ⊆ nodes, clique_nodes))
            if elem
                clique_nodes[idx] = nodes
                clique_structures[idx] = Mining.Model.Structures.Clique(ne(cG), nodes)
            end
        end
        if !any(map(elem -> nodes ⊆ elem, clique_nodes))
            if !isempty(nodes)
                c = Mining.Model.Structures.Clique(ne(cG), nodes)
                push!(clique_structures, c)
                push!(clique_nodes, Set(nodes))
            end
        end
    end
    return clique_structures
end


"""
Create *all* biclique or starclique candidates (as objects)
"""
function generate_all_biclique_structures(biclique_left_nodes, biclique_right_nodes, bicliques, G, bipartite_type="biclique")
    biclique_structures = []
    for (l, r, bG) in zip(biclique_left_nodes, biclique_right_nodes, bicliques)
        if !isempty(l)
            n_edges_in_left = ne(induced_subgraph(G, l)[1])
            n_edges_in_right = ne(induced_subgraph(G, r)[1])
            if bipartite_type == "biclique"
                b = Mining.Model.Structures.Biclique(ne(bG), l, r, n_edges_in_left, n_edges_in_right)
            else # assume "starclique"
                b = Mining.Model.Structures.Starclique(ne(bG), l, r, n_edges_in_left, n_edges_in_right)
            end
            push!(biclique_structures, b)
        end
    end
    return biclique_structures
end


"""
Model building entry function
"""
function compute_anytime_model(dataset::String, surprise_threshold::Float64, size_threshold::Float64, stop_threshold::Float64, stop_after_n_structures::Float64, sort_heuristic::String; overwrite::Bool=false, debug::Bool=false)
    dataset_core = splitdir(splitdir(dataset)[1])[end]*"/"*splitext(splitdir(dataset)[end])[1]
    logdir, logfile = splitdir("../results/$(dataset_core)_size-$(size_threshold)_max-$(stop_after_n_structures)_sort-$(sort_heuristic).log")
    # we don't want to accidentally overwrite logs when there is a corresponding model already
    if !overwrite && isfile(logdir * "/" * logfile[1:end-3] * "json")
        println("Skipping $(logfile) because there exists a corresponding model and overwrite is false...")
        return
    end
    mkpath(logdir)
    open("$(logdir)/$(logfile)", "w") do io
        M = Model.GraphModel(dataset, io)
        M.parameters["surprise_threshold"] = surprise_threshold
        M.parameters["size_threshold"] = size_threshold
        M.parameters["stop_after_n_structures"] = stop_after_n_structures
        if debug
            Model.log_progress!(M, "starting description length: $(M.start_description_length_maxent)\n")
            build_decomposition_based_model!(M, size_threshold, stop_threshold, stop_after_n_structures, sort_heuristic)
            Model.save_model(M, "../results"; overwrite=overwrite)
        else
            try
                Model.log_progress!(M, "starting description length: $(M.start_description_length_maxent)\n")
                build_decomposition_based_model!(M, size_threshold, stop_threshold, stop_after_n_structures, sort_heuristic)
            catch InterruptException
                println("\nCatching interrupt exception...")
            finally
                println("Saving state...")
                Model.save_model(M, "../results"; overwrite=overwrite)
                println("Done.")
            end
        end
    end
end

end