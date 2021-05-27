module NodeMap

using DataFrames
using CSV

using LightGraphs
import LightGraphs.loadgraph
import GraphIO.EdgeList.EdgeListFormat

function loadedgelist(io::IO, gname::String)
    srcs = Vector{String}()
    dsts = Vector{String}()
    while !eof(io)
        line = strip(chomp(readline(io)))
        if !startswith(line, "#") && (line != "")
            r = r"([\w-]+)[\s,]+([\w-]+)"
            src_s, dst_s = match(r, line).captures
            push!(srcs, src_s)
            push!(dsts, dst_s)
        end
    end
    vxset = unique(vcat(srcs, dsts))
    vxdict = Dict{String,Int}()
    for (v, k) in enumerate(vxset)
        vxdict[k] = v
    end

    n_v = length(vxset)
    g = LightGraphs.DiGraph(n_v)
    for (u, v) in zip(srcs, dsts)
        add_edge!(g, vxdict[u], vxdict[v])
    end
    return g, vxdict
end

loadgraph(io::IO, gname::String, ::EdgeListFormat) = loadedgelist(io, gname)

function create_node_map(graph_file)
    G, node_map = loadgraph(graph_file, EdgeListFormat())
    df = DataFrame(original_id = [keys(node_map)...], julia_id = [values(node_map)...])
    save_file = splitext(graph_file)[1] * "-nodemap.csv"
    CSV.write(save_file, df)
end

end