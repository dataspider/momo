using ArgParse
using Dates
using Distributed
@everywhere include("./mining.jl")
@everywhere using Distributed
@everywhere using Dates
@everywhere using .Mining

ccall(:jl_exit_on_sigint, Cvoid, (Cint,), 0)

function parse_cmd()
    s = ArgParseSettings()

    @add_arg_table! s begin
            "data_path"
                help = "data folder relative to data directory"
                required = true
            "--sort_heuristic"
                default = "plain" # alternative is "uncovered_nodes" but this doesn't work as well
                arg_type = String
            "--surprise_threshold"
                default = 10.
                arg_type = Float64
            "--size_threshold"
                default = 10.
                arg_type = Float64
            "--stop_threshold"
                default = Inf
                arg_type = Float64
            "--stop_after_n_structures"
                default = Inf
                arg_type = Float64
            "--overwrite", "-o"
                default = false
                arg_type = Bool
            "--single_process", "-s"
                default = false
                arg_type = Bool
            "--debug", "-d"
                default = false
                arg_type = Bool
    end
    return parse_args(s)
end

@everywhere function build_model_for_dataset(dataset::String, surprise_threshold::Float64, size_threshold::Float64, stop_threshold::Float64, stop_after_n_structures::Float64, sort_heuristic::String; overwrite::Bool=false, debug::Bool=false)
    println("$(now()) - Starting to compute model for $(dataset) with size threshold $(size_threshold), stop threshold $(stop_threshold), and at most $(stop_after_n_structures) structures...")
    flush(stdout)
    Mining.compute_anytime_model(dataset, surprise_threshold, size_threshold, stop_threshold, stop_after_n_structures, sort_heuristic; overwrite=overwrite, debug=debug)
end

function main()
    parsed_args = parse_cmd()
    dsbase = parsed_args["data_path"]
    datasets = [dsbase * "/" * x for x in readdir(dsbase) if (endswith(x, ".txt") || endswith(x, ".csv")) && !occursin("nodemap",x)]
    params = [parsed_args["surprise_threshold"], parsed_args["size_threshold"], parsed_args["stop_threshold"], parsed_args["stop_after_n_structures"], parsed_args["sort_heuristic"]]
    
    if parsed_args["single_process"]
        for dataset in datasets
            build_model_for_dataset(dataset, params...; overwrite=parsed_args["overwrite"], debug=parsed_args["debug"])
        end
    else
        pmap(x->build_model_for_dataset(x, params...; overwrite=parsed_args["overwrite"], debug=parsed_args["debug"]), datasets)
    end
end

main()

