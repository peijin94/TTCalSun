function parse_args(args::Vector{String})
    opts = Dict{String, Any}(
        "command" => "",
        "sources" => "",
        "column" => "CORRECTED_DATA",
        "maxiter" => 30,
        "tolerance" => 1e-4,
        "minuvw" => 10.0,
        "peeliter" => 3,
        "phase_only_maxiter" => 0,
        "timings" => false,
        "ms_files" => String[],
    )
    positional = String[]
    for arg in args
        if startswith(arg, "--column=")
            opts["column"] = split(arg, "=", limit=2)[2]
        elseif startswith(arg, "--maxiter=")
            opts["maxiter"] = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--tolerance=")
            opts["tolerance"] = parse(Float64, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--minuvw=")
            opts["minuvw"] = parse(Float64, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--peeliter=")
            opts["peeliter"] = parse(Int, split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--phase-only-maxiter=")
            opts["phase_only_maxiter"] = parse(Int, split(arg, "=", limit=2)[2])
        elseif arg == "--timings"
            opts["timings"] = true
        elseif arg == "--help" || arg == "-h"
            opts["help"] = true
        elseif startswith(arg, "--")
            error("unknown option: $arg")
        else
            push!(positional, arg)
        end
    end
    if get(opts, "help", false)
        return opts
    end
    length(positional) < 3 && error("usage: ttcalsun.jl <peel|shave|zest|prune> <sources.json> <ms1> [ms2...]")
    opts["command"] = positional[1]
    opts["sources"] = positional[2]
    opts["ms_files"] = positional[3:end]
    return opts
end

function run_cli(args::Vector{String}=ARGS)
    opts = parse_args(args)
    if get(opts, "help", false)
        println("usage: ttcalsun.jl <peel|shave|zest|prune> [--column=CORRECTED_DATA] [--maxiter=30] [--tolerance=1e-4] [--minuvw=10] [--peeliter=3] [--phase-only-maxiter=0] [--timings] <sources.json> <ms1> [ms2 ...]")
        println("threads=$(nthreads())")
        return 0
    end

    command = Symbol(opts["command"])
    sources = read_sources(opts["sources"])
    @printf("TTCalSun CPU mode: threads=%d, sources=%d, command=%s\n", nthreads(), length(sources), String(opts["command"]))
    for ms_path in opts["ms_files"]
        t0 = time()
        vis, meta, baseline_dict, nrows = read_ms(ms_path; column=opts["column"])
        result = if command == :peel
            peel!(vis, meta, sources; peeliter=opts["peeliter"], maxiter=opts["maxiter"], tolerance=opts["tolerance"], minuvw=opts["minuvw"], phase_only_maxiter=opts["phase_only_maxiter"], return_timing=opts["timings"])
        elseif command == :shave
            shave!(vis, meta, sources; peeliter=opts["peeliter"], maxiter=opts["maxiter"], tolerance=opts["tolerance"], minuvw=opts["minuvw"], phase_only_maxiter=opts["phase_only_maxiter"], return_timing=opts["timings"])
        elseif command == :zest
            zest!(vis, meta, sources; peeliter=opts["peeliter"], maxiter=opts["maxiter"], tolerance=opts["tolerance"], minuvw=opts["minuvw"], phase_only_maxiter=opts["phase_only_maxiter"], return_timing=opts["timings"])
        elseif command == :prune
            prune!(vis, meta, sources; peeliter=opts["peeliter"], maxiter=opts["maxiter"], tolerance=opts["tolerance"], minuvw=opts["minuvw"], phase_only_maxiter=opts["phase_only_maxiter"], return_timing=opts["timings"])
        else
            error("unsupported command: $(opts["command"])")
        end

        if opts["timings"]
            _, timing = result
            print_timing_summary(stdout, timing)
        end

        write_ms!(ms_path, vis, baseline_dict, nrows; column=opts["column"])
        @printf("Processed %s in %.2f s\n", ms_path, time() - t0)
    end
    return 0
end
