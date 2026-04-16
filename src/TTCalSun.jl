module TTCalSun

using Base.Threads: @threads, nthreads
using LinearAlgebra
using Printf
using PyCall

include("types.jl")
include("sources.jl")
include("msio.jl")
include("cpu.jl")
include("cli.jl")

export AbstractSource,
       PointSource,
       GaussianSource,
       MultiSource,
       PowerLawSpectrum,
       Visibilities,
       Calibration,
       Metadata,
       ProcessTiming,
       SourceTiming,
       read_sources,
       is_above_horizon,
       source_stokes_i,
       read_ms,
       write_ms!,
       peel!,
       zest!,
       shave!,
       prune!,
       print_timing_summary,
       run_cli

end
