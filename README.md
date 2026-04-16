# TTCalSun

CPU-only direction-dependent calibration for the local LWA workflow.

`TTCalSun` keeps the useful parts of `TTCal.jl` and `TTCalX`, but drops all CUDA codepaths:
- `TTCal.jl`: calibration / peeling structure
- `TTCalX`: Measurement Set I/O pattern and source JSON compatibility

This version is optimized for this workspace:
- Measurement Sets from `./MS`
- source models in the existing `TTCalX/sources.json` format
- `python-casacore` from `./lwa/bin/python`
- batch processing of multiple MS files in one Julia process

## Setup

Build `PyCall` against the local `lwa` Python once:

```bash
cd /fast/pipe2026solar/TTCalSun
export PYTHON=/fast/pipe2026solar/lwa/bin/python
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.build("PyCall")'
```

## Usage

Preferred entrypoint:

```bash
/fast/pipe2026solar/TTCalSun/bin/ttcalsun_cpu.sh \
  peel \
  /fast/pipe2026solar/TTCalX/sources.json \
  /path/to/file1.ms /path/to/file2.ms \
  --column=CORRECTED_DATA \
  --maxiter=30 \
  --tolerance=1e-4 \
  --minuvw=10 \
  --peeliter=3 \
  --timings
```

Available commands:
- `peel`
- `shave`
- `zest`
- `prune`

Batch mode is preferred. The same Julia process can process multiple MS files, so startup/JIT is only paid once.

`--timings` prints a per-source breakdown of:
- model generation
- model-square construction
- initial subtraction
- put-back time
- solve time
- measured-square time inside the solve
- solver iteration time
- subtraction after the updated solve

## Local Benchmark

On this machine, with `64` CPU threads and no GPU:

- `TTCalX` CPU-only `zest`, two-file batch (`73 + 78 MHz`): `232.50 s` wall, `85.33 s / MS` steady-state
- `TTCalSun` CPU-only `zest`, two-file batch (`73 + 78 MHz`): `31.58 s` wall
- `TTCalSun` per-file timings from that batch:
  - `73 MHz`: `16.56 s`
  - `78 MHz`: `13.40 s`

Single-file `73 MHz` run:

- `TTCalSun` wall time: `19.75 s`
- inner processing time reported by the CLI: `16.61 s`

## Notes

- This is intentionally CPU-only. There is no CUDA dependency in the package.
- The source reader accepts the current `TTCalX/sources.json` layout, including multi-component and Gaussian sources.
- The CLI reads and writes the specified MS column in place.
- The end-to-end peel benchmark driver for this workspace is in [tests/run_ttcalsun_peel_benchmark.sh](/fast/pipe2026solar/tests/run_ttcalsun_peel_benchmark.sh).
