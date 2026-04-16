const LIGHT_SPEED = 299792458.0

mutable struct SolveTiming
    phase_only_s::Float64
    phase_only_iters::Int
    phase_only_converged::Bool
    make_square_s::Float64
    iterate_s::Float64
    niters::Int
    converged::Bool
end

mutable struct SourceTiming
    name::String
    genvis_s::Float64
    model_square_s::Float64
    initial_subtract_s::Float64
    putsrc_s::Float64
    phase_only_s::Float64
    phase_only_iters::Int
    phase_only_converged::Int
    solve_total_s::Float64
    solve_make_square_s::Float64
    solve_iterate_s::Float64
    subsrc_s::Float64
    solves::Int
    niters::Int
    converged::Int
end

SourceTiming(source::AbstractSource) = SourceTiming(
    get_name(source),
    0.0,
    0.0,
    0.0,
    0,
    0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0,
    0,
    0,
)

mutable struct ProcessTiming
    mode::Symbol
    total_s::Float64
    source_timings::Vector{SourceTiming}
end

struct SourceWorkspace
    source::AbstractSource
    calibration::Calibration
    coherency::Visibilities
    corrupted::Visibilities
    model_sq::SquareVisibilities
end

struct SolveWorkspace
    measured_sq::SquareVisibilities
    step::Calibration
end

function baseline_keep_mask(meta::Metadata, minuvw::Float64)
    keep = trues(nbase(meta), nfreq(meta))
    minuvw <= 0 && return keep
    @threads for β in 1:nfreq(meta)
        λ = LIGHT_SPEED / meta.channels[β]
        threshold = minuvw * λ
        @inbounds for α in 1:nbase(meta)
            keep[α, β] = meta.baseline_lengths[α] >= threshold
        end
    end
    return keep
end

function make_square!(sq::SquareVisibilities, vis::Visibilities, meta::Metadata, keep::BitMatrix)
    fillzero!(sq)
    wideband = size(sq.xx, 3) == 1 && nfreq(vis) > 1
    if wideband
        @inbounds for β in 1:nfreq(vis), α in 1:nbase(meta)
            (vis.flags[α, β] || !keep[α, β]) && continue
            ant1 = meta.baselines[1, α]
            ant2 = meta.baselines[2, α]

            vxx = vis.xx[α, β]
            vxy = vis.xy[α, β]
            vyx = vis.yx[α, β]
            vyy = vis.yy[α, β]

            sq.xx[ant1, ant2, 1] += vxx
            sq.xy[ant1, ant2, 1] += vxy
            sq.yx[ant1, ant2, 1] += vyx
            sq.yy[ant1, ant2, 1] += vyy

            sq.xx[ant2, ant1, 1] += conj(vxx)
            sq.xy[ant2, ant1, 1] += conj(vyx)
            sq.yx[ant2, ant1, 1] += conj(vxy)
            sq.yy[ant2, ant1, 1] += conj(vyy)
        end
    else
        @threads for β in 1:nfreq(vis)
            @inbounds for α in 1:nbase(meta)
                (vis.flags[α, β] || !keep[α, β]) && continue
                ant1 = meta.baselines[1, α]
                ant2 = meta.baselines[2, α]

                vxx = vis.xx[α, β]
                vxy = vis.xy[α, β]
                vyx = vis.yx[α, β]
                vyy = vis.yy[α, β]

                sq.xx[ant1, ant2, β] += vxx
                sq.xy[ant1, ant2, β] += vxy
                sq.yx[ant1, ant2, β] += vyx
                sq.yy[ant1, ant2, β] += vyy

                sq.xx[ant2, ant1, β] += conj(vxx)
                sq.xy[ant2, ant1, β] += conj(vyx)
                sq.yx[ant2, ant1, β] += conj(vxy)
                sq.yy[ant2, ant1, β] += conj(vyy)
            end
        end
    end
    return sq
end

function subsrc!(vis::Visibilities, model::Visibilities)
    vis.xx .-= model.xx
    vis.xy .-= model.xy
    vis.yx .-= model.yx
    vis.yy .-= model.yy
    return vis
end

function putsrc!(vis::Visibilities, model::Visibilities)
    vis.xx .+= model.xx
    vis.xy .+= model.xy
    vis.yx .+= model.yx
    vis.yy .+= model.yy
    return vis
end

function corrupt!(vis::Visibilities, cal::Calibration, meta::Metadata)
    nfreq_cal = size(cal.xx, 2)
    @threads for β in 1:nfreq(vis)
        β_cal = nfreq_cal == 1 ? 1 : β
        @inbounds for α in 1:nbase(meta)
            ant1 = meta.baselines[1, α]
            ant2 = meta.baselines[2, α]

            if cal.is_diagonal
                j1_xx = cal.xx[ant1, β_cal]
                j1_yy = cal.yy[ant1, β_cal]
                j2_xx = cal.xx[ant2, β_cal]
                j2_yy = cal.yy[ant2, β_cal]

                vxx = vis.xx[α, β]
                vxy = vis.xy[α, β]
                vyx = vis.yx[α, β]
                vyy = vis.yy[α, β]

                vis.xx[α, β] = j1_xx * vxx * conj(j2_xx)
                vis.xy[α, β] = j1_xx * vxy * conj(j2_yy)
                vis.yx[α, β] = j1_yy * vyx * conj(j2_xx)
                vis.yy[α, β] = j1_yy * vyy * conj(j2_yy)
            else
                j1_xx = cal.xx[ant1, β_cal]
                j1_xy = cal.xy[ant1, β_cal]
                j1_yx = cal.yx[ant1, β_cal]
                j1_yy = cal.yy[ant1, β_cal]

                j2_xx = cal.xx[ant2, β_cal]
                j2_xy = cal.xy[ant2, β_cal]
                j2_yx = cal.yx[ant2, β_cal]
                j2_yy = cal.yy[ant2, β_cal]

                vxx = vis.xx[α, β]
                vxy = vis.xy[α, β]
                vyx = vis.yx[α, β]
                vyy = vis.yy[α, β]

                jv_xx = j1_xx * vxx + j1_xy * vyx
                jv_xy = j1_xx * vxy + j1_xy * vyy
                jv_yx = j1_yx * vxx + j1_yy * vyx
                jv_yy = j1_yx * vxy + j1_yy * vyy

                vis.xx[α, β] = jv_xx * conj(j2_xx) + jv_xy * conj(j2_xy)
                vis.xy[α, β] = jv_xx * conj(j2_yx) + jv_xy * conj(j2_yy)
                vis.yx[α, β] = jv_yx * conj(j2_xx) + jv_yy * conj(j2_xy)
                vis.yy[α, β] = jv_yx * conj(j2_yx) + jv_yy * conj(j2_yy)
            end
        end
    end
    return vis
end

function genvis!(vis::Visibilities, meta::Metadata, source::PointSource; phase_center_ra::Float64, phase_center_dec::Float64)
    l, m, n = source_direction_lmn(source, phase_center_ra, phase_center_dec)
    flux_xx = zeros(ComplexF64, nfreq(meta))
    flux_xy = zeros(ComplexF64, nfreq(meta))
    flux_yx = zeros(ComplexF64, nfreq(meta))
    flux_yy = zeros(ComplexF64, nfreq(meta))
    @inbounds for β in 1:nfreq(meta)
        I, Q, U, V = source.spectrum(meta.channels[β])
        flux_xx[β] = ComplexF64(I + Q)
        flux_xy[β] = ComplexF64(U + im * V)
        flux_yx[β] = ComplexF64(U - im * V)
        flux_yy[β] = ComplexF64(I - Q)
    end

    @threads for β in 1:nfreq(meta)
        λ = LIGHT_SPEED / meta.channels[β]
        fxx = flux_xx[β]
        fxy = flux_xy[β]
        fyx = flux_yx[β]
        fyy = flux_yy[β]
        @inbounds for α in 1:nbase(meta)
            u = meta.uvw[1, α]
            v = meta.uvw[2, α]
            w = meta.uvw[3, α]
            phase = -2π * (u * l + v * m + w * (n - 1.0)) / λ
            fringe = cis(phase)
            vis.xx[α, β] += fxx * fringe
            vis.xy[α, β] += fxy * fringe
            vis.yx[α, β] += fyx * fringe
            vis.yy[α, β] += fyy * fringe
        end
    end
    return vis
end

function genvis!(vis::Visibilities, meta::Metadata, source::GaussianSource; phase_center_ra::Float64, phase_center_dec::Float64)
    l, m, n = source_direction_lmn(source, phase_center_ra, phase_center_dec)
    flux_xx = zeros(ComplexF64, nfreq(meta))
    flux_xy = zeros(ComplexF64, nfreq(meta))
    flux_yx = zeros(ComplexF64, nfreq(meta))
    flux_yy = zeros(ComplexF64, nfreq(meta))
    @inbounds for β in 1:nfreq(meta)
        I, Q, U, V = source.spectrum(meta.channels[β])
        flux_xx[β] = ComplexF64(I + Q)
        flux_xy[β] = ComplexF64(U + im * V)
        flux_yx[β] = ComplexF64(U - im * V)
        flux_yy[β] = ComplexF64(I - Q)
    end

    dra = source.ra - phase_center_ra
    rotation = atan(sin(dra), cos(source.dec) * tan(phase_center_dec) - sin(source.dec) * cos(dra))
    pa = source.position_angle + rotation
    sigma_maj = source.major_fwhm / (2.0 * sqrt(2.0 * log(2.0)))
    sigma_min = source.minor_fwhm / (2.0 * sqrt(2.0 * log(2.0)))
    cos_pa = cos(pa)
    sin_pa = sin(pa)

    @threads for β in 1:nfreq(meta)
        λ = LIGHT_SPEED / meta.channels[β]
        fxx = flux_xx[β]
        fxy = flux_xy[β]
        fyx = flux_yx[β]
        fyy = flux_yy[β]
        @inbounds for α in 1:nbase(meta)
            u = meta.uvw[1, α]
            v = meta.uvw[2, α]
            w = meta.uvw[3, α]
            phase = -2π * (u * l + v * m + w * (n - 1.0)) / λ
            fringe = cis(phase)

            u_rot = u * cos_pa + v * sin_pa
            v_rot = -u * sin_pa + v * cos_pa
            u_lam = u_rot / λ
            v_lam = v_rot / λ
            envelope = exp(-2π^2 * (sigma_maj^2 * u_lam^2 + sigma_min^2 * v_lam^2))

            vis.xx[α, β] += fxx * fringe * envelope
            vis.xy[α, β] += fxy * fringe * envelope
            vis.yx[α, β] += fyx * fringe * envelope
            vis.yy[α, β] += fyy * fringe * envelope
        end
    end
    return vis
end

function genvis!(vis::Visibilities, meta::Metadata, source::MultiSource; phase_center_ra::Float64, phase_center_dec::Float64)
    for component in source.components
        genvis!(vis, meta, component; phase_center_ra=phase_center_ra, phase_center_dec=phase_center_dec)
    end
    return vis
end

function solve_step_diagonal!(step::Calibration, cal::Calibration, meas_sq::SquareVisibilities, model_sq::SquareVisibilities)
    na = nant(cal)
    nf = nfreq(cal)
    fill!(step.xx, 0.0 + 0.0im)
    fill!(step.yy, 0.0 + 0.0im)

    @threads for idx in 1:(na * nf)
        j = ((idx - 1) % na) + 1
        β = ((idx - 1) ÷ na) + 1

        num_xx = 0.0 + 0.0im
        num_yy = 0.0 + 0.0im
        den_xx = 0.0 + 0.0im
        den_yy = 0.0 + 0.0im

        @inbounds for i in 1:na
            gi_xx = cal.xx[i, β]
            gi_yy = cal.yy[i, β]

            m_xx = model_sq.xx[i, j, β]
            m_xy = model_sq.xy[i, j, β]
            m_yx = model_sq.yx[i, j, β]
            m_yy = model_sq.yy[i, j, β]

            v_xx = meas_sq.xx[i, j, β]
            v_xy = meas_sq.xy[i, j, β]
            v_yx = meas_sq.yx[i, j, β]
            v_yy = meas_sq.yy[i, j, β]

            gm_xx = gi_xx * m_xx
            gm_xy = gi_xx * m_xy
            gm_yx = gi_yy * m_yx
            gm_yy = gi_yy * m_yy

            num_xx += conj(gm_xx) * v_xx + conj(gm_yx) * v_yx
            num_yy += conj(gm_xy) * v_xy + conj(gm_yy) * v_yy
            den_xx += conj(gm_xx) * gm_xx + conj(gm_yx) * gm_yx
            den_yy += conj(gm_xy) * gm_xy + conj(gm_yy) * gm_yy
        end

        if abs(den_xx) > eps(Float64)
            step.xx[j, β] = conj(num_xx / den_xx) - cal.xx[j, β]
        end
        if abs(den_yy) > eps(Float64)
            step.yy[j, β] = conj(num_yy / den_yy) - cal.yy[j, β]
        end
    end
    return step
end

function solve_step_full!(step::Calibration, cal::Calibration, meas_sq::SquareVisibilities, model_sq::SquareVisibilities)
    na = nant(cal)
    nf = nfreq(cal)
    fill!(step.xx, 0.0 + 0.0im)
    fill!(step.xy, 0.0 + 0.0im)
    fill!(step.yx, 0.0 + 0.0im)
    fill!(step.yy, 0.0 + 0.0im)

    @threads for idx in 1:(na * nf)
        j = ((idx - 1) % na) + 1
        β = ((idx - 1) ÷ na) + 1

        num_xx = 0.0 + 0.0im
        num_xy = 0.0 + 0.0im
        num_yx = 0.0 + 0.0im
        num_yy = 0.0 + 0.0im
        den_xx = 0.0 + 0.0im
        den_xy = 0.0 + 0.0im
        den_yx = 0.0 + 0.0im
        den_yy = 0.0 + 0.0im

        @inbounds for i in 1:na
            gi_xx = cal.xx[i, β]
            gi_xy = cal.xy[i, β]
            gi_yx = cal.yx[i, β]
            gi_yy = cal.yy[i, β]

            m_xx = model_sq.xx[i, j, β]
            m_xy = model_sq.xy[i, j, β]
            m_yx = model_sq.yx[i, j, β]
            m_yy = model_sq.yy[i, j, β]

            v_xx = meas_sq.xx[i, j, β]
            v_xy = meas_sq.xy[i, j, β]
            v_yx = meas_sq.yx[i, j, β]
            v_yy = meas_sq.yy[i, j, β]

            gm_xx = gi_xx * m_xx + gi_xy * m_yx
            gm_xy = gi_xx * m_xy + gi_xy * m_yy
            gm_yx = gi_yx * m_xx + gi_yy * m_yx
            gm_yy = gi_yx * m_xy + gi_yy * m_yy

            num_xx += conj(gm_xx) * v_xx + conj(gm_yx) * v_yx
            num_xy += conj(gm_xx) * v_xy + conj(gm_yx) * v_yy
            num_yx += conj(gm_xy) * v_xx + conj(gm_yy) * v_yx
            num_yy += conj(gm_xy) * v_xy + conj(gm_yy) * v_yy

            den_xx += conj(gm_xx) * gm_xx + conj(gm_yx) * gm_yx
            den_xy += conj(gm_xx) * gm_xy + conj(gm_yx) * gm_yy
            den_yx += conj(gm_xy) * gm_xx + conj(gm_yy) * gm_yx
            den_yy += conj(gm_xy) * gm_xy + conj(gm_yy) * gm_yy
        end

        det_den = den_xx * den_yy - den_xy * den_yx
        if abs(det_den) > eps(Float64)
            inv_det = 1.0 / det_den
            inv_xx = den_yy * inv_det
            inv_xy = -den_xy * inv_det
            inv_yx = -den_yx * inv_det
            inv_yy = den_xx * inv_det

            sol_xx = inv_xx * num_xx + inv_xy * num_yx
            sol_xy = inv_xx * num_xy + inv_xy * num_yy
            sol_yx = inv_yx * num_xx + inv_yy * num_yx
            sol_yy = inv_yx * num_xy + inv_yy * num_yy

            step.xx[j, β] = conj(sol_xx) - cal.xx[j, β]
            step.xy[j, β] = conj(sol_yx) - cal.xy[j, β]
            step.yx[j, β] = conj(sol_xy) - cal.yx[j, β]
            step.yy[j, β] = conj(sol_yy) - cal.yy[j, β]
        end
    end
    return step
end

function step_norm(step::Calibration)
    if step.is_diagonal
        return sqrt(sum(abs2, step.xx) + sum(abs2, step.yy))
    end
    return sqrt(sum(abs2, step.xx) + sum(abs2, step.xy) + sum(abs2, step.yx) + sum(abs2, step.yy))
end

function gain_norm(cal::Calibration)
    if cal.is_diagonal
        return sqrt(sum(abs2, cal.xx) + sum(abs2, cal.yy))
    end
    return sqrt(sum(abs2, cal.xx) + sum(abs2, cal.xy) + sum(abs2, cal.yx) + sum(abs2, cal.yy))
end

function apply_step!(cal::Calibration, step::Calibration, damping::Float64)
    cal.xx .+= damping .* step.xx
    cal.yy .+= damping .* step.yy
    if !cal.is_diagonal
        cal.xy .+= damping .* step.xy
        cal.yx .+= damping .* step.yx
    end
    return cal
end

function project_phase_only!(cal::Calibration)
    @inbounds for i in eachindex(cal.xx)
        ax = abs(cal.xx[i])
        ay = abs(cal.yy[i])
        if ax > eps(Float64)
            cal.xx[i] /= ax
        else
            cal.xx[i] = 1.0 + 0.0im
        end
        if ay > eps(Float64)
            cal.yy[i] /= ay
        else
            cal.yy[i] = 1.0 + 0.0im
        end
    end
    if !cal.is_diagonal
        fill!(cal.xy, 0.0 + 0.0im)
        fill!(cal.yx, 0.0 + 0.0im)
    end
    return cal
end

function stefcal!(cal::Calibration, vis::Visibilities, model_sq::SquareVisibilities, meta::Metadata, solve::SolveWorkspace, keep::BitMatrix; maxiter::Int, tolerance::Float64, phase_only_maxiter::Int=0)
    t_make_square = time()
    make_square!(solve.measured_sq, vis, meta, keep)
    make_square_s = time() - t_make_square

    phase_only_s = 0.0
    phase_only_iters = 0
    phase_only_converged = false
    if phase_only_maxiter > 0 && cal.is_diagonal
        t_phase = time()
        for iter in 1:phase_only_maxiter
            phase_only_iters = iter
            solve_step_diagonal!(solve.step, cal, solve.measured_sq, model_sq)
            apply_step!(cal, solve.step, 0.5)
            project_phase_only!(cal)
            if step_norm(solve.step) < tolerance * gain_norm(cal)
                phase_only_converged = true
                break
            end
        end
        phase_only_s = time() - t_phase
    end

    converged = false
    niters = 0
    t_iterate = time()
    for iter in 1:maxiter
        niters = iter
        if cal.is_diagonal
            solve_step_diagonal!(solve.step, cal, solve.measured_sq, model_sq)
        else
            solve_step_full!(solve.step, cal, solve.measured_sq, model_sq)
        end
        if step_norm(solve.step) < tolerance * gain_norm(cal)
            converged = true
            apply_step!(cal, solve.step, 0.5)
            break
        end
        apply_step!(cal, solve.step, 0.5)
    end
    iterate_s = time() - t_iterate
    return SolveTiming(phase_only_s, phase_only_iters, phase_only_converged, make_square_s, iterate_s, niters, converged)
end

function mode_flags(mode::Symbol)
    if mode == :peel
        return true, false
    elseif mode == :shave
        return true, true
    elseif mode == :zest
        return false, false
    elseif mode == :prune
        return false, true
    end
    error("unsupported mode: $mode")
end

function prepare_workspaces(meta::Metadata, sources::Vector{<:AbstractSource}, mode::Symbol, keep::BitMatrix, source_timings::Vector{SourceTiming})
    diagonal, wideband = mode_flags(mode)
    cal_nfreq = wideband ? 1 : nfreq(meta)
    workspaces = SourceWorkspace[]
    for (idx, source) in enumerate(sources)
        coherency = Visibilities(nbase(meta), nfreq(meta))
        t_genvis = time()
        genvis!(coherency, meta, source; phase_center_ra=meta.phase_center_ra, phase_center_dec=meta.phase_center_dec)
        source_timings[idx].genvis_s += time() - t_genvis

        model_sq = SquareVisibilities(nant(meta), cal_nfreq)
        t_model_square = time()
        make_square!(model_sq, coherency, meta, keep)
        source_timings[idx].model_square_s += time() - t_model_square

        calibration = Calibration(nant(meta), cal_nfreq; diagonal=diagonal)
        corrupted = Visibilities(nbase(meta), nfreq(meta))
        push!(workspaces, SourceWorkspace(source, calibration, coherency, corrupted, model_sq))
    end
    return workspaces
end

function process_mode!(vis::Visibilities, meta::Metadata, sources::Vector{<:AbstractSource}, mode::Symbol; peeliter::Int, maxiter::Int, tolerance::Float64, minuvw::Float64, phase_only_maxiter::Int=0, return_timing::Bool=false)
    total_t0 = time()
    active_sources = AbstractSource[
        source for source in sources
        if is_above_horizon(source, meta.phase_center_ra, meta.phase_center_dec)
    ]
    reference_frequency = meta.channels[cld(length(meta.channels), 2)]
    sort!(
        active_sources;
        by=source -> source_stokes_i(source, reference_frequency),
        rev=true,
    )
    keep = baseline_keep_mask(meta, minuvw)
    source_timings = [SourceTiming(source) for source in active_sources]
    source_workspaces = prepare_workspaces(meta, active_sources, mode, keep, source_timings)
    diagonal, wideband = mode_flags(mode)
    solve = SolveWorkspace(
        SquareVisibilities(nant(meta), wideband ? 1 : nfreq(meta)),
        Calibration(nant(meta), wideband ? 1 : nfreq(meta); diagonal=diagonal),
    )

    for (idx, ws) in enumerate(source_workspaces)
        t_initial = time()
        copyto!(ws.corrupted, ws.coherency)
        corrupt!(ws.corrupted, ws.calibration, meta)
        subsrc!(vis, ws.corrupted)
        source_timings[idx].initial_subtract_s += time() - t_initial
    end

    for _ in 1:peeliter
        for (idx, ws) in enumerate(source_workspaces)
            t_put = time()
            copyto!(ws.corrupted, ws.coherency)
            corrupt!(ws.corrupted, ws.calibration, meta)
            putsrc!(vis, ws.corrupted)
            source_timings[idx].putsrc_s += time() - t_put

            solve_timing = stefcal!(ws.calibration, vis, ws.model_sq, meta, solve, keep; maxiter=maxiter, tolerance=tolerance, phase_only_maxiter=phase_only_maxiter)
            source_timings[idx].phase_only_s += solve_timing.phase_only_s
            source_timings[idx].phase_only_iters += solve_timing.phase_only_iters
            source_timings[idx].phase_only_converged += solve_timing.phase_only_converged ? 1 : 0
            source_timings[idx].solve_total_s += solve_timing.phase_only_s + solve_timing.make_square_s + solve_timing.iterate_s
            source_timings[idx].solve_make_square_s += solve_timing.make_square_s
            source_timings[idx].solve_iterate_s += solve_timing.iterate_s
            source_timings[idx].solves += 1
            source_timings[idx].niters += solve_timing.niters
            source_timings[idx].converged += solve_timing.converged ? 1 : 0

            t_sub = time()
            copyto!(ws.corrupted, ws.coherency)
            corrupt!(ws.corrupted, ws.calibration, meta)
            subsrc!(vis, ws.corrupted)
            source_timings[idx].subsrc_s += time() - t_sub
        end
    end

    timing = ProcessTiming(mode, time() - total_t0, source_timings)
    calibrations = [ws.calibration for ws in source_workspaces]
    return return_timing ? (calibrations, timing) : calibrations
end

peel!(vis::Visibilities, meta::Metadata, sources::Vector{<:AbstractSource}; peeliter::Int=3, maxiter::Int=20, tolerance::Float64=1e-3, minuvw::Float64=10.0, phase_only_maxiter::Int=0, return_timing::Bool=false) =
    process_mode!(vis, meta, sources, :peel; peeliter=peeliter, maxiter=maxiter, tolerance=tolerance, minuvw=minuvw, phase_only_maxiter=phase_only_maxiter, return_timing=return_timing)

shave!(vis::Visibilities, meta::Metadata, sources::Vector{<:AbstractSource}; peeliter::Int=3, maxiter::Int=20, tolerance::Float64=1e-3, minuvw::Float64=10.0, phase_only_maxiter::Int=0, return_timing::Bool=false) =
    process_mode!(vis, meta, sources, :shave; peeliter=peeliter, maxiter=maxiter, tolerance=tolerance, minuvw=minuvw, phase_only_maxiter=phase_only_maxiter, return_timing=return_timing)

zest!(vis::Visibilities, meta::Metadata, sources::Vector{<:AbstractSource}; peeliter::Int=3, maxiter::Int=20, tolerance::Float64=1e-3, minuvw::Float64=10.0, phase_only_maxiter::Int=0, return_timing::Bool=false) =
    process_mode!(vis, meta, sources, :zest; peeliter=peeliter, maxiter=maxiter, tolerance=tolerance, minuvw=minuvw, phase_only_maxiter=phase_only_maxiter, return_timing=return_timing)

prune!(vis::Visibilities, meta::Metadata, sources::Vector{<:AbstractSource}; peeliter::Int=3, maxiter::Int=20, tolerance::Float64=1e-3, minuvw::Float64=10.0, phase_only_maxiter::Int=0, return_timing::Bool=false) =
    process_mode!(vis, meta, sources, :prune; peeliter=peeliter, maxiter=maxiter, tolerance=tolerance, minuvw=minuvw, phase_only_maxiter=phase_only_maxiter, return_timing=return_timing)

function print_timing_summary(io::IO, timing::ProcessTiming)
    @printf(io, "TTCalSun timing summary (%s): total %.2f s\n", String(timing.mode), timing.total_s)
    for source in timing.source_timings
        @printf(
            io,
            "  %-8s genvis=%.2f model_sq=%.2f init_sub=%.2f put=%.2f solve=%.2f (phase=%.2f phase_iters=%d phase_conv=%d make_sq=%.2f iterate=%.2f solves=%d niters=%d converged=%d) sub=%.2f\n",
            source.name,
            source.genvis_s,
            source.model_square_s,
            source.initial_subtract_s,
            source.putsrc_s,
            source.solve_total_s,
            source.phase_only_s,
            source.phase_only_iters,
            source.phase_only_converged,
            source.solve_make_square_s,
            source.solve_iterate_s,
            source.solves,
            source.niters,
            source.converged,
            source.subsrc_s,
        )
    end
    return nothing
end
