const _tables = PyNULL()
const _np = PyNULL()

function init_pycasacore()
    if _tables == PyNULL()
        copy!(_tables, pyimport("casacore.tables"))
        copy!(_np, pyimport("numpy"))
    end
    return nothing
end

function read_ms(ms_path::AbstractString; column::AbstractString="DATA")
    init_pycasacore()
    ms = _tables.table(ms_path, readonly=true)

    raw_data = ms.getcol(column)
    raw_flags = ms.getcol("FLAG")
    ant1 = ms.getcol("ANTENNA1")
    ant2 = ms.getcol("ANTENNA2")
    raw_uvw = ms.getcol("UVW")
    nrows = length(ant1)

    spw = _tables.table(ms_path * "/SPECTRAL_WINDOW", readonly=true)
    chan_freq = spw.getcol("CHAN_FREQ")
    spw.close()
    channels = if ndims(chan_freq) == 1
        Float64.(collect(chan_freq))
    elseif size(chan_freq, 1) == 1
        vec(Float64.(chan_freq[1, :]))
    else
        vec(Float64.(chan_freq[:, 1]))
    end

    ant_table = _tables.table(ms_path * "/ANTENNA", readonly=true)
    positions = ant_table.getcol("POSITION")
    ant_table.close()
    antenna_positions = if ndims(positions) == 2 && size(positions, 1) == 3
        Float64.(positions)
    elseif ndims(positions) == 2
        Float64.(permutedims(positions))
    else
        reshape(Float64.(positions), 3, :)
    end

    field_table = _tables.table(ms_path * "/FIELD", readonly=true)
    phase_dir = field_table.getcol("PHASE_DIR")
    field_table.close()
    phase_center_ra, phase_center_dec = if ndims(phase_dir) == 3
        Float64(phase_dir[1, 1, 1]), Float64(phase_dir[1, 1, 2])
    elseif ndims(phase_dir) == 2 && size(phase_dir, 2) == 2
        Float64(phase_dir[1, 1]), Float64(phase_dir[1, 2])
    elseif ndims(phase_dir) == 2
        Float64(phase_dir[1, 1]), Float64(phase_dir[2, 1])
    elseif ndims(phase_dir) == 1 && length(phase_dir) >= 2
        Float64(phase_dir[1]), Float64(phase_dir[2])
    else
        0.0, π / 2.0
    end
    ms.close()

    uvw = if ndims(raw_uvw) == 2 && size(raw_uvw, 2) == 3
        Float64.(permutedims(raw_uvw))
    elseif ndims(raw_uvw) == 2
        Float64.(raw_uvw)
    else
        reshape(Float64.(raw_uvw), 3, nrows)
    end

    baseline_dict = Dict{Tuple{Int, Int}, Int}()
    baseline_row = Dict{Int, Int}()
    row_to_baseline = zeros(Int, nrows)
    for row in 1:nrows
        a1 = Int(ant1[row]) + 1
        a2 = Int(ant2[row]) + 1
        if a1 != a2
            key = (min(a1, a2), max(a1, a2))
            if !haskey(baseline_dict, key)
                idx = length(baseline_dict) + 1
                baseline_dict[key] = idx
                baseline_row[idx] = row
            end
            row_to_baseline[row] = baseline_dict[key]
        end
    end

    nbase_local = length(baseline_dict)
    baselines = zeros(Int32, 2, nbase_local)
    uvw_unique = zeros(Float64, 3, nbase_local)
    baseline_lengths = zeros(Float64, nbase_local)
    for ((a1, a2), idx) in baseline_dict
        baselines[1, idx] = a1
        baselines[2, idx] = a2
        uvw_unique[:, idx] = uvw[:, baseline_row[idx]]
        dx = antenna_positions[1, a1] - antenna_positions[1, a2]
        dy = antenna_positions[2, a1] - antenna_positions[2, a2]
        dz = antenna_positions[3, a1] - antenna_positions[3, a2]
        baseline_lengths[idx] = sqrt(dx * dx + dy * dy + dz * dz)
    end

    shape = size(raw_data)
    nfreq_local = if length(shape) == 3 && (shape[3] == 4 || shape[3] == 2)
        shape[2]
    elseif length(shape) == 3 && (shape[1] == 4 || shape[1] == 2)
        shape[2]
    else
        length(channels)
    end
    if length(channels) != nfreq_local
        channels = collect(range(channels[1], channels[end], length=nfreq_local))
    end

    vis = Visibilities(nbase_local, nfreq_local)
    if length(shape) == 3 && (shape[1] == 4 || shape[1] == 2)
        @inbounds for row in 1:nrows
            α = row_to_baseline[row]
            α == 0 && continue
            if shape[1] == 4
                for β in 1:nfreq_local
                    vis.xx[α, β] = raw_data[1, β, row]
                    vis.xy[α, β] = raw_data[2, β, row]
                    vis.yx[α, β] = raw_data[3, β, row]
                    vis.yy[α, β] = raw_data[4, β, row]
                    vis.flags[α, β] = raw_flags[1, β, row] | raw_flags[2, β, row] | raw_flags[3, β, row] | raw_flags[4, β, row]
                end
            else
                for β in 1:nfreq_local
                    vis.xx[α, β] = raw_data[1, β, row]
                    vis.yy[α, β] = raw_data[2, β, row]
                    vis.flags[α, β] = raw_flags[1, β, row] | raw_flags[2, β, row]
                end
            end
        end
    elseif length(shape) == 3 && (shape[3] == 4 || shape[3] == 2)
        @inbounds for row in 1:nrows
            α = row_to_baseline[row]
            α == 0 && continue
            if shape[3] == 4
                for β in 1:nfreq_local
                    vis.xx[α, β] = raw_data[row, β, 1]
                    vis.xy[α, β] = raw_data[row, β, 2]
                    vis.yx[α, β] = raw_data[row, β, 3]
                    vis.yy[α, β] = raw_data[row, β, 4]
                    vis.flags[α, β] = raw_flags[row, β, 1] | raw_flags[row, β, 2] | raw_flags[row, β, 3] | raw_flags[row, β, 4]
                end
            else
                for β in 1:nfreq_local
                    vis.xx[α, β] = raw_data[row, β, 1]
                    vis.yy[α, β] = raw_data[row, β, 2]
                    vis.flags[α, β] = raw_flags[row, β, 1] | raw_flags[row, β, 2]
                end
            end
        end
    else
        error("unsupported DATA shape $(size(raw_data))")
    end

    meta = Metadata(
        antenna_positions,
        baselines,
        channels,
        phase_center_ra,
        phase_center_dec,
        uvw_unique,
        baseline_lengths,
    )
    return vis, meta, baseline_dict, nrows
end

function write_ms!(ms_path::AbstractString, vis::Visibilities, baseline_dict::Dict{Tuple{Int, Int}, Int}, nrows::Int; column::AbstractString="DATA")
    init_pycasacore()
    ms = _tables.table(ms_path, readonly=false)
    ant1 = ms.getcol("ANTENNA1")
    ant2 = ms.getcol("ANTENNA2")
    existing = ms.getcol(column)
    flags = ms.getcol("FLAG")

    row_to_baseline = zeros(Int, nrows)
    for row in 1:nrows
        a1 = Int(ant1[row]) + 1
        a2 = Int(ant2[row]) + 1
        if a1 != a2
            row_to_baseline[row] = get(baseline_dict, (min(a1, a2), max(a1, a2)), 0)
        end
    end

    shape = size(existing)
    new_data = similar(existing)
    copyto!(new_data, existing)
    new_flags = similar(flags)
    copyto!(new_flags, flags)

    if shape[1] == 4 || shape[1] == 2
        @inbounds for row in 1:nrows
            α = row_to_baseline[row]
            α == 0 && continue
            if shape[1] == 4
                for β in 1:size(vis.xx, 2)
                    new_data[1, β, row] = vis.xx[α, β]
                    new_data[2, β, row] = vis.xy[α, β]
                    new_data[3, β, row] = vis.yx[α, β]
                    new_data[4, β, row] = vis.yy[α, β]
                    f = vis.flags[α, β]
                    new_flags[1, β, row] = f
                    new_flags[2, β, row] = f
                    new_flags[3, β, row] = f
                    new_flags[4, β, row] = f
                end
            else
                for β in 1:size(vis.xx, 2)
                    new_data[1, β, row] = vis.xx[α, β]
                    new_data[2, β, row] = vis.yy[α, β]
                    f = vis.flags[α, β]
                    new_flags[1, β, row] = f
                    new_flags[2, β, row] = f
                end
            end
        end
    elseif shape[3] == 4 || shape[3] == 2
        @inbounds for row in 1:nrows
            α = row_to_baseline[row]
            α == 0 && continue
            if shape[3] == 4
                for β in 1:size(vis.xx, 2)
                    new_data[row, β, 1] = vis.xx[α, β]
                    new_data[row, β, 2] = vis.xy[α, β]
                    new_data[row, β, 3] = vis.yx[α, β]
                    new_data[row, β, 4] = vis.yy[α, β]
                    f = vis.flags[α, β]
                    new_flags[row, β, 1] = f
                    new_flags[row, β, 2] = f
                    new_flags[row, β, 3] = f
                    new_flags[row, β, 4] = f
                end
            else
                for β in 1:size(vis.xx, 2)
                    new_data[row, β, 1] = vis.xx[α, β]
                    new_data[row, β, 2] = vis.yy[α, β]
                    f = vis.flags[α, β]
                    new_flags[row, β, 1] = f
                    new_flags[row, β, 2] = f
                end
            end
        end
    else
        ms.close()
        error("unsupported DATA shape $(size(existing))")
    end

    ms.putcol(column, _np.array(new_data, dtype=_np.complex128))
    ms.putcol("FLAG", _np.array(new_flags, dtype=_np.bool_))
    ms.flush()
    ms.close()
    return nothing
end
