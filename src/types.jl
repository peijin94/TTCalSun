struct Visibilities{T<:AbstractMatrix{ComplexF64}, F<:AbstractMatrix{Bool}}
    xx::T
    xy::T
    yx::T
    yy::T
    flags::F
end

function Visibilities(nbase::Int, nfreq::Int)
    Visibilities(
        zeros(ComplexF64, nbase, nfreq),
        zeros(ComplexF64, nbase, nfreq),
        zeros(ComplexF64, nbase, nfreq),
        zeros(ComplexF64, nbase, nfreq),
        falses(nbase, nfreq),
    )
end

Base.copy(vis::Visibilities) = Visibilities(
    copy(vis.xx),
    copy(vis.xy),
    copy(vis.yx),
    copy(vis.yy),
    copy(vis.flags),
)

nbase(vis::Visibilities) = size(vis.xx, 1)
nfreq(vis::Visibilities) = size(vis.xx, 2)

struct SquareVisibilities{T<:Array{ComplexF64, 3}}
    xx::T
    xy::T
    yx::T
    yy::T
end

function SquareVisibilities(nant::Int, nfreq::Int)
    SquareVisibilities(
        zeros(ComplexF64, nant, nant, nfreq),
        zeros(ComplexF64, nant, nant, nfreq),
        zeros(ComplexF64, nant, nant, nfreq),
        zeros(ComplexF64, nant, nant, nfreq),
    )
end

struct Calibration{T<:AbstractMatrix{ComplexF64}, F<:AbstractMatrix{Bool}}
    xx::T
    xy::T
    yx::T
    yy::T
    flags::F
    is_diagonal::Bool
end

function Calibration(nant::Int, nfreq::Int; diagonal::Bool)
    Calibration(
        ones(ComplexF64, nant, nfreq),
        zeros(ComplexF64, nant, nfreq),
        zeros(ComplexF64, nant, nfreq),
        ones(ComplexF64, nant, nfreq),
        falses(nant, nfreq),
        diagonal,
    )
end

Base.copy(cal::Calibration) = Calibration(
    copy(cal.xx),
    copy(cal.xy),
    copy(cal.yx),
    copy(cal.yy),
    copy(cal.flags),
    cal.is_diagonal,
)

nant(cal::Calibration) = size(cal.xx, 1)
nfreq(cal::Calibration) = size(cal.xx, 2)

struct Metadata
    antenna_positions::Matrix{Float64}
    baselines::Matrix{Int32}
    channels::Vector{Float64}
    phase_center_ra::Float64
    phase_center_dec::Float64
    uvw::Matrix{Float64}
    baseline_lengths::Vector{Float64}
end

nant(meta::Metadata) = size(meta.antenna_positions, 2)
nfreq(meta::Metadata) = length(meta.channels)
nbase(meta::Metadata) = size(meta.baselines, 2)

function fillzero!(vis::Visibilities)
    fill!(vis.xx, 0.0 + 0.0im)
    fill!(vis.xy, 0.0 + 0.0im)
    fill!(vis.yx, 0.0 + 0.0im)
    fill!(vis.yy, 0.0 + 0.0im)
    fill!(vis.flags, false)
    return vis
end

function Base.copyto!(dst::Visibilities, src::Visibilities)
    Base.copyto!(dst.xx, src.xx)
    Base.copyto!(dst.xy, src.xy)
    Base.copyto!(dst.yx, src.yx)
    Base.copyto!(dst.yy, src.yy)
    Base.copyto!(dst.flags, src.flags)
    return dst
end

function fillzero!(sq::SquareVisibilities)
    fill!(sq.xx, 0.0 + 0.0im)
    fill!(sq.xy, 0.0 + 0.0im)
    fill!(sq.yx, 0.0 + 0.0im)
    fill!(sq.yy, 0.0 + 0.0im)
    return sq
end
