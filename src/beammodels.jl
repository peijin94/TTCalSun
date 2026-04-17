# Adapted from TTCal.jl beam model definitions for CPU-only TTCalSun use.

abstract type AbstractBeamModel end

struct ConstantBeam <: AbstractBeamModel end

struct SineBeam <: AbstractBeamModel
    α::Float64
end

SineBeam() = SineBeam(1.6)

struct LWA178Beam <: AbstractBeamModel end

const Memo178Beam = LWA178Beam

beam_name(::ConstantBeam) = "constant"
beam_name(beam::SineBeam) = @sprintf("sine(%.3f)", beam.α)
beam_name(::LWA178Beam) = "lwa178"

beam_jones(::ConstantBeam, ν::Float64, az::Float64, el::Float64) = (
    1.0 + 0.0im,
    0.0 + 0.0im,
    0.0 + 0.0im,
    1.0 + 0.0im,
)

function beam_jones(beam::SineBeam, ν::Float64, az::Float64, el::Float64)
    gain = sin(el)^(beam.α / 2)
    value = ComplexF64(gain)
    return (value, 0.0 + 0.0im, 0.0 + 0.0im, value)
end

const _E178 = [-4.529931167425190e+01 -3.066691727279789e+01 +7.111192148086860e+01 +1.131338637814271e+01
               +1.723596273204143e+02 +1.372536555724785e+02 -2.664504470520252e+02 -3.493942140370373e+01
               -2.722311453669980e+02 -2.368121452949910e+02 +4.343953141004854e+02 +5.159758406047006e+01
               +2.408402047219155e+02 +2.302040670131884e+02 -3.993175413350526e+02 -4.152532958513045e+01
               -1.334589043702679e+02 -1.414814115813401e+02 +2.312015990805256e+02 +2.016753495756661e+01
               +4.917278320096442e+01 +5.820262041648866e+01 -8.952191930147437e+01 -6.170044495806915e+00
               -1.244683117802972e+01 -1.652458005311876e+01 +2.396407975769994e+01 +1.185153343629695e+00
               +2.195017889732819e+00 +3.281357788749097e+00 -4.500463893128127e+00 -1.311234364955677e-01
               -2.690339381925372e-01 -4.547844648614167e-01 +5.919912952153034e-01 +4.850030701166225e-03
               +2.243604641077215e-02 +4.311447324076417e-02 -5.345256010140458e-02 +7.056536178294128e-04
               -1.211643367267544e-03 -2.665185970775750e-03 +3.157455237232649e-03 -1.102532965104807e-04
               +3.810179416998095e-05 +9.681257383030349e-05 -1.099243217364324e-04 +6.250327364174989e-06
               -5.277176362640874e-07 -1.567746596003443e-06 +1.710529173897111e-06 -1.350506188926788e-07]

const _H178 = [+4.062920357822495e+02 +3.038713068453467e+01
               -1.706845366994521e+03 -1.337217221068207e+02
               +3.095438045596764e+03 +2.567017051330929e+02
               -3.164514198869798e+03 -2.757538828204631e+02
               +2.035008485840167e+03 +1.850998851042284e+02
               -8.713893456840954e+02 -8.227855478515090e+01
               +2.564477521133275e+02 +2.502479051323083e+01
               -5.262251671400085e+01 -5.287316195660030e+00
               +7.519512104996896e+00 +7.754610944702455e-01
               -7.337746902310935e-01 -7.744604652809532e-02
               +4.663414310713179e-02 +5.024254273102379e-03
               -1.740005497709271e-03 -1.908989166200470e-04
               +2.892116885178882e-05 +3.223985512686652e-06]

@eval function E178(ν::Float64, el::Float64)
    x = ν / 10e6
    θ = π / 2 - el
    α = @evalpoly(x, $(_E178[:, 1]...))
    β = @evalpoly(x, $(_E178[:, 2]...))
    γ = @evalpoly(x, $(_E178[:, 3]...))
    δ = @evalpoly(x, $(_E178[:, 4]...))
    return (1 - (θ / (π / 2))^α) * cos(θ)^β + γ * (θ / (π / 2)) * cos(θ)^δ
end

@eval function H178(ν::Float64, el::Float64)
    x = ν / 10e6
    θ = π / 2 - el
    α = @evalpoly(x, $(_H178[:, 1]...))
    β = @evalpoly(x, $(_H178[:, 2]...))
    return (1 - (θ / (π / 2))^α) * cos(θ)^β
end

function P178(ν::Float64, az::Float64, el::Float64)
    return sqrt((E178(ν, el) * cos(az))^2 + (H178(ν, el) * sin(az))^2)
end

function beam_jones(::LWA178Beam, ν::Float64, az::Float64, el::Float64)
    x = sqrt(max(P178(ν, az, el), 0.0))
    y = sqrt(max(P178(ν, az + π / 2, el), 0.0))
    return (ComplexF64(x), 0.0 + 0.0im, 0.0 + 0.0im, ComplexF64(y))
end

function apply_beam_to_coherency(
    beam::AbstractBeamModel,
    ν::Float64,
    az::Float64,
    el::Float64,
    xx::ComplexF64,
    xy::ComplexF64,
    yx::ComplexF64,
    yy::ComplexF64,
)
    jxx, jxy, jyx, jyy = beam_jones(beam, ν, az, el)

    jk_xx = jxx * xx + jxy * yx
    jk_xy = jxx * xy + jxy * yy
    jk_yx = jyx * xx + jyy * yx
    jk_yy = jyx * xy + jyy * yy

    return (
        jk_xx * conj(jxx) + jk_xy * conj(jxy),
        jk_xx * conj(jyx) + jk_xy * conj(jyy),
        jk_yx * conj(jxx) + jk_yy * conj(jxy),
        jk_yx * conj(jyx) + jk_yy * conj(jyy),
    )
end

function parse_beam_model(spec::AbstractString)
    value = lowercase(strip(spec))
    if value == "constant"
        return ConstantBeam()
    elseif value == "lwa178" || value == "memo178"
        return LWA178Beam()
    elseif value == "sine"
        return SineBeam()
    elseif startswith(value, "sine:")
        return SineBeam(parse(Float64, split(value, ":", limit=2)[2]))
    end
    error("unsupported beam model: $spec")
end
