const _py_json = PyNULL()

function _init_json()
    if _py_json == PyNULL()
        copy!(_py_json, pyimport("json"))
    end
    return _py_json
end

struct PowerLawSpectrum
    I::Float64
    Q::Float64
    U::Float64
    V::Float64
    reference_frequency::Float64
    index::Vector{Float64}
end

function (spectrum::PowerLawSpectrum)(frequency::Float64)
    ratio = log(frequency / spectrum.reference_frequency)
    log_flux = 0.0
    for (i, idx) in enumerate(spectrum.index)
        log_flux += idx * ratio^i
    end
    scale = exp(log_flux)
    return (
        spectrum.I * scale,
        spectrum.Q * scale,
        spectrum.U * scale,
        spectrum.V * scale,
    )
end

abstract type AbstractSource end

struct PointSource <: AbstractSource
    name::String
    ra::Float64
    dec::Float64
    spectrum::PowerLawSpectrum
end

struct GaussianSource <: AbstractSource
    name::String
    ra::Float64
    dec::Float64
    spectrum::PowerLawSpectrum
    major_fwhm::Float64
    minor_fwhm::Float64
    position_angle::Float64
end

struct MultiSource <: AbstractSource
    name::String
    components::Vector{AbstractSource}
end

get_name(source::AbstractSource) = source.name

function read_sources(filename::AbstractString)
    json = _init_json()
    raw = pycall(json.loads, PyAny, read(filename, String))
    return [construct_source(Dict(entry)) for entry in raw]
end

function construct_source(entry::Dict)
    name = String(get(entry, "name", ""))
    if haskey(entry, "components")
        components = AbstractSource[]
        for component in entry["components"]
            push!(components, construct_source(Dict(component)))
        end
        return MultiSource(name, components)
    end

    ra, dec = parse_source_direction(entry)
    spectrum = parse_source_spectrum(entry)
    if haskey(entry, "major-fwhm") && haskey(entry, "minor-fwhm") && haskey(entry, "position-angle")
        return GaussianSource(
            name,
            ra,
            dec,
            spectrum,
            deg2rad(Float64(entry["major-fwhm"]) / 3600.0),
            deg2rad(Float64(entry["minor-fwhm"]) / 3600.0),
            deg2rad(Float64(entry["position-angle"])),
        )
    end
    return PointSource(name, ra, dec, spectrum)
end

function parse_source_direction(entry::Dict)
    name = String(get(entry, "name", ""))
    if name == "Sun" || name == "Moon" || name == "Jupiter"
        error("moving sources are not supported in TTCalSun yet: $name")
    end
    if haskey(entry, "ra") && haskey(entry, "dec")
        return parse_angle(entry["ra"]), parse_angle(entry["dec"])
    end
    error("source $name is missing ra/dec")
end

function parse_source_spectrum(entry::Dict)
    return PowerLawSpectrum(
        Float64(entry["I"]),
        Float64(get(entry, "Q", 0.0)),
        Float64(get(entry, "U", 0.0)),
        Float64(get(entry, "V", 0.0)),
        Float64(entry["freq"]),
        Float64.(entry["index"]),
    )
end

function parse_angle(value)
    value isa Number && return Float64(value)
    s = String(value)

    m = match(r"([+-]?\d+(?:\.\d+)?)[hH](\d+(?:\.\d+)?)[mM](\d+(?:\.\d+)?)[sS]?", s)
    if m !== nothing
        hours = parse(Float64, m.captures[1])
        minutes = parse(Float64, m.captures[2])
        seconds = parse(Float64, m.captures[3])
        sign = hours >= 0 ? 1.0 : -1.0
        return sign * (abs(hours) + minutes / 60.0 + seconds / 3600.0) * (π / 12.0)
    end

    m = match(r"([+-]?\d+(?:\.\d+)?)[dD](\d+(?:\.\d+)?)[mM']?(\d+(?:\.\d+)?)[sS\"]?", s)
    if m !== nothing
        degrees = parse(Float64, m.captures[1])
        minutes = parse(Float64, m.captures[2])
        seconds = parse(Float64, m.captures[3])
        sign = startswith(strip(s), "-") ? -1.0 : 1.0
        return sign * (abs(degrees) + minutes / 60.0 + seconds / 3600.0) * (π / 180.0)
    end

    m = match(r"([+-]?\d+(?:\.\d+)?)", s)
    if m !== nothing
        return parse(Float64, m.captures[1]) * (π / 180.0)
    end

    error("cannot parse angle: $s")
end

function source_direction_lmn(source::Union{PointSource, GaussianSource}, phase_center_ra::Float64, phase_center_dec::Float64)
    dra = source.ra - phase_center_ra
    l = cos(source.dec) * sin(dra)
    m = sin(source.dec) * cos(phase_center_dec) - cos(source.dec) * sin(phase_center_dec) * cos(dra)
    n = sin(source.dec) * sin(phase_center_dec) + cos(source.dec) * cos(phase_center_dec) * cos(dra)
    return l, m, n
end

function source_direction_lmn(source::MultiSource, phase_center_ra::Float64, phase_center_dec::Float64)
    isempty(source.components) && return 0.0, 0.0, 1.0
    return source_direction_lmn(source.components[1], phase_center_ra, phase_center_dec)
end

function is_above_horizon(source::AbstractSource, phase_center_ra::Float64, phase_center_dec::Float64)
    _, _, n = source_direction_lmn(source, phase_center_ra, phase_center_dec)
    return n > 0.0
end

function source_stokes_i(source::Union{PointSource, GaussianSource}, frequency::Float64)
    I, _, _, _ = source.spectrum(frequency)
    return I
end

function source_stokes_i(source::MultiSource, frequency::Float64)
    total = 0.0
    for component in source.components
        total += source_stokes_i(component, frequency)
    end
    return total
end
