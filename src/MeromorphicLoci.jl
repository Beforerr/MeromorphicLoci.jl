module MeromorphicLoci

using LinearAlgebra

export findzeros, findpoles, findzeropoles, findloci, AAA, ArgumentPrinciple, Continuation, Muller
@static if VERSION >= v"1.11.0-DEV.469"
    eval(Meta.parse("public aaa, poles, Barycentric, nevals, muller, Survey, Branch, Sample"))
end

"""
    findzeros(f, zdomain, alg; polish = Muller(), mapper = map)

Zeros of `f(z)` in the box `zdomain = (zlo, zhi)` of opposite corners.
Returns a named tuple:

- `points`: sorted by real, then imaginary part.
- `windings`: multiplicity, negative for a pole. Points closer than `alg` resolves merge
  into one of their summed winding.
- `converged`: whether `polish` ([`Muller`](@ref)) met its tolerance. Only points of
  winding ±1 are polished; `polish = nothing` keeps the estimates of `alg`.
- `unresolved`: zeros or poles detected but not classified.
- `nevals`: evaluations of `f`.
- `fitconverged` ([`AAA`](@ref)) or `deficit` ([`ArgumentPrinciple`](@ref)).

Zeros and poles where `f` overflows to `Inf` or underflows to 0 go unfound: divide out
exponential growth, or shrink the box.
"""
findzeros(f, zdomain, alg; polish=Muller(), mapper=map) = _points(f, zdomain, alg, Val(1), polish, mapper)

"""
    findpoles(f, zdomain, alg; polish = Muller(), mapper = map)

Poles of `f(z)`; otherwise as [`findzeros`](@ref).
"""
findpoles(f, zdomain, alg; polish=Muller(), mapper=map) = _points(f, zdomain, alg, Val(-1), polish, mapper)

"""
    findzeropoles(f, zdomain, alg; polish = Muller(), mapper = map)

Zeros and poles of `f(z)` from one search, costing little more than either alone;
otherwise as [`findzeros`](@ref).
"""
findzeropoles(f, zdomain, alg; polish=Muller(), mapper=map) = _points(f, zdomain, alg, Val(0), polish, mapper)

# A locus estimate: centre `z`, winding `w` (0 if unknown), radius `r` holding no other locus.
struct Locus
    z::ComplexF64
    r::Float64
    w::Int
end

"""
    _enclose(f, zdomain, alg, ::Val{W}, mapper) -> (; loci::Vector{Locus}, nevals, info, pad)

Implemented by each point-finding algorithm `alg`.
"""
function _enclose end

# complex x/0 is NaN while 1/f → 0
_reciprocal(f) = z -> (w = f(z); isfinite(w) ? inv(w) : zero(w))

function _points(f, zdomain, alg, ::Val{W}, polish, mapper) where {W}
    (; loci, nevals, info, pad) = _enclose(f, zdomain, alg, Val(W), mapper)
    keep = zdomain isa Tuple{Number,Number} ? (z -> _inside(_zcorners(zdomain...), z, pad)) : Returns(true)
    found = @NamedTuple{z::ComplexF64, w::Int, converged::Bool}[]
    unresolved = ComplexF64[]
    for (; z, r, w) in loci
        W * w < 0 && continue
        if w == 0
            keep(z) && push!(unresolved, z)
            continue
        end
        converged = false
        if abs(w) == 1
            z, converged, ne = _polish(polish, w > 0 ? f : _reciprocal(f), z, r)
            nevals += ne
        end
        keep(z) && push!(found, (; z, w, converged))
    end
    sort!(found; by = t -> reim(t.z))
    sort!(unresolved; by = reim)
    points, windings, converged = [t.z for t in found], [t.w for t in found], [t.converged for t in found]
    return merge((; points, windings, converged, unresolved, nevals), info)
end

_is_valid(fz) = isfinite(fz) && !iszero(fz)

_zcorners(a, b) = (complex(min(real(a), real(b)), min(imag(a), imag(b))),
    complex(max(real(a), real(b)), max(imag(a), imag(b))))

_inside((lo, hi), z, pad) =
    real(lo) - pad <= real(z) <= real(hi) + pad && imag(lo) - pad <= imag(z) <= imag(hi) + pad

# Uncertainty of a locus estimate `z`: `scale`, but under half the way to any other locus,
# so neither the polish nor a circle around `z` reaches it.
_reach(z, others, scale) = min(scale, 0.4 * minimum(o -> abs(o - z), others; init = Inf))

include("muller.jl")
include("winding.jl")
include("survey.jl")
include("argument_principle.jl")
include("labels.jl")
include("aaa.jl")
include("continuation.jl")

end
