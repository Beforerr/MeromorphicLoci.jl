module MeromorphicLoci

using LinearAlgebra

export findzeros, findpoles, findzeropoles, findloci, AAA, ArgumentPrinciple, Muller
export Survey, Branch, Sample, winding
@static if VERSION >= v"1.11.0-DEV.469"
    eval(Meta.parse("public aaa, poles, Barycentric, nevals, muller"))
end

"""
    findzeros(f, zdomain, alg; polish = Muller())

Zeros of `f(z)` over `zdomain`, found by `alg` ([`ArgumentPrinciple`](@ref) or
[`AAA`](@ref)). The result is described in the README, [`MeromorphicLoci`](@ref).
"""
findzeros(f, zdomain, alg; polish=Muller()) = _points(f, zdomain, alg, Val(1), polish)

"""
    findpoles(f, zdomain, alg; polish = Muller())

Poles of `f(z)`; otherwise as [`findzeros`](@ref).
"""
findpoles(f, zdomain, alg; polish=Muller()) = _points(f, zdomain, alg, Val(-1), polish)

"""
    findzeropoles(f, zdomain, alg; polish = Muller())

Zeros and poles of `f(z)`; otherwise as [`findzeros`](@ref).
"""
findzeropoles(f, zdomain, alg; polish=Muller()) = _points(f, zdomain, alg, Val(0), polish)

# A locus estimate: centre `z`, winding `w` (0 if unknown), radius `r` holding no other locus.
struct Locus
    z::ComplexF64
    r::Float64
    w::Int
end

"""
    _enclose(f, zdomain, alg, ::Val{W}) -> (; loci::Vector{Locus}, nevals, info, pad)

Implemented by each point-finding algorithm `alg`.
"""
function _enclose end

# complex x/0 is NaN while 1/f → 0
_reciprocal(f) = z -> (w = f(z); isfinite(w) ? inv(w) : zero(w))

function _points(f, zdomain, alg, ::Val{W}, polish) where {W}
    (; loci, nevals, info, pad) = _enclose(f, zdomain, alg, Val(W))
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

"""
    findloci(f, zdomain, pdomain, alg)

Loci `z*(p)` of the zeros and poles of `f(p...)(z)` over `zdomain` as `p` ranges over
`pdomain`, by algorithm `alg` ([`ArgumentPrinciple`](@ref)), as a [`Survey`](@ref). The
result is described in the README, [`MeromorphicLoci`](@ref).
"""
findloci(f, zdomain, pdomain, alg) = _loci(f, zdomain, pdomain, alg)

_is_valid(fz) = isfinite(fz) && !iszero(fz)

_inside((lo, hi), z, pad) =
    real(lo) - pad <= real(z) <= real(hi) + pad && imag(lo) - pad <= imag(z) <= imag(hi) + pad

# Uncertainty of a locus estimate `z`: `scale`, but under half the way to any other locus,
# so neither the polish nor a circle around `z` reaches it.
_reach(z, others, scale) = min(scale, 0.4 * minimum(o -> abs(o - z), others; init = Inf))

include("muller.jl")
include("winding.jl")
include("argument_principle.jl")
include("link.jl")
include("aaa.jl")

end
