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

# A point-finding algorithm implements `_enclose(f, zdomain, alg, Val(W))`, returning
# `(; loci, nevals, info)`: each locus a centre `z`, a winding `w` (0 when its kind is
# unknown) and a radius `r` that holds no other locus; `info` holds extra result fields.
# Simple loci are polished within `r`.
function _points(f, zdomain, alg, ::Val{W}, polish) where {W}
    (; loci, nevals, info) = _enclose(f, zdomain, alg, Val(W))
    found = @NamedTuple{z::ComplexF64, w::Int, converged::Bool}[]
    unresolved = ComplexF64[]
    for (; z, r, w) in loci
        W * w < 0 && continue
        if w == 0
            push!(unresolved, z)
            continue
        end
        converged = false
        if abs(w) == 1
            z, converged, ne = _polish(polish, w > 0 ? f : inv ∘ f, z, r)
            nevals += ne
        end
        push!(found, (; z, w, converged))
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

# Uncertainty of a locus estimate `z`: `scale`, but under half the way to any other locus,
# so neither the polish nor a circle around `z` reaches it.
_reach(z, others, scale) = min(scale, 0.4 * minimum(o -> abs(o - z), others; init = Inf))

include("muller.jl")
include("winding.jl")
include("argument_principle.jl")
include("link.jl")
include("aaa.jl")

end
