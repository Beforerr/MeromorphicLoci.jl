module MeromorphicLoci

using LinearAlgebra

export findzeros, findpoles, findzeropoles, findloci, AAA, ArgumentPrinciple
export Survey, Branch, Sample, winding
@static if VERSION >= v"1.11.0-DEV.469"
    eval(Meta.parse("public aaa, poles, Barycentric, nevals, muller"))
end

"""
    findzeros(f, zdomain, alg)

Zeros of `f(z)` over `zdomain`, found by `alg` ([`ArgumentPrinciple`](@ref) or
[`AAA`](@ref)). The result is described in the README, [`MeromorphicLoci`](@ref).
"""
findzeros(f, zdomain, alg) = _points(f, zdomain, alg, Val(1))

"""
    findpoles(f, zdomain, alg)

Poles of `f(z)`; otherwise as [`findzeros`](@ref).
"""
findpoles(f, zdomain, alg) = _points(f, zdomain, alg, Val(-1))

"""
    findzeropoles(f, zdomain, alg)

Zeros and poles of `f(z)`; otherwise as [`findzeros`](@ref).
"""
findzeropoles(f, zdomain, alg) = _points(f, zdomain, alg, Val(0))

"""
    findloci(f, zdomain, pdomain, alg)

Loci `z*(p)` of the zeros and poles of `f(p...)(z)` over `zdomain` as `p` ranges over
`pdomain`, by algorithm `alg` ([`ArgumentPrinciple`](@ref)), as a [`Survey`](@ref). The
result is described in the README, [`MeromorphicLoci`](@ref).
"""
findloci(f, zdomain, pdomain, alg) = _loci(f, zdomain, pdomain, alg)

_is_valid(fz) = isfinite(fz) && !iszero(fz)

# Uncertainty of a locus estimate `z`: `scale`, but under half the way to any other locus,
# so a circle around `z` does not reach it.
_reach(z, others, scale) = min(scale, 0.4 * minimum(o -> abs(o - z), others; init = Inf))

include("winding.jl")
include("argument_principle.jl")
include("link.jl")
include("aaa.jl")
include("muller.jl")

end
