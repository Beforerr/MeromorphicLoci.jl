module MeromorphicLoci

export findzeros, findpoles, findzeropoles, findloci, ArgumentPrinciple
export Survey, Branch, Sample, winding
@static if VERSION >= v"1.11.0-DEV.469"
    eval(Meta.parse("public nevals"))
end

"""
    findzeros(f, zdomain, alg)

Zeros of `f(z)` over `zdomain`, found by `alg` ([`ArgumentPrinciple`](@ref)). The result
is described in the README, [`MeromorphicLoci`](@ref).
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

include("winding.jl")
include("argument_principle.jl")
include("link.jl")

end
