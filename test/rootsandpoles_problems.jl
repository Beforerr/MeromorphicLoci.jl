# The test problems of RootsAndPoles.jl (its test/*.jl at v3.0.1), with reference loci from
# GRPF (MATLAB) and the SA-GRPF paper, doi:10.1109/TMTT.2023.3238014. Windings are ours: the
# references list multiplicities only for some. Shared by test/test_rootsandpoles.jl and
# benchmark/rootsandpoles.jl.
using SpecialFunctions: besselj, besselj0, besselj1, bessely, bessely0, bessely1
using LinearAlgebra: det
using MeromorphicLoci: ArgumentPrinciple

inbox((lo, hi), z) = real(lo) < real(z) < real(hi) && imag(lo) < imag(z) < imag(hi)

simplefcn(z) = (z - 1) * (z - im)^2 * (z + 1)^3 / (z + im)

function defaultfcn(z)
    ϵᵣ, μᵣ, d = 5 - 2im, 1 - 2im, 1e-2
    k₀ = 2π * 1e9 / 3e8
    return ϵᵣ^2 * z^2 + z^2 * tan(z)^2 - ϵᵣ^2 * (k₀ * d)^2 * (ϵᵣ * μᵣ - 1)
end

const RF_ZA, RF_ZB, RF_ZC = 1 / 2 - sqrt(3) / 6 * im, sqrt(3) / 3 * im, -1 / 2 - sqrt(3) / 6 * im
rationalfcn(z, ϵ = 1e-2) = (z - RF_ZA) * (z - RF_ZB - ϵ) / (z - RF_ZC) / (z - RF_ZB + ϵ)

# Coaxial line partly filled with dielectric: TM/TE hybrid modes, m = 1
function complexmodes(z)
    z *= 10
    f, c = 5e9, 3e8
    μ₀, ϵ₀ = 4e-7π, 1e-9 / 36 / π
    a, b, ϵᵣ₁, ϵᵣ₂, m = 6.35e-3, 10e-3, 10, 1, 1
    ω = 2π * f
    k₀ = ω / c
    γ = complex(real(z), imag(z)) * k₀
    ϵ₁, ϵ₂, μ₁, μ₂ = ϵ₀ * ϵᵣ₁, ϵ₀ * ϵᵣ₂, μ₀, μ₀
    κ₁, κ₂ = sqrt(γ^2 + k₀^2 * ϵᵣ₁), sqrt(γ^2 + k₀^2 * ϵᵣ₂)
    η₁, η₂ = sqrt(μ₁ / ϵ₁), sqrt(μ₂ / ϵ₂)
    Ja1, Ja2, Ya2 = besselj1(κ₁ * a), besselj1(κ₂ * a), bessely1(κ₂ * a)
    Jb2, Yb2 = besselj1(κ₂ * b), bessely1(κ₂ * b)
    dJ(x) = (besselj0(x) - besselj(2, x)) / 2
    dY(x) = (bessely0(x) - bessely(2, x)) / 2
    dJa1, dJa2, dJb2, dYa2, dYb2 = dJ(κ₁ * a), dJ(κ₂ * a), dJ(κ₂ * b), dY(κ₂ * a), dY(κ₂ * b)
    W = [Ja1 0 -Ja2 -Ya2 0 0
         0 Ja1/η₁ 0 0 -Ja2/η₂ -Ya2/η₂
         γ*m*Ja1/(a*κ₁^2) -ω*μ₁*dJa1/(κ₁*η₁) -γ*m*Ja2/(a*κ₂^2) -γ*m*Ya2/(a*κ₂^2) ω*μ₂*dJa2/(κ₂*η₂) ω*μ₂*dYa2/(κ₂*η₂)
         -ω*ϵ₁*dJa1/κ₁ -m*γ*Ja1/(a*κ₁^2*η₁) ω*ϵ₂*dJa2/κ₂ ω*ϵ₂*dYa2/κ₂ m*γ*Ja2/(a*κ₂^2*η₂) m*γ*Ya2/(a*κ₂^2*η₂)
         0 0 Jb2 Yb2 0 0
         0 0 γ*m*Jb2/(b*κ₂^2) γ*m*Yb2/(b*κ₂^2) -ω*μ₂*dJb2/(κ₂*η₂) -ω*μ₂*dYb2/(κ₂*η₂)]
    return det(W)
end

# Lossy multilayered planar waveguide: the 2×2 determinant of its transfer matrix
function wvgd(z)
    ns, n1, nc = 0.065 - 4im, 1.5835, 1.0
    k₀d1 = 2π / 0.6328e-6 * 1.81e-6
    κ1, γs, γc = sqrt(n1^2 - z^2), sqrt(z^2 - ns^2), sqrt(z^2 - nc^2)
    m11, m12 = cos(κ1 * k₀d1), im / κ1 * sin(κ1 * k₀d1)
    m21, m22 = im * κ1 * sin(κ1 * k₀d1), cos(κ1 * k₀d1)
    return (-m21 + im * γc * m22) - im * γs * (-m11 + im * γc * m12)
end

# Graphene transmission line; the product over the four sheets of its two square roots is
# single-valued
function graphenefunction(z)
    f, c, μ₀ = 1e12, 299792458.0, 4π * 1e-7
    ϵ₀ = 1 / (μ₀ * c^2)
    e, kB, hk, vFe = 1.602176565e-19, 1.3806488e-23, 1.05457168e-34, 1e6
    muc, t, T, ϵᵣ₁, ϵᵣ₂ = 0.05e, 0.135e-12, 300, 1.0, 11.9
    ω = 2π * f
    k₀ = ω / c
    kᵣ₀ = -im * z * k₀
    Slo = -im * e^2 * kB * T * log(2 + 2cosh(muc / kB / T)) / (π * hk^2 * (ω - im / t))
    a = -3vFe^2 * Slo / (4(ω - im / t)^2)
    Y1 = ω * ϵᵣ₁ * ϵ₀ / sqrt(ϵᵣ₁ * k₀^2 - kᵣ₀^2)
    Y2 = ω * ϵᵣ₂ * ϵ₀ / sqrt(ϵᵣ₂ * k₀^2 - kᵣ₀^2)
    YS = Slo + (4 / 3) * a * kᵣ₀^2
    return (Y1 + Y2 + YS) * (-Y1 + Y2 + YS) * (Y1 - Y2 + YS) * (-Y1 - Y2 + YS)
end

# `loci` are (point, winding) pairs. `ap` is an `ArgumentPrinciple` and its box (`nothing`:
# `box`): `minsep` below the closest zero–pole spacing; for graphene, below its zeros 0.04
# from a double pole, affordable only in a box around them.
const RP_PROBLEMS = [
    (; name = "simple rational", f = simplefcn, box = (-2 - 2im, 2 + 2im), ap = (ArgumentPrinciple(zres = 0.02, minsep = 0.4), nothing),
        loci = [(-1.0 + 0im, 3), (1.0im, 2), (1.0 + 0im, 1), (-1.0im, -1)]),
    (; name = "default (tan)", f = defaultfcn, box = (-2 - 2im, 2 + 2im), ap = (ArgumentPrinciple(zres = 0.02, minsep = 0.4), nothing),
        loci = [(-1.624715288135189 + 0.182095877702038im, 1), (-1.520192978034417 - 0.173670452237129im, 1),
                (-0.515113098919392 + 0.507111597359180im, 1), (0.515113098795215 - 0.507111597284675im, 1),
                (1.520192978034417 + 0.173670452237129im, 1), (1.624715288135189 - 0.182095877702037im, 1),
                (-π / 2 + 0im, -2), (π / 2 + 0im, -2)]),
    (; name = "rational, a zero–pole pair 0.02 apart", f = rationalfcn, box = (-1 - 1im, 1 + 1im),
        ap = (ArgumentPrinciple(zres = 0.005, minsep = 0.0133), nothing),
        loci = [(RF_ZA, 1), (RF_ZB + 1e-2, 1), (RF_ZC, -1), (RF_ZB - 1e-2, -1)]),
    (; name = "complex modes", f = complexmodes, box = (-1 - 1im, 1 + 1im), ap = (ArgumentPrinciple(zres = 0.01, minsep = 0.2), nothing),
        loci = [vec([(s1 * 0.09664230246 + s2 * 0.06292339746im, 1) for s1 in (-1, 1), s2 in (-1, 1)]);
                vec([(s * x + 0im, 1) for x in (0.44442904311, 0.70377225022, 0.77502152220, 0.85611520391), s in (-1, 1)]);
                [(-0.1im, -2), (0.1im, -2)]]),
    (; name = "lossy multilayered waveguide", f = wvgd, box = (1 - 1im, 2.5 + 1im), ap = (ArgumentPrinciple(zres = 0.02, minsep = 0.4), nothing),
        loci = [(z, 1) for z in (1.09675254341 - 0.00019714688im, 1.24045447136 - 0.00013382215im,
                1.35314042918 - 0.00008613919im, 1.43979554425 - 0.00005200167im, 1.50416986640 - 0.00002802944im,
                1.54869224388 - 0.00001210101im, 1.57486304575 - 0.00000297462im)]),
    # its multiscale box, (-100 - 100im, 400 + 400im), defeats both point finders at an
    # affordable cost: see the README. Here the box around the origin, where zeros and
    # double poles crowd.
    (; name = "graphene transmission line", f = graphenefunction, box = (-5 - 5im, 5 + 5im),
        ap = (ArgumentPrinciple(zres = 0.004, minsep = 0.03), (-0.25 + 0.75im, 0.25 + 1.25im)),
        loci = [(-0.003206780 - 0.964810358im, 1), (0.003206780 + 0.964810358im, 1),
                (-0.004526720 + 0.955901829im, 1), (0.004526720 - 0.955901829im, 1),
                (-1.0im, -2), (1.0im, -2), (-3.449637662im, -2), (3.449637662im, -2)]),
]
