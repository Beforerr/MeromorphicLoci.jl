# Discrete fits of the same samples: MeromorphicLoci.aaa vs RationalFunctionApproximation
# (Barycentric, Thiele) vs BaryRational. Reports degree, off-grid fit error, pole error, time.
using MeromorphicLoci, Chairmarks, Logging
using MeromorphicLoci: aaa, poles
import RationalFunctionApproximation as RFA, BaryRational as BR
disable_logging(Logging.Info)

relerr(g, f, T) = maximum(abs(g(z) - f(z)) for z in T) / maximum(abs ∘ f, T)
poleerr(found, exact) = isempty(exact) ? NaN : maximum(p -> minimum(abs, found .- p), exact)

const tol = 1.0e-13
fitters = [
    "ours" => (y, z) -> (r = aaa(y, z; tol).fit; (length(r.nodes) - 1, r, () -> poles(r))),
    "RFA Barycentric" => (y, z) -> (a = RFA.approximate(y, z, RFA.Barycentric(); tol, max_degree = 150, max_iter = 151).fun; (RFA.degree(a), a, () -> RFA.poles(a))),
    "RFA Thiele" => (y, z) -> (a = RFA.approximate(y, z, RFA.Thiele(); tol, max_iter = 151).fun; (RFA.degree(a), a, () -> RFA.poles(a))),
    "BaryRational" => (y, z) -> (a = BR.aaa(z, y; tol, mmax = 151, clean = 0); (length(a.x) - 1, a, () -> BR.prz(a)[1])),
]

cases = [
    ("4 poles", z -> sum(1 / (z - p) for p in (0.3 + 0.6im, -1.1 + 0.5im, 2.0 - 0.3im, 0.1 - 0.7im)), [0.3 + 0.6im, -1.1 + 0.5im, 2.0 - 0.3im, 0.1 - 0.7im]),
    ("1/cos(3z²)+exp(2z)", z -> 1 / cos(3z^2) + exp(2z), [s * sqrt(complex((k + 0.5) * π / 3)) for k in -6:5 for s in (1, -1)]),
]

for (name, f, P) in cases, n in ((20, 16), (60, 40))
    Z = vec(range(-3, 3, n[1]) .+ im .* range(-1, 1, n[2])')
    T = range(-2.93, 2.91, 37) .+ im .* range(-0.97, 0.95, 23)'
    exact = filter(p -> abs(real(p)) < 2.5 && abs(imag(p)) < 0.8, P)
    y = f.(Z)
    println("== $name, $(n[1])×$(n[2]) samples")
    for (label, fit) in fitters
        deg, r, ps = fit(y, Z)
        t = @b fit($y, $Z)
        println(rpad(label, 16), " degree ", lpad(deg, 3), "  fit err ", round(relerr(r, f, T); sigdigits = 2),
            "  pole err ", round(poleerr(ps(), exact); sigdigits = 2), "  ", round(t.time * 1.0e3; sigdigits = 3), " ms")
    end
end
