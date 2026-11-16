# The test problems of RootsAndPoles.jl (test/rootsandpoles_problems.jl): evaluations,
# warm wall time and accuracy of AAA, ArgumentPrinciple and RootsAndPoles' SA-GRPF.
# Run: julia --project=benchmark benchmark/rootsandpoles.jl
using MeromorphicLoci
using RootsAndPoles
using Printf
include(joinpath(@__DIR__, "..", "test", "rootsandpoles_problems.jl"))

counting(f) = (n = Ref(0); (z -> (n[] += 1; f(z)), n))

# References matched by a point of their winding within 1e-6 relative, of the points
# returned, and the largest relative error of a matched point
function score((points, windings), loci)
    errs = map(loci) do (z, w)
        d = [windings[i] == w ? abs(points[i] - z) / max(1, abs(z)) : Inf for i in eachindex(points)]
        minimum(d; init = Inf)
    end
    matched = errs .< 1e-6
    return "$(count(matched))/$(length(loci)) of $(length(points))", maximum(errs[matched]; init = 0.0)
end

ours(alg) = (f, box) -> (r = findzeropoles(f, box, alg); (r.points, r.windings))
function sagrpf(f, (lo, hi))
    corners = [lo, complex(real(hi), imag(lo)), hi, complex(real(lo), imag(hi))]
    zr, zp, mr, mp = rootsandpoles(f, ComplexMesh(corners), ReturnMultiplicity())
    return [zr; zp], round.(Int, [mr; mp])
end

@printf("%-38s %-18s %8s %9s %12s %9s\n", "problem", "method", "evals", "wall (ms)", "found", "max err")
for p in RP_PROBLEMS
    apalg, apbox = p.ap
    for (label, run, box) in (("AAA()", ours(AAA()), p.box),
                              ("ArgumentPrinciple", ours(apalg), something(apbox, p.box)),
                              ("RootsAndPoles", sagrpf, p.box))
        g, n = counting(p.f)
        found = run(g, box)
        t = minimum(_ -> @elapsed(run(p.f, box)), 1:3)
        found, err = score(found, filter(l -> inbox(box, l[1]), p.loci))
        @printf("%-38s %-18s %8d %9.1f %12s %9.1e%s\n", p.name, label, n[], 1e3t, found, err,
            box == p.box ? "" : "  (box $box)")
    end
end
