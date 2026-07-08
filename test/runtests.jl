using MeromorphicLoci
using MeromorphicLoci: quadrant, nevals, _circle_winding
using Test

iszerobranch(b) = winding(b) > 0
ispolebranch(b) = winding(b) < 0

# half-open boundaries: 1:[0,π/2) 2:[π/2,π] 3:(-π,-π/2) 4:[-π/2,0)
@testset "quadrant logic" begin
    @test quadrant(1 + 0.1im) == 1
    @test quadrant(-0.1 + 1im) == 2
    @test quadrant(-1 - 0.1im) == 3
    @test quadrant(0.1 - 1im) == 4
    @test quadrant(1.0 + 0im) == quadrant(1.0 - 0im) == 1
    @test quadrant(1im) == 2
    @test quadrant(-1.0 + 0im) == 2
    @test quadrant(-1im) == 4
    @test quadrant(-1.0 - 0.0im) == quadrant(-1.0 + 0.0im) == 2
end

@testset "adaptive circle winding" begin
    include("test_winding.jl")
end

@testset "README example survives its own asserts" begin
    f(z, p) = (z - 0.80cispi(p)) * (z + 0.10 + 0.12im) /
              ((z - 0.804cispi(-p + 0.3)) * (z + 0.13 + 0.12im))
    box = ((-0.98 - 0.98im, 0.98 + 0.98im), (0.0, 1.0))
    s = survey(f, box; zres=0.02)
    @test sort([winding(b) for b in s]) == [-1, 1]        # graze escalated apart
    @test winding(only(survey(f, box; zres=0.02, zoom=0))) == 0
    @test length(survey(f, box; zres=0.02, minsep=0.3, zoom=0)) == 3
    @test sort([winding(b) for b in survey(f, box; zres=0.02, minsep=0.3)]) == [-1, -1, 1, 1]
end

@testset "each corner costs exactly one f call: neighbouring cells share corners and so share evaluations" begin
    calls = NTuple{3,Float64}[]
    lk = ReentrantLock()
    f(z, t) = (Base.@lock lk push!(calls, (real(z), imag(z), t)); z - (0.5 + 0.3t))
    box = ((0.0 - 0.3im, 1.0 + 0.3im), (0.0, 1.0))

    sc = scan(f, box; zres=0.02)
    @test length(calls) > 0
    @test allunique(calls)
    @test length(calls) == nevals(sc)

    # `keep` gets its own cache and the same guarantee
    kcalls = NTuple{3,Float64}[]
    empty!(calls)
    keep(z, t) = (Base.@lock lk push!(kcalls, (real(z), imag(z), t)); imag(z) > 0)
    sk = scan(f, box; zres=0.02, keep)
    @test allunique(kcalls)
end

@testset "show methods render" begin
    f(z, t) = z - (0.5 + 0.3t)
    sc = scan(f, ((0.0 - 0.3im, 1.0 + 0.3im), (0.0, 1.0)); zres=0.02)
    s = link(sc)
    @test occursin("Scan", repr(sc))
    @test occursin("Branch", repr(s[1])) && occursin("Zero", repr(s[1]))
    @test occursin("Survey", repr("text/plain", s))
    @test occursin("1 zero", repr("text/plain", s))
    @test occursin("Sample", repr(s[1][1]))
end

@testset "scan / link phases" begin
    f(z, t) = z - (0.5 + 0.3t)
    box = ((0.0 - 0.3im, 1.0 + 0.3im), (0.0, 1.0))

    sc = scan(f, box; zres=0.02)

    @test length(sc) == length(sc.cells) > 1
    @test eltype(sc) === Sample{1}

    s = link(sc)
    # linking partitions the scan's cells — none invented, none dropped
    @test sum(length(b.cells) for b in s) == length(sc.cells)
    # and costs only the classification evals on top of the scan's refinement
    @test nevals(s) ≥ nevals(sc) == nevals(scan(f, box; zres=0.02))

    # survey is exactly the two phases composed
    s2 = survey(f, box; zres=0.02)
    @test length(s2) == length(s) && nevals(s2) == nevals(s)
    @test [winding(b) for b in s2] == [winding(b) for b in s]

    # classification knobs thread through link
    s3 = link(sc; nreps=3)
    @test [winding(b) for b in s3] == [winding(b) for b in s]
end

@testset "zres must describe a reachable lattice" begin
    f(z) = z - 0.3
    box = ((-1.0 - 1.0im, 1.0 + 1.0im),)
    @test_throws ArgumentError survey(f, box; zres=0.0)
    @test_throws ArgumentError survey(f, box; zres=-1.0)
    @test_throws ArgumentError survey(f, box; zres=1e-12)   # more levels than Int32 indices span
    @test length(survey(f, box; zres=10.0)) == 1
end

@testset "Survey vector interface" begin
    f(z, t) = z - (0.5 + 0.3t)
    s = survey(f, ((0.0 - 0.3im, 1.0 + 0.3im), (0.0, 1.0)); zres=0.02)
    @test s isa AbstractVector{<:Branch}
    @test collect(s) == s.branches
    @test only(filter(iszerobranch, s)) === s[1]
    @test s[1] isa AbstractVector{<:Sample}
end

@testset "pole curve also fires (zero/pole agnostic)" begin
    # A moving simple pole: the winding criterion is blind to zero-vs-pole.
    f(z, t) = 1 / (z - (0.5 + 0.2t))
    box = ((0.0 - 0.3im, 1.0 + 0.3im), (0.0, 1.0))
    s = survey(f, box; zres=0.02)
    @test length(s) == 1
    @test winding(s[1]) < 0                    # opposite winding sign
    @test all(abs(smp.z - (0.5 + 0.2smp.p[1])) < 0.03 for smp in s[1])
end

@testset "two non-crossing branches stay separate" begin
    # Zeros at 0.3+0.5t and 1.5-0.5t: closest approach 0.4 (at t=1) ≫ zres, so
    # the components must not merge.
    r1(t) = 0.3 + 0.5t
    r2(t) = 1.5 - 0.5t
    f(z, t) = (z - r1(t)) * (z - r2(t))
    box = ((0.0 - 0.4im, 2.0 + 0.4im), (0.0, 1.0))
    s = survey(f, box; zres=0.02)
    @test length(s) == 2
    # each branch tracks exactly one of the two curves
    for b in s
        on1 = all(abs(smp.z - r1(smp.p[1])) < 0.04 for smp in b)
        on2 = all(abs(smp.z - r2(smp.p[1])) < 0.04 for smp in b)
        @test on1 ⊻ on2
    end
end

@testset "`minsep` sets the base mesh, hence the evaluation floor" begin
    f(z, t) = z - (0.5 + 0.3t)
    box = ((0.0 - 0.3im, 1.0 + 0.3im), (0.0, 1.0))
    ev(minsep) = survey(f, box; zres=0.02, minsep).nevals
    # Uniform base mesh over a 3-D box ⇒ halving `minsep` costs sharply more evaluations.
    # Measured ~3.7×, not the naive 2^3: neighbouring base cells share corners and the
    # cache charges each only once. This cost floor is why `minsep` is the knob to reach
    # for last, despite being the only one that recovers a locus a base cell never
    # fired on.
    @test ev(0.4) < ev(0.2) < ev(0.1)
    @test ev(0.1) / ev(0.2) > 3
end

@testset "escalation splits a fused zero/pole pair below zres" begin
    box = ((0.0 - 0.4im, 2.0 + 0.4im), (0.0, 1.0))
    pair(sep) = (z, t) -> (z - (1.0 - sep / 2)) / (z - (1.0 + sep / 2))

    # Separated at cell scale (0.08 ≳ 2 cells): two branches, each settled by its
    # own enclosing faces — correct labels even though the 2zres = 0.1 vote
    # circles would swallow the partner.
    wide = survey(pair(0.08), box; zres=0.05)
    @test sort([sign(winding(b)) for b in wide]) == [-1, 1]

    # Fused (closer than a cell): the loci merge into one branch whose faces are
    # contaminated and whose vote circles enclose both partners — no reshuffling
    # of the zres-level corners can label it. `link` escalates the branch below
    # zres locally until the pair falls apart into clean sub-branches.
    sc = scan(pair(0.02), box; zres=0.05)
    tight = link(sc)
    @test sort([winding(b) for b in tight]) == [-1, 1]
    for b in tight   # escalated samples still sit on their locus
        z = winding(b) > 0 ? 0.99 : 1.01
        @test all(abs(smp.z - z) < 0.05 for smp in b)
    end

    # zoom=0 disables escalation: one merged branch, winding 0 — honest, where a
    # Bool `isroot` would have silently mislabelled the zero as a pole.
    @test winding(only(link(sc; zoom=0))) == 0
end

@testset "escalation stays local to the ambiguous branch (spending only a fraction of the scan's budget)" begin
    f(z, p) = (z - 0.80cispi(p)) * (z + 0.11 + 0.12im) /
              ((z - 0.55cispi(-p + 0.3)) * (z + 0.13 + 0.12im))
    box = ((-0.98 - 0.98im, 0.98 + 0.98im), (0.0, 1.0))
    sc = scan(f, box; zres=0.02, minsep=0.12)
    n0 = nevals(sc)
    s = link(sc)
    @test sort([sign(winding(b)) for b in s]) == [-1, -1, 1, 1]
    @test nevals(s) - n0 < n0 ÷ 4
end

@testset "circle fallback: multiplicity, and veto near the box edge" begin
    # A double zero defeats the per-face winding (every face reads ambiguous), so
    # classification must come from the circle fallback — with true multiplicity.
    # Same-sign votes are trusted as-is: no escalation triggers.
    sc = scan(z -> (z - 0.2)^2, ((-1.0 - 1.0im, 1.0 + 1.0im),); zres=0.02)
    s = link(sc)
    @test length(s) == 1
    @test winding(s[1]) == 2                   # a zero of multiplicity 2
    @test nevals(s) == nevals(link(sc; zoom=0))   # cheap path, zoom untouched

    # Near the box edge every top-level vote circle (radius 2zres) leaves the box,
    # so the vote is withheld rather than evaluating f outside. Escalation halves
    # the circles until they fit and recovers the label.
    esc = scan(z -> (z - 0.985)^2, ((-1.0 - 1.0im, 1.0 + 1.0im),); zres=0.02)
    @test winding(only(link(esc))) == 2
    @test winding(only(link(esc; zoom=0))) == 0   # withheld without escalation
end

@testset "m=0: pure 2-D box finds and classifies point zeros/poles" begin
    # Two zeros (0.3, −0.4i) and one pole (0.6) in a plain complex box, no parameters.
    f(z) = (z - 0.3) * (z + 0.4im) / (z - 0.6)
    s = survey(f, ((-1.0 - 1.0im, 1.0 + 1.0im),); zres=0.01)
    @test length(s) == 3
    zs, ps = filter(iszerobranch, s), filter(ispolebranch, s)
    @test length(zs) == 2
    @test length(ps) == 1
    ctr(b) = sum(smp -> smp.z, b) / length(b)
    @test any(abs(ctr(b) - 0.3) < 0.02 for b in zs)
    @test any(abs(ctr(b) + 0.4im) < 0.02 for b in zs)
    @test abs(ctr(ps[1]) - 0.6) < 0.02
    @test all(b[1].p === () for b in s)                  # m=0 ⇒ empty params
end

@testset "4-D survey: zero surface z*(p₁,p₂) over a 2-parameter box" begin
    # One zero locus z*(p₁,p₂) = 0.5 + 0.3p₁ − 0.2p₂ threading a 4-D box.
    f(z, p1, p2) = z - (0.5 + 0.3p1 - 0.2p2)
    box = ((0.0 - 0.3im, 1.0 + 0.3im), (0.0, 1.0), (0.0, 1.0))
    s = survey(f, box; zres=0.03)
    @test length(s) == 1
    @test winding(s[1]) > 0
    pts = s[1]
    @test all(abs(smp.z - (0.5 + 0.3smp.p[1] - 0.2smp.p[2])) < 0.05 for smp in pts)
    @test any(smp.p[2] > 0.7 for smp in pts) && any(smp.p[2] < 0.3 for smp in pts)  # spans p₂
end

@testset "keep vetoes a region: unevaluated, unrefined, unreported" begin
    # Two zero curves; `keep` admits only the upper one. The lower is as
    # winding-worthy as the upper, so only the veto can suppress it — the case of an
    # `f` that is an approximation with a limited domain of validity.
    f(z, t) = (z - (0.5 + 0.2im + 0.1t)) * (z - (0.5 - 0.2im - 0.1t))
    box = ((0.0 - 0.4im, 1.0 + 0.4im), (0.0, 1.0))

    s = survey(f, box; zres=0.02)
    @test length(s) == 2

    # `f` is called from multiple threads, so the probe log needs a lock
    probed = Float64[]
    lk = ReentrantLock()
    g(z, t) = (Base.@lock lk push!(probed, imag(z)); f(z, t))
    sk = survey(g, box; zres=0.02, keep=(z, t) -> imag(z) > 0)
    @test length(sk) == 1
    @test all(imag(smp.z) > 0 for smp in sk[1])
    @test sk.nevals < s.nevals
    # cells wholly below the cut are never evaluated; those straddling it are kept,
    # so probes reach one base cell down but not the vetoed curve at Im z ≈ -0.2
    @test minimum(probed) > -0.2
end
