using MeromorphicLoci
using MeromorphicLoci: quadrant, nevals, scan, link
using Test

iszerobranch(b) = winding(b) > 0

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

@testset "winding" begin
    include("test_winding.jl")
end

@testset "AAA" begin
    include("test_aaa.jl")
end

@testset "muller" begin
    include("test_muller.jl")
end

@testset "point finders polish neighbouring loci apart, to machine precision" begin
    a = 0.3 + 0.2im
    f(z) = (z - a) * (z - a - 5e-4im) / (z + 0.5)
    box = (-1 - 1im, 1 + 1im)
    for alg in (AAA(), ArgumentPrinciple(zres=1e-4))
        r = findzeropoles(f, box, alg)
        @test r.windings == [-1, 1, 1] && isapprox(r.points[1], -0.5; atol = 1e-14)
        # equal real parts: their reim order is rounding
        @test sort(r.points[2:3]; by=imag) ≈ [a, a + 5e-4im] atol = 1e-14
        @test all(r.converged) && isempty(r.unresolved)
        raw = findzeropoles(f, box, alg; polish=nothing)
        @test raw.nevals < r.nevals && !any(raw.converged)
    end
end

@testset "README example survives its own asserts" begin
    f(p) = z -> (z - 0.80cispi(p)) * (z + 0.10 + 0.12im) /
                ((z - 0.804cispi(-p + 0.3)) * (z + 0.13 + 0.12im))
    zbox, pbox = (-0.98 - 0.98im, 0.98 + 0.98im), ((0.0, 1.0),)
    loci(; kw...) = findloci(f, zbox, pbox, ArgumentPrinciple(; zres=0.02, kw...))
    @test sort(winding.(loci())) == [-1, 1]        # graze escalated apart
    @test winding(only(loci(zoom=0))) == 0
    @test length(loci(minsep=0.3, zoom=0)) == 3
    @test sort(winding.(loci(minsep=0.3))) == [-1, -1, 1, 1]

end

@testset "each corner costs exactly one f call: neighbouring cells share corners and so share evaluations" begin
    calls = NTuple{3,Float64}[]
    lk = ReentrantLock()
    nsetup = Threads.Atomic{Int}(0)
    f(t) = (Threads.atomic_add!(nsetup, 1);
    z -> (Base.@lock lk push!(calls, (real(z), imag(z), t)); z - (0.5 + 0.3t)))
    box = ((0.0 - 0.3im, 1.0 + 0.3im), (0.0, 1.0))

    sc = scan(f, box; zres=0.02)
    @test length(calls) > 0
    @test allunique(calls)
    @test length(calls) == nevals(sc)
    # the per-p setup is shared by the z's evaluated together at that p
    @test nsetup[] < length(calls) ÷ 10
end

@testset "scan / link phases" begin
    f(t) = z -> z - (0.5 + 0.3t)
    box = ((0.0 - 0.3im, 1.0 + 0.3im), (0.0, 1.0))

    sc = scan(f, box; zres=0.02)

    @test length(sc) > 1 && !isempty(sc.cells)
    @test eltype(sc) === Sample{1}

    s = link(sc)
    # linking partitions the scan's cells — none invented, none dropped
    @test sum(length(b.cells) for b in s) == length(sc.cells)
    # and costs only the classification evals on top of the scan's refinement
    @test nevals(s) ≥ nevals(sc) == nevals(scan(f, box; zres=0.02))

    # findloci is exactly the two phases composed
    s2 = findloci(f, box[1], Base.tail(box), ArgumentPrinciple(zres=0.02))
    @test length(s2) == length(s) && nevals(s2) == nevals(s)
    @test [winding(b) for b in s2] == [winding(b) for b in s]
end

@testset "zres must describe a reachable lattice" begin
    f(z) = z - 0.3
    zeros_at(zres) = findzeros(f, (-1.0 - 1.0im, 1.0 + 1.0im), ArgumentPrinciple(; zres)).points
    @test_throws ArgumentError zeros_at(0.0)
    @test_throws ArgumentError zeros_at(-1.0)
    @test_throws ArgumentError zeros_at(1e-12)   # more levels than Int32 indices span
    @test length(zeros_at(10.0)) == 1
end

@testset "Survey vector interface" begin
    f(t) = z -> z - (0.5 + 0.3t)
    s = findloci(f, (0.0 - 0.3im, 1.0 + 0.3im), ((0.0, 1.0),), ArgumentPrinciple(zres=0.02))
    @test s isa AbstractVector{<:Branch}
    @test collect(s) == s.branches
    @test only(filter(iszerobranch, s)) === s[1]
    @test s[1] isa AbstractVector{<:Sample}
end

@testset "pole curve also fires (zero/pole agnostic)" begin
    # A moving simple pole: the winding criterion is blind to zero-vs-pole.
    f(t) = z -> 1 / (z - (0.5 + 0.2t))
    s = findloci(f, (0.0 - 0.3im, 1.0 + 0.3im), ((0.0, 1.0),), ArgumentPrinciple(zres=0.02))
    @test length(s) == 1
    @test winding(s[1]) < 0                    # opposite winding sign
    @test all(abs(smp.z - (0.5 + 0.2smp.p[1])) < 0.03 for smp in s[1])
end

@testset "two non-crossing branches stay separate" begin
    # Zeros at 0.3+0.5t and 1.5-0.5t: closest approach 0.4 (at t=1) ≫ zres, so
    # the components must not merge.
    r1(t) = 0.3 + 0.5t
    r2(t) = 1.5 - 0.5t
    f(t) = z -> (z - r1(t)) * (z - r2(t))
    s = findloci(f, (0.0 - 0.4im, 2.0 + 0.4im), ((0.0, 1.0),), ArgumentPrinciple(zres=0.02))
    @test length(s) == 2
    # each branch tracks exactly one of the two curves
    for b in s
        on1 = all(abs(smp.z - r1(smp.p[1])) < 0.04 for smp in b)
        on2 = all(abs(smp.z - r2(smp.p[1])) < 0.04 for smp in b)
        @test on1 ⊻ on2
    end
end

@testset "`minsep` sets the base mesh, hence the evaluation floor" begin
    f(t) = z -> z - (0.5 + 0.3t)
    ev(minsep) = nevals(findloci(f, (0.0 - 0.3im, 1.0 + 0.3im), ((0.0, 1.0),),
        ArgumentPrinciple(; zres=0.02, minsep)))
    # A uniform z mesh at every knot: halving `minsep` multiplies the floor by up to 4.
    @test ev(0.4) < ev(0.2) < ev(0.1)
end

@testset "escalation splits a fused zero/pole pair below zres" begin
    box = ((0.0 - 0.4im, 2.0 + 0.4im), (0.0, 1.0))
    pair(sep) = t -> z -> (z - (1.0 - sep / 2)) / (z - (1.0 + sep / 2))

    # Separated at cell scale (0.08 ≳ 2 cells): two branches, each voted on its own
    # outline, which stops short of the partner's cells.
    wide = link(scan(pair(0.08), box; zres=0.05))
    @test sort([sign(winding(b)) for b in wide]) == [-1, 1]

    # Fused (closer than a cell): the loci merge into one branch whose outline
    # encloses both partners and votes 0 — no reshuffling of the zres-level corners
    # can label it. `link` escalates the branch below zres locally until the pair
    # falls apart into clean sub-branches.
    sc = scan(pair(0.02), box; zres=0.05)
    tight = link(sc)
    @test sort([winding(b) for b in tight]) == [-1, 1]
    for b in tight   # escalated samples still sit on their locus
        z = winding(b) > 0 ? 0.99 : 1.01
        @test all(abs(smp.z - z) < 0.05 for smp in b)
    end

    # zoom=0 disables escalation: one merged branch, winding 0.
    @test winding(only(link(sc; zoom=0))) == 0
end

@testset "escalation stays local to the ambiguous branch (spending only a fraction of the scan's budget)" begin
    f(p) = z -> (z - 0.80cispi(p)) * (z + 0.11 + 0.12im) /
                ((z - 0.55cispi(-p + 0.3)) * (z + 0.13 + 0.12im))
    box = ((-0.98 - 0.98im, 0.98 + 0.98im), (0.0, 1.0))
    sc = scan(f, box; zres=0.02, minsep=0.12)
    n0 = nevals(sc)
    s = link(sc)
    @test sort([sign(winding(b)) for b in s]) == [-1, -1, 1, 1]
    @test nevals(s) - n0 < n0 ÷ 4
end

@testset "a vote whose outline would leave the box is withheld" begin
    # Within a column of the box edge the outline would cross it, so the vote is
    # withheld rather than evaluating f outside. Escalation narrows the columns until
    # it fits and recovers the label.
    esc = scan(Returns(z -> (z - 0.997)^2), ((-1.0 - 1.0im, 1.0 + 1.0im),); zres=0.02)
    @test winding(only(link(esc))) == 2
    @test winding(only(link(esc; zoom=0))) == 0   # withheld without escalation
end

@testset "m=0: distinct points with their windings" begin
    # Zeros 0.3 and −0.4i (double), pole 0.6.
    f(z) = (z - 0.3) * (z + 0.4im)^2 / (z - 0.6)
    box, alg = (-1.0 - 1.0im, 1.0 + 1.0im), ArgumentPrinciple(zres=0.01)
    r = findzeropoles(f, box, alg)
    @test all(abs.(r.points .- [-0.4im, 0.3, 0.6]) .< 0.01) && r.windings == [2, 1, -1]
    @test r.converged == [false, true, true] && isempty(r.unresolved)
    @test findzeros(f, box, alg).points == r.points[1:2] && findpoles(f, box, alg).points == r.points[3:3]
    @test findzeropoles(f, (1.0 - 1.0im, -1.0 + 1.0im), alg) == r
    # a fused pair without escalation is reported
    fused = findzeropoles(z -> (z - 0.5) / (z - 0.505), box, ArgumentPrinciple(zres=0.02, zoom=0))
    @test isempty(fused.points) && length(fused.unresolved) == 1

    # multiple zeros off the lattice and centred in a face, where a 4-corner face walk
    # aliases the multiplicity
    h = 2 / 2^8                  # zres-level column side for this box at zres = 0.02
    for (k, z0) in ((2, 0.1234 + 0.0567im), (3, complex(-1 + 100.5h, -1 + 131.5h)))
        g(z) = (z - z0)^k
        r = findzeropoles(g, box, ArgumentPrinciple(zres=0.02))
        @test r.windings == [k]
        @test findzeros(g, box, ArgumentPrinciple(zres=0.02)).points == r.points
    end
    # simple zeros a column apart: neither's vote counts the other
    @test length(findzeros(z -> (z - 0.3) * (z - 0.32), box, ArgumentPrinciple(zres=0.02)).points) == 2
    # simple zeros fused into one branch stay at its mean: a polish would put both on one
    r = findzeros(z -> (z - 0.3) * (z - 0.31), box, ArgumentPrinciple(zres=0.02))
    @test r.windings == [2] && r.converged == [false] && abs(only(r.points) - 0.3) < 0.02 && abs(only(r.points) - 0.31) < 0.02
end

@testset "4-D: zero surface z*(p₁,p₂) over a 2-parameter box" begin
    # One zero locus z*(p₁,p₂) = 0.5 + 0.3p₁ − 0.2p₂ threading a 4-D box.
    f(p1, p2) = z -> z - (0.5 + 0.3p1 - 0.2p2)
    s = findloci(f, (0.0 - 0.3im, 1.0 + 0.3im), ((0.0, 1.0), (0.0, 1.0)), ArgumentPrinciple(zres=0.03))
    @test length(s) == 1
    @test winding(s[1]) > 0
    pts = s[1]
    @test all(abs(smp.z - (0.5 + 0.3smp.p[1] - 0.2smp.p[2])) < 0.05 for smp in pts)
    @test any(smp.p[2] > 0.7 for smp in pts) && any(smp.p[2] < 0.3 for smp in pts)  # spans p₂
end

@testset "parameter resolution follows the loci, not zres" begin
    zbox, alg = (0.0 - 0.4im, 2.0 + 0.4im), ArgumentPrinciple(zres=0.02)
    # A locus slower than a cell per knot step is sampled exactly at the knots, uneven or not.
    knots = [0.0, 0.1, 0.5, 1.0]
    b = only(findloci(t -> z -> z - (1.0 + 0.01t), zbox, (knots,), alg))
    @test sort(unique(s.p[1] for s in b)) == knots

    # A fast one (1.5 per unit p, ~100 cells) with only the end knots still comes out as
    # one branch, sampled at least every cell along z.
    s = findloci(t -> z -> z - (0.2 + 1.5t), zbox, ([0.0, 1.0],), alg)
    b = only(s)
    @test winding(b) == 1
    zs = sort([real(x.z) for x in b])
    @test maximum(diff(zs)) < 0.02 && zs[1] < 0.22 && zs[end] > 1.68

    @test_throws ArgumentError findloci(t -> z -> z, zbox, ([0.0, 0.0, 1.0],), alg)
    @test_throws ArgumentError findloci(t -> z -> z, zbox, ([0.5],), alg)
end
