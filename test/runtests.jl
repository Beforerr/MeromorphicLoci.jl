using MeromorphicLoci
using MeromorphicLoci: quadrant, nevals
using Test

tmap(g, xs) = fetch.([Threads.@spawn g(x) for x in xs])

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

@testset "continuation" begin
    include("test_continuation.jl")
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

@testset "nevals counts every call of f" begin
    for alg in (AAA(), ArgumentPrinciple(zres=0.02)), polish in (Muller(), nothing)
        n = Threads.Atomic{Int}(0)
        f(z) = (Threads.atomic_add!(n, 1); (z - 0.3) * (z + 0.4im)^2 / (z - 0.6))
        r = findzeropoles(f, (-1 - 1im, 1 + 1im), alg; polish)
        @test r.nevals == n[]
        @test findzeropoles(f, (-1 - 1im, 1 + 1im), alg; polish, mapper=tmap) == r
    end
    # the search, the certificate and the votes share their points: none is evaluated twice
    seen = ComplexF64[]
    findzeropoles(z -> (push!(seen, z); (z - 0.3) * (z + 0.4im)^2 / (z - 0.6)), (-1 - 1im, 1 + 1im),
        ArgumentPrinciple(zres=0.02); polish=nothing)
    @test allunique(seen)
end

@testset "zres must describe a reachable lattice" begin
    f(z) = z - 0.3
    zeros_at(zres) = findzeros(f, (-1.0 - 1.0im, 1.0 + 1.0im), ArgumentPrinciple(; zres)).points
    @test_throws ArgumentError zeros_at(0.0)
    @test_throws ArgumentError zeros_at(-1.0)
    @test_throws ArgumentError zeros_at(1e-12)   # more levels than Int32 indices span
    @test_throws ArgumentError findzeros(f, (-1.0 + 0.0im, 1.0 + 0.0im), ArgumentPrinciple(zres=0.1))
    @test length(zeros_at(10.0)) == 1
end

@testset "a locus beside another is not misread away" begin
    # Each pair's partner turns an edge of the cell holding the other past ±π, so the
    # corner walk alone reads that cell empty.
    box, alg = (-1 - 1im, 1 + 1im), ArgumentPrinciple(zres=0.02)
    for (sg, a, b) in ((1, 0.182 + 0.075im, 0.367 + 0.069im), (-1, -0.032 - 0.174im, 0.242 - 0.16im),
                       (1, 0.118 - 0.07im, 0.371 - 0.022im))
        r = findzeropoles(z -> ((z - a) * (z - b))^sg, box, alg; polish=nothing)
        @test r.windings == [sg, sg] && all(z -> minimum(abs.(r.points .- z)) < 0.02, (a, b))
    end
end

@testset "loci a coarse cell aliases away are recovered where the box's count exceeds the found" begin
    # eight zeros on a circle 0.1 across, inside one root cell: its corners turn 8 times
    for c in (0.3 + 0.2im, 0.123 - 0.456im)
        r = findzeros(z -> (z - c)^8 - 0.1^8, (-1 - 1im, 1 + 1im), ArgumentPrinciple(zres=0.02))
        @test r.windings == fill(1, 8) && r.deficit == 0 && r.nevals < 3500
    end
    # a zero on the box edge leaves the box's count unreadable: unchecked, not passed
    g8(z) = ((z - 0.3 - 0.2im)^8 - 0.1^8) * (z - 1 - 0.1234im)
    @test findzeros(g8, (-1 - 1im, 1 + 1im), ArgumentPrinciple(zres=0.02)).deficit === nothing
    # a 4-fold zero on a lattice point: only cells beside it fire, so the cells' count
    # passes but labelling drops them; the deficit counts the labels, not the cells
    @test findzeros(z -> (z - (1 + 1im) / 128)^4, (-1 - 1im, 1 + 1im), ArgumentPrinciple(zres=0.02)).deficit == 4
end

@testset "`minsep` sets the base mesh, hence the evaluation floor" begin
    ev(minsep) = findzeros(z -> z - 0.65, (0.0 - 0.3im, 1.0 + 0.3im), ArgumentPrinciple(; zres=0.02, minsep)).nevals
    # A uniform mesh: halving `minsep` multiplies the floor by up to 4.
    @test ev(0.4) < ev(0.2) < ev(0.1)
    # below `zres` it still sets the mesh: a zero and a pole 0.015 apart are told apart
    r = findzeropoles(z -> (z - 0.5) / (z - 0.515), (0.0 - 0.3im, 1.0 + 0.3im), ArgumentPrinciple(zres=0.05, minsep=0.01, zoom=0))
    @test r.windings == [1, -1]
end

@testset "escalation splits a zero–pole pair fused below zres" begin
    pair(z) = (z - 0.99) / (z - 1.01)
    box = (0.0 - 0.4im, 2.0 + 0.4im)
    r = findzeropoles(pair, box, ArgumentPrinciple(zres=0.05); polish=nothing)
    @test r.windings == [1, -1] && all(abs.(r.points .- [0.99, 1.01]) .< 0.05)
    # zoom = 0 disables escalation: the pair is reported, its kind unknown
    fused = findzeropoles(pair, box, ArgumentPrinciple(zres=0.05, zoom=0))
    @test isempty(fused.points) && length(fused.unresolved) == 1
end

@testset "a locus within a cell of the box edge is labelled, its outline clipped to the box" begin
    box = (-1 + 0im, 1 + 1im)
    for (d, zoom) in ((1e-6, 4), (1e-3, 4), (1e-3, 0))
        outside = Threads.Atomic{Int}(0)
        f(z) = (-1 <= real(z) <= 1 && 0 <= imag(z) <= 1 || Threads.atomic_add!(outside, 1); z - (0.3 + d * im))
        r = findzeros(f, box, ArgumentPrinciple(; zres=0.02, zoom); polish=nothing)
        @test r.windings == [1] && isempty(r.unresolved) && abs(only(r.points) - 0.3) < 0.02
        @test outside[] == 0
    end
end

@testset "distinct points with their windings" begin
    # Zeros 0.3 and −0.4i (double), pole 0.6.
    f(z) = (z - 0.3) * (z + 0.4im)^2 / (z - 0.6)
    box, alg = (-1.0 - 1.0im, 1.0 + 1.0im), ArgumentPrinciple(zres=0.01)
    r = findzeropoles(f, box, alg)
    @test all(abs.(r.points .- [-0.4im, 0.3, 0.6]) .< 0.01) && r.windings == [2, 1, -1]
    @test r.converged == [false, true, true] && isempty(r.unresolved)
    @test findzeros(f, box, alg).points == r.points[1:2] && findpoles(f, box, alg).points == r.points[3:3]
    @test findzeropoles(f, (1.0 - 1.0im, -1.0 + 1.0im), alg) == r
    # multiple zeros off the lattice and centred in a cell, where a 4-corner walk
    # aliases the multiplicity
    h = 2 / 2^8                  # zres-level cell side for this box at zres = 0.02
    for (k, z0) in ((2, 0.1234 + 0.0567im), (3, complex(-1 + 100.5h, -1 + 131.5h)))
        g(z) = (z - z0)^k
        r = findzeropoles(g, box, ArgumentPrinciple(zres=0.02))
        @test r.windings == [k]
        @test findzeros(g, box, ArgumentPrinciple(zres=0.02)).points == r.points
    end
    # cells of a 5:1 box stay square: a double zero is discovered, not walked past
    @test findzeropoles(z -> (z + 0.21 - 0.13im)^2, (-1.0 - 0.2im, 1.0 + 0.2im),
        ArgumentPrinciple(zres=0.02)).windings == [2]
    # and so do those of a strip thinner than zres
    @test findzeropoles(z -> (z - 3.3 - 0.05im)^2 * (z - 7.1 - 0.03im), (0.0 + 0.0im, 10.0 + 0.1im),
        ArgumentPrinciple(zres=0.2)).windings == [2, 1]
    # simple zeros a cell apart: neither's vote counts the other
    @test length(findzeros(z -> (z - 0.3) * (z - 0.32), box, ArgumentPrinciple(zres=0.02)).points) == 2
    # simple zeros fused into one branch stay at its mean: a polish would put both on one
    r = findzeros(z -> (z - 0.3) * (z - 0.31), box, ArgumentPrinciple(zres=0.02))
    @test r.windings == [2] && r.converged == [false] && abs(only(r.points) - 0.3) < 0.02 && abs(only(r.points) - 0.31) < 0.02
end
