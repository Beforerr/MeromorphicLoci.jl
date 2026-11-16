using MeromorphicLoci: nevals, Branch, Sample

winding(b) = b.winding

box = (-1 - 1im, 1 + 1im)
prange(b) = extrema(x -> x.p[1], b)
# each branch's samples lie on one of the analytic `loci(p)`, a different one per branch
function onloci(s, loci; atol=1e-8)
    ids = [findfirst(k -> all(x -> abs(x.z - loci(x.p[1])[k]) < atol, b), eachindex(loci(0.0))) for b in s]
    return !any(isnothing, ids) && allunique(ids)
end

@testset "a zero grazing a pole 0.004 apart: one branch each, on the loci" begin
    f(p) = z -> (z - 0.8cispi(p)) / (z - 0.804cispi(-p + 0.3))
    loci(p) = (0.8cispi(p), 0.804cispi(-p + 0.3))
    for (alg, cap) in ((AAA(), 4500), (Continuation(ArgumentPrinciple(zres=0.02, minsep=0.18)), 5500))
        s = findloci(f, box, ((0.0, 1.0),), alg)
        @test winding.(s) == [-1, 1] && onloci(s, loci)
        @test all(b -> prange(b) == (0.0, 1.0) && isempty(b.uncertified), s) && isempty(s.deficit)
        @test nevals(s) < cap
    end
    # knots only discover: the two ends find both loci; a threaded mapper changes nothing
    s = findloci(f, box, ([0.0, 1.0],), AAA())
    @test winding.(s) == [-1, 1] && onloci(s, loci) && nevals(s) < 1550
    t = findloci(f, box, ([0.0, 1.0],), AAA(); mapper=tmap)
    @test nevals(t) == nevals(s) && [b.samples for b in t] == [b.samples for b in s]
    @test_throws ArgumentError findloci(f, box, (), AAA())
end

@testset "a zero passing through a pole: both continuous, uncertified where they meet" begin
    g(p) = z -> (z - (-0.5 + p) - 0.05im) / (z - 0.1 - 0.05im)
    s = findloci(g, box, ((0.0, 1.0),), AAA())
    @test sort(winding.(s)) == [-1, 1] && onloci(s, p -> (-0.5 + p + 0.05im, 0.1 + 0.05im))
    @test all(s) do b
        (p1,), (p2,) = only(b.uncertified)
        prange(b) == (0.0, 1.0) && p1 < 0.6 < p2 && p2 - p1 < 1e-3
    end
    @test occursin("2 branches uncertified", summary(s)) && nevals(s) < 4350
end

@testset "zeros of one kind: grazing ones keep their identity, a branch point ends them" begin
    s = findloci(p -> z -> (z - 0.8cispi(p)) * (z - 0.804cispi(-p + 0.3)), box, ((0.0, 1.0),), AAA())
    @test winding.(s) == [1, 1] && onloci(s, p -> (0.8cispi(p), 0.804cispi(-p + 0.3)))
    @test all(b -> prange(b) == (0.0, 1.0) && isempty(b.uncertified), s) && nevals(s) < 4600
    # z = 0.1i ± √(0.5(p − 0.45)): which zero continues which past p = 0.45 has no answer
    s = findloci(p -> z -> (z - 0.1im)^2 - 0.5(p - 0.45), box, ((0.0, 1.0),), AAA())
    @test length(s) == 4 && all(==(1), winding.(s))
    @test all(s) do b
        lo, hi = prange(b)
        q = lo == 0 ? hi : lo
        (lo == 0 || hi == 1) && b.uncertified == [((q,), (q,))] && abs(q - 0.45) < 1e-6
    end
end

@testset "an avoided crossing: each zero keeps to its own locus, flagged by nothing" begin
    # z = ±√(a² + g²/4), a = 0.8(p − p0) + i d0/2. The lines ±a pass d0 apart; coupled by
    # g > d0, the zeros exchange them, each keeping the sign of Re z, never closer than
    # √(g² − d0²): 10 floors, then wider than `zres`. A step over the crossing lands
    # each track on the other's locus with no error in either.
    cross(p0, d0, g) = p -> z -> z^2 - (0.8(p - p0) + im * d0 / 2)^2 - g^2 / 4
    for (f, knots) in ((cross(0.34, 1e-3, 1.5e-3), ((0.0, 1.0),)), (cross(0.565, 0.03, 0.05), ([0.0, 1.0],)))
        s = findloci(f, box, knots, AAA())
        @test length(s) == 2 && isempty(s.deficit)
        @test all(b -> prange(b) == (0.0, 1.0) && isempty(b.uncertified) && allequal(x -> real(x.z) > 0, b), s)
    end
    # The same exchange seen moving with one line: a zero at rest until the other arrives,
    # stops in its place and sends it on. The other is in the box at neither knot, so no
    # seed tells of it, and the root found at rest after the crossing is not the one before.
    d0, g = 1e-3, 1.5e-3
    a(p) = 0.8(p - 0.5) + im * d0 / 2
    s = findloci(p -> z -> (z - 0.8(p - 0.5))^2 - a(p)^2 - g^2 / 4, 0.3 .* box, ([0.0, 1.0],), AAA())
    @test onloci(s, p -> 0.8(p - 0.5) .+ (1, -1) .* sqrt(a(p)^2 + g^2 / 4)) && all(b -> isempty(b.uncertified), s)
    @test [round.(prange(b); digits=2) for b in s] == [(0.0, 0.69), (0.31, 1.0)]
    # over two parameters an exchange along some grid lines would join the two sheets
    s = findloci((p, q) -> cross(0.37 + 0.2q, 1e-3, 2e-3)(p), box, ([0.0, 0.5, 1.0], [0.0, 0.5, 1.0]), AAA())
    @test [b.degree for b in s] == [1, 1] && all(b -> allequal(x -> real(x.z) > 0, b), s)
end

@testset "multiple loci: a double zero is one winding-2 branch, a fused pair two branches" begin
    s = findloci(p -> z -> (z - 0.5cispi(p))^2 / (z + 0.3), box, ((0.0, 1.0),), AAA())
    @test winding.(s) == [-1, 2] && onloci(s, p -> (-0.3, 0.5cispi(p)))
    @test all(b -> prange(b) == (0.0, 1.0) && isempty(b.uncertified), s) && nevals(s) < 4850
    s = findloci(p -> z -> (z + 0.3) / (z - 0.5cispi(p))^2, box, ((0.0, 1.0),), AAA())
    @test winding.(s) == [1, -2] && all(b -> prange(b) == (0.0, 1.0) && isempty(b.uncertified), s)
    # whatever the scale of the box
    s = findloci(p -> z -> (1e6z - 0.5cispi(p))^2 / (1e6z + 0.3), 1e-6 .* box, ((0.0, 1.0),), AAA())
    @test winding.(s) == [-1, 2] && all(b -> prange(b) == (0.0, 1.0) && isempty(b.uncertified), s)
    # the seeder reads the pair, 0.01 apart, as one double zero at both knots
    pair(p) = z -> (z - 0.3 - 0.2p) * (z - 0.31 - 0.2p)
    seed = ArgumentPrinciple(zres=0.02, minsep=0.18)
    @test findzeros(pair(0.0), box, seed).windings == findzeros(pair(0.75), box, seed).windings == [2]
    s = findloci(pair, box, ([0.0, 0.75],), Continuation(seed))
    @test winding.(s) == [1, 1] && onloci(s, p -> (0.3 + 0.2p, 0.31 + 0.2p))
    # closer than the floor `zfloor`, a pair is one double locus
    s = findloci(p -> z -> (z - 0.5cispi(p)) * (z - 0.5cispi(p) - 2e-5), box, ((0.0, 1.0),), AAA())
    @test winding.(s) == [2] && isempty(only(s).uncertified)
end

@testset "a knot where two loci coincide: they end on the double locus found there" begin
    # ω = ±k from k = 0
    s = findloci(k -> z -> (z^2 - k^2) * (z - 1.5 - 0.5im), (-2 - 1im, 2 + 1im), ((0.0, 1.0),), AAA())
    @test sort(winding.(s)) == [1, 1, 1, 2] && isempty(s.deficit) && nevals(s) < 4700
    pair = filter(b -> abs(last(b).z) ≈ 1, s)
    @test length(pair) == 2 && all(b -> prange(b) == (0.0, 1.0) && all(x -> abs(abs(x.z) - x.p[1]) < 1e-6, b), pair)
    @test all(b -> all(u -> u[2][1] < 1e-3, b.uncertified), pair)
    double = only(filter(b -> winding(b) == 2, s))
    @test length(double) == 1 && double.uncertified == [((0.0,), (0.0,))]
end

@testset "loci entering between knots are found; f is evaluated outside the box only by the seeder" begin
    zbox = (-2 - 12im, 1 + 12im)
    outside = Threads.Atomic{Int}(0)
    inside(z) = -2 <= real(z) <= 1 && -12 <= imag(z) <= 12
    f(τ) = z -> (inside(z) || Threads.atomic_add!(outside, 1); z + 1 + 2exp(-z * τ))
    s = findloci(f, zbox, ((0.5, 3.0),), AAA())
    n = outside[]
    @test length(s) == 12 && all(==(1), winding.(s)) && all(b -> isempty(b.uncertified), s)
    @test count(b -> prange(b) == (0.5, 3.0), s) == 2 && all(b -> prange(b)[2] == 3.0, s)
    @test all(b -> all(x -> abs(f(x.p[1])(x.z)) < 1e-9, b), s) && nevals(s) < 6900
    foreach(τ -> findzeropoles(f(τ), zbox, AAA()), range(0.5, 3.0, 8))   # the default knots
    @test outside[] == 2n > 0
end

@testset "a locus the seeder misses at a knot is carried there by its track" begin
    # the pair is closer than the seeder's `minsep` at p = q only, where its windings cancel
    seed = ArgumentPrinciple(zres=0.02, minsep=0.18)
    z0 = 0.3 + 0.2im
    ks = [0.0, 0.5, 1.0]
    for q in (0.0, 0.5)    # carried from the middle knot into both ends, and through the middle one
        f(p) = z -> (z - z0) / (z - z0 - 0.03 - 0.4sinpi(p - q)^2)
        @test [length(findzeropoles(f(k), box, seed).points) for k in ks] == 2 .* (ks .!= q .&& ks .!= q + 1)
        s = findloci(f, box, (ks,), seed)
        @test winding.(s) == [1, -1] && all(b -> isempty(b.uncertified) && issubset(ks, [x.p[1] for x in b]), s)
    end
    # missed along the whole grid line q = 0.5: the lines across carry it into its nodes,
    # and the line is followed again from them
    h(p, q) = z -> (z - z0) * (z + 0.5 - 0.3p) / (z - z0 - 0.03 - 1.6(q - 0.5)^2)
    @test all(p -> length(findzeropoles(h(p, 0.5), box, seed).points) == 1, ks)
    s = findloci(h, box, (ks, ks), seed)
    @test sort(winding.(s)) == [-1, 1, 1] && all(b -> b.degree == 1 && isempty(b.uncertified), s)
    @test all(b -> any(x -> x.p[2] == 0.5 && !(x.p[1] in ks), b), s)
end

@testset "a seed that cannot be followed is a branch of one sample" begin
    # zoom = 0 leaves the pair, 0.005 apart, unresolved at each knot
    s = findloci(p -> z -> (z - 0.5) / (z - 0.505), box, ([0.0, 1.0],),
        Continuation(ArgumentPrinciple(zres=0.02, minsep=0.18, zoom=0)))
    @test length(s) == 2 && all(b -> winding(b) == 0 && length(b) == 1, s)
    @test [b.uncertified for b in s] == [[((0.0,), (0.0,))], [((1.0,), (1.0,))]]
end

@testset "the count at each knot flags a seeder that came up short" begin
    cs = -3.0:0.15:3.0
    f(t) = z -> prod(c -> z - (c + 2t) * cis(0.05), cs)
    zbox = (-1 - 0.2im, 1 + 0.2im)
    inbox(p) = count(c -> abs((c + 2p) * cos(0.05)) <= 1, cs)
    held(s, p) = count(b -> any(x -> x.p[1] == p, b), s)
    s = findloci(f, zbox, ([0.0, 1.0],), AAA(nsamples=30))    # too few samples for the 14 zeros at p = 1
    @test s.deficit == [(1.0,) => inbox(1.0) - held(s, 1.0)] && last(only(s.deficit)) > 0
    @test occursin("count off at 1 knot", summary(s))
    t = findloci(f, zbox, ([0.0, 1.0],), Continuation(AAA(nsamples=30); certify=false))
    @test isempty(t.deficit) && length(t) == length(s) && nevals(t) < nevals(s)
    @test isempty(findloci(f, zbox, ([0.0, 1.0],), AAA()).deficit)
end

@testset "a zero closing on a pole is followed while they are apart at all" begin
    # Whistler dispersion: the zero closes on the cyclotron pole as 0.5/k², to 0.008 at k = 8
    w0 = 0.5 - 0.01im
    s = findloci(k -> ω -> k^2 - ω^2 + ω / (ω - w0), (-1.0 - 0.5im, 3.0 + 0.5im), ((0.5, 8.0),), AAA())
    pair = filter(b -> abs(last(b).z - w0) < 0.05, s)
    @test sort(winding.(pair)) == [-1, 1] && all(b -> prange(b) == (0.5, 8.0) && isempty(b.uncertified), pair)
    @test nevals(s) < 4900
end

@testset "a train of loci through the box is one branch per locus" begin
    cs = -3.0:0.15:3.0
    s = findloci(t -> z -> prod(c -> z - (c + 2t) * cis(0.05), cs), (-1 - 0.2im, 1 + 0.2im), ((0.0, 1.0),), AAA())
    @test length(s) == 27 && all(==(1), winding.(s)) && nevals(s) < 9100
end

@testset "knots: any sorted real vector, each one sampled" begin
    zbox = (0.0 - 0.4im, 2.0 + 0.4im)
    for knots in ([0.0, 0.1, 0.5, 1.0], range(0, 1, 3))
        b = only(findloci(t -> z -> z - (1.0 + 0.01t), zbox, (knots,), AAA()))
        @test issubset(knots, [x.p[1] for x in b]) && b isa AbstractVector{<:Sample}
    end
    @test_throws ArgumentError findloci(t -> z -> z, zbox, ([0.0, 0.0, 1.0],), AAA())
    @test_throws ArgumentError findloci(t -> z -> z, zbox, ([0.5],), AAA())
    @test_throws ArgumentError findloci(t -> z -> z, zbox, ([0.0, 1.0],), Continuation(AAA(); zres=0.01, zfloor=0.02))
    # a point finder stands for Continuation seeded by it
    s = findloci(t -> z -> z - (1.0 + 0.01t), zbox, ([0.0, 1.0],), ArgumentPrinciple(zres=0.02, minsep=0.18))
    @test s isa AbstractVector{<:Branch} && winding(only(s)) == 1
end

@testset "two parameters: the loci are followed along the lines of the knot grid" begin
    # the pole surface passes 0.004–0.05 from the zero surface near p = 0.15
    h(p, q) = z -> (z - 0.8cispi(p)) / (z - (0.804 + 0.05q) * cispi(-p + 0.3))
    s = findloci(h, box, ([0.0, 0.5, 1.0], [0.0, 0.5, 1.0]), AAA())
    @test winding.(s) == [-1, 1] && all(b -> b.degree == 1 && isempty(b.uncertified), s) && isempty(s.deficit)
    @test all(x -> abs(x.z - (0.804 + 0.05x.p[2]) * cispi(-x.p[1] + 0.3)) < 1e-8, s[1])
    @test all(x -> abs(x.z - 0.8cispi(x.p[1])) < 1e-8, s[2]) && nevals(s) < 5900
    @test all(b -> extrema(x -> x.p[1], b) == extrema(x -> x.p[2], b) == (0.0, 1.0), s)
    # two zeros that swap around a point of the parameter plane are one branch of degree 2
    ep(p, q) = z -> z^2 - (p + im * q - (0.43 + 0.37im))
    b = only(findloci(ep, box, ((0.0, 1.0), (0.0, 1.0)), AAA()))
    @test b.degree == 2 && winding(b) == 1 && all(x -> abs(x.z^2 - (x.p[1] + im * x.p[2] - (0.43 + 0.37im))) < 1e-8, b)
    # and two branches over a range that does not hold the point
    @test [b.degree for b in findloci(ep, box, ((0.0, 0.4), (0.0, 0.3)), AAA())] == [1, 1]
end
