using MeromorphicLoci: aaa, poles

# Multiplicity found within `tol` of `z`, however the finder split it.
mult_near(r, z, tol) = sum(abs, r.windings[abs.(r.points .- z) .< tol]; init = 0)

@testset "aaa recovers a rational's poles at its exact degree" begin
    ps = [0.21 + 0.1im, -0.5im, -0.7 + 0.05im, 0.9 - 0.15im]
    f(x) = 3 + sum(c / (x - p) for (c, p) in zip((1, 2im, -0.5, 1.0e-3), ps))
    # 51 samples: an odd count; 400: the fit must stop far short of its sample budget
    for z in (complex.(range(-1, 1; length = 51)), [complex(x, y) for x in range(-1, 1; length = 20) for y in range(-0.2, 0.2; length = 20)])
        (; fit, converged) = aaa(f.(z), z)
        @test converged
        @test length(fit.nodes) == length(ps) + 1
        @test fit.(fit.nodes) == f.(fit.nodes)
        @test fit(0.3 + 0.05im) ≈ f(0.3 + 0.05im) rtol = 1.0e-10
        found = poles(fit)
        @test all(p -> minimum(abs, found .- p) < 1.0e-9, ps)
    end
end

@testset "findzeros skips non-finite samples" begin
    zs = [1.0, -1.0, 0.5im]
    f(z) = prod(z .- zs) / (z - 2)
    # the 13×5 grid holds z = 2 exactly, where f is infinite
    Z = range(-3, 3, 13) .+ im .* range(-1, 1, 5)'
    r = findzeros(f, Z, AAA())
    @test r.fitconverged
    @test all(z -> minimum(abs, r.points .- z) < 1.0e-8, zs)
    r = findzeros(z -> Inf, Z, AAA())
    @test isempty(r.points) && !r.fitconverged
end

@testset "findzeros on a filtered grid and on boundary points" begin
    zs = [1.0, -1.0, 0.5im]
    f(z) = prod(z .- zs) / (z - 2)
    for Z in (filter(z -> abs(z) < 1.5, range(-1.5, 1.5, 30) .+ im .* range(-1.5, 1.5, 30)'), 1.5 .* cispi.(2 .* (0:99) ./ 100))
        r = findzeros(f, Z, AAA())
        @test all(z -> minimum(abs, r.points .- z) < 1.0e-8, zs)
        # sample points bound no region: the pole beyond them is reported
        @test only(findpoles(f, Z, AAA()).points) ≈ 2 atol = 1.0e-8
    end
end

@testset "a box zdomain reports the loci inside it" begin
    r = findzeropoles(z -> (z - 0.5) * (z + 2) / (z - 1im), (-1 - 0.5im, 1 + 0.5im), AAA())
    @test only(r.points) ≈ 0.5 && r.windings == [1]
end

@testset "findpoles and findzeropoles" begin
    zs, ps = [1.0, -1.0, 0.5im], [2.0, -0.3im]
    f(z) = prod(z .- zs) / prod(z .- ps)
    Z = range(-3, 3, 20) .+ im .* range(-1, 1, 16)'
    near(found, exact) = all(x -> minimum(abs, found .- x) < 1.0e-8, exact)
    l = findzeropoles(f, Z, AAA())
    @test l.fitconverged && near(l.points[l.windings .> 0], zs) && near(l.points[l.windings .< 0], ps)
    @test l.points[l.windings .> 0] == findzeros(f, Z, AAA()).points
    @test near(findpoles(f, Z, AAA()).points, ps)
end

@testset "real samples: real zeros, off-axis poles, one result type" begin
    f(x) = (x^2 - 0.25) / (x^2 + 0.25)
    l = findzeropoles(f, range(-1, 1, 40), AAA())
    near(found, exact) = length(found) == length(exact) && all(x -> minimum(abs, found .- x) < 1.0e-8, exact)
    @test l.fitconverged && near(l.points[l.windings .> 0], [0.5, -0.5]) && near(l.points[l.windings .< 0], [0.5im, -0.5im])
    @inferred findzeros(f, range(-1, 1, 40), AAA())
end

@testset "noisy multiple zeros keep their multiplicity" begin
    # Noise splits each into candidates whose own circles miss it; relative and absolute noise.
    f(z) = (z - 0.4 - 0.2im)^2 * (z + 1.1)^3
    n(z) = cis(1.0e15 * real(z) + 3.0e14 * imag(z))
    for g in (z -> f(z) * (1 + 1.0e-8 * n(z)), z -> f(z) + 1.0e-7 * n(z))
        r = findzeros(g, (-3 - 1im, 3 + 1im), AAA())
        @test mult_near(r, 0.4 + 0.2im, 1.0e-3) == 2
        @test mult_near(r, -1.1, 1.0e-3) == 3
        @test sum(r.windings) == 5 && isempty(r.unresolved)
    end
    # a double zero closer to the box edge than its candidates' split
    a = 3 - 1.0e-9 + 0.25im
    @test mult_near(findzeros(z -> (z - a)^2 * (z + 1), (-3 - 1im, 3 + 1im), AAA()), a, 1.0e-6) == 2
end

@testset "a high-order zero keeps its order" begin
    # a symmetry-forced z^k factor: the circle is centred on the zero, where |f| is constant
    for k in (6, 7)
        @test findzeros(z -> z^k * exp(z), (-3 - 1im, 3 + 1im), AAA()).windings == [k]
    end
end
