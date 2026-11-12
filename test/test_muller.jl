using MeromorphicLoci: muller

@testset "simple zero" begin
    r = muller(z -> (z - (1 + 2im)) * (z + 3) * exp(z), 0.8 + 1.7im; h = 1.0e-3)
    @test r.converged && r.root ≈ 1 + 2im && r.nevals < 20
end

@testset "seed exactly on a double zero, where the Muller step is 0/0" begin
    @test muller(z -> (z - 0.25)^2, 0.25; h = 0.25) == (; root = 0.25 + 0im, converged = true, nevals = 3)
end

@testset "non-finite f" begin
    @test isnan(muller(x -> ComplexF64(Inf, Inf), 1.0, 1.1, 1.2).root)
    @test isnan(muller(x -> ComplexF64(NaN, NaN), 1.0, 1.1, 1.2).root)

    hits = Ref(0)
    fovf(x) = (v = exp(10x) - 1; isfinite(v) || (hits[] += 1); v)
    @test abs(muller(fovf, -2.0, -1.9, -1.8).root) < 1.0e-10 && hits[] > 0
end

@testset "decades of |f| across the seeds collapse the step, but are no zero" begin
    slope(x) = exp(1.0e7 * x^2)
    @test !any(muller(slope, 0.002 * cis(θ); h = 2.0e-6).converged for θ in 0:0.2:6.2)
end

@testset "absolute noise floor" begin
    noise(x) = 1.0e-7 * cis(1.0e15 * real(x) + 3.0e14 * imag(x))
    for x0 in (0.0, 0.2im, -0.3)
        @test !muller(z -> exp(-5z) + noise(z), x0; h = 1.0e-3).converged
    end
    # Stops early at the floor; a `tol` within the noise converges.
    r0 = 0.3 - 0.1im
    noisy(x) = (x - r0) * (1 + x) + noise(x)
    for seed in (r0 + 0.01, r0 - 0.02im, r0 - 0.005)
        r = muller(noisy, seed; h = 1.0e-3)
        @test r.nevals < 40 && abs(r.root - r0) < 1.0e-6
        r = muller(noisy, seed; h = 1.0e-3, tol = 1.0e-5)
        @test r.converged && abs(r.root - r0) < 1.0e-5
    end
end

@testset "estimate near the rim of its disk" begin
    grow(z) = exp(6z) * (z^3 - 1) + 0.1
    ζ = -0.3936314947062384
    for θ in (0.0, 3.0, 6.0)
        r = muller(grow, ζ + 0.32 * cis(θ); h = 0.09, radius = 0.36)
        @test r.converged && abs(r.root - ζ) < 1.0e-12
    end
end
