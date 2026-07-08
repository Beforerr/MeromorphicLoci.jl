using MeromorphicLoci: _circle_winding

@testset "adaptive circle winding" begin
    # range tracks the seed=8 default: ±12 sits on the π/2 acceptance boundary
    # and 13 ≡ -3 (mod 2seed) aliases undetectably; raising `seed` is the only cure.
    for m in -11:11
        @test _circle_winding(z -> z^m, 0.0 + 0.0im, 0.1)[1] == m
    end
    @test_broken _circle_winding(z -> z^13, 0.0 + 0.0im, 0.1)[1] == 13
    @test _circle_winding(z -> z^13, 0.0 + 0.0im, 0.1; seed=16)[1] == 13

    @test _circle_winding(z -> (z - 1e-3) * (z + 1e-3), 0.0im, 0.1)[1] == 2
    # Unusable samples yield `nothing` instead of an InexactError crash.
    @test _circle_winding(z -> NaN, 0.0im, 0.1)[1] === nothing
end