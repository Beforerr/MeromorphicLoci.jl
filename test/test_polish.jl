# `_reach` sets r ≤ 0.4 × the distance to every other estimate; a locus nearer went unfound.
@testset "a multiple point beside an unfound locus is polished to its mean, or unconverged" begin
    z0, r = 0.3 + 0.1im, 0.1
    for R in (2.5, 1.6, 1.3, 1.15), off in (0.02, 0.3), p in (1, -1)
        f = z -> (z - z0)^2 * (z - z0 - R * r)^p
        z, ok, _ = MeromorphicLoci._polish(Polish(), f, z0 + off * r, r, 2, map)
        @test ok ? abs(z - z0) < 1e-14 : R < 1.6 && z == z0 + off * r
    end
end
