using MeromorphicLoci: _polygon_winding

# Outline of a 3×3 block of unit columns, one edge per column side, as `_outline` builds it.
corners = (0, 3, 3 + 3im, 3im)
zs = ComplexF64[corners[s÷3+1] + (s % 3) * (corners[mod1(s ÷ 3 + 2, 4)] - corners[s÷3+1]) / 3 for s in 0:11]
edges = [(i, mod1(i + 1, 12)) for i in 1:12]

# Loci centred, and hugging an edge or a vertex, where a phase-only bisection aliases.
for z0 in (1.5 + 1.5im, 1.5 + 1e-4im, 1.0001 + 1e-4im, 1e-3 + 2.9im), k in -13:13
    @test _polygon_winding(z -> (z - z0)^k * cis(z), zs, edges)[1] == k
end
@test _polygon_winding(z -> NaN, zs, edges)[1] === nothing
