using MeromorphicLoci: _polygon_winding, _circle_winding

# Outline of a 3×3 block of unit columns, one edge per column side, as `_outline` builds it.
corners = (0, 3, 3 + 3im, 3im)
zs = ComplexF64[corners[s÷3+1] + (s % 3) * (corners[mod1(s ÷ 3 + 2, 4)] - corners[s÷3+1]) / 3 for s in 0:11]
edges = [(i, mod1(i + 1, 12)) for i in 1:12]

# Loci centred, and hugging an edge or a vertex, where a phase-only bisection aliases.
for z0 in (1.5 + 1.5im, 1.5 + 1e-4im, 1.0001 + 1e-4im, 1e-3 + 2.9im), k in -13:13
    @test _polygon_winding(z -> (z - z0)^k * cis(z), zs, edges)[1] == k
end
@test _polygon_winding(z -> NaN, zs, edges)[1] === nothing

# A centred locus keeps |f| constant on the circle, so the seed alone sets the range.
for m in -11:11
    @test _circle_winding(z -> z^m, 0.0im, 0.1)[1] == m
end
@test _circle_winding(z -> NaN, 0.0im, 0.1)[1] === nothing
# Noise swamping f on the circle: unresolved, within the budget.
w, ne = _circle_winding(z -> cis(1.0e15 * real(z) + 3.0e14 * imag(z)), 0.0im, 0.1)
@test w === nothing && ne <= 128
