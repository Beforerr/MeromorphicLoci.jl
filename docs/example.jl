using MeromorphicLoci, CairoMakie

# helices that graze at p ≈ 0.15 (0.004 apart) + a fixed zero–pole pair 0.03 apart
f(p) = z -> (z - 0.80 * cispi(p)) * (z + 0.10 + 0.12im) /
            ((z - 0.804 * cispi(-p + 0.3)) * (z + 0.13 + 0.12im))
zbox, pbox = (-0.98 - 0.98im, 0.98 + 0.98im), ((0.0, 1.0),)
s = findloci(f, zbox, pbox, Continuation(AAA()))

const JULIA_GREEN = "#389826"
const JULIA_RED = "#CB3C33"
color(b) = b.winding > 0 ? JULIA_GREEN : JULIA_RED

fig = Figure(; size=(700, 650), backgroundcolor=:white)
ax = Axis3(fig[1, 1];
    xlabel="Re z", ylabel="Im z", zlabel="p",
    title="$(length(s)) branches, $(s.nevals) evaluations",
    azimuth=0.55π, elevation=0.16π, perspectiveness=0.35,
)
for b in s
    xyz = ([real(x.z) for x in b], [imag(x.z) for x in b], [x.p[1] for x in b])
    scatterlines!(ax, xyz...; color=color(b), markersize=5)
end
Legend(fig[0, 1],
    [LineElement(; color=c) for c in (JULIA_GREEN, JULIA_RED)], ["zero", "pole"];
    orientation=:horizontal, framevisible=false)

out = joinpath(@__DIR__, "src", "assets", "survey.png")
save(out, fig; px_per_unit=2)
println("wrote ", out)
