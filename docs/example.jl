# Visualization: the zero/pole loci z*(p) as (Re z, Im z, p), one panel per
# README config — escalation (zoom) and aliasing (minsep) behaviors.
using MeromorphicLoci, CairoMakie

# helices that graze at p ≈ 0.15 (0.004 apart) + a fixed zero–pole pair 0.03 apart
f(z, p) = (z - 0.80 * cispi(p)) * (z + 0.10 + 0.12im) /
          ((z - 0.804 * cispi(-p + 0.3)) * (z + 0.13 + 0.12im))
box = ((-0.98 - 0.98im, 0.98 + 0.98im), (0.0, 1.0))

configs = [
    "zoom = 0" => (; zoom=0),                      # graze fuses helices: 1 branch, winding 0
    "default" => (;),                              # escalation splits them: 2 clean branches
    "minsep = 0.3, zoom = 0" => (; minsep=0.3, zoom=0),  # pair discovered, helices fused: 3
    "minsep = 0.3" => (; minsep=0.3),              # everything found and labeled: 4
]

const JULIA_GREEN = "#389826"
const JULIA_RED = "#CB3C33"
const JULIA_PURPLE = "#9558B2"
color(b) = winding(b) > 0 ? JULIA_GREEN : winding(b) < 0 ? JULIA_RED : JULIA_PURPLE

fig = Figure(; size=(1200, 1100), backgroundcolor=:white)
for (i, (label, kw)) in enumerate(configs)
    s = survey(f, box; zres=0.02, kw...)
    ax = Axis3(fig[fld1(i, 2), mod1(i, 2)];
        xlabel="Re z", ylabel="Im z", zlabel="p",
        title="$label — $(length(s)) branches, $(s.nevals) evals",
        azimuth=0.55π, elevation=0.16π, perspectiveness=0.35,
    )
    for b in s
        xyz = ([real(x.z) for x in b], [imag(x.z) for x in b], [x.p[1] for x in b])
        # a fused branch holds several loci, so a p-sorted polyline scribbles
        # across all of them — scatter it instead
        winding(b) == 0 ? scatter!(ax, xyz...; color=color(b), markersize=2) :
        lines!(ax, xyz...; color=color(b))
    end
end

Legend(fig[0, :],
    [LineElement(; color=c) for c in (JULIA_GREEN, JULIA_RED, JULIA_PURPLE)],
    ["zero", "pole", "indeterminate"];
    orientation=:horizontal, framevisible=false)

out = joinpath(@__DIR__, "survey.png")
save(out, fig; px_per_unit=2)
println("wrote ", out)
