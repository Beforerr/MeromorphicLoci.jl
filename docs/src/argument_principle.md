# How `ArgumentPrinciple` works

`ArgumentPrinciple` finds the zeros and poles (**loci**) of `f` in a box in four steps:

1. **Search**: refine a mesh of cells toward loci, reading only the phase of `f` at cell corners.
2. **Certificate**: count the box's loci by its winding; search again where the cells fall short.
3. **Label**: group the cells and count each group by the argument principle.
4. **Escalate**: refine groups whose label is unclear.

## Cells

- The box is tiled by near-square root cells of diagonal at most `minsep`.
- Cells halve down to `zres`, and up to `zoom` levels further when escalating.
- Cells stay near-square whatever the box's shape: four corners of a long, thin cell miss
  phase turns.

## Search

Cells are refined breadth-first. A cell **fires**, and is split further, when either check
says it may hold a locus:

- **Walk**: step the quadrant of `f` around its four corners. A nonzero winding fires
  (+4 for a simple zero, −4 for a simple pole), and so does a jump across two quadrants,
  whose direction is ambiguous.
- **Guard**: the corners, plus the midpoints of edges that change quadrant, show at least
  3 quadrants. Then both Re f = 0 and Im f = 0 cross the cell, and they meet only at loci.

Firing cells split down to `zres`. The rest go dark.

```@raw html
<details><summary>Code</summary>
```

```@example ap
using MeromorphicLoci, CairoMakie, Base64
using MeromorphicLoci: Lattice, Memo, Cell, Tree,
    _zcorners, _refine!, _certify!, _vote, _group, _outline, _mixed, _deepen!,
    _polygon_winding, _cellwinding, _qstep, quadrant, _darkleaves

const GREEN, RED, GREY = "#389826", "#CB3C33", "#777777"
const QCOLORS = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7"]   # quadrants 1–4
const BOX = (-1 - 1im, 1 + 1im)
signcolor(w) = isnothing(w) || w == 0 ? GREY : w > 0 ? GREEN : RED
signed(w) = isnothing(w) ? "unusable" : w > 0 ? "+$w" : w < 0 ? "−$(-w)" : "0"

# Shown at its logical width, so text keeps its size, and rendered at twice that. Inlined
# as a data URI: a relative src breaks under prettyurls, which nest each page a level down.
CairoMakie.activate!(px_per_unit=2)
function embed(fig, name)
    png = base64encode(io -> show(io, MIME"image/png"(), fig))
    HTML("<img src=\"data:image/png;base64,$png\" width=\"$(size(fig.scene)[1])\" alt=\"$name\">")
end

rect(lat, c::Cell) = Rect2f(reim(lat[c.base])..., c.step[1] * lat.zh[1], c.step[2] * lat.zh[2])
grid(lo, hi, n) = (range(real(lo), real(hi), n), range(imag(lo), imag(hi), n))
function phase!(ax, f, lo, hi; n=400)
    xs, ys = grid(lo, hi, n)
    heatmap!(ax, xs, ys, [angle(f(complex(x, y))) for x in xs, y in ys];
        colormap=:phase, colorrange=(-π, π), alpha=0.55)
    limits!(ax, real(lo), real(hi), imag(lo), imag(hi))
end
function loci!(ax, zeros, poles; markersize=13)
    scatter!(ax, reim.(zeros); marker=:circle, color=:white, strokecolor=:black, strokewidth=1.5, markersize)
    scatter!(ax, reim.(poles); marker=:xcross, color=:black, markersize)
end
# The mesh the search has reached: its dark cells, its evaluations, and a marker on each firing cell.
function mesh!(ax, memo, cells, dark)
    lat = memo.lat
    poly!(ax, rect.(Ref(lat), dark); color=:transparent, strokecolor=(:black, 0.45), strokewidth=0.5)
    scatter!(ax, [Point2f(reim(lat[ix])...) for ix in keys(memo.slot)]; color=:black, markersize=2)
    scatter!(ax, [Point2f(reim(lat[c.base.+c.step./2])...) for c in cells]; marker=:rect, color=:black, strokecolor=:white, strokewidth=1, markersize=9)
end
# The search as `_enclose` runs it: the tree, and its firing cells at `zres`.
function search(f, box, alg)
    zbox = _zcorners(box...)
    minsep = Float64(something(alg.minsep, max(abs(zbox[2] - zbox[1]) / 8, alg.zres)))
    lat = Lattice(zbox, alg.zres, minsep)
    memo = Memo(f, lat, map)
    st = (Int32(1) << lat.kz, Int32(1) << lat.kz)
    roots = [Cell(memo, (Int32(i - 1) * st[1], Int32(j - 1) * st[2]), st) for i in 1:(lat.zn[1]>>lat.kz), j in 1:(lat.zn[2]>>lat.kz)]
    tree = Tree(memo, roots)
    return tree, _refine!(tree, vec(roots), alg.zres)
end
# Outline of a group in the z-plane, its vote, the points the vote evaluated, and whether it is mixed.
function outline(tree, cells)
    lat = tree.memo.lat
    vs, edges = _outline(lat, cells)
    pts = NTuple{2,Float64}[]
    _polygon_winding(p -> (append!(pts, p); tree.memo(p)), vs, edges)
    segs = [Point2f(reim(lat[vs[i]])...) for e in edges for i in e]
    return (; segs, w=_vote(tree, cells), mixed=_mixed(tree.memo, cells), pts=[Point2f(reim(lat[p])...) for p in pts])
end

f = z -> (z - 0.3) * (z + 0.4im)^2 / (z - 0.6)
tree, cells = search(f, BOX, ArgumentPrinciple(zres=0.05))
fig = Figure(; size=(420, 400))
ax = Axis(fig[1, 1]; aspect=DataAspect(), xlabel="Re z", ylabel="Im z")
phase!(ax, f, BOX...)
mesh!(ax, tree.memo, cells, _darkleaves(tree.roots))
nothing # hide
```

```@raw html
</details>
```

```@example ap
embed(fig, "search") # hide
```

*The search for (z − 0.3)(z + 0.4i)²/(z − 0.6) at `zres = 0.05`, over the phase of `f`
(hue). Dots are evaluations, squares the cells firing at `zres`.*

The guard is needed because the walk reads each edge as the shortest step: a +270° turn
reads −90°. A neighbouring locus can bend an edge that far (middle panel below). Neither
check sees a 4-fold zero centred in a cell (right): each edge turns a full 2π, so all
corners share one quadrant.

```@raw html
<details><summary>Code</summary>
```

```@example ap
# One cell over the quadrant of f. Corners are in `_quadrants` order.
function cell!(ax, f, lo, hi; zeros=ComplexF64[], view=(lo, hi))
    h = hi - lo
    pad = 0.3 * max(real(view[2] - view[1]), imag(view[2] - view[1]))
    vlo, vhi = view[1] - pad * (1 + im), view[2] + pad * (1 + im)
    xs, ys = grid(vlo, vhi, 500)
    heatmap!(ax, xs, ys, [quadrant(f(complex(x, y))) for x in xs, y in ys];
        colormap=cgrad(QCOLORS, 4; categorical=true), colorrange=(0.5, 4.5), alpha=0.4)
    limits!(ax, real(vlo), real(vhi), imag(vlo), imag(vhi))
    corners = [lo, lo + real(h), lo + imag(h) * im, hi]
    q = quadrant.(f.(corners))
    lines!(ax, reim.(corners[[1, 2, 4, 3, 1]]); color=:black, linewidth=2)
    w = _cellwinding(q)
    walkfires = isnothing(w) || w != 0
    seen = Set(q)
    guarded = !walkfires && length(seen) < 3   # the only cells whose midpoints are read
    for (i, j) in ((1, 2), (2, 4), (4, 3), (3, 1))
        m = (corners[i] + corners[j]) / 2
        out = (m - (lo + hi) / 2) / abs(m - (lo + hi) / 2)
        step = _qstep(q[i], q[j])
        text!(ax, reim(m + 0.12 * abs(h) * out)...; text=isnothing(step) ? "±2" : signed(step),
            align=(:center, :center), fontsize=15, font=:bold)
        guarded && q[i] != q[j] || continue
        push!(seen, quadrant(f(m)))
        scatter!(ax, reim(m)...; marker=:diamond, color=QCOLORS[quadrant(f(m))], strokecolor=:black, strokewidth=1.5, markersize=15)
    end
    scatter!(ax, reim.(corners); color=QCOLORS[q], strokecolor=:black, strokewidth=1.5, markersize=17)
    loci!(ax, zeros, ComplexF64[])
    ax.title = walkfires ? "walk $(isnothing(w) ? "ambiguous" : signed(w)): fires" :
               length(seen) ≥ 3 ? "walk 0, $(length(seen)) quadrants: guard fires" : "walk 0, one quadrant: dark"
end

fig = Figure(; size=(760, 300))
axs = [Axis(fig[1, i]; aspect=DataAspect()) for i in 1:3]
hidedecorations!.(axs)
a = 0.37 + 0.62im
cell!(axs[1], z -> (z - a) * cis(0.4), 0im, 1 + 1im; zeros=[a])
# a pair of test/runtests.jl, and the lattice cell of `BOX` that holds `a`
a, b = 0.367 + 0.069im, 0.182 + 0.075im
cell!(axs[2], z -> (z - a) * (z - b), 0.25 + 0im, 0.375 + 0.125im; zeros=[a, b], view=(real(b) + 0im, 0.375 + 0.125im))
cell!(axs[3], z -> (z - 0.5 - 0.5im)^4 * cis(0.5), 0im, 1 + 1im; zeros=[0.5 + 0.5im])
Legend(fig[2, :], [PolyElement(; color=(c, 0.6)) for c in QCOLORS], ["quadrant $i" for i in 1:4];
    orientation=:horizontal, framevisible=false)
nothing # hide
```

```@raw html
</details>
```

```@example ap
embed(fig, "walk_guard") # hide
```

*Corners (circles) and guard midpoints (diamonds) over the quadrant of `f`, with each
edge's quadrant step. Left: one simple zero. Middle: one of two zeros (○). Right: a 4-fold
zero.*

## Certificate

The winding of `f` around the box counts its loci, zeros minus poles.

- The search passes if the firing cells, grown by one cell, wind the same.
- Otherwise each dark cell is checked by its own winding. It holds no firing cell, so a
  nonzero winding is a missed locus: the cell is searched again, and fires if it still
  winds at `zres`.
- What the labels (below) leave of the box's count is reported as `deficit` (see
  [`ArgumentPrinciple`](@ref)).

```@raw html
<details><summary>Code</summary>
```

```@example ap
# Eight zeros inside one root cell, before and after the count is checked.
c = 0.3 + 0.2im
f = z -> (z - c)^8 - 0.1^8
alg = ArgumentPrinciple(zres=0.02, minsep=0.36)
tree, cells = search(f, BOX, alg)
fig = Figure(; size=(720, 380))
for (i, title) in enumerate(("search", "after the certificate"))
    i == 2 && _certify!(tree, cells, alg.zres)
    local ax = Axis(fig[1, i]; aspect=DataAspect(), xlabel="Re z", ylabel="Im z", title="$title: $(length(cells)) firing cells")
    phase!(ax, f, BOX...)
    mesh!(ax, tree.memo, cells, _darkleaves(tree.roots))
end
nothing # hide
```

```@raw html
</details>
```

```@example ap
embed(fig, "certificate") # hide
```

*Eight zeros on a circle 0.1 across, inside one root cell (`minsep = 0.36`). The search
finds some; the box counts 8, so the dark cells that wind are searched again.*

## Labels

- Firing cells at `zres` that share an edge or a corner form a **group**.
- A group's label is the winding of `f` around its **outline**: its cells grown by one
  cell and clipped to the box. Separate groups are at least two cells apart, so an outline
  holds no other group's cells.
- Outlines are integrated adaptively: each edge is bisected until every step turns the
  phase by less than π/2 and changes |f| by less than 2×.
- A group whose cells' walks show both signs may hold a zero and a pole that a net vote
  hides, e.g. a simple zero beside a double pole votes −1.

```@raw html
<details><summary>Code</summary>
```

```@example ap
zs, ps = [0.10 + 0.07im, -0.2 - 0.15im], [0.27 + 0.12im]
f = z -> (z - zs[1]) * (z - zs[2])^2 / (z - ps[1])
tree, cells = search(f, BOX, ArgumentPrinciple(zres=0.05))
fig = Figure(; size=(500, 380))
ax = Axis(fig[1, 1]; aspect=DataAspect(), xlabel="Re z", ylabel="Im z")
phase!(ax, f, -0.42 - 0.36im, 0.5 + 0.32im)
for g in _group(cells)
    o = outline(tree, g)
    poly!(ax, rect.(Ref(tree.memo.lat), g); color=(signcolor(o.w), 0.55), strokecolor=:black, strokewidth=0.5)
    linesegments!(ax, o.segs; color=:black, linewidth=2.5)
    scatter!(ax, o.pts; color=:white, strokecolor=:black, strokewidth=1, markersize=6)
    text!(ax, minimum(first, o.segs), maximum(last, o.segs); text="w = $(signed(o.w))", fontsize=16, font=:bold, align=(:left, :bottom))
end
loci!(ax, zs, ps; markersize=10)
nothing # hide
```

```@raw html
</details>
```

```@example ap
embed(fig, "labels") # hide
```

*Groups of firing cells (filled), their outlines, and the points each vote evaluated, for
a simple zero, a simple pole and a double zero at `zres = 0.05`.*

## Escalation

- A group whose cells show both signs, votes 0, or cannot be read is refined one level
  below `zres`, regrouped and labelled again, up to `zoom` levels.
- A group that votes 0 or cannot be read, and has nothing firing one level finer, is
  dropped: its cells fired on a misread phase (a locus just outside the box, a phase
  turning faster than a cell), not on a locus.
- A group still unlabelled at the last level is reported `unresolved`: a zero and a pole
  closer than a last-level cell.
- Each labelled group is reported at the mean of its cell centres, then polished.

```@raw html
<details><summary>Code</summary>
```

```@example ap
# One panel per level until every group is labelled.
rows = [
    (; f=z -> (z - 0.985 - 0.011im) / (z - 1.012 - 0.006im), box=(0 - 0.4im, 2 + 0.4im), zres=0.05,
        zeros=[0.985 + 0.011im], poles=[1.012 + 0.006im], view=(0.94 - 0.04im, 1.06 + 0.06im)),
]
fig = Figure(; size=(580, 250))
for (r, (; f, box, zres, zeros, poles, view)) in enumerate(rows)
    local tree, cells = search(f, box, ArgumentPrinciple(; zres))
    undecided(g) = (o = outline(tree, g); o.mixed || something(o.w, 0) == 0)
    levels = [_group(cells)]
    while any(undecided, levels[end]) && length(levels) < 3
        push!(levels, reduce(vcat, [undecided(g) ? _group(_deepen!(tree, g)) : [g] for g in levels[end]]))
    end
    for (i, groups) in enumerate(levels)
        local ax = Axis(fig[r, i]; aspect=DataAspect(), title=i == 1 ? "zres" : "zres / $(2^(i - 1))")
        hidedecorations!(ax)
        phase!(ax, f, view...)
        for g in groups
            o = outline(tree, g)
            poly!(ax, rect.(Ref(tree.memo.lat), g); color=(signcolor(o.w), 0.6), strokecolor=:black, strokewidth=0.5)
            linesegments!(ax, o.segs; color=:black, linewidth=2.5)
            left = something(o.w, 0) ≥ 0    # opposite corners, so neighbouring labels stay apart
            text!(ax, (left ? minimum : maximum)(first, o.segs), maximum(last, o.segs); text="w = $(signed(o.w))" * (o.mixed ? ", mixed" : ""),
                fontsize=15, font=:bold, align=(left ? :left : :right, :bottom))
        end
        loci!(ax, zeros, poles; markersize=10)
    end
end
nothing # hide
```

```@raw html
</details>
```

```@example ap
embed(fig, "escalation") # hide
```

*A zero (○) and a pole (×) 0.027 apart at `zres = 0.05` form one group voting 0; two levels
finer they split. Cells are green, red or grey as their group votes positive, negative or
0.*

## What can go wrong

| Case | Effect | Recovered by |
|---|---|---|
| Zero–pole pair much closer than `minsep` | cancels on the walk and in the box's count; its phase fades within a few pair-widths | nothing |
| Same-sign loci closer than `minsep`, or phase turning > π/2 between corners | walk reads 0 | certificate |
| k-fold locus centred in a cell | walk misreads it (triple −4, 4-fold 0) at every scale, since (z − z0)^k is scale-free | certificate |
| Simple locus with a second one hugging the same edge | can slip past walk and guard | certificate, or fused into its partner's label |
| k-fold locus on a lattice point | only cells beside it fire; outlines through it can't be read and it is dropped | flagged by `deficit` |
