"""
    ArgumentPrinciple(; zres, minsep, zoom = 4)

Zeros and poles of `f(z)` in the box `zdomain = (zlo, zhi)` by
adaptive refinement of cells where the phase of `f` winds.

- `zres` — resolution of the points before the polish, and the scale at which zeros and
  poles are told apart: those closer than about `zres` fuse into one of their net winding,
  e.g. two simple zeros read as one double zero.
- `minsep` — the closest zero–pole spacing the search must see: a closer pair can cancel
  and vanish, unflagged by `deficit`. Zeros (or poles) closer than `minsep` are still
  found, at extra cost. Sets the cost floor, about `2·area/minsep²` evaluations.
- `zoom` — how many times, at most 8, to halve below `zres` to separate a zero and a pole
  of equal order; a pair closer than `zres/2^zoom` is reported `unresolved`. Of unequal
  orders, they stay fused as under `zres`.

Besides the common fields, returns `deficit`: the box's count of zeros minus poles, less
the windings found (of both kinds, as `findzeropoles` reports them). Nonzero means some
were missed, often beside a found zero or pole of high order: divide known factors out of
`f`. `nothing` means the winding around the box can't be read: `f` is non-finite on its
edge, or a zero or pole sits there (enlarge the box slightly).
"""
struct ArgumentPrinciple
    zres::Float64
    minsep::Float64
    zoom::Int
    function ArgumentPrinciple(zres, minsep, zoom)
        zres > 0 || throw(ArgumentError("zres must be positive, got $zres"))
        minsep > 0 || throw(ArgumentError("minsep must be positive, got $minsep"))
        0 ≤ zoom ≤ _HEADROOM || throw(ArgumentError("zoom must be in 0:$_HEADROOM, got $zoom"))
        return new(zres, minsep, zoom)
    end
end
ArgumentPrinciple(; zres, minsep, zoom=4) = ArgumentPrinciple(Float64(zres), Float64(minsep), zoom)

# The search: `f` on the lattice, the cells firing at `zres`, and the dark cells.
function _search(f, zdomain, alg::ArgumentPrinciple, mapper)
    zbox = _zcorners(zdomain...)
    lat = Lattice(zbox, alg.zres, alg.minsep)
    memo = Memo(f, lat, mapper)
    roots = _cells(memo, (Int32(0), Int32(0)), (Int32(1) << lat.kz, Int32(1) << lat.kz), lat.zn .>> lat.kz)
    dark = Cell[]
    cells = _refine!(memo, roots, alg.zres, dark)
    return memo, cells, dark
end

# One locus per labelled group of cells, at the mean of its cell centres. The mean is
# within a cell of the locus, but an escalated group can sit nearer another than that,
# hence the reach.
function _enclose(f, zdomain, alg::ArgumentPrinciple, ::Val, mapper)
    memo, cells, dark = _search(f, zdomain, alg, mapper)
    (; lat) = memo
    wbox = _certify!(memo, cells, dark, alg.zres)
    groups = _label(memo, cells; alg.zoom)
    deficit = isnothing(wbox) ? nothing : wbox - sum(last, groups; init=0)
    while something(deficit, 0) != 0 && !isempty(cells) && _beside!(memo, cells, dark, alg.zres)
        groups = _label(memo, cells; alg.zoom)
        deficit = wbox - sum(last, groups; init=0)
    end
    zs = ComplexF64[sum(c -> _zcenter(lat, c), cs) / length(cs) for (cs, _) in groups]
    loci = map(enumerate(groups)) do (i, (_, w))
        others = (zs[j] for j in eachindex(zs) if j != i)
        Locus(zs[i], _reach(zs[i], others, alg.zres), w)
    end
    return (; loci, nevals=nevals(memo), info=(; deficit), pad=0.0)
end

# Argument-principle quadrant, CCW half-open: 1:[0,π/2) 2:[π/2,π] 3:(-π,-π/2) 4:[-π/2,0).
# f = 0 or NaN bins arbitrarily (q2/q4) — a singularity pinned to a corner is detected by
# the neighbouring cells and labelled by the outline vote (`_vote`)
@inline function quadrant(z)
    re, im = reim(z)
    return im ≥ 0 ? (re > 0 ? 1 : 2) : (re < 0 ? 3 : 4)
end

# Corners are addressed by integer lattice index, not raw coordinates, so cells meeting at
# a corner share one evaluation. The box is tiled by near-square root cells (`_zroots`),
# each halved `kz` times.
struct Lattice
    zo::NTuple{2,Float64}
    zh::NTuple{2,Float64}
    zn::NTuple{2,Int32}             # index extent per axis, root cells × 2^kz
    kz::Int
end

# Headroom below the `zres` level so escalation can index sub-`zres` corners
const _HEADROOM = 8

# Root cells: the base mesh, refined to diagonal `dbase` (`minsep`) before any winding is
# read, tiled near-square with the fewest cells across the short axis.
function _zroots(w, dbase)
    L, S = max(w...), min(w...)
    for ns in ceil(Int, min(S / dbase, 2.0^30)):(1<<30)
        nl = round(Int, ns * L / S)
        hypot(L / nl, S / ns) ≤ dbase && return w[1] ≥ w[2] ? (nl, ns) : (ns, nl)
    end
    throw(ArgumentError("minsep=$dbase is too fine for this box: its lattice exceeds Int32 indices"))
end

function Lattice(zbox, zres, minsep)
    zlo, zhi = zbox
    w = (real(zhi) - real(zlo), imag(zhi) - imag(zlo))
    all(>(0), w) || throw(ArgumentError("the z-box must have positive width and height, got $zbox"))
    nz = _zroots(w, minsep)
    lz = max(0, ceil(Int, log2(hypot((w ./ nz)...) / zres)))
    kmax = 30 - ceil(Int, log2(max(nz...)))
    lz ≤ kmax || throw(ArgumentError("zres=$zres is too fine for this box: its lattice exceeds Int32 indices"))
    kz = min(lz + _HEADROOM, kmax)
    return Lattice((real(zlo), imag(zlo)), w ./ (nz .* 2.0^kz), Int32.(nz) .<< kz, kz)
end

@inline Base.getindex(lat::Lattice, ix) = complex(lat.zo[1] + ix[1] * lat.zh[1], lat.zo[2] + ix[2] * lat.zh[2])

# A box of the lattice: `base` index of its low corner, `step` its index extent per axis,
# `corners` the memo slots of its corners (0,0), (1,0), (0,1), (1,1).
struct Cell
    base::NTuple{2,Int32}
    step::NTuple{2,Int32}
    corners::NTuple{4,Int32}
end

function _cells(memo, base, step, n)
    at(i, j) = (base[1] + Int32(i) * step[1], base[2] + Int32(j) * step[2])
    s = [_slot!(memo, at(i, j)) for i in 0:n[1], j in 0:n[2]]
    return vec([Cell(at(i, j), step, (s[i+1, j+1], s[i+2, j+1], s[i+1, j+2], s[i+2, j+2])) for i in 0:n[1]-1, j in 0:n[2]-1])
end

_split(memo, c::Cell) = _cells(memo, c.base, c.step .>> 1, (2, 2))

_zcenter(lat::Lattice, c::Cell) = complex(
    lat.zo[1] + (c.base[1] + c.step[1] / 2) * lat.zh[1],
    lat.zo[2] + (c.base[2] + c.step[2] / 2) * lat.zh[2])

@inline _zdiag(lat::Lattice, c::Cell) = hypot(c.step[1] * lat.zh[1], c.step[2] * lat.zh[2])

# `f` on the lattice, shared by the search, the certificate and the votes. `_slot!` claims a
# point's slot in `vals` and queues it; `_fill!` evaluates the queue as one batch.
struct Memo{F,P}
    f::F
    lat::Lattice
    slot::Dict{NTuple{2,Int32},Int32}
    vals::Vector{ComplexF64}
    queue::Vector{NTuple{2,Int32}}
    off::Base.RefValue{Int}    # evaluations off the lattice
    mapper::P
end
Memo(f, lat, mapper) = Memo(f, lat, Dict{NTuple{2,Int32},Int32}(), ComplexF64[], NTuple{2,Int32}[], Ref(0), mapper)

nevals(m::Memo) = length(m.vals) + m.off[]

function _slot!(m::Memo, ix)
    s = get!(m.slot, ix, Int32(length(m.slot) + 1))
    s > length(m.vals) + length(m.queue) && push!(m.queue, ix)
    return s
end

function _fill!(m::Memo)
    append!(m.vals, m.mapper(m.f, [m.lat[ix] for ix in m.queue]))
    empty!(m.queue)
    return
end

# Quadrant of `f` at a slot, 0 where it awaits `_fill!`; `_quadrant!` claims the slot of a
# lattice point first.
@inline _quadrant(m::Memo, s::Int32) = s ≤ length(m.vals) ? quadrant(@inbounds m.vals[s]) : 0
@inline _quadrant!(m::Memo, ix) = _quadrant(m, _slot!(m, ix))

# `f` at outline points in lattice coordinates: memoised on the lattice, evaluated directly
# where a bisection passes below it.
function (m::Memo)(pts)
    slots = [all(isinteger, x) ? _slot!(m, Int32.(x)) : Int32(0) for x in pts]
    _fill!(m)
    vals = [s == 0 ? complex(NaN) : m.vals[s] for s in slots]
    off = findall(iszero, slots)
    vals[off] = m.mapper(m.f, [m.lat[pts[i]] for i in off])
    m.off[] += length(off)
    return vals
end

@inline _quadrants(memo, c::Cell) = map(s -> _quadrant(memo, s), c.corners)

# Signed phase winding of a cell's corner quadrants (in `_quadrants` order) walked CCW: +4
# for an enclosed zero, −4 for a pole, `nothing` if a step jumps two quadrants (±π,
# direction ambiguous).
@inline function _cellwinding(q)
    w = 0
    for (a, b) in ((q[1], q[2]), (q[2], q[4]), (q[4], q[3]), (q[3], q[1]))
        w += @something(_qstep(a, b), return nothing)
    end
    return w
end

# Quadrant step from `a` to `b`, `nothing` for an ambiguous ±π.
@inline function _qstep(a, b)
    d = mod(b - a, 4)
    return d == 2 ? nothing : d == 3 ? -1 : d
end

# Whether a cell fires: the walk around its corners winds or is ambiguous, or the guard
# (docs/src/argument_principle.md, Search) finds ≥ 3 quadrants at its corners and at the
# midpoints of its quadrant-crossing edges. A midpoint is a lattice point, reread by the
# halves if the cell splits. `nothing` while a point it needs is unevaluated.
function _fires(memo, c::Cell)
    q = _quadrants(memo, c)
    any(iszero, q) && return nothing
    w = _cellwinding(q)
    (isnothing(w) || w != 0) && return true
    m = 0
    for x in q
        m |= 1 << x
    end
    count_ones(m) ≥ 3 && return true
    all(>(1), c.step) || return false
    h = c.step .>> 1
    # (corner, corner, midpoint in half-steps) for each edge
    for (k1, k2, x, y) in ((1, 2, 1, 0), (2, 4, 2, 1), (4, 3, 1, 2), (3, 1, 0, 1))
        q[k1] == q[k2] && continue
        m |= 1 << _quadrant!(memo, (c.base[1] + Int32(x) * h[1], c.base[2] + Int32(y) * h[2]))
    end
    isodd(m) && return nothing    # bit 0: a midpoint unevaluated
    return count_ones(m) ≥ 3
end

# `_fires` of each cell, evaluating the points they read in batches.
function _fires(memo, cells::Vector{Cell})
    fires = falses(length(cells))
    pending = collect(eachindex(cells))
    while !isempty(pending)
        filter!(pending) do i
            f = _fires(memo, cells[i])
            isnothing(f) && return true
            fires[i] = f
            return false
        end
        _fill!(memo)
    end
    return fires
end

# Refine `frontier` where it fires, down to diagonal `zres`, and return the firing cells
# there; cells that do not fire are pushed to `dark`. Breadth-first: each round's
# evaluations are known up front and batched.
function _refine!(memo, frontier, zres, dark=Cell[])
    out = Cell[]
    while !isempty(frontier)
        next = Cell[]
        for (c, fires) in zip(frontier, _fires(memo, frontier))
            if !fires
                push!(dark, c)
            elseif _zdiag(memo.lat, c) > zres
                append!(next, _split(memo, c))
            else
                push!(out, c)
            end
        end
        frontier = next
    end
    return out
end

# Certificate (docs/src/argument_principle.md): the box's winding, `nothing` if unreadable.
# Where the found cells' vote differs from it, dark cells that wind are searched again.
function _certify!(memo, cells, dark, zres)
    (; lat) = memo
    wbox = _polygon_winding(memo, _boxpoly(lat, (Int32(1) << lat.kz, Int32(1) << lat.kz))...)
    isnothing(wbox) && return nothing
    wf = isempty(cells) ? 0 : _vote(memo, cells)
    (isnothing(wf) || wf == wbox) && return wbox
    _recover!(memo, cells, dark, zres)
    return wbox
end

# Search again the dark cells that wind, by their own winding; `dark` keeps those that don't.
function _recover!(memo, cells, dark, zres)
    (; lat) = memo
    suspects = copy(dark)
    empty!(dark)
    while !isempty(suspects)
        frontier = Cell[]
        for c in suspects
            w = _polygon_winding(memo, _boxpoly(c)...)
            if w == 0
                push!(dark, c)
            elseif _zdiag(lat, c) > zres
                append!(frontier, _split(memo, c))
            else
                push!(isnothing(w) ? dark : cells, c)
            end
        end
        suspects = Cell[]
        append!(cells, _refine!(memo, frontier, zres, suspects))
    end
    return dark
end

# Search again the dark cells within a cell of `cells`. A coarse one can wind 0 holding a
# counted locus and a missed partner, so it is split; one at `zres` is read by its winding,
# as a k-fold locus centred in it fires neither check. Returns whether any cell was split
# or fired, since a split half can still hold both.
function _beside!(memo, cells, dark, zres)
    (; lat) = memo
    n = length(cells)
    st = cells[1].step
    near = map(dark) do d
        any(c -> all(i -> c.base[i] - st[i] ≤ d.base[i] + d.step[i] && d.base[i] ≤ c.base[i] + 2st[i], 1:2), cells)
    end
    coarse = [c for c in dark[near] if _zdiag(lat, c) > zres]
    suspects = [c for c in dark[near] if _zdiag(lat, c) ≤ zres]
    deleteat!(dark, near)
    isempty(coarse) || append!(cells, _refine!(memo, reduce(vcat, [_split(memo, c) for c in coarse]), zres, suspects))
    append!(dark, _recover!(memo, cells, suspects, zres))
    return !isempty(coarse) || length(cells) > n
end

# The outline of a box of `n` cells, each `st` lattice steps, through every cell corner.
function _boxpoly(st, n; base=(0, 0))
    nx, ny = n
    ring = [(0:nx-1) .=> 0; nx .=> (0:ny-1); (nx:-1:1) .=> ny; 0 .=> (ny:-1:1)]
    return [(Float64(base[1] + x * st[1]), Float64(base[2] + y * st[2])) for (x, y) in ring], [(i, mod1(i + 1, length(ring))) for i in eachindex(ring)]
end
_boxpoly(lat::Lattice, st) = _boxpoly(st, (Int(lat.zn[1] ÷ st[1]), Int(lat.zn[2] ÷ st[2])))
_boxpoly(c::Cell) = _boxpoly(c.step, (1, 1); c.base)
