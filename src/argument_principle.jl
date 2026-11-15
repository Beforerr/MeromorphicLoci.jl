"""
    ArgumentPrinciple(; zres, minsep = nothing, zoom = 4)

Zeros and poles of `f(z)` in the box `zdomain = (zlo, zhi)` by
adaptive refinement of cells where the phase of `f` winds.

- `zres` — resolution of the points before the polish, and the scale at which loci are
  told apart: loci closer than about `zres` fuse into one of their net winding, e.g. two
  simple zeros read as one double zero.
- `minsep` — cells are refined to this size before any is discarded. A zero and a pole
  closer than `minsep` can cancel inside one cell and vanish silently, so keep it below
  their closest spacing; below `zres`, it is the resolution in place of `zres`. Loci of
  one kind closer than `minsep` can escape the search; the box's count recovers them at
  extra cost (`deficit` below). `minsep` also sets the cost floor, about
  `2·area/minsep²` evaluations. `nothing` means the box diagonal / 8, or `zres` if that
  is coarser.
- `zoom` — halvings below `zres`, at most 8, to separate a fused zero–pole pair. A pair of
  equal orders closer than `zres/2^zoom` is reported `unresolved`; one of different orders
  closer than about `zres` stays fused, as one locus of their net winding.

Besides the common fields, returns `deficit`: zeros minus poles in the box that the loci
found (both kinds, as `findzeropoles` reports them) leave unaccounted. Nonzero means
missed loci, often a zero or pole of order ≥ 4: divide known factors out of `f`.
`nothing` means the winding around the box can't be read: `f` is non-finite on its edge,
or a locus sits there (pad the box). A zero–pole pair cancels in the count, so `deficit`
can't see one missed.
"""
struct ArgumentPrinciple
    zres::Float64
    minsep::Union{Nothing,Float64}
    zoom::Int
end
ArgumentPrinciple(; zres, minsep=nothing, zoom=4) = ArgumentPrinciple(Float64(zres), minsep, zoom)

# One locus per labelled group of cells, at the mean of its cell centres. The mean is
# within a cell of the locus, but an escalated group can sit nearer another than that,
# hence the reach.
function _enclose(f, zdomain, alg::ArgumentPrinciple, ::Val, mapper)
    zbox = _zcorners(zdomain...)
    minsep = Float64(something(alg.minsep, max(abs(zbox[2] - zbox[1]) / 8, alg.zres)))
    lat = Lattice(zbox, alg.zres, minsep)
    memo = Memo(f, lat, mapper)
    st = (Int32(1) << lat.kz, Int32(1) << lat.kz)
    roots = [Cell(memo, (Int32(i - 1) * st[1], Int32(j - 1) * st[2]), st) for i in 1:(lat.zn[1]>>lat.kz), j in 1:(lat.zn[2]>>lat.kz)]
    tree = Tree(memo, roots)
    cells = _refine!(tree, vec(roots), alg.zres)
    wbox = _certify!(tree, cells, alg.zres)
    groups = _label(tree, cells; alg.zoom)
    deficit = isnothing(wbox) ? nothing : wbox - sum(last, groups; init=0)
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
    zres > 0 || throw(ArgumentError("zres must be positive, got $zres"))
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

# A classified cell is dark, or fires on its face or only on the guard (`_classify`).
const _PENDING, _DARK, _FACE, _GUARD, _SPLIT = Int8(0), Int8(1), Int8(2), Int8(3), Int8(4)

# A box of the lattice: `base` index of its low corner, `step` its index extent per axis,
# `corners` the memo slots of its corners (0,0), (1,0), (0,1), (1,1).
mutable struct Cell
    base::NTuple{2,Int32}
    step::NTuple{2,Int32}
    corners::NTuple{4,Int32}
    state::Int8
    children::Vector{Cell}
end
function Cell(memo, base, step)
    corners = ntuple(i -> _slot!(memo, (base[1] + Int32((i - 1) & 1) * step[1], base[2] + Int32((i - 1) >> 1) * step[2])), 4)
    return Cell(base, step, corners, _PENDING, Cell[])
end

function split!(memo, c::Cell)
    h = c.step .>> 1
    for j in 0:1, i in 0:1
        push!(c.children, Cell(memo, (c.base[1] + Int32(i) * h[1], c.base[2] + Int32(j) * h[2]), h))
    end
    c.state = _SPLIT
    return c.children
end

_firing(c::Cell) = c.state == _FACE || c.state == _GUARD

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

# Quadrant of `f` at a slot, or at a lattice point; 0 where it awaits `_fill!`.
@inline _quadrant(m::Memo, s::Int32) = s ≤ length(m.vals) ? quadrant(@inbounds m.vals[s]) : 0
@inline Base.getindex(m::Memo, ix) = _quadrant(m, _slot!(m, ix))

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

# Signed phase winding of the quadrants at points walked CCW: +4 for an enclosed zero, −4
# for a pole, `nothing` if any step jumps two quadrants (an unresolved ±π step whose
# direction is ambiguous — the cell must be refined).
@inline function _face_winding(q::Tuple)
    w = 0
    n = length(q)
    @inbounds for i in 1:n
        w += @something(_qstep(q[i], q[i%n+1]), return nothing)
    end
    return w
end
# Around a cell, from its corner quadrants in `_quadrants` order.
@inline _cellwinding(q) = _face_winding((q[1], q[2], q[4], q[3]))

# Quadrant step from `a` to `b`, `nothing` for an ambiguous ±π.
@inline function _qstep(a, b)
    d = mod(b - a, 4)
    return d == 2 ? nothing : d == 3 ? -1 : d
end

_fw_fires(fw) = isnothing(fw) || fw != 0

# A cell fires on its face (`_FACE`) if the walk around its corners winds or is ambiguous.
#
# A cell whose face is clean still fires on the guard (`_GUARD`, docs/src/argument_principle.md,
# Search) if the quadrants at its corners and at the midpoints of its quadrant-crossing
# edges number ≥ 3. A midpoint is a lattice point, reread by the children if the cell
# splits.
#
# Returns `_PENDING` while a point it needs reads 0 (unevaluated).
function _classify(memo, c::Cell)
    q = _quadrants(memo, c)
    any(iszero, q) && return _PENDING
    _fw_fires(_cellwinding(q)) && return _FACE
    m = 0
    for x in q
        m |= 1 << x
    end
    count_ones(m) ≥ 3 && return _GUARD
    all(>(1), c.step) || return _DARK
    h = c.step .>> 1
    # (corner, corner, midpoint in half-steps) for each edge
    for (k1, k2, x, y) in ((1, 2, 1, 0), (2, 4, 2, 1), (4, 3, 1, 2), (3, 1, 0, 1))
        q[k1] == q[k2] && continue
        m |= 1 << memo[(c.base[1] + Int32(x) * h[1], c.base[2] + Int32(y) * h[2])]
    end
    isodd(m) && return _PENDING
    return count_ones(m) ≥ 3 ? _GUARD : _DARK
end

# Set each cell's state from `_classify`, evaluating the points it reads in batches.
function _classify!(memo, cells)
    while !isempty(cells)
        cells = filter(cells) do c
            c.state = _classify(memo, c)
            c.state == _PENDING
        end
        _fill!(memo)
    end
end

# The refinement state the labelling continues from: escalation refines further.
struct Tree{M}
    memo::M
    roots::Matrix{Cell}
end

# Refine `frontier` where it fires, down to diagonal `zres`, and return the firing cells
# there. Breadth-first: each round's evaluations are known up front and batched.
function _refine!(tree::Tree, frontier, zres)
    out = Cell[]
    while !isempty(frontier)
        _classify!(tree.memo, frontier)
        next = Cell[]
        for c in frontier
            _firing(c) || continue
            _zdiag(tree.memo.lat, c) > zres ? append!(next, split!(tree.memo, c)) : push!(out, c)
        end
        frontier = next
    end
    return out
end

# Certificate (docs/src/argument_principle.md): the box's winding, `nothing` if unreadable.
# Where the found cells' vote differs from it, dark leaves that wind are searched again.
function _certify!(tree::Tree, cells, zres)
    (; memo) = tree
    (; lat) = memo
    wbox = _polygon_winding(memo, _boxpoly(lat, (Int32(1) << lat.kz, Int32(1) << lat.kz))...)
    isnothing(wbox) && return nothing
    wf = isempty(cells) ? 0 : _vote(tree, cells)
    (isnothing(wf) || wf == wbox) && return wbox
    suspects = _darkleaves(tree.roots)
    while !isempty(suspects)
        frontier = Cell[]
        for c in suspects
            w = _polygon_winding(memo, _boxpoly(c)...)
            w == 0 && continue
            if _zdiag(lat, c) > zres
                append!(frontier, split!(memo, c))
            elseif !isnothing(w)
                c.state = _FACE
                push!(cells, c)
            end
        end
        append!(cells, _refine!(tree, frontier, zres))
        suspects = _darkleaves(frontier)
    end
    return wbox
end

# The outline of a box of `n` cells, each `st` lattice steps, through every cell corner.
function _boxpoly(st, n; base=(0, 0))
    nx, ny = n
    ring = [(0:nx-1) .=> 0; nx .=> (0:ny-1); (nx:-1:1) .=> ny; 0 .=> (ny:-1:1)]
    return [(Float64(base[1] + x * st[1]), Float64(base[2] + y * st[2])) for (x, y) in ring], [(i, mod1(i + 1, length(ring))) for i in eachindex(ring)]
end
_boxpoly(lat::Lattice, st) = _boxpoly(st, (Int(lat.zn[1] ÷ st[1]), Int(lat.zn[2] ÷ st[2])))
_boxpoly(c::Cell) = _boxpoly(c.step, (1, 1); c.base)

function _darkleaves(cells)
    out = Cell[]
    stack = vec(collect(cells))
    while !isempty(stack)
        c = pop!(stack)
        isempty(c.children) ? (c.state == _DARK && push!(out, c)) : append!(stack, c.children)
    end
    return out
end
