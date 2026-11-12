"""
    ArgumentPrinciple(; zres, minsep=nothing, zoom=4)

Adaptive refinement of `(z, p...)` near the zero/pole loci by phase winding.
`zdomain = (zlo, zhi)` gives opposite complex corners of the z-rectangle. `pdomain` gives
one entry per parameter: an interval `(lo, hi)`, sampled at 9 equispaced knots, or a
sorted vector of knots.

- `zres` — z-plane resolution of the reported loci, and the scale at which they are told
  apart: loci closer than about `zres` fuse into one branch of their net winding. Two
  simple zeros read as one double zero; a zero–pole pair reads 0 until `zoom` separates it.
- `minsep` — cells are refined to this size before any is discarded. Distinct loci closer
  than `minsep` can cancel inside one cell and vanish silently, so keep it below their
  closest spacing. It also sets the cost floor, about `2·area/minsep²` evaluations per
  parameter knot. `nothing` means the z-diagonal / 8.
- `zoom` — halvings below `zres` a branch of indeterminate winding is refined to separate
  a fused zero–pole pair. Winding 0 past that means a pair closer than `zres/2^zoom`, or
  loci crossing at some `p`.

Loci are discovered at the parameter knots; one existing only between two knots can be
missed. Between knots `p` is refined only as finely as the loci move, so evaluations
concentrate where loci meet or a zero–pole pair annihilates.

What else it can miss or leave indeterminate:
- a zero or pole of order ≥ 3 — occasionally at order 3, often from order 4; divide known
  factors out of `f`;
- a locus within about a cell of the box edge, so pad the box;
- loci touching another: a zero of `f` along all of `p` (a symmetry-forced factor `z^k`,
  say) fuses with every locus crossing it; divide such factors out of `f` first.

`f` is evaluated across all Julia threads: `f(p...)` and the functions it returns must be
thread-safe. How the search works: `docs/argument_principle.md`.
"""
struct ArgumentPrinciple
    zres::Float64
    minsep::Union{Nothing,Float64}
    zoom::Int
end
ArgumentPrinciple(; zres, minsep=nothing, zoom=4) =
    ArgumentPrinciple(Float64(zres), minsep, zoom)

_loci(f, zdomain, pdomain::Tuple, alg::ArgumentPrinciple) =
    link(scan(f, (zdomain, pdomain...); alg.zres, alg.minsep); alg.zoom)

# One locus per branch, at its mean sample. The mean is within a cell of the locus, but an
# escalated branch can sit nearer another than that, hence the reach.
function _enclose(f, zdomain, alg::ArgumentPrinciple, ::Val)
    s = _loci(Returns(f), zdomain, (), alg)
    zs = ComplexF64[sum(x -> x.z, b) / length(b) for b in s]
    loci = map(enumerate(s)) do (i, b)
        others = (zs[j] for j in eachindex(zs) if j != i)
        (z=zs[i], r=_reach(zs[i], others, alg.zres), w=winding(b))
    end
    return (; loci, nevals=nevals(s), info=(;))
end

# Argument-principle quadrant, CCW half-open: 1:[0,π/2) 2:[π/2,π] 3:(-π,-π/2) 4:[-π/2,0).
# f = 0 or NaN bins arbitrarily (q2/q4) — pinned singularities are
# detected by neighboring faces and labelled by the outline vote (`_votes`)
@inline function quadrant(z)
    re, im = reim(z)
    return im ≥ 0 ? (re > 0 ? 1 : 2) : (re < 0 ? 3 : 4)
end

# Corners are addressed by integer lattice index, not raw coordinates, so cells meeting at
# a corner share one evaluation. The z axes halve the box `kz` times. A parameter axis
# halves each knot interval `kp` times: knots are lattice points and splits stay dyadic
# within a knot interval, whatever the knot spacing.
struct Lattice{M}
    zo::NTuple{2,Float64}
    zh::NTuple{2,Float64}
    zn::Int32                       # z index extent, 2^kz
    knots::NTuple{M,Vector{Float64}}
    kp::Int
end

# Headroom below the `zres` level so `link`'s escalation can index sub-`zres` corners
const _HEADROOM = 8

function Lattice(zbox, knots::NTuple{M,Vector{Float64}}, zres) where {M}
    zres > 0 || throw(ArgumentError("zres must be positive, got $zres"))
    zlo, zhi = zbox
    w = (real(zhi) - real(zlo), imag(zhi) - imag(zlo))
    lz = max(0, ceil(Int, log2(hypot(w...) / zres)))
    lz ≤ 30 || throw(ArgumentError("zres=$zres is too fine for this box: it needs $lz refinement levels"))
    kz = min(lz + _HEADROOM, 30)
    nslab = maximum(length, knots; init=2) - 1
    kp = min(kz, 30 - ceil(Int, log2(nslab)))
    kp ≥ 1 || throw(ArgumentError("too many knots: $nslab intervals"))
    return Lattice{M}((real(zlo), imag(zlo)), w ./ 2.0^kz, Int32(1) << kz, knots, kp)
end

_knots(p::Tuple{Real,Real}) = _knots(collect(range(Float64(p[1]), Float64(p[2]); length=_PKNOTS)))
function _knots(t::AbstractVector)
    length(t) ≥ 2 && all(isfinite, t) && issorted(t; lt=≤) ||
        throw(ArgumentError("parameter knots must be ≥ 2 finite, strictly increasing values, got $t"))
    return t
end
_knots(x) = throw(ArgumentError(
    "each pdomain entry must be an interval (lo, hi) or a knot vector, got $(repr(x)); one parameter is ((lo, hi),)"))

# Knots an interval `(lo, hi)` is discovered at. A locus living only between two knots
# can be missed, as a pair closer than `minsep` can in z.
const _PKNOTS = 9

@inline function _knotcoord(t, i, kp)
    k, j = divrem(Int(i), 1 << kp)
    return j == 0 ? t[k+1] : t[k+1] + j * (t[k+2] - t[k+1]) / (1 << kp)
end

@inline Base.getindex(lat::Lattice{M}, ix) where {M} = (
    lat.zo[1] + ix[1] * lat.zh[1], lat.zo[2] + ix[2] * lat.zh[2],
    ntuple(d -> _knotcoord(lat.knots[d], ix[d+2], lat.kp), Val(M))...)

const _PENDING, _DARK, _FIRING, _SPLIT = Int8(0), Int8(1), Int8(2), Int8(3)

# A box of the lattice: `base` index of its low corner, `step` its index extent per axis.
# Axes 1–2 are Re z, Im z, the rest parameters. Children split any subset of axes.
mutable struct Cell{N}
    base::NTuple{N,Int32}
    step::NTuple{N,Int32}
    state::Int8
    children::Vector{Cell{N}}
end
Cell(base::NTuple{N,Int32}, step::NTuple{N,Int32}) where {N} = Cell{N}(base, step, _PENDING, Cell{N}[])

nparams(::Cell{N}) where {N} = N - 2

const _ZMASK = 0b11

# Halve the axes whose bit is set in `mask`.
function split!(c::Cell{N}, mask) where {N}
    h = ntuple(d -> (mask >> (d - 1)) & 1 == 1 ? c.step[d] >> 1 : c.step[d], Val(N))
    for k in 0:((1<<N)-1)
        k & ~mask == 0 || continue
        push!(c.children, Cell(ntuple(d -> c.base[d] + Int32((k >> (d - 1)) & 1) * h[d], Val(N)), h))
    end
    c.state = _SPLIT
    return c.children
end

# Corner `i` (bit d → axis d) of a cell.
@inline _corner(c::Cell{N}, i) where {N} =
    ntuple(d -> @inbounds(c.base[d] + Int32((i >> (d - 1)) & 1) * c.step[d]), Val(N))

# Parameter-lattice index of p-corner `j` (bit d → parameter d).
@inline _pcorner(c::Cell{N}, j) where {N} =
    ntuple(d -> @inbounds(c.base[d+2] + Int32((j >> (d - 1)) & 1) * c.step[d+2]), Val(N - 2))

# Adapts a lattice point (x, y, p...) to the curried f(p...)(z::Complex).
struct FunctionWrapper{F}
    f::F
end
(g::FunctionWrapper)(pt) = _setup(g, pt)(complex(pt[1], pt[2]))
_setup(g::FunctionWrapper, pt) = g.f(_ptail(pt)...)
_ptail(pt) = Base.tail(Base.tail(pt))

# Memo of `post ∘ f` on the shared corner lattice to avoid redundant evaluations
struct CornerCache{V,N,F,P,L<:Lattice}
    f::F
    post::P
    lat::L
    d::Dict{NTuple{N,Int32},V}
end
CornerCache(f, post, lat::Lattice{M}) where {M} =
    CornerCache(f, post, lat, Dict{NTuple{M + 2,Int32},Int8}())

_qval(c::CornerCache, ix) = Int8(c.post(c.f(c.lat[ix])))
@inline Base.getindex(c::CornerCache, ix) = get!(() -> _qval(c, ix), c.d, ix)
nevals(c::CornerCache) = length(c.d)

# All 2^N corner quadrants of a cell, gathered once
@inline _quadrants(cache, cell::Cell{N}) where {N} =
    ntuple(i -> cache[_corner(cell, i - 1)], Val(1 << N))

# Signed phase winding of the quadrants at points walked CCW around a z-plane face: +4
# for an enclosed zero, −4 for a pole, `nothing` if any step jumps two quadrants (an
# unresolved ±π step whose direction is ambiguous — the cell must be refined).
@inline _face_winding(qa, qb, qc, qd) = _face_winding((qa, qb, qc, qd))
@inline function _face_winding(q::Tuple)
    w = 0
    n = length(q)
    @inbounds for i in 1:n
        d = mod(q[i%n+1] - q[i], 4)
        d == 2 && return nothing
        w += d == 3 ? -1 : d
    end
    return w
end

_fw_fires(fw) = isnothing(fw) || fw != 0

# Winding of the z-face with low corner `zb`, extent `zs`, at parameter index `pix`.
@inline function _facewinding(cache, zb, zs, pix)
    q(a, b) = cache[(zb[1] + Int32(a) * zs[1], zb[2] + Int32(b) * zs[2], pix...)]
    return _face_winding(q(0, 0), q(1, 0), q(1, 1), q(0, 1))
end

_zbase(c::Cell) = (c.base[1], c.base[2])
_zstep(c::Cell) = (c.step[1], c.step[2])
_facewinding(cache, c::Cell, j) = _facewinding(cache, _zbase(c), _zstep(c), _pcorner(c, j))

# The z-plane faces span axes 1,2; the m = N−2 parameter axes are each pinned to their
# low or high end, giving 2^m faces (one per corner of the parameter subcube). A locus
# piercing the cell winds at least one of them. Corner indices for a face with
# parameter-combo `j`: base = j<<2, then (0,1,3,2) = CCW in (axis1, axis2).
@inline function _fires(cache, cell)
    q = _quadrants(cache, cell)
    for j in 0:((1<<nparams(cell))-1)
        b = j << 2
        _fw_fires(@inbounds _face_winding(q[b+1], q[b+2], q[b+4], q[b+3])) && return true
    end
    # Grazing guard: a locus with steep |dz/dp| can cross the cell through its
    # parameter side faces, leaving every z-face winding clean. Corners spanning
    # ≥3 quadrants still betray it — the Re f = 0 and Im f = 0 sheets both cross
    # the cell, and they intersect only on loci. Costs no extra evaluations (the
    # face loop above already touched every corner), only extra refinement.
    m = 0
    for x in q
        m |= 1 << x
    end
    return count_ones(m) ≥ 3
end

_zcenter(lat::Lattice, c::Cell) = complex(
    lat.zo[1] + (c.base[1] + c.step[1] / 2) * lat.zh[1],
    lat.zo[2] + (c.base[2] + c.step[2] / 2) * lat.zh[2])

# Coordinates are linear within a knot interval, which a cell never straddles.
_pcenter(lat::Lattice{M}, c::Cell) where {M} = ntuple(Val(M)) do d
    t = lat.knots[d]
    (_knotcoord(t, c.base[d+2], lat.kp) + _knotcoord(t, c.base[d+2] + c.step[d+2], lat.kp)) / 2
end

_pcoords(lat::Lattice{M}, pix) where {M} = ntuple(d -> _knotcoord(lat.knots[d], pix[d], lat.kp), Val(M))

@inline _zdiag(lat::Lattice, c::Cell) = hypot(c.step[1] * lat.zh[1], c.step[2] * lat.zh[2])

# Reserve-then-fill: `_reserve!` claims a key with `_UNSET` before the threaded
# fill, so a corner shared by several cells is enqueued exactly once.
const _UNSET = Int8(0)

function _reserve!(pts, dict, ix)
    n = length(dict)
    get!(dict, ix, _UNSET)
    length(dict) > n && push!(pts, ix)
end

# Evaluate a batch of corners, calling `f(p...)` once per distinct parameter point so its setup is shared by every z there.
function _parfill!(c::CornerCache{V,N}, pts) where {V,N}
    slot = Vector{Int}(undef, length(pts))
    rep = eltype(pts)[]
    at = Dict{NTuple{N - 2,Int32},Int}()
    for (i, ix) in pairs(pts)
        slot[i] = get!(at, _ptail(ix), length(rep) + 1)
        slot[i] > length(rep) && push!(rep, ix)
    end
    gs = Vector{Any}(undef, length(rep))
    Threads.@threads for j in eachindex(gs, rep)
        gs[j] = _setup(c.f, c.lat[rep[j]])
    end
    vals = Vector{V}(undef, length(pts))
    _fill!(vals, c, identity.(gs), slot, pts)    # narrowed eltype: concrete calls below
    for (pt, v) in zip(pts, vals)
        c.d[pt] = v
    end
end

function _fill!(vals, c, gs, slot, pts)
    Threads.@threads for i in eachindex(vals, pts)
        pt = c.lat[pts[i]]
        vals[i] = c.post(gs[slot[i]](complex(pt[1], pt[2])))
    end
end

# The refinement state `link` continues from: escalation refines further.
struct Tree{M,N,C<:CornerCache,A<:AbstractArray{Cell{N},M}}
    cache::C
    lat::Lattice{M}
    slabs::A        # one root cell per knot interval (tensor over parameters)
    minsep::Float64
end

# Parameter faces of `c` that lose a locus: `(j2, d)` for each face at p-corner `j` that
# fires while its neighbour `j2` along parameter `d` does not. The locus left the cell's
# z-column over that parameter step; `_connect!` looks for it next door.
function _lost(cache, c::Cell{N}) where {N}
    M = N - 2
    out = Tuple{Int,Int}[]
    M == 0 && return out
    fires = ntuple(j -> _fw_fires(_facewinding(cache, c, j - 1)), Val(1 << M))
    for j in 0:((1<<M)-1), d in 1:M
        j2 = j ⊻ (1 << (d - 1))
        fires[j+1] && !fires[j2+1] && push!(out, (j2, d))
    end
    return out
end

_grazes(cache, c::Cell{N}) where {N} = N > 2 && !any(j -> _fw_fires(_facewinding(cache, c, j)), 0:((1<<(N-2))-1))

_inbox(nb, zs, zn) = 0 ≤ nb[1] && nb[1] + zs[1] ≤ zn && 0 ≤ nb[2] && nb[2] + zs[2] ≤ zn

_ring(c::Cell) = ((c.base[1] + Int32(dx) * c.step[1], c.base[2] + Int32(dy) * c.step[2])
                  for dx in -1:1 for dy in -1:1 if (dx, dy) != (0, 0))

function _prefetch!(tree::Tree, cells)
    (; cache, lat, minsep) = tree
    pts = keytype(cache.d)[]
    for c in cells
        _zdiag(lat, c) > minsep && continue
        for i in 0:((1<<length(c.base))-1)
            _reserve!(pts, cache.d, _corner(c, i))
        end
    end
    _parfill!(cache, pts)
end

function _prefetch_ring!(tree::Tree, cells)
    (; cache, lat) = tree
    pts = keytype(cache.d)[]
    for c in cells, (j2, _) in _lost(cache, c)
        pix, zs = _pcorner(c, j2), _zstep(c)
        for nb in _ring(c)
            _inbox(nb, zs, lat.zn) || continue
            for a in 0:1, b in 0:1
                _reserve!(pts, cache.d, (nb[1] + Int32(a) * zs[1], nb[2] + Int32(b) * zs[2], pix...))
            end
        end
    end
    _parfill!(cache, pts)
end

# Parameter axes (as a split mask) to refine `c` along: those along which `c` loses a
# locus that no z-neighbour picks up at the same parameter step — the locus moved more
# than a cell, and adjacency could no longer link it. Neighbours that do pick it up are
# pushed to `handoffs`; leaving through the box edge is not a loss.
#
# A cell firing only on the grazing guard (no z-face fires) has a locus cross it between
# faces, or pass nearby: all parameters split, so the band of such cells along a fast
# locus narrows with z.
function _connect!(handoffs, cache, lat, c::Cell{N}) where {N}
    M = N - 2
    M == 0 && return 0
    zs = _zstep(c)
    _grazes(cache, c) && return sum(d -> c.step[d+2] > 1 ? 1 << (d + 1) : 0, 1:M)
    mask = 0
    for (j2, d) in _lost(cache, c)
        ok = false
        pix = _pcorner(c, j2)
        for nb in _ring(c)
            if !_inbox(nb, zs, lat.zn)
                ok = true
            elseif _fw_fires(_facewinding(cache, nb, zs, pix))
                ok = true
                push!(handoffs, (c, nb, j2))
            end
        end
        ok || c.step[d+2] ≤ 1 || (mask |= 1 << (d + 1))
    end
    return mask
end

_contains(c::Cell{N}, q) where {N} = all(d -> 2Int(c.base[d]) < q[d] < 2(Int(c.base[d]) + c.step[d]), 1:N)

# Leaf containing `q` (doubled index coordinates, odd so never on a boundary)
function _locate(tree::Tree{M}, q) where {M}
    c = tree.slabs[ntuple(d -> q[d+2] >> (tree.lat.kp + 1) + 1, Val(M))...]
    while !isempty(c.children)
        c = c.children[findfirst(ch -> _contains(ch, q), c.children)]
    end
    return c
end

# Make sure a cell at `c`'s z-level has the face at p-corner `j2` where neighbour column
# `nb` took over the locus: refine the leaf there toward it — a dark one, or a firing one
# that only grazes the locus within a coarser parameter step. Firing leaves of another
# level belong to other groups (escalation) and are left alone. Returns `false` if the
# leaf is still pending, to retry once it is classified.
function _handoff!(frontier, tree::Tree, (c, nb, j2)::Tuple{Cell{N},Any,Int}) where {N}
    q = ntuple(Val(N)) do d
        d ≤ 2 && return 2Int(nb[d]) + c.step[d]
        bit = (j2 >> (d - 3)) & 1
        return 2(Int(c.base[d]) + bit * c.step[d]) + (bit == 1 ? -1 : 1)
    end
    leaf = _locate(tree, q)
    leaf.state == _PENDING && return false
    leaf.state == _DARK || (leaf.state == _FIRING && leaf.step[1] == c.step[1]) || return true
    while true
        mask = leaf.step[1] > c.step[1] ? _ZMASK : 0
        for d in 3:N
            (leaf.base[d] < c.base[d] || leaf.base[d] + leaf.step[d] > c.base[d] + c.step[d]) &&
                (mask |= 1 << (d - 1))
        end
        mask == 0 && break
        next = leaf
        for k in split!(leaf, mask)
            _contains(k, q) ? (next = k) : push!(frontier, k)
        end
        leaf = next
    end
    leaf.state == _PENDING && push!(frontier, leaf)
    return true
end

# Refine `frontier` to z-diagonal `zres` and return the firing cells there. z splits where
# a cell fires (or is coarser than `minsep`); a parameter splits, at any level, only where
# a locus moves more than a cell across the cell's parameter step (`_connect!`), so the
# parameter resolution follows the loci's speed rather than `zres`. Breadth-first: each
# round's evaluations are known up front and batched.
function _refine!(tree::Tree{M,N}, frontier::Vector{Cell{N}}, zres) where {M,N}
    (; cache, lat, minsep) = tree
    out = Cell{N}[]
    next = Cell{N}[]
    firing = Cell{N}[]
    handoffs = Tuple{Cell{N},NTuple{2,Int32},Int}[]
    while !isempty(frontier)
        while !isempty(frontier)
            _prefetch!(tree, frontier)
            empty!(firing)
            for c in frontier
                zd = _zdiag(lat, c)
                if zd > zres && zd > minsep
                    append!(next, split!(c, _ZMASK))
                elseif !_fires(cache, c)
                    c.state = _DARK
                else
                    push!(firing, c)
                end
            end
            _prefetch_ring!(tree, firing)
            for c in firing
                coarse = _zdiag(lat, c) > zres
                # At `zres` the handoffs below already chain the face-holding cells
                # along every locus; a grazed cell would only add an unlabelled island
                # beside a fast one. Splitting it instead would chase a locus clipping
                # its corner through ever-thinner slices.
                if !coarse && _grazes(cache, c)
                    c.state = _DARK
                    continue
                end
                n = length(handoffs)
                mask = _connect!(handoffs, cache, lat, c) | (coarse ? _ZMASK : 0)
                if mask == 0
                    c.state = _FIRING
                    push!(out, c)
                else
                    resize!(handoffs, n)
                    append!(next, split!(c, mask))
                end
            end
            frontier, next = next, frontier
            empty!(next)
        end
        filter!(h -> !_handoff!(frontier, tree, h), handoffs)
    end
    return filter!(c -> c.state == _FIRING, out)
end

"""
    Sample{M}

A point on a locus: parameters `p` where the winding was measured, and the center `z`
of the z-face it was found in.
"""
struct Sample{M}
    z::ComplexF64
    p::NTuple{M,Float64}
end

Base.length(::Sample{M}) where {M} = M + 1
Base.iterate(s::Sample, i=1) = i > length(s) ? nothing : (i == 1 ? s.z : s.p[i-1], i + 1)

# One sample per firing z-face of `cells`, shared faces once; cell centers when none fires
# (no parameters, or a locus pinned to the lattice that every face reads ambiguous).
function _samples(tree::Tree{M}, cells) where {M}
    (; cache, lat) = tree
    out = Sample{M}[]
    if M > 0
        seen = Set{NTuple{M + 2,Int32}}()
        for c in cells, j in 0:((1<<M)-1)
            _fw_fires(_facewinding(cache, c, j)) || continue
            pix = _pcorner(c, j)
            key = (_zbase(c)..., pix...)
            key in seen && continue
            push!(seen, key)
            push!(out, Sample{M}(_zcenter(lat, c), _pcoords(lat, pix)))
        end
    end
    isempty(out) && (out = [Sample{M}(_zcenter(lat, c), _pcenter(lat, c)) for c in cells])
    return sort!(out; by=s -> s.p)
end

"""
    Branch{M} <: AbstractVector{Sample{M}}

One locus `z*(p)` over `M` parameters in a [`Survey`](@ref), as its [`Sample`](@ref)s;
described in the README, [`MeromorphicLoci`](@ref).
"""
struct Branch{M,C} <: AbstractVector{Sample{M}}
    samples::Vector{Sample{M}}
    cells::Vector{C}
    winding::Int
end
Branch(tree::Tree{M}, cells::Vector{C}, winding) where {M,C} =
    Branch{M,C}(_samples(tree, cells), cells, winding)

Base.size(b::Branch) = size(b.samples)
Base.getindex(b::Branch, i::Int) = b.samples[i]

"""
    winding(branch)::Int

Net zeros − poles the branch encloses, by phase winding.
"""
winding(b::Branch) = b.winding

_kind(b::Branch) = b.winding > 0 ? "Zero" : b.winding < 0 ? "Pole" : "Indeterminate"

Base.summary(io::IO, b::Branch{M}) where {M} =
    print(io, "Branch{", M, "}(", winding(b), " ", _kind(b), ", ", length(b), " samples)")
Base.show(io::IO, b::Branch) = summary(io, b)

# The refined tree and its candidate cells; `link` evaluates `f` further and splits
# ambiguous cells.
struct Scan{M,C,T<:Tree{M}} <: AbstractVector{Sample{M}}
    samples::Vector{Sample{M}}
    cells::Vector{C}
    tree::T
    zres::Float64
end

Base.size(s::Scan) = size(s.samples)
Base.getindex(s::Scan, i::Int) = s.samples[i]

nevals(s::Scan) = nevals(s.tree.cache)

Base.summary(io::IO, s::Scan{M}) where {M} =
    print(io, "Scan{", M, "}(", length(s.cells), " candidate cells, ", nevals(s), " evals)")
Base.show(io::IO, s::Scan) = summary(io, s)

_zcorners(a, b) = (complex(min(real(a), real(b)), min(imag(a), imag(b))),
    complex(max(real(a), real(b)), max(imag(a), imag(b))))

# Refine over the curried `f(p...)(z)` and collect the candidate cells.
# `box = ((zlo, zhi), p...)`, each `p` an interval `(lo, hi)` or knot vector.
function scan(f, box; zres, minsep=nothing)
    zbox = _zcorners(box[1]...)
    minsep = Float64(something(minsep, abs(zbox[2] - zbox[1]) / 8))
    knots = map(_knots, Base.tail(box))
    M = length(knots)
    lat = Lattice(zbox, knots, zres)
    cache = CornerCache(FunctionWrapper(f), quadrant, lat)
    ps = Int32(1) << lat.kp
    slabs = map(CartesianIndices(map(t -> length(t) - 1, knots))) do k
        Cell((Int32(0), Int32(0), ntuple(d -> Int32(k[d] - 1) * ps, Val(M))...),
            (lat.zn, lat.zn, ntuple(_ -> ps, Val(M))...))
    end
    tree = Tree(cache, lat, slabs, minsep)
    cells = _refine!(tree, vec(collect(slabs)), Float64(zres))
    return Scan(_samples(tree, cells), cells, tree, Float64(zres))
end
