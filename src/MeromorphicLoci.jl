module MeromorphicLoci

using RegionTrees: Cell, split!, isleaf, children
using StaticArrays

export survey, scan, link, Survey, Scan, Branch, Sample, winding

include("winding.jl")

"""
    survey(f, box; zres, minsep=abs(zmax-zmin)/8, keep=nothing, nreps=5, zoom=4) -> Survey

Find every zero and pole locus of `f(z::Complex, p...)` over `box`: [`scan`](@ref)
for candidate cells, then [`link`](@ref) them into classified [`Branch`](@ref)es.
"""
survey(f, box; nreps=5, zoom=4, kw...) = link(scan(f, box; kw...); nreps, zoom)


# Argument-principle quadrant, CCW half-open: 1:[0,π/2) 2:[π/2,π] 3:(-π,-π/2) 4:[-π/2,0).
# f = 0 or NaN bins arbitrarily (q2/q4) — pinned singularities are
# detected by neighboring faces and classified by the circle fallback
@inline function quadrant(z)
    re, im = reim(z)
    return im ≥ 0 ? (re > 0 ? 1 : 2) : (re < 0 ? 3 : 4)
end

# RegionTrees halves every axis each level, and stops refinement at the first level whose
# z-diagonal fits `zres`, so every corner lands exactly on the lattice of that deepest level.
# Addressing corners by integer lattice index not raw coordinates so two cells meeting
# at a corner only need to agree to within rounding instead of bit-for-bit.
struct Lattice{N,T}
    o0::SVector{N,T}
    h::SVector{N,T}
    invh::SVector{N,T}
end

# Headroom below the `zres` level so `link`'s escalation can index sub-`zres` corners
const _HEADROOM = 8

function Lattice(origin::SVector{N,T}, widths, zres) where {N,T}
    zres > 0 || throw(ArgumentError("zres must be positive, got $zres"))
    lmax = max(0, ceil(Int, log2(hypot(widths[1], widths[2]) / zres)))
    lmax ≤ 30 || throw(ArgumentError("zres=$zres is too fine for this box: it needs $lmax refinement levels"))
    h = widths ./ 2.0^min(lmax + _HEADROOM, 30)
    return Lattice{N,T}(origin, h, 1 ./ h)
end

@inline Base.getindex(lat::Lattice{N}, ix) where {N} =
    ntuple(d -> @inbounds(lat.o0[d] + ix[d] * lat.h[d]), Val(N))

# Lattice index of a cell's origin, and the index stride between its corners.
@inline function _cellindex(lat::Lattice{N}, cell::Cell{<:Any,N}) where {N}
    o, w = cell.boundary.origin, cell.boundary.widths
    base = ntuple(d -> @inbounds(round(Int32, (o[d] - lat.o0[d]) * lat.invh[d])), Val(N))
    step = ntuple(d -> @inbounds(round(Int32, w[d] * lat.invh[d])), Val(N))
    return base, step
end

# Corner `i` (bit d → axis d) of the cell at `base`.
@inline _cornerindex(base::NTuple{N}, step, i) where {N} =
    ntuple(d -> @inbounds(base[d] + Int32((i >> (d - 1)) & 1) * step[d]), Val(N))

# Adapts a lattice point (x, y, p...) to the user's f(z::Complex, p...) signature.
struct FunctionWrapper{F}
    f::F
end
(g::FunctionWrapper)(pt) = g.f(complex(pt[1], pt[2]), Base.tail(Base.tail(pt))...)
_unwrap(g::FunctionWrapper) = g.f

# Memo of `post ∘ f` on the shared corner lattice to avoid redundant evaluations
struct CornerCache{V,N,F,P,L<:Lattice{N}}
    f::F
    post::P
    lat::L
    d::Dict{NTuple{N,Int32},V}
end
CornerCache(f, post, lat::Lattice{N}) where {N} =
    CornerCache(f, post, lat, Dict{NTuple{N,Int32},Int8}())

_qval(c::CornerCache, ix) = Int8(c.post(c.f(c.lat[ix])))
@inline Base.getindex(c::CornerCache, ix) = get!(() -> _qval(c, ix), c.d, ix)
nevals(c::CornerCache) = length(c.d)

# All 2^N corner quadrants of a cell, gathered once
@inline function _quadrants(cache, cell::Cell{<:Any,N}) where {N}
    base, step = _cellindex(cache.lat, cell)
    return ntuple(i -> cache[_cornerindex(base, step, i - 1)], Val(1 << N))
end

# Signed phase winding of 4 quadrants walked CCW around a z-plane face: +4 for an
# enclosed zero, −4 for a pole, `nothing` if any edge jumps two quadrants (an
# unresolved ±π step whose direction is ambiguous — the cell must be refined).
@inline function _face_winding(qa, qb, qc, qd)
    w = 0
    q = (qa, qb, qc, qd)
    @inbounds for i in 1:4
        d = mod(q[i%4+1] - q[i], 4)
        d == 2 && return nothing
        w += d == 3 ? -1 : d
    end
    return w
end

# The z-plane 2-faces span axes 1,2; the remaining m = N−2 axes are each pinned to
# 0/1, giving 2^m faces (one per corner of the parameter subcube). A singularity
# locus piercing the cell winds at least one of them. Corner indices for a face
# with parameter-combo `j`: base = j<<2, then (0,1,3,2) = CCW in (axis1, axis2).
@inline function _fires(cache, cell)
    q = _quadrants(cache, cell)
    for j in 0:((1<<nparams(cell))-1)
        b = j << 2
        fw = @inbounds _face_winding(q[b+1], q[b+2], q[b+4], q[b+3])
        (isnothing(fw) || fw != 0) && return true
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

# z-plane center and parameter coordinates of a cell.
_zcenter(c) = complex(
    c.boundary.origin[1] + c.boundary.widths[1] / 2,
    c.boundary.origin[2] + c.boundary.widths[2] / 2
)

nparams(::Cell{<:Any,N}) where {N} = N - 2

_pcenter(c::Cell) = ntuple(
    d -> c.boundary.origin[d+2] + c.boundary.widths[d+2] / 2, Val(nparams(c))
)

# A cell survives only if `keep` holds somewhere on it. Testing every corner so
#  a cell straddling the boundary of the kept region is always refined.
_kept(::Nothing, cell) = true      # no veto requested
function _kept(c::CornerCache{<:Any,N}, cell) where {N}
    base, step = _cellindex(c.lat, cell)
    return any(i -> c[_cornerindex(base, step, i)] > 0, 0:((1<<N)-1))
end

@inline _zdiag(cell) = (w=cell.boundary.widths; hypot(w[1], w[2]))
@inline _atres(cell, zres) = _zdiag(cell) ≤ zres

function _needs_refinement(cell, cache, kcache, zres, minsep)
    _atres(cell, zres) && return false
    _kept(kcache, cell) || return false
    # A base mesh fine enough that each cell brackets ≤1 branch — otherwise two
    # loci in one coarse cell alias the winding to zero and it never fires.
    _zdiag(cell) > minsep && return true
    return _fires(cache, cell)
end

# Reserve-then-fill: `_corners!` claims a key with `_UNSET` before the threaded
# fill, so a corner shared by several frontier cells is enqueued exactly once.
const _UNSET = Int8(0)

function _corners!(pts, dict, lat, cell::Cell{<:Any,N}) where {N}
    base, step = _cellindex(lat, cell)
    for i in 0:((1<<N)-1)
        ix = _cornerindex(base, step, i)
        n = length(dict)
        get!(dict, ix, _UNSET)
        length(dict) > n && push!(pts, ix)
    end
end

function _parfill!(g, dict, pts, vals)
    resize!(vals, length(pts))
    Threads.@threads for i in eachindex(vals, pts)
        vals[i] = g(pts[i])
    end
    for (pt, v) in zip(pts, vals)
        dict[pt] = v
    end
end

# Fill both caches with every corner the frontier's tests will read — including
# at-`zres` cells, whose corners the final `_candidates!` scan reads.
function _prefetch!(frontier, cache, kcache, minsep, pts, vals)
    if !isnothing(kcache)
        empty!(pts)
        for c in frontier
            _corners!(pts, kcache.d, kcache.lat, c)
        end
        _parfill!(ix -> _qval(kcache, ix), kcache.d, pts, vals)
    end
    empty!(pts)
    for c in frontier
        (_zdiag(c) > minsep || !_kept(kcache, c)) && continue
        _corners!(pts, cache.d, cache.lat, c)
    end
    _parfill!(ix -> _qval(cache, ix), cache.d, pts, vals)
end

# Breadth-first refinement: each level's `f`/`keep` evaluations are known up
# front, so `_prefetch!` batches them; the per-cell decisions then replay from cache.
function _refine!(root, cache, kcache, zres, minsep)
    frontier = [root]
    next = empty(frontier)
    pts = keytype(cache.d)[]
    vals = valtype(cache.d)[]
    while !isempty(frontier)
        _prefetch!(frontier, cache, kcache, minsep, pts, vals)
        for cell in frontier
            _needs_refinement(cell, cache, kcache, zres, minsep) || continue
            split!(cell)
            append!(next, children(cell))
        end
        frontier, next = next, frontier
        empty!(next)
    end
end

# Candidate leaves: fully refined (coarser leaves were already vetoed by `_kept`
# or `_fires` during refinement) and still firing at the final resolution.
function _candidates!(out, cell, cache, kcache, zres)
    if isleaf(cell)
        _atres(cell, zres) && _kept(kcache, cell) && _fires(cache, cell) && push!(out, cell)
    else
        for child in children(cell)
            _candidates!(out, child, cache, kcache, zres)
        end
    end
    return out
end

"""
    Sample{M}

A candidate-cell center: parameters `p`, z-plane center `z`.
"""
struct Sample{M}
    z::ComplexF64
    p::NTuple{M,Float64}
end

Base.length(::Sample{M}) where {M} = M + 1
Base.iterate(s::Sample, i=1) = i > length(s) ? nothing : (i == 1 ? s.z : s.p[i-1], i + 1)

_samples(cells, ::Val{M}) where {M} =
    sort!([Sample{M}(_zcenter(c), _pcenter(c)) for c in cells]; by=s -> s.p)

"""
    Branch{M} <: AbstractVector{Sample{M}}

One connected component of candidate cells — a single locus `z*(p)` over `M`
parameters, as its [`Sample`](@ref)s sorted by parameter. [`winding`](@ref)
classifies it zero-vs-pole by sign: `+` zero, `−` pole, `0` fused/indeterminate.
"""
struct Branch{M,C} <: AbstractVector{Sample{M}}
    samples::Vector{Sample{M}}
    cells::Vector{C}
    winding::Int
end
Branch{M}(cells::Vector{C}, winding) where {M,C} =
    Branch{M,C}(_samples(cells, Val(M)), cells, winding)

Base.size(b::Branch) = size(b.samples)
Base.getindex(b::Branch, i::Int) = b.samples[i]

"""
    winding(branch)::Int

Net zeros − poles the branch encloses, by phase winding: `±m` for a
multiplicity-`m` zero or pole, `0` for a fused pair.
"""
winding(b::Branch) = b.winding

# Display label from the winding sign — internal, no longer a public enum.
_kind(b::Branch) = b.winding > 0 ? "Zero" : b.winding < 0 ? "Pole" : "Indeterminate"

Base.summary(io::IO, b::Branch{M}) where {M} =
    print(io, "Branch{", M, "}(", winding(b), " ", _kind(b), ", ", length(b), " samples)")
Base.show(io::IO, b::Branch) = summary(io, b)

"""
    Scan{M} <: AbstractVector{Sample{M}}

Result of [`scan`](@ref): every candidate cell of the refined tree as an [`Sample`](@ref).

`cells` are the candidate RegionTrees leaves, `tree` the raw root.
"""
struct Scan{M,C,A,K,R} <: AbstractVector{Sample{M}}
    samples::Vector{Sample{M}}
    cells::Vector{C}
    cache::A
    kcache::K
    tree::R
    zres::Float64
end
Scan{M}(samples, cells::Vector{C}, cache::A, kcache::K, tree::R, zres) where {M,C,A,K,R} =
    Scan{M,C,A,K,R}(samples, cells, cache, kcache, tree, zres)

Base.size(s::Scan) = size(s.samples)
Base.getindex(s::Scan, i::Int) = s.samples[i]

nevals(s::Scan) = nevals(s.cache)

Base.summary(io::IO, s::Scan{M}) where {M} =
    print(io, "Scan{", M, "}(", length(s.cells), " candidate cells, ", nevals(s), " evals)")
Base.show(io::IO, s::Scan) = summary(io, s)


"""
    scan(f, box; zres, minsep=abs(zmax-zmin)/8, keep=nothing)::Scan

Refine a region tree over `f(z::Complex, p...)` and return the candidate cells.

`box` gives one `(lo, hi)` per argument of `f`; the leading z entry is opposite
complex corners of the z-rectangle. Cells refine until their z-diagonal is ≤ `zres`
near every locus, where the winding criterion fires on zeros and poles alike.

- `zres` — z-plane resolution of the reported samples.
- `minsep` — detection guarantee: the cell size the tree refines to *before*
  winding may reject a cell, so loci pairs farther apart than `minsep` cannot
  alias to winding 0 inside one cell and vanish.
- `keep(z, p...)` — domain veto; cells with no kept corner are never evaluated,
  refined, or reported.

See the README's *Knobs* section for how `zres` and `minsep` trade off.
"""
function scan(f, box; zres, minsep=abs(box[1][2] - box[1][1]) / 8, keep=nothing)
    zlo, zhi = box[1]                  # complex corners of the z-rectangle
    ps = Base.tail(box)               # one (lo, hi) real interval per parameter
    N = length(box) + 1
    origin = SVector{N,Float64}(real(zlo), imag(zlo), map(first, ps)...)
    widths = SVector{N,Float64}(real(zhi) - real(zlo), imag(zhi) - imag(zlo),
        map(p -> p[2] - p[1], ps)...)
    root = Cell(origin, widths)
    lat = Lattice(origin, widths, zres)
    kcache = isnothing(keep) ?
             nothing : CornerCache(FunctionWrapper(keep), k -> Int8(k ? 1 : -1), lat)
    cache = CornerCache(FunctionWrapper(f), quadrant, lat)
    return _scan(root, cache, kcache, Float64(zres), Float64(minsep))
end

function _scan(root::Cell{<:Any,N}, cache, kcache, zres, minsep) where {N}
    _refine!(root, cache, kcache, zres, minsep)
    leaves = _candidates!(typeof(root)[], root, cache, kcache, zres)
    return Scan{N - 2}(_samples(leaves, Val(N - 2)), leaves, cache, kcache, root, zres)
end

include("link.jl")

end
