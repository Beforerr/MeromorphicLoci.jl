"""
    Survey{M} <: AbstractVector{Branch}

The [`Branch`](@ref)es a [`survey`](@ref) over `M` parameters found.
"""
struct Survey{M,B} <: AbstractVector{B}
    branches::Vector{B}
    nevals::Int
end


Survey{M}(branches::Vector{B}, nevals) where {M,B} = Survey{M,B}(branches, nevals)

"""
    link(scan; nreps=5, zoom=4) -> Survey

Group a [`Scan`](@ref)'s candidate cells into connected [`Branch`](@ref)es — 
one per locus `z*(p)` — classified by phase winding.
"""
function link(s::Scan{M}; nreps=5, zoom=4) where {M}
    C = eltype(s.cells)
    branches = Branch{M,C}[]
    extra = 0
    inside = _insider(s)
    for g in (isempty(s.cells) ? Vector{C}[] : _group(s.cells, s.cache.lat))
        subs, ne = _resolve(s.cache, s.kcache, g, 2 * s.zres, inside; nreps, zoom)
        extra += ne
        for (cs, w) in subs
            push!(branches, Branch{M}(cs, w))
        end
    end
    return Survey{M}(branches, nevals(s) + extra)
end

# Connected components of a nonempty set of same-size cells (union-find). Full
# corner/edge/face adjacency is exactly "origin index differs by one stride on
# every axis" — O(n·3^N) through a lookup table.
function _group(cells::Vector{C}, lat) where {C<:Cell{<:Any,N}} where {N}
    n = length(cells)
    _, stride = _cellindex(lat, cells[1])
    at = Dict{NTuple{N,Int32},Int}(_cellindex(lat, cells[i])[1] => i for i in 1:n)

    parent = collect(1:n)
    function find(i)
        while parent[i] != i
            parent[i] = parent[parent[i]]
            i = parent[i]
        end
        return i
    end
    for i in 1:n
        k, _ = _cellindex(lat, cells[i])
        for off in CartesianIndices(ntuple(_ -> -1:1, Val(N)))
            j = get(at, k .+ Tuple(off) .* stride, 0)
            j == 0 && continue
            parent[find(i)] = find(j)
        end
    end
    groups = Dict{Int,Vector{C}}()
    for i in 1:n
        push!(get!(() -> C[], groups, find(i)), cells[i])
    end
    return collect(values(groups))
end

function _insider(s::Scan)
    o, w = s.tree.boundary.origin, s.tree.boundary.widths
    keep = isnothing(s.kcache) ? nothing : _unwrap(s.kcache.f)
    return function (z, p)
        o[1] ≤ real(z) ≤ o[1] + w[1] && o[2] ≤ imag(z) ≤ o[2] + w[2] || return false
        return isnothing(keep) || keep(z, p...)
    end
end

# Label one connected group of same-size cells; returns `cells => winding` pairs
# and the evaluations spent.
#
# Cheap paths first. (1) The cached fully-enclosing faces — a 4-corner walk
# resolves at most one turn, so multiplicity > 1 loci reads ambiguous; trusted only
# when every enclosure agrees in sign. (2) Circle votes (`r = 2zres`), trusted only when 
# every usable vote agrees in sign (covers multiplicity and corner-pinned loci).
#
# Otherwise the label is ambiguous — a fused zero/pole pair, or nothing usable.
# Escalate: split the cells one level, regroup the children that still fire, 
# recurse with the vote radius halved and one less `zoom`. 
# A zero-pole pair a few sub-cells apart falls into cleanly labeled sub-branches; 
# one that stays connected to the floor (genuinely fused, or loci crossing at some p) 
# keeps winding 0. Cost is local to the ambiguous group.
function _resolve(cache, kcache, cells, r, inside; nreps, zoom)
    npos = nneg = 0
    for c in cells, fw in _enclosing_faces(cache, c)
        fw > 0 && (npos += 1)
        fw < 0 && (nneg += 1)
    end
    (npos == 0) ⊻ (nneg == 0) && return [cells => (npos > 0 ? 1 : -1)], 0
    ws, ne = _votes(cache, cells, r, inside; nreps)
    if !isempty(ws) && (first(ws) > 0 || last(ws) < 0)
        return [cells => ws[(length(ws)+1)÷2]], ne
    end
    kids = zoom > 0 ? _deepen!(cache, kcache, cells) : empty(cells)
    isempty(kids) && return [cells => 0], ne
    out = Pair{typeof(cells),Int}[]
    for g in _group(kids, cache.lat)
        subs, n = _resolve(cache, kcache, g, r / 2, inside; nreps, zoom=zoom - 1)
        ne += n
        append!(out, subs)
    end
    # Still one connected piece — same locus set, possibly a sharper label; keep
    # the original coarser cells so samples stay at the requested resolution.
    length(out) == 1 && return [cells => out[1].second], ne
    return out, ne
end

# Sign of each **fully-enclosing** (±4) z-plane face winding of a cell.
# Partial/ambiguous faces (a singularity sitting on a
# face boundary) give 0, so a locus pinned to the corner lattice — which no single
# face encloses — reports nothing and leaves the vote to the circle fallback.
@inline function _enclosing_faces(cache, cell)
    q = _quadrants(cache, cell)
    return ntuple(Val(1 << nparams(cell))) do j
        b = (j - 1) << 2
        fw = @inbounds _face_winding(q[b+1], q[b+2], q[b+4], q[b+3])
        fw == 4 ? 1 : fw == -4 ? -1 : 0
    end
end

# Split each cell one level and keep the children that still fire; 
# empty at the lattice floor, where child corners have no index left
function _deepen!(cache, kcache, cells)
    _, step = _cellindex(cache.lat, cells[1])
    step[1] ≤ 1 && return empty(cells)
    kids = empty(cells)
    for c in cells
        isleaf(c) && split!(c)
        for ch in children(c)
            _kept(kcache, ch) && _fires(cache, ch) && push!(kids, ch)
        end
    end
    return kids
end

# Sorted usable circle windings of `nreps` representatives spread along the group,
# and the evaluations they cost. `r` must be wide enough that a merged opposite
# partner lands *inside* the circle and cancels the vote, rather than hiding
# outside it while every pick lands on one side — hence 2×`zres` at the top level,
# not 2× the cell diagonal. When `f` *at* a corner lattice point has a singularity
# (0 or NaN), no face walk can be trusted to close, so this is also the only
# usable signal there. A representative whose circle leaves `inside` is skipped.
function _votes(cache, cells, r, inside; nreps)
    f = _unwrap(cache.f)
    n = length(cells)
    sorted = sort(cells; by=c -> (_pcenter(c), real(_zcenter(c))))
    picks = unique(round.(Int, range(1, n; length=min(nreps, n))))
    ws = Int[]
    ne = 0
    for i in picks
        c = sorted[i]
        z0, p = _zcenter(c), _pcenter(c)
        w, nw = _circle_winding(z0, r) do z
            inside(z, p) ? f(z, p...) : NaN
        end
        ne += nw
        isnothing(w) || push!(ws, w)
    end
    return sort!(ws), ne
end


nevals(s::Survey) = s.nevals
Base.size(s::Survey) = size(s.branches)
Base.getindex(s::Survey, i::Int) = s.branches[i]
Base.IndexStyle(::Type{<:Survey}) = IndexLinear()

function Base.summary(io::IO, s::Survey{M}) where {M}
    nz = count(b -> b.winding > 0, s.branches)
    np = count(b -> b.winding < 0, s.branches)
    print(io, "Survey{", M, "} with ", length(s.branches), " branches (",
        nz, " zero, ", np, " pole")
    (i = count(b -> b.winding == 0, s.branches)) > 0 && print(io, ", ", i, " indeterminate")
    return print(io, "), ", s.nevals, " evals")
end