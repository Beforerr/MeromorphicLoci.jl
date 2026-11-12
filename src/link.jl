"""
    Survey{M} <: AbstractVector{Branch}

The [`Branch`](@ref)es an [`ArgumentPrinciple`](@ref) search over `M` parameters
found.
"""
struct Survey{M,B} <: AbstractVector{B}
    branches::Vector{B}
    nevals::Int
end


Survey{M}(branches::Vector{B}, nevals) where {M,B} = Survey{M,B}(branches, nevals)

# Group a `Scan`'s candidate cells into connected branches — one per locus
# `z*(p)` — classified by phase winding.
function link(s::Scan{M}; zoom=4) where {M}
    C = eltype(s.cells)
    branches = Branch{M,C}[]
    extra = 0
    for g in (isempty(s.cells) ? Vector{C}[] : _group(s.cells))
        subs, ne = _resolve(s.tree, g; zoom)
        extra += ne
        for (cs, w) in subs
            push!(branches, Branch(s.tree, cs, w))
        end
    end
    return Survey{M}(branches, nevals(s) + extra)
end

# Connected components (union-find) of a nonempty set of cells sharing one z size;
# their parameter extents vary. Adjacent: z columns within one stride on both z axes,
# parameter ranges touching or overlapping on every parameter axis.
function _group(cells::Vector{C}) where {N,C<:Cell{N}}
    n = length(cells)
    zs = _zstep(cells[1])
    bycol = Dict{NTuple{2,Int32},Vector{Int}}()
    for (i, c) in pairs(cells)
        push!(get!(Vector{Int}, bycol, _zbase(c)), i)
    end
    parent = collect(1:n)
    function find(i)
        while parent[i] != i
            parent[i] = parent[parent[i]]
            i = parent[i]
        end
        return i
    end
    touch(a, b) = all(d -> a.base[d] ≤ b.base[d] + b.step[d] && b.base[d] ≤ a.base[d] + a.step[d], 3:N)
    for (i, c) in pairs(cells), dx in -1:1, dy in -1:1
        col = (c.base[1] + Int32(dx) * zs[1], c.base[2] + Int32(dy) * zs[2])
        for j in get(bycol, col, ())
            touch(c, cells[j]) && (parent[find(i)] = find(j))
        end
    end
    groups = Dict{Int,Vector{C}}()
    for i in 1:n
        push!(get!(() -> C[], groups, find(i)), cells[i])
    end
    return collect(values(groups))
end

# Label one connected group of same-size cells; returns `cells => winding` pairs
# and the evaluations spent.
#
# Labels come from `_votes`, not from faces: a 4-corner face walk aliases multiplicity
# (a triple zero centred in a face turns 270° per edge and reads as a pole), so a ±4
# face only says "a locus is here". Faces enclosing both signs do flag a zero and pole
# fused in the group, which a net vote could hide (a simple zero beside a double pole).
#
# Such a group escalates before voting; so does one whose votes are not unanimous in
# sign or all unusable. Escalation splits the cells one z level, regroups the children
# that still fire, and recurses with one less `zoom`. A zero-pole pair a few sub-cells
# apart falls into cleanly labeled sub-branches; one that stays connected to the floor
# (genuinely fused, or loci crossing at some p) keeps winding 0.
function _resolve(tree, cells; zoom)
    ws, ne = zoom > 0 && _mixed(tree.cache, cells) ? (Int[], 0) : _votes(tree, cells)
    if !isempty(ws) && (first(ws) > 0 || last(ws) < 0)
        return [cells => ws[(length(ws)+1)÷2]], ne
    end
    kids = zoom > 0 ? _deepen!(tree, cells) : empty(cells)
    isempty(kids) && return [cells => 0], ne
    out = Pair{typeof(cells),Int}[]
    for g in _group(kids)
        subs, n = _resolve(tree, g; zoom=zoom - 1)
        ne += n
        append!(out, subs)
    end
    # Still one connected piece — same locus set, possibly a sharper label; keep
    # the original coarser cells so samples stay at the requested resolution.
    length(out) == 1 && return [cells => out[1].second], ne
    return out, ne
end

function _mixed(cache, cells)
    pos = neg = false
    for c in cells, fw in _enclosing_faces(cache, c)
        pos |= fw > 0
        neg |= fw < 0
    end
    return pos && neg
end

# Sign of each **fully-enclosing** (±4) z-plane face winding of a cell; 0 for the others.
@inline function _enclosing_faces(cache, cell)
    q = _quadrants(cache, cell)
    return ntuple(Val(1 << nparams(cell))) do j
        b = (j - 1) << 2
        fw = @inbounds _face_winding(q[b+1], q[b+2], q[b+4], q[b+3])
        fw == 4 ? 1 : fw == -4 ? -1 : 0
    end
end

# Refine the cells one z level, by the same rules as the scan (so a moving locus stays
# linked); empty at the lattice floor, where child corners have no index left.
function _deepen!(tree, cells)
    cells[1].step[1] ≤ 1 && return empty(cells)
    kids = empty(cells)
    for c in cells
        append!(kids, split!(c, _ZMASK))
    end
    return _refine!(tree, kids, _zdiag(tree.lat, kids[1]))
end

# Slices a group is voted at. The outline's net winding is the same at every slice
# unless a locus enters or leaves the box between them, so extra slices mostly
# re-measure it; more than one, so a slice whose outline leaves the box or meets an
# unusable sample still leaves a vote.
const _NVOTES = 3

# Sorted usable windings of the group's outline (`_outline`) at up to `_NVOTES` parameter
# slices spread along it — slices where one of its faces fires — and the evaluations
# they cost.
function _votes(tree, cells)
    (; cache, lat) = tree
    f = cache.f.f
    M = nparams(cells[1])
    pixs = [_pcorner(c, j) for c in cells for j in 0:((1<<M)-1) if _fw_fires(_facewinding(cache, c, j))]
    isempty(pixs) && (pixs = [_pcorner(c, j) for c in cells for j in 0:((1<<M)-1)])
    unique!(sort!(pixs))
    n = length(pixs)
    ws = Int[]
    ne = 0
    for i in unique(round.(Int, range(1, n; length=min(_NVOTES, n))))
        poly = _outline(lat, cells, pixs[i])
        isnothing(poly) && continue
        w, nw = _polygon_winding(f(_pcoords(lat, pixs[i])...), poly...)
        ne += nw
        isnothing(w) || push!(ws, w)
    end
    return sort!(ws), ne
end

# The z-columns of `cells` covering parameter slice `pix`, grown by one column, as polygon
# vertices and edges with the region on their left; `nothing` if it leaves the box. It
# holds the group's every locus at the slice — a fused opposite partner cancels — and
# no other group's: their columns sit ≥ 2 columns away.
function _outline(lat, cells, pix)
    st = _zstep(cells[1])
    cols = Set{NTuple{2,Int}}()
    for c in cells
        all(d -> c.base[d+2] ≤ pix[d] ≤ c.base[d+2] + c.step[d+2], eachindex(pix)) || continue
        x, y = c.base[1] ÷ st[1], c.base[2] ÷ st[2]
        for dx in -1:1, dy in -1:1
            push!(cols, (x + dx, y + dy))
        end
    end
    nx, ny = lat.zn ÷ st[1], lat.zn ÷ st[2]
    all(((x, y),) -> 0 ≤ x < nx && 0 ≤ y < ny, cols) || return nothing
    at = Dict{NTuple{2,Int},Int}()
    zs = ComplexF64[]
    vertex(x, y) = get!(at, (x, y)) do
        push!(zs, complex(lat.zo[1] + x * st[1] * lat.zh[1], lat.zo[2] + y * st[2] * lat.zh[2]))
        length(zs)
    end
    edges = NTuple{2,Int}[]
    for (x, y) in cols
        (x, y - 1) in cols || push!(edges, (vertex(x, y), vertex(x + 1, y)))
        (x + 1, y) in cols || push!(edges, (vertex(x + 1, y), vertex(x + 1, y + 1)))
        (x, y + 1) in cols || push!(edges, (vertex(x + 1, y + 1), vertex(x, y + 1)))
        (x - 1, y) in cols || push!(edges, (vertex(x, y + 1), vertex(x, y)))
    end
    return zs, edges
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