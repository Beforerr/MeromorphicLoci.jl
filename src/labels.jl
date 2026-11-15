# Group the firing cells into loci and label each by winding: `cells => winding` pairs.
function _label(tree::Tree, cells; zoom)
    isempty(cells) && return Pair{Vector{Cell},Int}[]
    return _resolve(tree, _group(cells); zoom)
end

# Connected components (union-find) of a nonempty set of cells of one size, adjacent
# where they share an edge or a corner.
function _group(cells::Vector{Cell})
    parent = collect(eachindex(cells))
    function find(i)
        while parent[i] != i
            parent[i] = parent[parent[i]]
            i = parent[i]
        end
        return i
    end
    st = cells[1].step
    at = Dict(c.base => i for (i, c) in pairs(cells))
    for (i, c) in pairs(cells), dx in -1:1, dy in -1:1
        k = get(at, (c.base[1] + Int32(dx) * st[1], c.base[2] + Int32(dy) * st[2]), i)
        find(i) != find(k) && (parent[find(i)] = find(k))
    end
    groups = Dict{Int,Vector{Cell}}()
    for i in eachindex(cells)
        push!(get!(() -> Cell[], groups, find(i)), cells[i])
    end
    return collect(values(groups))
end

# Label connected groups of same-size cells by `_vote`, escalating as in
# docs/src/argument_principle.md (Labels, Escalation); returns `cells => winding`
# pairs, winding 0 for a group left unlabelled.
function _resolve(tree, groups; zoom)
    out = Pair{Vector{Cell},Int}[]
    # (cells, the coarsest ancestor cells of the same locus set): while escalation leaves a
    # group in one piece, the label sharpens but its cells stay at the requested resolution.
    todo = [(g, g) for g in groups]
    for level in zoom:-1:0
        deepen = Tuple{Vector{Cell},Vector{Cell},Bool}[]
        for (cells, top) in todo
            w, lo, hi = _outlined(tree, cells)
            # at the last level, or at the lattice floor, the net vote stands
            last = level == 0 || cells[1].step[1] ≤ 1
            mixed = !last && !isnothing(w) && _mixed(tree.memo, cells, lo, hi)
            if !isnothing(w) && w != 0 && !mixed
                push!(out, top => w)
            elseif last
                push!(out, top => 0)
            else
                push!(deepen, (cells, top, !mixed))
            end
        end
        empty!(todo)
        for (cells, top, uncounted) in deepen
            kids = _deepen!(tree, cells)
            if isempty(kids)
                # nothing fires a level finer and no outline counts a locus: dropped
                uncounted || push!(out, top => 0)
                continue
            end
            gs = _group(kids)
            length(gs) == 1 ? push!(todo, (gs[1], top)) : append!(todo, ((g, g) for g in gs))
        end
        isempty(todo) && break
    end
    return out
end

# Whether the group holds both zero and pole, maximum modulus principle.
function _mixed(memo, cells, lo, hi)
    zero = pole = false
    for c in cells, s in c.corners
        a = abs2(memo.vals[s])
        zero |= a < lo
        pole |= a > hi
    end
    return zero && pole
end

# Refine the cells one level, above the lattice floor.
function _deepen!(tree, cells)
    kids = Cell[]
    for c in cells
        append!(kids, split!(tree.memo, c))
    end
    return _refine!(tree, kids, _zdiag(tree.memo.lat, kids[1]))
end

# Winding of `f` around the group's outline, `nothing` if unusable, and the
# extremes of |f|² sampled on it.
function _outlined(tree, cells)
    (; memo) = tree
    seen = ComplexF64[]
    w = _polygon_winding(_outline(memo.lat, cells)...) do pts
        vals = memo(pts)
        append!(seen, vals)
        return vals
    end
    lo, hi = extrema(abs2, seen)
    return w, lo, hi
end
_vote(tree, cells) = first(_outlined(tree, cells))

# The cells grown by one cell all round and clipped to the box, as polygon vertices
# and edges with the region on their left (docs/src/argument_principle.md, Labels).
function _outline(lat, cells)
    st = cells[1].step
    nx, ny = lat.zn[1] ÷ st[1], lat.zn[2] ÷ st[2]
    block = Set{NTuple{2,Int}}()
    for c in cells
        x, y = c.base[1] ÷ st[1], c.base[2] ÷ st[2]
        for dx in -1:1, dy in -1:1
            0 ≤ x + dx < nx && 0 ≤ y + dy < ny && push!(block, (x + dx, y + dy))
        end
    end
    return _blockpoly(st, block)
end

# Boundary of a set of cells `st` wide, indexed in units of `st`, as polygon vertices and
# edges with the region on their left.
function _blockpoly(st, block)
    at = Dict{NTuple{2,Int},Int}()
    vs = NTuple{2,Float64}[]
    vertex(x, y) = get!(at, (x, y)) do
        push!(vs, (Float64(x * st[1]), Float64(y * st[2])))
        length(vs)
    end
    edges = NTuple{2,Int}[]
    for (x, y) in block
        (x, y - 1) in block || push!(edges, (vertex(x, y), vertex(x + 1, y)))
        (x + 1, y) in block || push!(edges, (vertex(x + 1, y), vertex(x + 1, y + 1)))
        (x, y + 1) in block || push!(edges, (vertex(x + 1, y + 1), vertex(x, y + 1)))
        (x - 1, y) in block || push!(edges, (vertex(x, y + 1), vertex(x, y)))
    end
    return vs, edges
end
