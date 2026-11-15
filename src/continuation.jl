"""
    Continuation(seed; zres = nothing, zfloor = nothing, certify = true)

Loci over parameters for [`findloci`](@ref): the point finder `seed` finds them at the knots, and [`muller`](@ref) follows them
between. [How it works and what it can miss](@ref continuation).

- `zres` — the scale of a step's correction. `nothing` takes `seed.zres` of an
  `ArgumentPrinciple` seed, and a quarter of the sample spacing of an `AAA` seed.
- `zfloor` — loci closer than this are not told apart (`uncertified` of a
  [`Branch`](@ref)); two closer throughout are one locus.
  At most `zres`; `nothing` takes `zres/256`.
- `certify` — count the zeros minus poles in the box at each knot against the branches
  (`deficit` of the [`Survey`](@ref)). Always on for an `ArgumentPrinciple` seed.
"""
struct Continuation{S}
    seed::S
    zres::Union{Nothing,Float64}
    zfloor::Union{Nothing,Float64}
    certify::Bool
end
Continuation(seed; zres=nothing, zfloor=nothing, certify=true) =
    Continuation(seed, _float(zres), _float(zfloor), certify)
_float(x) = isnothing(x) ? nothing : Float64(x)

_zres(seed::ArgumentPrinciple, zbox) = seed.zres
_zres(seed::AAA, (lo, hi)) = sqrt(prod(reim(hi - lo)) / seed.nsamples) / 4
_zres(seed, zbox) = throw(ArgumentError("Continuation needs `zres` with a $(typeof(seed)) seed"))

"""
    findloci(f, zdomain, pdomain, alg; mapper = map) :: Survey

Loci `z*(p)` of the zeros and poles of `f(p...)(z)` in the box `zdomain = (zlo, zhi)` as `p` ranges over `pdomain`.

`pdomain` holds one entry per parameter: an interval `(lo, hi)`, searched at equispaced
knots (9 for one parameter, 5 each for more), or a sorted vector of knots.
"""
findloci(f, zdomain, pdomain, alg::Union{AAA,ArgumentPrinciple}; mapper=map) =
    findloci(f, zdomain, pdomain, Continuation(alg); mapper)

# Equispaced knots of an interval `(lo, hi)`: over one parameter, and per parameter over more.
const _PKNOTS = (9, 5)

_knots(p::Tuple{Real,Real}, n) = _knots(collect(range(Float64(p[1]), Float64(p[2]); length=n)), n)
function _knots(t::AbstractVector, _)
    length(t) ≥ 2 && all(isfinite, t) && issorted(t; lt=≤) ||
        throw(ArgumentError("parameter knots must be ≥ 2 finite, strictly increasing values, got $t"))
    return collect(Float64, t)
end
_knots(x, _) = throw(ArgumentError(
    "each pdomain entry must be an interval (lo, hi) or a knot vector, got $(repr(x)); one parameter is ((lo, hi),)"))

# A corrector within `tol` of its prediction takes another locus for its own only within
# `_CAPTURE·tol` of it: a resolution `zres` allows corrections up to `zres/_CAPTURE`.
const _CAPTURE = 2.5

# A step is accepted if `log d` of the new root's germ (`_polishin`) is within this of its
# extrapolation. An exchange of two loci turns `d` by a half turn, a meeting by a quarter.
const _TURN = 0.5

struct _Config
    box::NTuple{2,ComplexF64}
    tol::Float64      # `zres/_CAPTURE`: the largest correction a step is accepted with
    floor::Float64    # `zfloor`: loci closer than this are not told apart
    mapper::Any
end

# `_Config` and one line followed. `f` is untyped so that only the corrector and the
# windings compile anew for each `f`.
struct _Tracker
    f::Any            # `p -> f(p)` along the line followed
    box::NTuple{2,ComplexF64}
    tol::Float64
    floor::Float64
    K::Vector{Float64}
    gs::Vector{Any}   # `f(K[k])`
    mapper::Any
end

function _fun(c::_Tracker, p)
    k = searchsortedfirst(c.K, p)
    return k <= length(c.K) && c.K[k] == p ? c.gs[k] : c.f(p)
end

# `f` inside the box, NaN outside (the corrector backs off); counts evaluations.
struct _Boxed{F}
    f::F
    box::NTuple{2,ComplexF64}
    n::Base.RefValue{Int}
end
(g::_Boxed)(z) = _inside(g.box, z, 0.0) ? (g.n[] += 1; g.f(z)) : complex(NaN, NaN)

# A locus of multiplicity `m` as a simple zero: of `g`, or of `g/g′`.
_simple(g, m, h) = m == 1 ? g : z -> g(z) / ((g(z + h) - g(z - h)) / 2h)

# eps^(1/3) of the box: balances truncation against rounding.
_fdstep(c) = 1.0e-5 * abs(c.box[2] - c.box[1])

# What a locus of `kind` is followed as: a function with a simple zero on it.
_target(c, g, kind, n) = _simple(_Boxed(kind > 0 ? g : _reciprocal(g), c.box, n), abs(kind), _fdstep(c))

# Muller's method on `target` from `z`, its starting triangle turned into the box. Also the
# germ of the root `z*`, from the quadratic through the triangle: `d = target′(z*)` and the
# radius `R = |2target′/target″|` within which no other locus lies. With
# `target = (z − z*)·q`, `d` is `q(z*)`: the separations of `z*` from every other locus
# (those of the opposite kind inverted) times a factor that varies slowly, and `1/R` the
# sum of their inverses. Whatever comes near, seeded or not, shows in both; loci placed
# evenly around `z*` cancel in `R`.
function _polishin(box, target, z, tol)
    lo, hi = box
    h = tol / 4
    a = real(z) - h < real(lo) ? z + h : z - h
    b = imag(z) + h > imag(hi) ? z - h * im : z + h * im
    fa, fb, fz = target(a), target(b), target(z)
    seen(x) = x == z ? fz : x == a ? fa : x == b ? fb : target(x)
    r = _muller(seen, a, b, z, z, 2tol)
    root = ComplexF64(r.root)
    slope = (fz - fb) / (z - b)
    curv = (slope - (fb - fa) / (b - a)) / (z - a)
    d = slope + curv * (z - b) + 2curv * (root - z)
    R = abs(d) / abs(curv)
    return root, r.converged, ComplexF64(d), isnan(R) ? 0.0 : R
end

# All `_advance!` asks of a corrector: the root, and its germ `d`, `R` (`_polishin`), which
# derivatives of `g` would give exactly. Return types asserted: `g` is not inferred.
function _correct(c, g, kind, z, tol)::Tuple{ComplexF64,Bool,Int,ComplexF64,Float64}
    n = Ref(0)
    root, ok, d, R = _polishin(c.box, _target(c, g, kind, n), z, tol)
    return root, ok, n[], d, R
end

# The germ at a node, read at the tolerance its radius allows.
function _germ(c, g, kind, z)
    tol, ne = c.tol, 0
    while true
        _, _, n, d, R = _correct(c, g, kind, z, tol)
        ne += n
        tl = _tolerance(c, R)
        tl >= tol / 2 && return d, R, ne
        tol = tl
    end
end

# The corrector's tolerance beside a locus `R` away: its starting triangle, `tol/4` wide,
# is then small against `R`.
_tolerance(c, R) = clamp(R / 8, c.floor / _CAPTURE, c.tol)

function _winds(c, g, z, r, kind)::Tuple{Bool,Int}
    n = Ref(0)
    w, _ = _circle_winding(_Boxed(g, c.box, n), z, r)
    return w == kind, n[]
end

mutable struct _Node
    z::ComplexF64
    w::Int
    knot::Int
    hist::Vector{Tuple{Float64,ComplexF64}}   # ≤ 2 points beside the knot, nearest last: a predictor's history
    dp::Float64
    d::ComplexF64     # germ (`_polishin`), `R` NaN until measured
    R::Float64
    s::ComplexF64     # d(log d)/dp of the track that arrived
end

mutable struct _Track
    kind::Int
    ps::Vector{Float64}
    zs::Vector{ComplexF64}
    hist::Vector{Tuple{Float64,ComplexF64}}   # points of the locus before `ps[1]`
    target::Float64
    from::Int          # node it starts at
    to::Int            # node it is linked to at `target`, 0 if none
    status::Symbol     # :arrived at `target`, :exit through the box edge, :stuck
    dp::Float64        # next step
    nevals::Int
    d::ComplexF64      # germ at the last sample
    R::Float64
    s::ComplexF64      # d(log d)/dp there
    hot::Int           # last sample before a stretch within the floor of another locus, 0 outside one
    dhot::ComplexF64   # `d` there
    flags::Vector{NTuple{2,Float64}}   # stretches left uncertified
end
_Track(kind, p, z, hist, target, from) =
    _Track(kind, [p], [z], hist, target, from, 0, :stuck, 0.0, 0, NaN, NaN, 0.0, 0, NaN, [])

function _truncate!(t::_Track, n)
    resize!(t.ps, n)
    resize!(t.zs, n)
    t.status = :stuck
    t.hot = 0
    return t
end

_forward(t::_Track) = t.target > t.ps[1]

# Extrapolate the locus to `p` through its last ≤ 3 points; also the polynomial's degree.
function _predict(t::_Track, p)
    n = length(t.ps)
    m = min(3, n + length(t.hist))
    pt(i) = i <= n ? (t.ps[n+1-i], t.zs[n+1-i]) : t.hist[end-(i-n)+1]
    z = zero(ComplexF64)
    for i in 1:m
        q, zi = pt(i)
        l = 1.0
        for j in 1:m
            j == i || (l *= (p - pt(j)[1]) / (q - pt(j)[1]))
        end
        z += l * zi
    end
    return z, m - 1
end

# Follow `t` from its last point to `p1`. The state followed is the germ `(z, d)`: a step
# is accepted if the root is within the tolerance of its prediction and `log d` within
# `_TURN` of its own, linear in `p`. A step onto another locus, or over an avoided
# crossing where the loci themselves exchange directions, lands on a root that predicts
# well with `d` turned; steps therefore shorten against the distance to any locus closing
# in, as `1/|d(log d)/dp|` is the parameter distance to a meeting with it.
#
# Within the floor of another locus (`R < floor`) steps are accepted on the root alone.
# Out the far side, `d` has turned by about a half turn if the two passed or crossed,
# and the stretch is flagged; by about a quarter if they met and turned (a branch point
# in `p`), where which outgoing locus continues the track has no answer: it ends there.
function _advance!(c, t::_Track, p1, dp)
    L = abs(p1 - t.ps[end])
    dp > 0 || (dp = L / 16)
    dir = sign(p1 - t.ps[end])
    function ended(status)
        if t.hot > 0
            pa = t.ps[t.hot]
            status === :stuck ? _truncate!(t, t.hot) : push!(t.flags, minmax(pa, t.ps[end]))
        end
        return status
    end
    while true
        pc, zc = t.ps[end], t.zs[end]
        left = abs(p1 - pc)
        left == 0 && return ended(:arrived)
        tol = _tolerance(c, t.R)
        inner(z) = _inside(c.box, z, -tol / 4)
        step = min(dp, left)
        p = step >= left * (1 - 1.0e-12) ? p1 : pc + dir * step
        zp, ord = _predict(t, p)
        if !_inside(c.box, zp, 0.0)
            # predicted to leave: step only as far as the prediction stays `tol/4` inside,
            # where the corrector's triangle fits
            lo, hi = 0.0, 1.0
            for _ in 1:30
                mid = (lo + hi) / 2
                inner(_predict(t, pc + mid * (p - pc))[1]) ? (lo = mid) : (hi = mid)
            end
            p = pc + (inner(zc) ? lo : 0.0) * (p - pc)
            zp, ord = _predict(t, p)
            abs(zp - zc) < tol / 4 && return ended(:exit)
            step = abs(p - pc)
        end
        g = _fun(c, p)
        z, ok, ne, d, R = _correct(c, g, t.kind, zp, tol)
        t.nevals += ne
        # a multiple locus is followed on its root alone: its `d` is `1/m` whatever
        # lies near, and the loci it stands for lie within the floor of it
        abs(t.kind) == 1 || ((d, R) = (t.d, t.R))
        err = abs(z - zp)
        turn = log(d / (t.d * exp(t.s * (p - pc))))
        ok &= err <= tol && isfinite(d) && (t.hot > 0 || abs(turn) <= _TURN)
        if ok
            push!(t.ps, p)
            push!(t.zs, z)
            grow = _grow(tol, err, ord)
            if t.hot == 0 && R < c.floor
                t.hot, t.dhot = length(t.ps) - 1, t.d
            elseif t.hot > 0 && R >= c.floor
                pa = t.ps[t.hot]
                if π / 4 < abs(angle(d / t.dhot)) < 3π / 4
                    _truncate!(t, t.hot)
                    return :stuck
                end
                push!(t.flags, minmax(pa, p))
                t.hot, t.s = 0, 0.0
            elseif t.hot == 0
                t.s += turn / (p - pc)
                grow = min(grow, _grow(_TURN, abs(turn), 1))
            end
            t.d, t.R = d, R
            # a step cut short to land on `p1` says nothing against the longer one
            t.dp = dp = max(step * grow, p == p1 ? dp : 0.0)
        else
            # beside the edge, a locus the corrector cannot follow is leaving
            step < 1.0e-10L && return ended(inner(zc) ? :stuck : :exit)
            dp = step / 3
        end
    end
end

_grow(tol, err, ord) = clamp(0.8 * (tol / max(err, 1.0e-3tol))^(1 / (ord + 1)), 0.3, 2.0)

function _launch(c, nodes, ids, p0, p1)
    ts = c.mapper(eachindex(ids)) do i
        n = nodes[ids[i]]
        # history from the far side of the knot only
        hist = all(q -> (q[1] - p0) * (p1 - p0) < 0, n.hist) ? n.hist : empty(n.hist)
        t = _Track(n.w, p0, n.z, hist, p1, ids[i])
        isnan(n.R) && ((n.d, n.R, t.nevals) = abs(n.w) == 1 ? _germ(c, _fun(c, p0), n.w, n.z) : (1.0, Inf, 0))
        t.d, t.R, t.s = n.d, n.R, n.s
        t.R < c.floor && ((t.hot, t.dhot) = (1, t.d))
        dp = isempty(hist) ? abs(p1 - p0) / 16 : n.dp
        if abs(n.w) == 1 || !isempty(hist)
            t.status = _advance!(c, t, p1, dp)
        else
            # A multiple seed no multiple track arrived on may be simple loci that only
            # coincide at the knot: if they part a step off it, they are followed from
            # the knots where they are apart, not from this seed.
            t.status = _advance!(c, t, p0 + (p1 - p0) / 16, dp)
            ok, ne = _winds(c, _fun(c, t.ps[end]), t.zs[end], c.floor, t.kind)
            t.nevals += ne
            ok || _truncate!(t, 1)
            ok && t.status === :arrived && (t.status = _advance!(c, t, p1, t.dp))
        end
        return t
    end
    return collect(_Track, ts)
end

struct _State
    c::_Tracker
    nodes::Vector{_Node}
    at::Vector{Vector{Int}}        # nodes of each knot
    spans::Vector{Vector{_Track}}  # tracks of each knot interval
    nseeds::Int                    # the leading nodes it was given; the rest, tracks carried in
end

# The node of knot `k` that a track of `kind` arriving at `z` continues, 0 if none. A
# multiple node stands for every locus within the floor of it.
function _match(st::_State, k, kind, z)
    i = findfirst(i -> _same(st.c, st.nodes[i], kind, z), st.at[k])
    return isnothing(i) ? 0 : st.at[k][i]
end

_same(c, o, w, z) = o.w == w && abs(o.z - z) <= (abs(w) > 1 ? c.floor : c.floor / 20)

function _addnode!(st::_State, k, z, w)
    push!(st.nodes, _Node(z, w, k, [], 0.0, NaN, NaN, 0.0))
    push!(st.at[k], length(st.nodes))
    return length(st.nodes)
end

function _remember!(n::_Node, t::_Track)
    m = length(t.ps)
    n.hist = [(t.ps[i], t.zs[i]) for i in max(1, m - 2):m-1]
    n.dp, n.d, n.R, n.s = t.dp, t.d, t.R, t.s
    return
end

_flag!(t::_Track, p) = push!(t.flags, (p, p))

# A track arriving at knot `k` on no node of its kind: a locus the seeder missed there,
# which becomes a node. Not so if it arrives on a node of higher multiplicity, as one of the
# loci that coincide there, nor if it is multiple and no longer winds as its kind (it is
# then uncertified throughout). Returns the node, 0 if none.
function _newnode!(st::_State, t::_Track, k)
    if any(j -> (n = st.nodes[j]; sign(n.w) == sign(t.kind) && abs(n.w) > abs(t.kind) && abs(n.z - t.zs[end]) <= st.c.floor), st.at[k])
        _flag!(t, t.ps[end])
        return 0
    end
    if abs(t.kind) > 1
        ok, ne = _winds(st.c, st.c.gs[k], t.zs[end], st.c.floor, t.kind)
        t.nevals += ne
        if !ok
            t.status = :stuck
            push!(t.flags, minmax(t.ps[1], t.ps[end]))
            return 0
        end
    end
    return _addnode!(st, k, t.zs[end], t.kind)
end

# Follow the loci of interval `k`: forward from the nodes of knot `k`, then backward from
# each node of knot `k+1` no forward track reached. A backward arrival on no node of knot
# `k` becomes one, and is followed on backward through the earlier intervals
# (`forward = false`: only the nodes new at knot `k+1` are followed).
function _span!(st::_State, k; forward=true)
    (; c, nodes, at) = st
    K = c.K
    Ts = st.spans[k]
    if forward
        append!(Ts, _launch(c, nodes, at[k], K[k], K[k+1]))
        arrived = Dict{Int,_Track}()
        for t in Ts
            t.status === :arrived || continue
            j = _match(st, k + 1, t.kind, t.zs[end])
            j == 0 && (j = _newnode!(st, t, k + 1))
            j == 0 && continue
            if haskey(arrived, j)
                # two tracks on one node met there
                _flag!(t, K[k+1])
                _flag!(arrived[j], K[k+1])
                continue
            end
            arrived[j] = t
            t.to = j
            _remember!(nodes[j], t)
        end
    end
    reached = Set(_forward(t) ? t.to : t.from for t in Ts)
    new = _launch(c, nodes, filter(!in(reached), at[k+1]), K[k+1], K[k])
    append!(Ts, new)
    carried = false
    for t in new
        m = min(3, length(t.ps))
        nodes[t.from].hist = [(t.ps[i], t.zs[i]) for i in m:-1:2]
        nodes[t.from].dp = m > 1 ? abs(t.ps[2] - t.ps[1]) : 0.0
        t.status === :arrived || continue
        i = _match(st, k, t.kind, t.zs[end])
        if i == 0
            t.to = _newnode!(st, t, k)
            t.to == 0 && continue
            _remember!(nodes[t.to], t)
            carried = true
            continue
        end
        # the node's own track followed this locus the other way: the two are one branch
        # if it ended stuck on the way; ending elsewhere, one of them changed locus
        others = filter(u -> u !== t && (_forward(u) ? u.from : u.to) == i, Ts)
        if all(u -> u.status === :stuck, others)
            t.to = i
        else
            foreach(u -> _flag!(u, K[k]), (t, others...))
        end
    end
    carried && k > 1 && _span!(st, k - 1; forward=false)
    return
end

# The loci a seed of winding `w` near `z` stands for, as `(; z, w)`, and the winding left
# unaccounted for. A seed of |w| ≥ 2 may be simple loci closer than the seeder resolves,
# of which `g/g′` would follow one: each is polished with those found divided out, and
# counted by the winding within the floor of it.
function _members(c, g, z, w)
    n = Ref(0)
    gb = _Boxed(g, c.box, n)
    gt = w > 0 ? gb : _Boxed(_reciprocal(g), c.box, n)
    found = @NamedTuple{z::ComplexF64, w::Int}[]
    left = w
    while left != 0
        rest(x) = gt(x) / prod(m -> (x - m.z)^abs(m.w), found; init=one(x))
        root, ok = _polishin(c.box, _simple(rest, abs(left), _fdstep(c)), z, c.tol)
        ok || break
        wv, _ = _circle_winding(gb, root, c.floor)
        wv isa Int && sign(wv) == sign(w) && abs(wv) <= abs(left) || break
        push!(found, (; z=root, w=wv))
        left -= wv
    end
    return found, left, n[]
end

# The loci of `g = f(p)` at a knot: `found` to follow, `odd` that cannot be (kind
# unreadable, or left over by `_members`), and the seeder's count of zeros minus poles
# (`nothing` if unreadable, `missing` if it does not count).
function _seeds(c, seed, g)
    r = _points(g, c.box, seed, Val(0), Muller(), c.mapper)
    count = !haskey(r, :deficit) ? missing : isnothing(r.deficit) ? nothing : r.deficit + sum(r.windings; init=0)
    ne = r.nevals
    found = @NamedTuple{z::ComplexF64, w::Int}[]
    # `AAA` reports points just outside the box, where nothing is followed nor counted
    odd = [(; z, w=0) for z in r.unresolved if _inside(c.box, z, 0.0)]
    for (z, w, ok) in zip(r.points, r.windings, r.converged)
        _inside(c.box, z, 0.0) || continue
        if ok && abs(w) == 1
            push!(found, (; z, w))
        else
            ms, left, n = _members(c, g, z, w)
            ne += n
            append!(found, ms)
            left == 0 || push!(odd, (; z, w=left))
        end
    end
    return (; found, odd, nevals=ne, count)
end

# Zeros minus poles of `g` in the box by the winding of its outline. A phase that turns
# many times between outline samples reads as a wrong count, not as unreadable, so a count
# differing from `found` is taken again twice as finely, down to `zres`.
function _boxcount(c, g, found)
    lo, hi = c.box
    w = reim(hi - lo)
    n = _zroots(w, hypot(w...) / 8)
    ne = 0
    while true
        vs, edges = _boxpoly(w ./ n, n)
        # rounding can put the far corners an ulp outside
        zs = [complex(clamp(real(lo) + x, real(lo), real(hi)), clamp(imag(lo) + y, imag(lo), imag(hi))) for (x, y) in vs]
        count = _polygon_winding(zs -> (ne += length(zs); c.mapper(g, zs)), zs, edges)
        (isnothing(count) || count == found || hypot((w ./ n)...) <= _CAPTURE * c.tol) && return count, ne
        n = 2 .* n
    end
end

function findloci(f, zdomain, pdomain, alg::Continuation; mapper=map)
    pdomain isa Tuple && !isempty(pdomain) || throw(ArgumentError(
        "pdomain holds one entry per parameter, got $(repr(pdomain)); one parameter is ((lo, hi),)"))
    zdomain isa Tuple{Number,Number} || throw(ArgumentError(
        "Continuation needs a box (zlo, zhi) for zdomain, got $(repr(zdomain))"))
    box = ComplexF64.(_zcorners(zdomain...))
    zres = something(alg.zres, _zres(alg.seed, box))
    zres > 0 || throw(ArgumentError("zres must be positive, got $zres"))
    zfloor = something(alg.zfloor, zres / 256)
    0 < zfloor <= zres || throw(ArgumentError("zfloor must be in (0, zres = $zres], got $zfloor"))
    Ks = map(p -> _knots(p, _PKNOTS[min(length(pdomain), 2)]), pdomain)
    c = _Config(box, zres / _CAPTURE, zfloor, mapper)
    return _survey(f, Ks, c, alg.seed, zres, alg.certify)
end

# Seed the nodes of the knot grid, follow each grid line from the loci at its nodes
# (`_follow`), and join node loci linked by tracks into branches.
function _survey(f, Ks::NTuple{M,Vector{Float64}}, c::_Config, seed, zres, certify) where {M}
    G = CartesianIndices(map(length, Ks))
    coords(I) = ntuple(d -> Ks[d][I[d]], Val(M))
    seeded = c.mapper(vec(G)) do I
        g = f(coords(I)...)
        return g, _seeds(c, seed, g)
    end
    gs = reshape(Any[first(x) for x in seeded], size(G))
    seeds = reshape(Any[last(x) for x in seeded], size(G))
    # the loci at each node: its seeds, then those a track carried in
    loci = [copy(s.found) for s in seeds]
    lines = [(d, J) for d in 1:M for J in CartesianIndices(ntuple(e -> e == d ? (1:1) : axes(G, e), Val(M)))]
    node(d, J, k) = CartesianIndex(Base.setindex(Tuple(J), k, d))
    states = Vector{_State}(undef, length(lines))
    ids = Vector{Vector{Int}}(undef, length(lines))    # each line node's index among its grid node's loci
    ne = sum(s -> s.nevals, seeds)
    # A line is followed again when another carries a new locus into one of its nodes: a
    # locus hops one line a pass, across the grid in at most this many.
    dirty = trues(length(lines))
    for _ in 1:sum(length, Ks)
        todo = findall(dirty)
        isempty(todo) && break
        dirty .= false
        states[todo] = c.mapper(todo) do l
            d, J = lines[l]
            ks = eachindex(Ks[d])
            along(t) = f(Base.setindex(coords(J), t, d)...)
            return _follow(c, along, Ks[d], [gs[node(d, J, k)] for k in ks], [loci[node(d, J, k)] for k in ks])
        end
        for l in todo
            d, J = lines[l]
            st = states[l]
            ne += sum(sp -> sum(t -> t.nevals, sp; init=0), st.spans)
            id = ids[l] = zeros(Int, length(st.nodes))
            seen = zeros(Int, length(Ks[d]))
            for (i, n) in pairs(st.nodes)
                at = loci[node(d, J, n.knot)]
                if i <= st.nseeds
                    id[i] = (seen[n.knot] += 1)
                    continue
                end
                # carried in by a track: known at the node from another line, or new there
                j = findfirst(o -> _same(c, o, n.w, n.z), at)
                if isnothing(j)
                    push!(at, (; z=n.z, w=n.w))
                    j = length(at)
                    for e in 1:M
                        e == d || (dirty[findfirst(==((e, node(e, node(d, J, n.knot), 1))), lines)] = true)
                    end
                end
                id[i] = j
            end
        end
    end

    # one branch per set of node loci linked by tracks
    first_ = cumsum([0; length.(vec(loci))])
    gid(I, j) = first_[LinearIndices(G)[I]] + j
    parent = collect(1:first_[end])
    find(i) = parent[i] == i ? i : (parent[i] = find(parent[i]))
    of(l, i) = ((d, J) = lines[l]; gid(node(d, J, states[l].nodes[i].knot), ids[l][i]))
    for (l, st) in pairs(states), sp in st.spans, t in sp
        t.to != 0 && (parent[find(of(l, t.from))] = find(of(l, t.to)))
    end
    samples = Dict{Int,Vector{Sample{M}}}()
    unc = Dict{Int,Vector{NTuple{2,NTuple{M,Float64}}}}()
    winds = Dict{Int,Int}()
    degree = Dict{Int,Int}()
    for I in G
        here = Dict{Int,Int}()
        for (j, (; z, w)) in pairs(loci[I])
            r = find(gid(I, j))
            push!(get!(Vector{Sample{M}}, samples, r), Sample{M}(z, coords(I)))
            winds[r] = w
            here[r] = get(here, r, 0) + 1
        end
        mergewith!(max, degree, here)
    end
    for (l, st) in pairs(states)
        d, J = lines[l]
        on(t) = Base.setindex(coords(J), t, d)
        flags = Dict{Int,Vector{NTuple{2,Float64}}}()
        flag(t, p1, p2) = push!(get!(Vector{NTuple{2,Float64}}, flags, find(of(l, t.from))), (p1, p2))
        for sp in st.spans, t in sp
            for i in 2:length(t.ps)-(t.to != 0)
                push!(samples[find(of(l, t.from))], Sample{M}(t.zs[i], on(t.ps[i])))
            end
            t.status === :stuck && flag(t, t.ps[end], t.ps[end])
            foreach(fl -> flag(t, fl...), t.flags)
        end
        for (r, fl) in flags
            append!(get!(Vector{NTuple{2,NTuple{M,Float64}}}, unc, r), ((on(p1), on(p2)) for (p1, p2) in _union(fl)))
        end
    end
    branches = [Branch{M}(sort!(xs; by=x -> x.p), winds[r], sort!(get(unc, r, NTuple{2,NTuple{M,Float64}}[])), degree[r])
                for (r, xs) in samples]

    # seeds not followed, with no followed locus of their node within `zres`
    loose = [filter(o -> all(l -> abs(l.z - o.z) >= zres, loci[i]), seeds[i].odd) for i in eachindex(seeds)]
    for i in eachindex(loose), (; z, w) in loose[i]
        p = coords(G[i])
        push!(branches, Branch{M}([Sample{M}(z, p)], w, [(p, p)], 1))
    end
    sort!(branches; by=b -> (first(b).p, reim(first(b).z)))

    # zeros minus poles in the box at each node: `missing` where not counted
    held = [sum(o -> o.w, loci[i]; init=0) + sum(o -> o.w, loose[i]; init=0) for i in eachindex(loci)]
    counts = c.mapper(eachindex(seeds)) do i
        return ismissing(seeds[i].count) && certify ? _boxcount(c, gs[i], held[i]) : (seeds[i].count, 0)
    end
    # a line still to be followed again leaves its nodes unchecked
    stale = falses(size(G))
    for l in findall(dirty), k in eachindex(Ks[lines[l][1]])
        stale[node(lines[l]..., k)] = true
    end
    deficit = Pair{NTuple{M,Float64},Union{Nothing,Int}}[]
    for (i, (n, _)) in pairs(counts)
        stale[i] && (n = nothing)
        ismissing(n) || n == held[i] || push!(deficit, coords(G[i]) => (isnothing(n) ? nothing : n - held[i]))
    end
    return Survey{M}(branches, ne + sum(last, counts), deficit)
end

# Follow the loci `found` at each knot `K` of one parameter, `along(t)` being `f` there and
# `gs` its value at the knots.
function _follow(c::_Config, along, K, gs, found)
    c = _Tracker(along, c.box, c.tol, c.floor, K, gs, c.mapper)
    st = _State(c, _Node[], [Int[] for _ in K], [_Track[] for _ in 2:length(K)], sum(length, found))
    for (k, s) in pairs(found), (; z, w) in s
        _addnode!(st, k, z, w)
    end
    for k in 1:length(K)-1
        _span!(st, k)
    end
    return st
end

# Sorted union of closed intervals.
function _union(I)
    out = NTuple{2,Float64}[]
    for (a, b) in sort(I)
        isempty(out) || a > out[end][2] ? push!(out, (a, b)) : (out[end] = (out[end][1], max(out[end][2], b)))
    end
    return out
end
