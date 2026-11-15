# Evaluations per edge: a clean k-fold locus hugging an edge costs O(k·log(edge/distance));
# noise swamping `f` exhausts it.
const _EDGE_BUDGET = 1024

# Winding of `f` around closed polygons (`edges` index pairs into the vertices `vs`);
# `nothing` on an unusable sample (0/NaN/Inf), an unresolved phase jump, or an edge past
# `_EDGE_BUDGET`.
function _polygon_winding(F, vs, edges)
    fv = F(vs)
    all(_is_valid, fv) || return nothing
    lv = [log(complex(x)) for x in fv]
    # (edge, the stretch [a, b] of it, log f at its ends)
    segs = [(k, 0.0, 1.0, lv[i], lv[j]) for (k, (i, j)) in pairs(edges)]
    spent = zeros(Int, length(edges))
    acc = 0.0
    depth = 0
    while !isempty(segs)
        depth > 40 && return nothing
        fm = F(map(segs) do (k, a, b)
            i, j = edges[k]
            vs[i] .+ ((a + b) / 2) .* (vs[j] .- vs[i])
        end)
        next = empty(segs)
        for ((k, a, b, la, lb), x) in zip(segs, fm)
            _is_valid(x) || return nothing
            lm = log(complex(x))
            Δ = _settled_turn(la, lm, lb)
            if isnothing(Δ)
                (spent[k] += 2) < _EDGE_BUDGET || return nothing
                m = (a + b) / 2
                push!(next, (k, a, m, la, lm), (k, m, b, lm, lb))
            else
                acc += Δ
            end
        end
        segs = next
        depth += 1
    end
    return round(Int, acc / 2π)
end

# Winding of `f` on the circle |z − z0| = r from `seed` equispaced samples refined within `budget` evaluations.
function _circle_winding(f, z0, r; seed=8, budget=16seed)
    g(t) = f(z0 + r * cispi(2t))
    ls = Vector{ComplexF64}(undef, seed)
    for k in 1:seed
        gk = g((k - 1) / seed)
        _is_valid(gk) || return nothing, k
        ls[k] = log(complex(gk))
    end
    ne, acc = seed, 0.0
    for k in 1:seed
        Δ, n = _wind_refine(g, (k - 1) / seed, k / seed, ls[k], ls[mod1(k + 1, seed)], 0, budget - ne)
        ne += n
        isnan(Δ) && return nothing, ne
        acc += Δ
    end
    return round(Int, acc / 2π), ne
end

# Lifted phase increment of `g` over [a, b] from `log g` at the ends, and evaluations
# spent; NaN past `budget`, which bounds the bisection where noise swamps `g`.
function _wind_refine(g, a, b, la, lb, depth, budget)
    budget < 1 && return NaN, 0
    m = (a + b) / 2
    gm = g(m)
    _is_valid(gm) || return NaN, 1
    lm = log(complex(gm))
    Δ = _settled_turn(la, lm, lb)
    isnothing(Δ) || return Δ, 1
    depth ≥ 40 && return NaN, 1
    Δ1, n1 = _wind_refine(g, a, m, la, lm, depth + 1, budget - 1)
    isnan(Δ1) && return NaN, n1 + 1
    Δ2, n2 = _wind_refine(g, m, b, lm, lb, depth + 1, budget - 1 - n1)
    return Δ1 + Δ2, n1 + n2 + 1
end

# Lifted phase increment over a stretch from `log g` at its ends and middle; `nothing` if its halves need bisecting.
function _settled_turn(la, lm, lb)
    dl, dr = _dlog(la, lm), _dlog(lm, lb)
    return _small(dl) && _small(dr) ? imag(dl) + imag(dr) : nothing
end

_dlog(la, lb) = complex(real(lb) - real(la), rem2pi(imag(lb) - imag(la), RoundNearest))
_small(d) = abs(imag(d)) < π / 2 && abs(real(d)) < log(2)
