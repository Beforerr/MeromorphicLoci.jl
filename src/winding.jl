# Winding of `f` around closed polygons — `edges` index pairs into the vertices `zs` — and
# the evaluations spent; `nothing` marks an unusable sample (0/NaN/Inf) or an
# unresolvable phase jump.
function _polygon_winding(f, zs, edges)
    lf = Vector{ComplexF64}(undef, length(zs))
    Threads.@threads for i in eachindex(zs, lf)
        fz = f(zs[i])
        lf[i] = _is_valid(fz) ? log(complex(fz)) : NaN
    end
    ne = length(zs)
    any(isnan, lf) && return nothing, ne
    Δ = Vector{Float64}(undef, length(edges))
    n = Vector{Int}(undef, length(edges))
    Threads.@threads for k in eachindex(edges, Δ, n)
        i, j = edges[k]
        za, zb = zs[i], zs[j]
        Δ[k], n[k] = _wind_refine(t -> f(za + t * (zb - za)), 0.0, 1.0, lf[i], lf[j])
    end
    ne += sum(n; init=0)
    any(isnan, Δ) && return nothing, ne
    return round(Int, sum(Δ; init=0.0) / 2π), ne
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

# Lifted phase increment of `g` over [a, b], from `log g` at the ends, and evaluations spent.
# Bisect until each half turns the phase by < π/2 and scales |g| by < 2. The phase test
# alone aliases a half that turns by a full extra turn — a multiple locus hugging the
# path; such a locus also swings |g| across the half, which the modulus test catches.
# Where evaluation noise swamps g, every half fails the tests at every depth, so the
# bisection would branch everywhere: past `budget` evaluations the increment is NaN.
function _wind_refine(g, a, b, la, lb, depth=0, budget=typemax(Int))
    budget < 1 && return NaN, 0
    m = (a + b) / 2
    gm = g(m)
    _is_valid(gm) || return NaN, 1
    lm = log(complex(gm))
    dl, dr = _dlog(la, lm), _dlog(lm, lb)
    _small(dl) && _small(dr) && return imag(dl) + imag(dr), 1
    depth ≥ 40 && return NaN, 1
    Δ1, n1 = _wind_refine(g, a, m, la, lm, depth + 1, budget - 1)
    isnan(Δ1) && return NaN, n1 + 1
    Δ2, n2 = _wind_refine(g, m, b, lm, lb, depth + 1, budget - 1 - n1)
    return Δ1 + Δ2, n1 + n2 + 1
end

_dlog(la, lb) = complex(real(lb) - real(la), rem2pi(imag(lb) - imag(la), RoundNearest))
_small(d) = abs(imag(d)) < π / 2 && abs(real(d)) < log(2)
