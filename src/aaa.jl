"""
    AAA(; nsamples = 320, tol = 1e-13, max_degree = 150, stagnation = 5)

Zeros and poles of `f(z)` from a rational fit to its samples (the AAA algorithm).
`zdomain` is either a box `(zlo, zhi)` of opposite complex
corners, sampled evenly at about `nsamples` points, with loci reported inside it; or an
array of sample points, with loci reported wherever the fit finds them.

- `nsamples` — the fit finds at most about half as many zeros and poles as it has samples;
  raise it for a box that holds many.
- `tol` — the fit stops once it matches the samples to `tol` relative to their largest
  magnitude, so a locus whose effect on the samples is smaller goes unfound. Noise in `f`
  above `tol` stalls the fit short (`fitconverged = false`) with spurious and split
  candidates: the confirming windings drop the first and count the multiplicity of the
  second, but loci buried in the noise (multiple zeros first) stay lost.
- `max_degree`, `stagnation` — the fit stops at this degree, or after this many steps
  without improving.

`unresolved` lists candidates that `f` was too noisy or undefined around to confirm. `f` can
return `NaN` where its value cannot be trusted; those samples are skipped. A zero and a pole
of `f` much closer together than the samples are not reported.
"""
Base.@kwdef struct AAA
    nsamples::Int = 320
    tol::Float64 = 1.0e-13
    max_degree::Int = 150
    stagnation::Int = 5
end

_loci(f, zdomain, pdomain, ::AAA) =
    throw(ArgumentError("AAA finds loci at fixed parameters only: use findzeros, findpoles or findzeropoles, or ArgumentPrinciple for findloci"))

# Samples and a filter `inside(z, pad = 0)` on loci. Sample points bound no region.
_region(z::AbstractArray, m) = vec(z), Returns(true)
function _region(zbox::Tuple{Number, Number}, m)
    lo, hi = _zcorners(zbox...)
    w, h = reim(hi - lo)
    w > 0 && h > 0 || throw(ArgumentError("the box must have positive width and height, got $zbox"))
    # Near-square cells: the spacing caps the confirming circles, which must span a noisy
    # multiple locus's split candidates. ≥ 4 rows: with 2, noisy fits lost loci.
    ny = max(4, round(Int, sqrt(m * h / w)))
    nx = max(4, round(Int, m / ny))
    Z = vec(range(real(lo), real(hi), nx) .+ im .* range(imag(lo), imag(hi), ny)')
    inside(z, pad = 0) = real(lo) - pad <= real(z) <= real(hi) + pad && imag(lo) - pad <= imag(z) <= imag(hi) + pad
    return Z, inside
end

# Median distance from a sample to its nearest neighbour.
function _spacing(Z)
    d = [minimum(j -> j == i ? Inf : abs(Z[j] - Z[i]), eachindex(Z); init = Inf) for i in eachindex(Z)]
    return partialsort!(d, (length(d) + 1) ÷ 2)
end

# Zeros of f are the poles of the fit to 1/f, and its poles that fit's zeros (findpoles fits f
# itself). Each candidate is confirmed by the winding of f on a circle around it, kept clear of
# the other candidates (see `_circle`) so that it holds one candidate's locus and the polish
# cannot reach a neighbour's.
#
# Candidates within spacing/1000 are one locus. Noise splits a multiple locus wider, and each
# fragment's circle misses it or drowns in noise, so leftovers are retried as clusters (up to
# spacing/2) on one circle enclosing them; f unusable there gives one `unresolved` per cluster.
function _enclose(f, zdomain, alg::AAA, ::Val{W}) where {W}
    Z, inside = _region(zdomain, alg.nsamples)
    spacing = _spacing(Z)
    (; fit, nevals, converged) = _fit(W == -1 ? f : inv ∘ f, Z, alg)
    zc, pc = W == -1 ? (_zeros(fit), poles(fit)) : (poles(fit), _zeros(fit))
    # A zero and a pole as close as one locus's candidates cancel: noise fills the fit with such
    # Froissart doublets, which f does not show, so they neither count nor bound a circle.
    lone(zs, others) = filter(z -> all(o -> abs(o - z) >= spacing / 1000, others), zs)
    zc, pc = lone(zc, pc), lone(pc, zc)
    cands = [zc; pc]
    loci = @NamedTuple{z::ComplexF64, r::Float64, w::Int}[]
    n = nevals
    mean(g) = sum(g) / length(g)
    for (sgn, zs) in ((1, zc), (-1, pc))
        W == -sgn && continue
        zs = filter(z -> inside(z, spacing / 2), zs)
        groups = [zs[i] for i in _clusters(zs, spacing / 1000)]
        state = fill(:left, length(groups))
        for (i, group) in enumerate(groups)
            z, r = _circle(group, cands, spacing / 2)
            (r > 0 && inside(z)) || continue
            w, ne = _circle_winding(f, z, r)
            n += ne
            w === nothing && (state[i] = :unresolved; continue)
            sign(w) == sgn || continue
            push!(loci, (; z, r, w))
            state[i] = :done
        end
        # Retry the leftovers, dropping the outermost group while other candidates crowd the circle.
        left = findall(!=(:done), state)
        for cl in _clusters(map(i -> mean(groups[i]), left), spacing / 2)
            ids = left[cl]
            while length(ids) >= 2
                z, r = _circle(reduce(vcat, groups[ids]), cands, spacing / 2)
                if r == 0
                    deleteat!(ids, argmax(k -> abs(mean(groups[ids[k]]) - z), eachindex(ids)))
                    continue
                end
                inside(z) || break
                w, ne = _circle_winding(f, z, r)
                n += ne
                w === nothing || sign(w) == sgn || break
                push!(loci, (; z, r, w = something(w, 0)))
                state[ids] .= :done
                break
            end
            unres = filter(i -> state[i] == :unresolved, left[cl])
            isempty(unres) || push!(loci, (z = mean(reduce(vcat, groups[unres])), r = 0.0, w = 0))
        end
    end
    return (; loci, nevals = n, info = (; fitconverged = converged))
end

# Circle around the candidates `pts`: centred at their mean, as wide as the other candidates
# allow (see `_reach`) up to `cap`; radius 0 when that leaves too little room around `pts`.
function _circle(pts, cands, cap)
    z = sum(pts) / length(pts)
    r = _reach(z, Iterators.filter(∉(pts), cands), cap)
    return z, r < 1.5 * maximum(p -> abs(p - z), pts) ? zero(r) : r
end

# Single-linkage groups, as index vectors, of points closer than `δ`.
function _clusters(zs, δ)
    groups = Vector{Int}[]
    for (i, z) in enumerate(zs)
        hit = findall(g -> any(j -> abs(zs[j] - z) < δ, g), groups)
        push!(groups, reduce(vcat, groups[hit]; init = [i]))
        deleteat!(groups, hit)
    end
    return groups
end

# A zero or non-finite value is dropped: fitting it would spend the fit's degree on its error.
function _fit(g, Z, alg::AAA)
    y = map(g, Z)
    ok = map(_is_valid, y)
    all(ok) || ((Z, y) = (Z[ok], y[ok]))
    (; fit, converged) = aaa(y, Z; alg.tol, alg.max_degree, alg.stagnation)
    return (; fit, nevals = length(ok), converged)
end

struct Barycentric{Z <: AbstractVector, V <: AbstractVector, W <: AbstractVector}
    nodes::Z
    values::V
    weights::W
end

function (r::Barycentric)(x)
    (; nodes, values, weights) = r
    T = promote_type(typeof(x), eltype(nodes), eltype(values), eltype(weights))
    num = den = zero(T)
    for k in eachindex(nodes, values, weights)
        δ = x - nodes[k]
        iszero(δ) && return T(values[k])
        c = weights[k] / δ
        num += c * values[k]
        den += c
    end
    return num / den
end

# Poles of a barycentric rational = finite eigenvalues of the arrowhead pencil (Klein 2012);
# the infinite ones are the pencil's own, not the fit's. Zero weights leave the pencil singular,
# so drop those nodes first. Its zeros: the same pencil with the numerator's weights.
poles(r::Barycentric) = _pencil_roots(r, r.weights)
_zeros(r::Barycentric) = _pencil_roots(r, r.weights .* r.values)

function _pencil_roots(r, a)
    keep = map(!iszero, r.weights)
    z, a = r.nodes[keep], a[keep]
    m = length(a)
    T = real(eltype(z))
    B = diagm([zero(T); ones(T, m)])
    E = [zero(T) transpose(a); ones(T, m) diagm(z)]
    return filter!(isfinite, complex(eigvals(E, B)))
end

# The Loewner matrix changes by one column and one zeroed row (the new node's) per step, so
# its thin QR is updated rather than refactored and the weights come from inverse iteration
# on R: O(m·n) per step where a fresh SVD costs O(m·n²). Derivation: docs/aaa.md.
function aaa(y, z; tol = 1.0e-13, max_degree = 150, stagnation = 5)
    m = length(z)
    T = promote_type(eltype(z), eltype(y))
    m == 0 && return (; fit = Barycentric(z[1:0], y[1:0], T[]), converged = false)
    fmax = maximum(abs, y)
    # A node's own Loewner row is dropped, so the fit can never outrank the surviving samples.
    N = min(max_degree + 1, m >> 1)
    C = Matrix{T}(undef, m, N)      # Cauchy 1/(zᵢ − σₖ)
    Q = Matrix{T}(undef, m, N)      # Loewner (yᵢ − fσₖ)·C over the test rows = Q·R
    R = zeros(T, N, N)
    a = Vector{T}(undef, m); u = Vector{T}(undef, m); h = Vector{T}(undef, N)
    num = Vector{T}(undef, m); den = Vector{T}(undef, m)
    test = trues(m)
    _, i0 = findmin(abs, y)         # first node: same seed as `approximate`
    σ = [i0]; test[i0] = false
    errs = Float64[]
    ws = Vector{Vector{T}}()        # weights of every step; step s interpolates at σ[1:s]
    for n in 1:N
        j = σ[n]
        _qr_zero_row!(Q, R, j, n - 1, u, h)
        @inbounds for i in 1:m
            δ = z[i] - z[j]
            C[i, n] = iszero(δ) ? 1 / eps() : 1 / δ   # δ=0 only on duplicate samples
            a[i] = test[i] ? (y[i] - y[j]) * C[i, n] : zero(T)
        end
        _qr_append!(Q, R, a, n - 1, h)
        w = _nullvector!(Vector{T}(undef, n), R)

        Cn = view(C, :, 1:n)
        hn = view(h, 1:n)
        hn .= w .* view(y, σ)
        mul!(num, Cn, hn)
        mul!(den, Cn, w)
        err, k = -Inf, 0
        @inbounds for i in 1:m
            test[i] || continue
            e = abs(y[i] - num[i] / den[i])
            e > err && ((err, k) = (e, i))
        end
        push!(errs, err)
        push!(ws, w)

        s = _quitting_check(errs, stagnation, tol, fmax, N)
        # s > 0, unconverged: fall back on the best fit
        if s != 0
            b = s < 0 ? n : s
            return (; fit = Barycentric(z[σ[1:b]], y[σ[1:b]], ws[b]), converged = s < 0)
        end
        push!(σ, k); test[k] = false
    end
    # A single sample: nothing to fit.
    return (; fit = Barycentric(z[σ], y[σ], ones(T, length(σ))), converged = false)
end

# Thin QR `A = Q[:, 1:n]·R[1:n, 1:n]`, updated in place as `A` gains a column or has a row
# zeroed. `h` (length ≥ n) and `u` (length m) are scratch.

# Append column `a` (destroyed): classical Gram–Schmidt, repeated while the remainder
# keeps collapsing — one pass leaves it unorthogonal exactly when `a` is nearly dependent.
function _qr_append!(Q, R, a, n, h)
    Qn = view(Q, :, 1:n); r = view(R, 1:n, n + 1); hn = view(h, 1:n)
    fill!(r, 0)
    nrm = norm(a)
    for _ in 1:3
        mul!(hn, Qn', a)
        mul!(a, Qn, hn, -1, 1)
        r .+= hn
        prev, nrm = nrm, norm(a)
        nrm > prev / 2 && break
    end
    R[n + 1, n + 1] = nrm
    q = view(Q, :, n + 1)
    iszero(nrm) ? fill!(q, 0) : (q .= a ./ nrm)
    return
end

# Zero row `p` of `A`: extend Q by the part `u` of e_p outside its range, then rotate columns
# n…1 against `u` until row p of Q vanishes. The same rotations keep R triangular and gather
# the removed row in a discarded extra row of R (`h`).
function _qr_zero_row!(Q, R, p, n, u, h)
    n == 0 && return
    Qn = view(Q, :, 1:n); hn = view(h, 1:n)
    fill!(u, 0); u[p] = 1
    for _ in 1:2
        mul!(hn, Qn', u)
        mul!(u, Qn, hn, -1, 1)
    end
    τ = norm(u)
    iszero(τ) || (u ./= τ)
    fill!(hn, 0)
    b = real(u[p])
    @inbounds for j in n:-1:1
        q = Q[p, j]
        ρ = sqrt(abs2(q) + b^2)
        iszero(ρ) && continue
        c, s = b / ρ, q / ρ
        for i in axes(Q, 1)
            qi, ui = Q[i, j], u[i]
            Q[i, j] = c * qi - s * ui
            u[i] = conj(s) * qi + c * ui
        end
        Q[p, j] = 0
        for l in j:n
            rj, rh = R[j, l], h[l]
            R[j, l] = c * rj - conj(s) * rh
            h[l] = s * rj + c * rh
        end
        b = ρ
    end
    return
end

# Right singular vector of the smallest singular value of the upper-triangular `R`, by
# inverse iteration on RᴴR. A few sweeps reach a vector whose `‖R·w‖` is the smallest
# singular value to within a small factor, which is all the fit asks of its weights; the
# floored diagonal keeps an exactly singular `R` (a reproduced rational) solvable.
function _nullvector!(w, R)
    n = length(w)
    T = real(eltype(R))
    tiny = max(eps(T) * maximum(i -> abs(R[i, i]), 1:n), sqrt(floatmin(T)))
    d(i) = (x = R[i, i]; abs(x) < tiny ? oftype(x, tiny) : x)
    fill!(w, 1)
    @inbounds for _ in 1:4
        for i in 1:n                 # Rᴴx = w
            s = w[i]
            for l in 1:(i - 1)
                s -= conj(R[l, i]) * w[l]
            end
            w[i] = s / conj(d(i))
        end
        normalize!(w)
        for i in n:-1:1              # Rw = x
            s = w[i]
            for l in (i + 1):n
                s -= R[i, l] * w[l]
            end
            w[i] = s / d(i)
        end
        normalize!(w)
    end
    return w
end

# `RationalFunctionApproximation.quitting_check`, reproduced so the fit stops exactly where
# `approximate` did: -1 converged (keep the last fit), 0 continue, n rewind to fit n.
function _quitting_check(errs, stagnation, tol, fmax, max_iter)
    n = length(errs)
    errs[end] <= tol * fmax && return -1
    stagnant = false
    min_err, min_k = findmin(errs)
    if n >= stagnation + 5 && min_err < fmax / 100
        plateau = 5min_err
        stagnant = all(plateau < e for e in last(errs, stagnation)) || (min_k < n - 2stagnation)
    end
    return (n >= max_iter || stagnant) ? argmin(errs) : 0
end
