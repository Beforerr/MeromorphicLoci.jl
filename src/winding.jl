_is_valid(fz) = isfinite(fz) && !iszero(fz)

# Winding of `f` on the circle |z − z0| = r and the number of evaluations.
# `seed` sets the Nyquist floor: refinement recovers |winding| up to ~3seed/2,
# but a winding ≡ j (mod 2seed) with |j| < seed/2 aliases to j undetectably;
# raise `seed` when higher multiplicities are expected. `nothing` marks an
# unusable sample (0/NaN/Inf) or an unresolvable phase jump on the circle.
function _circle_winding(f, z0, r; seed=8)
    f0 = f(z0 + r)
    ne = 1
    _is_valid(f0) || return nothing, ne
    φ0 = φprev = angle(f0)
    acc = 0.0
    for k in 1:seed
        if k < seed
            fz = f(z0 + r * cispi(2k / seed))
            ne += 1
            _is_valid(fz) || return nothing, ne
            φ = angle(fz)
        else
            φ = φ0  # close the loop on the θ=0 sample; z0 + r*cispi(2) == z0 + r
        end
        Δ, n = _wind_refine(f, z0, r, (k - 1) / seed, k / seed, φprev, φ)
        ne += n
        isnan(Δ) && return nothing, ne
        acc += Δ
        φprev = φ
    end
    return round(Int, acc / 2π), ne
end

# Lifted phase increment of `f` over the arc 2π*[a, b], and evaluations spent.
# Bisect until both half-steps are < π/2: a 2× Nyquist margin
function _wind_refine(f, z0, r, a, b, φa, φb, depth=0)
    m = (a + b) / 2
    fm = f(z0 + r * cispi(2m))
    _is_valid(fm) || return NaN, 1
    φm = angle(fm)
    Δl = rem2pi(φm - φa, RoundNearest)
    Δr = rem2pi(φb - φm, RoundNearest)
    abs(Δl) < π / 2 && abs(Δr) < π / 2 && return Δl + Δr, 1
    depth ≥ 40 && return NaN, 1
    Δ1, n1 = _wind_refine(f, z0, r, a, m, φa, φm, depth + 1)
    isnan(Δ1) && return NaN, n1 + 1
    Δ2, n2 = _wind_refine(f, z0, r, m, b, φm, φb, depth + 1)
    return Δ1 + Δ2, n1 + n2 + 1
end
