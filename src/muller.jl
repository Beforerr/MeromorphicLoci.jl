"""
    muller(f, x0; h, radius = Inf, tol = 1e-10, maxiter = 100)
    muller(f, x0, x1, x2; tol = 1e-10, maxiter = 100)

Muller's method for a zero of `f`, from a guess `x0` uncertain by about `h` or from three
starting points. Returns `(; root, converged, nevals)`: `root` is the smallest-`|f|` point
evaluated (`NaN` if none was finite), a zero when `converged`, i.e. a step fell below
`tol·max(|root|, h)` (`h` the starting points' spread). Steps stay within `radius` of `x0`;
keep `h` and `radius` below the distance to any other zero, or it may converge there.

Where noise or rounding in `f` limits accuracy above `tol`, as at a multiple zero, `root`
lands about that limit away, typically unconverged.
"""
muller(f, x0; h, radius = Inf, kw...) = _muller(f, x0 - h, x0 + h * im, x0, x0, radius; kw...)
muller(f, x0, x1, x2; kw...) = _muller(f, x0, x1, x2, x2, Inf; kw...)

function _muller(f, x0, x1, x2, center, radius; tol = 1.0e-10, maxiter = 100)
    x0, x1, x2 = promote(complex(x0), complex(x1), complex(x2))
    f0, f1, f2 = f(x0), f(x1), f(x2)
    a0, a1, a2 = abs(f0), abs(f1), abs(f2)
    nevals = 3
    # Step scale and |f| scale from the middle value, so a start on a zero
    # does not set it.
    width = max(abs(x1 - x0), abs(x2 - x1), abs(x2 - x0))
    mid = max(min(a0, a1), min(max(a0, a1), a2))
    best, abest = typeof(x2)(NaN, NaN), oftype(a2, Inf)
    for (x, a) in ((x0, a0), (x1, a1), (x2, a2))
        a < abest && ((best, abest) = (x, a))
    end
    iszero(abest) && return (; root = best, converged = true, nevals)
    stall = 0
    for _ in 1:maxiter
        h1 = x1 - x0
        h2 = x2 - x1
        (h1 == 0 || h2 == 0 || h1 + h2 == 0) && return (; root = best, converged = false, nevals)
        delta1 = (f1 - f0) / h1
        delta2 = (f2 - f1) / h2
        a = (delta2 - delta1) / (h2 + h1)
        b = a * h2 + delta2
        c = f2

        disc = sqrt(b^2 - 4a * c)
        # The larger of |b ± disc|: |b + disc|² − |b − disc|² = 4 Re(b̄·disc)
        denom = real(conj(b) * disc) > 0 ? b + disc : b - disc
        denom == 0 && return (; root = best, converged = false, nevals)

        dx = -2c / denom
        step = abs(dx)
        while abs(x2 + dx - center) > radius
            dx /= 2
        end
        x3 = x2 + dx
        f3 = f(x3)
        a3 = abs(f3)
        nevals += 1
        # Back off only from huge f, which overflows the divided differences. Backing off from
        # any rise in |f| would, on a noise floor, bunch the iterates and fake convergence.
        flimit = 1.0e6 * max(a0, a1, a2)
        for _ in 1:12
            (isfinite(f3) && a3 <= flimit) && break
            dx /= 2
            x3 = x2 + dx
            f3 = f(x3)
            a3 = abs(f3)
            nevals += 1
        end
        isfinite(f3) || return (; root = best, converged = false, nevals)
        a3 < abest ? ((best, abest, stall) = (x3, a3, 0)) : (stall += 1)
        # A small step alone is no zero: decay toward tiny |f|, 
        # or seeds decades apart in |f|, shrink it too. 
        # Require |f| no larger than a zero `xtol` away would give at
        # the starting slope mid/width.
        xtol = tol * max(abs(x3), width)
        if iszero(f3) || (step <= xtol && a3 * width <= mid * xtol)
            return (; root = x3, converged = true, nevals)
        end
        # Steps without improving on the best |f|, seeds included, mean a noise floor or
        # decay to one, which shrinks steps as a zero does: give up with `best`.
        stall >= 5 && return (; root = best, converged = false, nevals)

        x0, x1, x2 = x1, x2, x3
        f0, f1, f2 = f1, f2, f3
        a0, a1, a2 = a1, a2, a3
    end
    return (; root = best, converged = false, nevals)
end

"""
    Muller(; tol = 1e-10, maxiter = 100)

Polish for the simple zeros and poles using [`muller`](@ref) from each estimate.
"""
Base.@kwdef struct Muller
    tol::Float64 = 1.0e-10
    maxiter::Int = 100
end

# Unconverged, `root` still has no larger |g| than `z`, one of the seeds.
function _polish(alg::Muller, g, z, radius)
    r = muller(g, z; h = radius / 4, radius, alg.tol, alg.maxiter)
    return (isnan(r.root) ? z : r.root), r.converged, r.nevals
end
_polish(::Nothing, g, z, radius) = (z, false, 0)
