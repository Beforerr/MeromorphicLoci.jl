"""
    Polish(; tol = 1e-10, maxiter = 100)

Refines each point to `tol·max(|z|, r)`, `r` the resolution of its estimate. A simple
zero or pole by [`muller`](@ref) on `f` or `1/f`, within `maxiter` steps. A point of
winding `w`, `|w| > 1`, to the winding-weighted mean of the loci it holds (poles count
negatively), `(2πi w)⁻¹ ∮ z f′/f dz`, from 32 to 256 evaluations on a circle: an exact
`|w|`-fold zero or pole to near machine precision, a cluster to its mean.

Where `f` is computed less accurately than `tol` demands (e.g. by quadrature), points land
about its error away, unconverged: loosen `tol` to match.
"""
Base.@kwdef struct Polish
    tol::Float64 = 1.0e-10
    maxiter::Int = 100
end

# -> (z, converged, nevals)
function _polish(alg::Polish, f, z, r, w, mapper)
    r > 0 || return z, false, 0
    abs(w) == 1 || return _moment_mean(f, z, r, w, mapper; alg.tol)
    m = muller(w > 0 ? f : _reciprocal(f), z; h = r / 4, radius = r, alg.tol, alg.maxiter)
    return (isnan(m.root) ? z : m.root), m.converged, m.nevals
end
_polish(::Nothing, f, z, r, w, mapper) = (z, false, 0)

# Mean z₀ = z + δ of the loci of total winding `w` inside the circle z + r·e^{iθ}: with
# f = h·∏(ζ − ζⱼ)^{wⱼ}, h zero-free inside, log f − iwθ has e^{−iθ} Fourier coefficient
# −wδ/r. Trapezoid aliasing from `n` samples decays as (r/d)^n, d the distance to the
# nearest locus outside, perhaps unfound; the estimate from every other sample aliases at
# n/2, so their agreement bounds the error. The circle stays put: recentring on z₀ barely
# speeds convergence ((δ/r)^n is negligible by n = 128 for |δ| < r/2) and can drag it
# onto an unfound locus.
function _moment_mean(f, z, r, w, mapper; tol, nmax = 256)
    n = 32
    θ = (0:n-1) .* (2π / n)
    fs = mapper(f, z .+ r .* cis.(θ))
    ne = n
    while true
        all(_is_valid, fs) || return z, false, ne
        lifted = _lift(log.(complex.(fs)))
        if !isnothing(lifted)
            φ, turns = lifted
            turns == w || return z, false, ne
            g = (φ .- im * w .* θ) .* cis.(θ)
            δ = -r / w * sum(g) / n
            δhalf = -r / w * sum(@view g[1:2:end]) / (n ÷ 2)
            abs(δ - δhalf) <= tol * max(abs(z + δ), r) && return z + δ, true, ne
        end
        2n > nmax && return z, false, ne
        θmid = θ .+ π / n
        fmid = mapper(f, z .+ r .* cis.(θmid))
        ne += n
        θ, fs = vec(permutedims([θ θmid])), vec(permutedims([fs fmid]))
        n *= 2
    end
end

# `log f` on a closed curve with its phase lifted from sample to sample, and the turns it
# makes; `nothing` if a step turns by π/2 or more.
function _lift(l)
    φ = copy(l)
    for k in 2:length(l)
        d = _dlog(l[k-1], l[k])
        abs(imag(d)) < π / 2 || return nothing
        φ[k] = φ[k-1] + d
    end
    d = _dlog(l[end], l[1])
    abs(imag(d)) < π / 2 || return nothing
    return φ, round(Int, imag(φ[end] + d - φ[1]) / 2π)
end
