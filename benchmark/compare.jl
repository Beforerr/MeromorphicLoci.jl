# One (ω,k) survey vs the naive per-k GRPF slicing, on a cold electron–proton parallel
# dispersion (Stix R/L transverse branches). Both cover the same box in (ω, k): the survey
# finds the roots at a few k and follows them in between, while per-slice GRPF re-pays a
# full 2-D mesh at every k. Run: julia --project=benchmark benchmark/compare.jl
using MeromorphicLoci
using RootsAndPoles
using Printf

# --- cold e-p parallel dispersion, Ω_ref = |Ωe| = 1 --------------------------
const Ωe, Ωi = -1.0, 1 / 1836.15
const P2e = 10.0
const P2i = P2e / 1836.15
R(ω) = 1 - P2e / (ω * (ω - Ωe)) - P2i / (ω * (ω - Ωi))
L(ω) = 1 - P2e / (ω * (ω + Ωe)) - P2i / (ω * (ω + Ωi))
f(ω, k) = (R(ω) - (k / ω)^2) * (L(ω) - (k / ω)^2)   # transverse: R=n² and L=n²

const BOX = ((0.05 - 0.15im, 4.05 + 0.15im), (1.0, 3.0))   # (ω, k)
const ZTOL = 0.02
const NSLICE = 128

# --- one survey: roots found at the knots in k and followed between them -----
function run_survey()
    s = findloci(k -> ω -> f(ω, k), BOX[1], Base.tail(BOX), Continuation(AAA(); zres=ZTOL))
    return (; nbranch=length(s), nzero=count(b -> b.winding > 0, s),
        nevals=s.nevals)
end

# --- naive baseline: independent 2-D GRPF at each k --------------------------
function run_slices(; nslice=NSLICE)
    ll, ur = BOX[1]
    ks = range(BOX[2][1], BOX[2][2]; length=nslice)
    evals = 0
    nroot = 0
    coords = RootsAndPoles.rectangulardomain(ll, ur, ZTOL)
    mesh = ComplexMesh(coords)
    for k in ks
        calls = Ref(0)
        fk = ω -> (calls[] += 1; f(ω, k))
        roots, _ = rootsandpoles(fk, mesh, FinderParams(tol=ZTOL);
            refinement_mode=:regular)
        nroot += length(roots)
        evals += calls[]
    end
    return (; nevals=evals, nroot_total=nroot, nslice)
end

# --- warm up, then time -------------------------------------------------------
run_survey()
run_slices()

t_sv = @elapsed sv = run_survey()
t_sl = @elapsed sl = run_slices()

@printf("\n%-18s %12s %12s %10s\n", "method", "evals", "wall (s)", "branches")
@printf("%-18s %12d %12.3f %10d  (%d root-curves)\n",
    "findloci", sv.nevals, t_sv, sv.nbranch, sv.nzero)
@printf("%-18s %12d %12.3f %10s  (%d roots over %d slices)\n",
    "per-k GRPF", sl.nevals, t_sl, "—", sl.nroot_total, sl.nslice)
@printf("\nspeedup: %.1f× evals, %.1f× wall\n",
    sl.nevals / sv.nevals, t_sl / t_sv)
