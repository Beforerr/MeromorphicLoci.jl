# MeromorphicLoci.jl

Find the zeros and poles of a [meromorphic](#branch-cuts) function of complex `z` inside a region — 
at fixed parameters, or as loci `z*(p)` traced over `m ≥ 0` parameters `p`.

```julia
findzeros(f, zdomain, alg)       # (; points, windings, converged, unresolved, nevals, …)
findpoles(f, zdomain, alg)       # same for poles
findzeropoles(f, zdomain, alg)   # both, told apart by the sign of `windings`
findloci(f, zdomain, pdomain, alg)   # loci z*(p)
```

The point finders take `f(z)` and return

- `points` — each distinct zero or pole found, sorted by real then imaginary part;
- `windings` — the winding of `f` around each: its multiplicity, negative for a pole, or
  the count of simple ones closer together than the algorithm resolves;
- `converged` — whether each was polished to the tolerance of the keyword `polish`
  (`Muller()`; `nothing` keeps the algorithm's own estimates). Only points of winding ±1
  are polished: a polish would put fused simple ones on one, and rounding in `f` limits
  its accuracy at a multiple one;
- `unresolved` — places where a zero or pole was detected but its kind could not be told;
- `nevals` — the evaluations of `f` spent;

and fields particular to the algorithm. `findzeropoles` costs about as much as `findzeros`
or `findpoles` alone.

`findloci` takes `f(p...)` returning the function of `z`, and
returns a `Survey` of `Branch`es (below).

Algorithms:

- **`AAA(; nsamples, tol, …)`** — a point finder: rational fit of samples of `f(z)`.
- **`ArgumentPrinciple(; zres, …)`** — a point finder: adaptive cells over `z`, refined
  only near the zeros and poles by Cauchy's argument principle.
- **`Continuation(seed; …)`** — loci over parameters: the point finder `seed` finds them
  at a few parameter knots, and each is followed in between. `findloci` with a point
  finder stands for `Continuation` seeded by it.

## Loci over parameters

`m = 1` parameter → loci are curves `z*(p)`; `m = 2` → surfaces `z*(p₁,p₂)`.

```julia
using MeromorphicLoci

# a zero and a pole winding around helices that graze at p ≈ 0.15 (0.004 apart),
# plus a fixed zero–pole pair 0.03 apart
f(p) = z -> (z - 0.80 * cispi(p)) * (z + 0.10 + 0.12im) /
            ((z - 0.804 * cispi(-p + 0.3)) * (z + 0.13 + 0.12im))
zbox = (-0.98 - 0.98im, 0.98 + 0.98im)  # opposite corners of the z-rectangle
pbox = ((0.0, 1.0),)                    # one (lo, hi) or knot vector per parameter
s = findloci(f, zbox, pbox, AAA())      # Continuation(AAA())
@assert length(s) == 4 && isempty(s.deficit) && all(b -> isempty(b.uncertified), s)
zero_branches = [b for b in s if b.winding > 0]
```

A `Branch` is one locus `z*(p)`: its `Sample`s `(z, p...)` sorted by parameter, each a
polished point of the locus. `b.winding` classifies it — `> 0` zero, `< 0` pole, `0` a
lone sample of unreadable kind — with magnitude the multiplicity. `b.uncertified` lists the stretches where the branch came
too close to another locus to be told from it, and `s.deficit` the knots where the box
holds more zeros minus poles than the branches there.

What is found is what the seeder finds at some knot:

```julia
# the fixed pair is closer than `minsep`: its windings cancel in every cell, and in the count
@assert length(findloci(f, zbox, pbox, ArgumentPrinciple(zres = 0.02))) == 2
```

![zero and pole loci in 3D](docs/src/assets/survey.png)

How it works, what is certified and what it can miss:
[continuation.md](docs/src/continuation.md).
[benchmark/compare.jl](./benchmark/compare.jl) follows the roots of a dispersion relation
over `k` with one survey, and with 128 independent slices by
[RootsAndPoles.jl](https://github.com/fgasdia/RootsAndPoles.jl) (SA-GRPF by Delaunay
triangulation, at fixed parameters).

## Argument principle

```julia
r = findzeropoles(z -> (z - 0.3) * (z + 0.4im)^2 / (z - 0.6), (-1 - 1im, 1 + 1im),
                  ArgumentPrinciple(zres = 0.01))
@assert r.windings == [2, 1, -1] && isempty(r.unresolved)   # −0.4i, 0.3, 0.6
```

How it works: [argument_principle.md](docs/src/argument_principle.md).

## AAA

```julia
f(z) = (z^2 - 1) * (z - 0.5im) / (z - 2)
r = findzeros(f, (-3 - 1im, 3 + 1im), AAA())   # sampled evenly, about 320 points
Z = filter(z -> abs(z) < 1.5, range(-1.5, 1.5, 30) .+ im .* range(-1.5, 1.5, 30)')
r = findzeros(f, Z, AAA())                     # any sample points; loci wherever the fit puts them
```

`findzeros` fits `1/f` by AAA, reads candidate zeros off its poles (and poles of `f` off
its zeros), and confirms each candidate by the winding of `f` on a small circle around it,
which rejects spurious candidates and counts multiplicity (see the `AAA` docstring). `fitconverged = false` means the fit
stopped short of `tol` and may have missed loci.
See [aaa.md](docs/src/aaa.md) for the fit, its cost and what it can miss;
`benchmark/aaa.jl` compares it against RationalFunctionApproximation.jl and
BaryRational.jl.

## Branch cuts

A cut in the region is not detected as such. `ArgumentPrinciple` reports it as one
`unresolved` (winding-0) locus when its phase jump skips a quadrant (`sqrt`), and
misses it silently otherwise (`cbrt`), along with any zero at its branch point.
`AAA` fills it with a string of interleaved zeros and poles.
