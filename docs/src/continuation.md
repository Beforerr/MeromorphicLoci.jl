# [How `Continuation` follows loci](@id continuation)

Write $f(z; p)$ for `f(p)(z)`, $f' = \partial_z f$, $B$ for the box in $z$, and
$p_1 < \dots < p_K$ for the **knots**. `Continuation(seed)` finds the zeros and poles of
$f(\cdot\,; p_k)$ in $B$ with a point finder, the **seeder**, and follows each from knot to knot.

## Following one locus

A simple zero $z^*(p)$ continues analytically while $f'(z^*; p) \ne 0$: until it meets
another zero or a pole, or leaves $B$. A pole is followed as a zero of $1/f$, and a
locus of constant multiplicity $m \ge 2$ as a simple zero of $f/f'$.

Each step extrapolates the last samples (degree $\le 2$ in $p$) to the next $p$ and
corrects with Muller's method. A track ends at the next knot, on $\partial B$, or
**stuck**: no step is accepted however short.

## What is certified

A `Branch` asserts that its samples lie on one locus. Write $z_j$ for the other zeros and
poles of $f$ in $B$, of winding $w_j$. Near a simple zero $z^*$,

$$f(z) = (z - z^*)\,q(z), \qquad d = f'(z^*) = q(z^*) = c \prod_j (z^* - z_j)^{w_j},$$

with $c$ free of zeros and poles in $B$, varying slowly. A track follows the **germ**
$(z^*, d)$: when $z^*$ and a simple $z_j$ exchange places, $z^* - z_j$ changes
sign and $\arg d$ turns by $\pi$. A step to $p$ is accepted if

$$|z^*(p) - \hat z(p)| \le \tau, \qquad |\log d(p) - \log \hat d(p)| \le \tfrac12,$$

with $\hat z$, $\hat d$ extrapolated from earlier samples ($\log \hat d$ linear in
$p$) and $\tau$ set by `zres`. As another locus approaches, $\log d$ turns faster and
steps shorten, and $\tau$ shrinks against the distance to the nearest other locus,

$$R = \left|\frac{2f'}{f''}\right|_{z^*} = \left|\frac{q}{q'}\right|_{z^*}.$$

Below the **floor** $\delta$, `zfloor`, nothing is told apart.
Across a stretch with $R < \delta$, the turn $\Delta = \arg(d_\text{out}/d_\text{in})$
tells what happened:

- $\Delta \approx 0$ or $\pi$: the loci passed or crossed. The track is followed
  through, and the branch lists the stretch in `uncertified`.
- $\Delta \approx \pm\pi/2$: two loci of a kind met at a branch point,
  $z_\pm \approx a \pm b\sqrt{p - p_0}$, and which outgoing locus continues which incoming
  one has no answer. The track ends there, uncertified. Expected where real roots of a
  real $f$ collide and leave as a complex pair.

At each knot, arriving tracks are matched to the seeds. A locus found at any knot is
followed through the knots where the seeder missed it, and one that entered $B$ between
knots is followed back from where it was found.

These checks rest on estimates, not theorems: $d$ and $R$ are leading-order readings,
and an exchange among three loci equally close could leave $d$ unturned. The step
control follows a priori path tracking for polynomial homotopies [^TVV]; $1/R$ is the
leading term of Smale's $\gamma$ in certified tracking [^BL].

[^TVV]: S. Telen, M. Van Barel, J. Verschelde, *A robust numerical path tracking algorithm
    for polynomial homotopy continuation*, SIAM J. Sci. Comput. 42 (2020),
    [arXiv:1909.04984](https://arxiv.org/abs/1909.04984).
[^BL]: C. Beltrán, A. Leykin, *Certified numerical homotopy tracking*, Exp. Math. 21 (2012),
    69–83, [doi:10.1080/10586458.2011.606184](https://doi.org/10.1080/10586458.2011.606184).

## The count

Tracks keep to their own loci beside an unseeded one, but do not report it. What notices
a seeder that came up short is the argument principle at each knot,

$$N(p_k) = \frac{1}{2\pi i} \oint_{\partial B} \frac{f'}{f}\, dz = \#\text{zeros} - \#\text{poles},$$

counted with multiplicity, against $\sum_b w_b$ over the branches at $p_k$. The
[`Survey`](@ref MeromorphicLoci.Survey)'s `deficit` lists the knots where they differ, and those where $N$
cannot be read: a locus on $\partial B$, or $f$ unusable along it. $N$ cannot see a
zero and a pole missing together, nor anything between knots. An `ArgumentPrinciple`
seeder computes $N$ anyway; for others it costs one sampling of $\partial B$ per knot,
skipped with `certify = false`.

## Several parameters

The knots span a grid $K_1 \times \dots \times K_M$. Every node is seeded, and every grid
line, one $p_i$ varying with the others at knots, is followed as above. A locus found at
any node is followed along every line through the nodes it reaches.

A **branch** is a set of node loci connected by tracks along any lines. Two zeros
generically coincide only at isolated points of a parameter plane, and continuation around
such a point swaps them: tracks around a grid cell holding one connect two loci of one
node, and the branch has `degree` 2. A branch of degree 1 is a single-valued $z^*(p)$.

For real $f$ of real parameters, meetings lie along whole curves, where real roots
collide and leave as a complex pair: the branches on either side end on the curve,
uncertified at each grid line that crosses it.

## What it can miss

- A locus in $B$ at no knot: one living only within $(p_k, p_{k+1})$, or entering
  after $p_{K-1}$ and gone by $p_K$.
- A locus the seeder misses at every knot it is in $B$ at. One found at any knot is
  followed to all the others it reaches.
- Loci of noisy $f$: the corrector needs $f$ accurate to about $10^{-10}$ relative,
  and to $10^{-10m}$ at a locus of multiplicity $m$. Tracks then end stuck: keep known
  multiple factors factored. At a branch cut tracks end stuck too.
- $f$ that jumps in $p$, in scale or phase (a piecewise definition, a normalisation
  that switches): $d$ does not continue across the jump, and tracks end stuck there. A
  factor that varies rapidly in $p$ only costs steps: divide it out.
- A knot placed exactly where loci meet cuts their branches there, uncertified. Move the
  knot.
- Over several parameters, what only the inside of a grid cell shows: loci in $B$ only
  there, loci connected only there (separate branches), and two swap points in one cell,
  whose swaps undo each other (two branches of degree 1). A finer grid resolves each.

## Cost

Seeding usually dominates, so the knots are the main cost. They serve only to discover
loci: where none enters $B$ between them, $\{p_\text{lo}, p_\text{hi}\}$ is enough.
Close approaches between loci cost the most steps. Accuracy does not depend on `zres`, and cost only weakly.

`f(p)` is called once per knot and once per step attempted, so setup work inside `f(p)`
is paid per step and per locus.
