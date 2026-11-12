# The AAA fit

`AAA()` finds zeros of $f$ from samples, without seeds. It samples $f$ on a grid
$z_1,\dots,z_m$, fits a rational function to $y_i = 1/f(z_i)$, and reads the zeros of $f$
off as the poles of the fit. This page describes the fit and what it costs.

## The AAA iteration

AAA [1] writes the fit in barycentric form over $n$ *support points*
$\sigma_1,\dots,\sigma_n$ chosen among the samples,

```math
r(z) = \frac{\sum_{k=1}^{n} \dfrac{w_k\, y_{\sigma_k}}{z - z_{\sigma_k}}}
            {\sum_{k=1}^{n} \dfrac{w_k}{z - z_{\sigma_k}}} ,
```

a rational function of degree $n-1$ that interpolates $y$ at the support points for any
nonzero weights $w$. The weights are free to make $r$ match the other samples. Multiplying
$y_i - r(z_i)$ by the denominator makes the mismatch linear in $w$: it is $(Lw)_i$ with the
Loewner matrix

```math
L_{ik} = \frac{y_i - y_{\sigma_k}}{z_i - z_{\sigma_k}}, \qquad i \notin \sigma ,
```

so $w$ is taken as the unit vector minimizing $\lVert Lw\rVert$, the right singular vector
of the smallest singular value. The iteration is greedy: the sample where $|y_i - r(z_i)|$
is largest becomes the next support point, until the mismatch falls below
`tol`$\cdot\max|y|$. The poles are the zeros of the denominator, the finite eigenvalues of
an $(n+1)\times(n+1)$ arrowhead pencil.

Computing $w$ from a fresh SVD at step $n$ costs $O(m n^2)$, so a fit that ends at $N$
support points costs $O(m N^3)$. A box holding dozens of zeros needs $N \sim 100$, where the
SVDs can outweigh the evaluations of $f$.

## Updating instead of refactoring

Between consecutive steps $L$ changes in two small ways. Choosing sample $p$ as the new
support point

1. removes row $p$, which is a rank-one change, and
2. appends one column, the Loewner column of $p$.

A QR factorization can follow both changes in $O(mn)$ [2],
[3, §6.5]; Lietaert et al. [7] update AAA's Loewner QR the same way. Keep $L = QR$ with $Q \in \mathbb{C}^{m\times n}$ having
orthonormal columns and $R$ upper triangular, and let removed rows stay in place as zero
rows, so no index bookkeeping is needed.

**Zeroing row $p$.** Let $q = Q^{H} e_p$ and $u = (e_p - Qq)/\tau$ with
$\tau = \lVert e_p - Qq\rVert$. Then $[\,Q\ \ u\,]$ has orthonormal columns, its range
contains $e_p$, and its row $p$ is $[\,q^{H}\ \ \tau\,]$. Apply $n$ Givens rotations, each
mixing column $j$ with the last column for $j = n,\dots,1$, chosen to move row $p$ onto
$[\,0\ \cdots\ 0\ \ 1\,]$. Afterwards the last column is $e_p$ itself, since it has unit
norm and a unit entry, and the first $n$ columns $Q'$ vanish in row $p$. The same rotations
applied to the rows of $\begin{bmatrix} R \\ 0 \end{bmatrix}$ give
$\begin{bmatrix} R' \\ v^{T} \end{bmatrix}$ with $R'$ still upper triangular, because the
order $j = n,\dots,1$ only fills the extra row. So

```math
L = Q'R' + e_p v^{T} ,
```

where $v^T$ is the old row $p$ of $L$, and $Q'R'$ is $L$ with that row zeroed. Each
rotation touches two columns of length $m$ and two rows of length at most $n$: $O(mn)$ in
total.

**Appending a column.** Orthogonalize the new column $a$ against $Q$ by classical
Gram–Schmidt, $r = Q^{H}a$, $a \leftarrow a - Qr$; the coefficients $r$ and the remaining
norm form the new column of $R$. Near convergence the new column is almost in the range of
$Q$, which is exactly when one pass loses orthogonality, so the pass is repeated while the
remainder keeps collapsing; two passes are enough [4]. Cost $O(mn)$.

## Weights by inverse iteration

The smallest right singular vector of $L$ is that of $R$, and for a triangular matrix it
comes cheaply from inverse iteration on $R^{H}R$ [5]: each sweep
$w \leftarrow R^{-1}R^{-H}w$, then normalize, is two triangular solves, $O(n^2)$. A sweep
damps every other singular direction by $(\sigma_n/\sigma_j)^2$. Convergence *to the
singular vector* is slow when $\sigma_{n-1}\approx\sigma_n$, but the fit does not need the
vector: it needs $\lVert Rw\rVert$ small, and directions that survive the damping have
$\sigma_j \approx \sigma_n$ by construction. A fixed handful of sweeps therefore gives
weights as good as the SVD's.

The two triangular solves work on $R$ and $R^{H}$ separately rather than on the Gram matrix
$L^{H}L$. That matters: a converged fit has $\sigma_n/\sigma_1 \sim$ `tol` $= 10^{-13}$,
and factoring $L^{H}L$ would square it past machine precision. Updating a Cholesky factor of
the Gram matrix is the route taken by FastAAA [6], which suits looser
tolerances than a zero search can afford.

## Scaling

With $m$ samples and a fit ending at $N$ support points:

| stage | refactoring every step | updating |
|:--|:--|:--|
| Loewner factorization | $O(mN^3)$ | $O(mN^2)$ |
| weights | $O(N^4)$ | $O(N^3)$ |
| mismatch at the samples | $O(mN^2)$ | $O(mN^2)$ |
| poles, once per fit | $O(N^3)$ | $O(N^3)$ |

The fit drops from $O(mN^3)$ to $O(mN^2)$, a factor that grows with the number of roots in
the box. What remains at large $N$ is the pole
eigenproblem.

## What the fit can see

A zero $z_0$ of $f$ is a pole of $1/f$ with residue $1/f'(z_0)$. The fit resolves it only if
that residue stands above `tol`$\cdot\max|1/f|$ at the samples. Zeros where $|f'|$ is
exponentially large are found or missed depending on details of the sampling, so a zero set
is reproducible only above this threshold.

## References

1. Y. Nakatsukasa, O. Sète, L. N. Trefethen, *The AAA algorithm for
  rational approximation*, SIAM J. Sci. Comput. **40**, A1494 (2018).
2. J. W. Daniel, W. B. Gragg, L. Kaufman, G. W. Stewart, *Reorthogonalization
  and stable algorithms for updating the Gram–Schmidt QR factorization*, Math. Comp.
  **30**, 772 (1976).
3. G. H. Golub, C. F. Van Loan, *Matrix Computations*, 4th ed., Johns Hopkins
  University Press (2013).
4. L. Giraud, J. Langou, M. Rozložník, J. van den Eshof, *Rounding error
  analysis of the classical Gram–Schmidt orthogonalization process*, Numer. Math. **101**,
  87 (2005).
5. I. C. F. Ipsen, *Computing an eigenvector with inverse iteration*, SIAM Rev.
  **39**, 254 (1997).
6. A. Hochman, *FastAAA: A fast rational-function fitter*, IEEE EPEPS (2017).
7. P. Lietaert, J. Pérez, B. Vandereycken, K. Meerbergen, *Automatic rational
  approximation and linearization of nonlinear eigenvalue problems*, IMA J. Numer. Anal.
  **42**, 1087 (2022).
