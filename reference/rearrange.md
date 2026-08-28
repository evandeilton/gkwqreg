# Monotone rearrangement of crossing quantile curves

Restores monotonicity in the quantile level by sorting each row of the
fitted quantile matrix. This is the rearrangement operator of
Chernozhukov, Fernandez-Val and Galichon (2010), which repairs a
non-monotone estimate of a monotone curve without refitting anything and
without ever making the estimate worse.

## Usage

``` r
rearrange(object, newdata = NULL, ...)
```

## Arguments

- object:

  A `"gkwqregs"` container, or a plain list of `"gkwqreg"` fits at
  different levels. A single `"gkwqreg"` fit is accepted but has nothing
  to repair: a message says so and its quantiles are returned unchanged.

- newdata:

  Optional data frame on which to evaluate the quantiles before
  rearranging. Passed straight to
  [`check_crossing()`](https://evandeilton.github.io/gkwqreg/reference/check_crossing.md),
  so the same rules apply.

- ...:

  Passed to
  [`check_crossing()`](https://evandeilton.github.io/gkwqreg/reference/check_crossing.md);
  in particular `taus` and `tol`.

## Value

A numeric matrix of rearranged quantiles with the same dimensions and
dimnames as the input quantile matrix: one row per evaluation
observation, one column per level, columns named by level. The matrix
carries an attribute `"crossing"` holding the `"gkwq_crossing"` object
computed *before* rearrangement, so that what was repaired remains
recoverable through `attr(x, "crossing")`.

## Details

Fix a row \\i\\ and let \\\hat{Q}\_{i1}, \ldots, \hat{Q}\_{im}\\ be its
fitted quantiles at the increasing levels \\\tau_1 \< \cdots \<
\tau_m\\. The rearranged curve is simply the sorted vector,

\$\$\hat{Q}^{\*}\_{ij} \\=\\ \hat{Q}\_{i(j)},\$\$

the \\j\\th order statistic of the row. Equivalently, and this is the
way to see why it is the right operation, \\\hat{Q}^{\*}\\ is the
quantile function of the random variable \\\hat{Q}(U)\\ obtained by
feeding a uniform level \\U\\ through the unsorted curve: rearrangement
replaces a non-monotone curve by the monotone curve with the same
distribution of values.

## Why sorting is a legitimate repair, not a cosmetic one

Sorting looks like tidying up an embarrassment, and it is natural to
suspect it of hiding a problem rather than solving one. It does not, and
the reason is a theorem rather than a convention.

The estimand \\Q_0(\cdot \mid x)\\ is a genuinely non-decreasing
function of \\\tau\\. Chernozhukov, Fernandez-Val and Galichon (2010)
show that the rearrangement operator is a contraction towards the set of
non-decreasing functions: for every \\p \ge 1\\,

\$\$\left\\ \hat{Q}^{\*} - Q_0 \right\\\_p \\\le\\ \left\\ \hat{Q} - Q_0
\right\\\_p,\$\$

and, under their mild regularity conditions, strictly so whenever
\\\hat{Q}\\ actually crosses. Sorting therefore *weakly dominates* the
unsorted estimate in every \\L_p\\ distance simultaneously, whatever the
true curve happens to be. The intuition is direct: a decreasing stretch
of the estimate cannot be tracking a non-decreasing target, so
exchanging the two values must move both closer to it. Nothing here is
asymptotic, and nothing depends on the estimator that produced
\\\hat{Q}\\.

Two limits are worth stating plainly. First, the guarantee is about the
*quantile curve*, not about the coefficients: after rearranging you hold
a matrix of fitted values, not a model, and the coefficients of the fits
that produced it are unchanged and still non-monotone. Second, monotone
is not the same as correct. Rearrangement removes an internal
contradiction; it does not repair a misspecified family, and it cannot
be read as evidence that the fit is adequate. If a great many rows need
repair, that is a diagnostic about the model, and a single parametric
fit whose implied quantiles cannot cross at all (see
[`check_crossing()`](https://evandeilton.github.io/gkwqreg/reference/check_crossing.md))
is usually the better answer than a large correction.

## References

Chernozhukov, V., Fernandez-Val, I. and Galichon, A. (2010). Quantile
and probability curves without crossing. *Econometrica* **78**,
1093-1125.

## See also

[`check_crossing()`](https://evandeilton.github.io/gkwqreg/reference/check_crossing.md)
for the diagnostic that motivates this repair.

## Examples

``` r
set.seed(1)
n <- 120
x <- runif(n, -1, 1)
mu <- plogis(0.2 + 1.2 * x)
y <- rbeta(n, 4 * mu, 4 * (1 - mu))
d <- data.frame(y = y, x = x)

fits <- gkwqreg(y ~ x | x, data = d, tau = seq(0.05, 0.95, by = 0.05),
                family = "kw")

## Evaluate on a grid that extrapolates well beyond the observed x, which is
## where independently fitted levels are least constrained.
grid <- data.frame(x = seq(-3, 3, length.out = 201))
R <- rearrange(fits, newdata = grid)

## What was repaired is kept with the result.
cr <- attr(R, "crossing")
c(rows_repaired = cr$n_crossing, of = nrow(R), worst = cr$worst)
#> rows_repaired            of         worst 
#>  4.400000e+01  2.010000e+02  6.146179e-03 
#> rows_repaired            of         worst
#>  4.400000e+01  2.010000e+02  6.146179e-03

## Monotone in tau afterwards, row by row, by construction.
all(apply(R, 1, function(r) all(diff(r) >= 0)))
#> [1] TRUE
#> TRUE

## The repair is local: only the rows that crossed were touched.
Q <- cr$Q
identical(which(rowSums(abs(R - Q)) > 0), cr$which)
#> [1] TRUE
#> TRUE

## And it is small -- the largest fitted value moved by sorting.
max(abs(R - Q))
#> [1] 0.0100315
#> 0.0100315

## A single fit has nothing to rearrange and says so.
fit <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
Q1 <- rearrange(fit, taus = c(0.1, 0.5, 0.9))
#> a single fit has no crossing to rearrange: its quantiles are those of one distribution function. Returning them unchanged.
#> a single fit has no crossing to rearrange: its quantiles are those of one
#> distribution function. Returning them unchanged.
head(round(Q1, 4), 3)
#>         0.1    0.5    0.9
#> [1,] 0.1478 0.3856 0.6565
#> [2,] 0.1788 0.4584 0.7488
#> [3,] 0.2439 0.5977 0.8851
```
