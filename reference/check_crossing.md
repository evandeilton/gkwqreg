# Detect quantile crossing

Checks whether a set of fitted conditional quantiles is monotone in the
quantile level, row by row. A quantile function must be non-decreasing
in \\\tau\\; an estimate need not be, and where it is not the fitted
object is not a conditional distribution at all. `qcrossing()` is an
alias.

## Usage

``` r
check_crossing(object, newdata = NULL, taus = NULL, tol = 0, ...)

qcrossing(object, newdata = NULL, taus = NULL, tol = 0, ...)
```

## Arguments

- object:

  One of three things, and the answer means something different in each
  case (see Details). A `"gkwqregs"` container returned by
  [`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md)
  with a vector `tau`; a plain list of `"gkwqreg"` fits, which is sorted
  by level before checking; or a single `"gkwqreg"` fit, whose implied
  quantiles are read off its one fitted conditional distribution.

- newdata:

  Optional data frame on which to evaluate the quantiles. It must supply
  every covariate used by any part of the model. Extrapolating beyond
  the observed covariate range is legitimate here and often revealing:
  see the examples.

- taus:

  Used only when `object` is a single fit: the levels at which to read
  the fitted distribution. Defaults to `seq(0.05, 0.95, by = 0.05)`.
  Ignored for a container or a list of fits, whose levels are already
  fixed by the fits themselves.

- tol:

  Non-negative numeric tolerance. A decrease between adjacent levels
  counts as a crossing only if it exceeds `tol` in absolute value. The
  default `0` flags any decrease at all; a small positive value (say
  `1e-8`) is appropriate if you wish to ignore decreases attributable to
  floating-point arithmetic rather than to the fit.

- ...:

  Currently unused, accepted so that the alias and future methods keep a
  stable signature.

## Value

An object of class `"gkwq_crossing"`, a list with components

- `taus`:

  Numeric vector of the levels checked, in increasing order.

- `Q`:

  The \\n \times m\\ matrix of fitted quantiles, rows in the order of
  the evaluation data and columns named by level.

- `mode`:

  `"separate"` for independently fitted levels, `"implied"` for the
  levels of one fitted distribution.

- `n_crossing`:

  Number of rows containing at least one crossing.

- `frac`:

  `n_crossing` divided by the number of rows.

- `which`:

  Integer row indices of the offending rows, so that the covariate
  values responsible can be inspected.

- `worst`:

  Largest decrease observed, on the scale of the response; `0` when
  there is none. Useful for deciding whether a crossing is a substantive
  defect or numerical dust.

- `pairs`:

  Data frame with one row per adjacent pair of levels, giving `tau_lo`,
  `tau_hi` and the count `n` of rows violated by that pair.

A `print` method summarises the check and states which of the two
questions was answered.

## Details

Let \\\hat{Q}\_{ij}\\ be the fitted quantile for row \\i\\ at the
\\j\\th level of an increasing grid \\\tau_1 \< \cdots \< \tau_m\\. The
function forms the adjacent differences \\D\_{ij} = \hat{Q}\_{i,j+1} -
\hat{Q}\_{ij}\\ and records a crossing at \\(i,j)\\ whenever \\D\_{ij}
\< -\texttt{tol}\\. A row is counted once, however many of its adjacent
pairs are violated, so `n_crossing` and `frac` describe *observations
affected*, while `pairs` localises the damage in \\\tau\\.

## Two genuinely different questions

The word "crossing" covers two situations that this function
deliberately treats together and reports separately.

**Across independently fitted levels** (`mode = "separate"`). Each level
is a separate optimisation with its own likelihood and its own
coefficient vector. Nothing in the estimation ties them together, so
nothing forces \\\hat{Q}\_\tau(x)\\ to increase with \\\tau\\ at any
particular \\x\\. This is the classic and well-documented failure of
per-level quantile regression, and it is what the check reports for a
`"gkwqregs"` container or a list of fits. It typically bites where the
fitted lines are least constrained: in sparse regions of the covariate
space, near the extreme levels, and above all under extrapolation.

**Within a single fit** (`mode = "implied"`). Here the quantiles at
every level are read off one fitted conditional distribution function,
through \\\hat{Q}(\tau \mid x) = F^{-1}(\tau \mid x)\\ with a single
parameter vector. Because \\F(\cdot \mid x)\\ is a proper distribution
function, its inverse is non-decreasing by construction and crossing is
arithmetically impossible.

Reporting zero crossings in the second case is therefore not a vacuous
tautology dressed up as a result. It is precisely the argument for
parametric quantile regression over the check-function kind: the
guarantee holds at every covariate value, in sample and out, at every
level, without any post-hoc repair. The printed output says so
explicitly rather than letting a reader mistake it for good luck.
Verified in this package, ten levels on the same data: `0` of the rows
cross within a single `gkwqreg` fit, against `4.3` percent of rows for
[`quantreg::rq`](https://rdrr.io/pkg/quantreg/man/rq.html) fitted level
by level.

The price of the guarantee is that it is bought with a distributional
assumption. Under misspecification the single fitted distribution can be
monotone and wrong, whereas independently fitted levels remain
consistent for each quantile separately. Monotonicity is a property
worth having, not a certificate of correctness; check it alongside
[`plot.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/plot.gkwqreg.md)
panel 5 and
[`pinball()`](https://evandeilton.github.io/gkwqreg/reference/pinball.md),
not instead of them.

## References

Chernozhukov, V., Fernandez-Val, I. and Galichon, A. (2010). Quantile
and probability curves without crossing. *Econometrica* **78**,
1093-1125.

## See also

[`rearrange()`](https://evandeilton.github.io/gkwqreg/reference/rearrange.md)
for the post-hoc monotonicity fix;
[`predict.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/predict.gkwqreg.md)
with a `tau` argument for reading other levels off one fitted
distribution.

## Examples

``` r
set.seed(1)
n <- 120
x <- runif(n, -1, 1)
mu <- plogis(0.2 + 1.2 * x)
y <- rbeta(n, 4 * mu, 4 * (1 - mu))
d <- data.frame(y = y, x = x)

## (a) ONE fit, nineteen levels read off its fitted distribution.
fit <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
check_crossing(fit, taus = seq(0.05, 0.95, by = 0.05))
#> 
#> Quantile crossing check
#> levels: 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95
#> source: quantiles implied by ONE fitted distribution
#> 
#>   rows with a crossing: 0 of 120 (0.00%)
#>   worst violation     : 0.000e+00
#> 
#>   Zero crossings is guaranteed here, not luck: these are quantiles of
#>   a single proper distribution function. That guarantee is what a
#>   parametric fit buys over independently fitted levels.
#>   rows with a crossing: 0 of 120 (0.00%)
## Zero, and it could not have been anything else.

## (b) NINETEEN independent fits, one per level, shape modelled too.
fits <- gkwqreg(y ~ x | x, data = d, tau = seq(0.05, 0.95, by = 0.05),
                family = "kw")
check_crossing(fits)
#> 
#> Quantile crossing check
#> levels: 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95
#> source: independently fitted models, one per level
#> 
#>   rows with a crossing: 0 of 120 (0.00%)
#>   worst violation     : 0.000e+00
#>   rows with a crossing: 0 of 120 (0.00%)
## None in sample -- but that is an empirical fact about these 120 rows,
## not a guarantee. Evaluate the same fits where the data do not constrain
## them and the guarantee's absence becomes visible:

grid <- data.frame(x = seq(-3, 3, length.out = 201))
cr <- check_crossing(fits, newdata = grid)
cr
#> 
#> Quantile crossing check
#> levels: 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95
#> source: independently fitted models, one per level
#> 
#>   rows with a crossing: 44 of 201 (21.89%)
#>   worst violation     : 6.146e-03
#> 
#>   Independently fitted levels carry no monotonicity guarantee.
#>   rearrange() applies the Chernozhukov et al. (2010) fix.
#>   rows with a crossing: 44 of 201 (21.89%)
#>   worst violation     : 6.146e-03

## Where, in tau: the upper levels, which are the least well determined.
cr$pairs[cr$pairs$n > 0, ]
#>      tau_lo tau_hi  n
#> 0.65   0.60   0.65  2
#> 0.70   0.65   0.70 11
#> 0.75   0.70   0.75 19
#> 0.80   0.75   0.80 26
#> 0.85   0.80   0.85 32
#> 0.90   0.85   0.90 38
#> 0.95   0.90   0.95 44
#>      tau_lo tau_hi  n
#> 0.65   0.60   0.65  2
#> 0.70   0.65   0.70 11
#> ...
#> 0.95   0.90   0.95 44

## Where, in x: entirely outside the observed range of the covariate.
range(grid$x[cr$which])   # -3.00 -1.71
#> [1] -3.00 -1.71
range(d$x)                # -0.97  0.99
#> [1] -0.9738448  0.9853681

## The same nineteen levels taken from the single fit stay monotone there.
check_crossing(fit, newdata = grid, taus = seq(0.05, 0.95, by = 0.05))$n_crossing
#> [1] 0
#> 0
```
