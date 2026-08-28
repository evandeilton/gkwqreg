# Predictions from a fixed-level quantile regression fit

Evaluates a fitted
[`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md)
model at the estimation data or at `newdata`. The first part of the
model *is* the conditional quantile, so the default prediction is the
conditional \\\tau\\-quantile itself and never a conditional mean.
Densities, distribution functions and moments are all derived from the
same fitted conditional distribution.

## Usage

``` r
# S3 method for class 'gkwqreg'
predict(
  object,
  newdata = NULL,
  type = c("quantile", "response", "link", "parameter", "mu", "density", "probability",
    "terms", "mean", "variance"),
  at = NULL,
  tau = NULL,
  elementwise = NULL,
  na.action = stats::na.pass,
  ...
)
```

## Arguments

- object:

  A `"gkwqreg"` fit, as returned by
  [`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md).

- newdata:

  Optional data frame in which to evaluate the model. It must contain
  every variable used by *every* part of the formula, not only the
  quantile part. Factor levels and contrasts recorded at fit time are
  reapplied, so a factor whose levels are incompletely represented in
  `newdata` still produces columns consistent with the fit. Defaults to
  the estimation data.

- type:

  The quantity to predict; see the table below. `"response"` is an alias
  for `"quantile"`, retained so that code written against other
  regression classes keeps running – but note that here the "response
  scale" prediction is a quantile, not a mean.

- at:

  Numeric vector of points in `(0,1)` at which the density or the
  distribution function is evaluated. Required for `type = "density"`
  and `type = "probability"`, and ignored otherwise.

- tau:

  Optional numeric vector of quantile levels in `(0,1)`, used only by
  `type = "quantile"` and `type = "response"`. **These are further
  quantiles of the same fitted conditional distribution, not a refit**;
  see the section below. `NULL` (the default) returns the level the
  model was fitted at.

- elementwise:

  Logical, governing how `at` is combined with the rows. `TRUE` pairs
  the `i`th element of `at` with the `i`th row, recycling `at` if
  necessary, and returns a vector; `FALSE` crosses every element of `at`
  with every row and returns a matrix. The default is `TRUE` when `at`
  has exactly one element per row and there is more than one row, and
  `FALSE` otherwise.

- na.action:

  How to treat missing values in `newdata`. The default,
  [`stats::na.pass()`](https://rdrr.io/r/stats/na.fail.html), retains
  the rows and propagates `NA` into the predictions, so the result stays
  aligned with the rows of `newdata`.

- ...:

  Currently unused; present for consistency with the generic.

## Value

The shape depends on `type`, and for `"quantile"` also on `tau`.

- `"quantile"`, `"response"`, `"mu"`: a numeric vector with one element
  per row of `newdata`, or per observation used in the fit. If `tau` is
  supplied to `"quantile"`/`"response"`, an `n` by `length(tau)` matrix
  whose columns are named `"tau=<level>"`.

- `"link"`: a data frame with `n` rows and one column per model part,
  named after the parts (`mu` first), holding the linear predictors.

- `"parameter"`: a data frame with `n` rows and the five columns
  `alpha`, `beta`, `gamma`, `delta`, `lambda`, in `gkwdist` order.

- `"mean"`, `"variance"`: numeric vectors of length `n`.

- `"density"`, `"probability"`: an `n` by `length(at)` matrix whose
  columns are named after `at`, or a numeric vector of length `n` when
  `elementwise = TRUE`.

- `"terms"`: a named list with one matrix per model part. Each matrix
  has `n` rows and one column per non-intercept term of that part; a
  part with no covariates contributes a matrix with zero columns.

Applying [`predict()`](https://rdrr.io/r/stats/predict.html) to the
`"gkwqregs"` container returned by a vector-valued `tau` predicts from
each fit in turn, and binds the results column-wise, named by level,
when every level returned a plain vector of the same length.

## Details

Write \\g\\ for the link of the quantile part (`"logit"` by default) and
\\o_i\\ for its offset. The default prediction is \$\$\widehat{Q}(\tau
\mid x_i) \\=\\ \hat\mu\_\tau(x_i) \\=\\
g^{-1}\\\left(x_i^{\top}\hat\beta\_\mu + o_i\right).\$\$

Every other `type` is obtained by first reconstructing, row by row, the
full parameter vector \\(\alpha, \beta, \gamma, \delta, \lambda)\\ of
the Generalized Kumaraswamy distribution: the non-anchored parameters
come from their own linear predictors and links, and the anchored
parameter is solved for so that the resulting distribution has
\\\hat\mu\_\tau(x_i)\\ as its exact \\\tau\\-quantile. A useful
consequence is that \$\$F\\\left(\widehat{Q}(\tau \mid x_i) \mid
x_i\right) = \tau\$\$ holds to machine precision for every row, which is
the cheapest available check that a prediction pipeline has not silently
drifted onto some other scale (see the examples).

## Values of `type`

|  |  |  |
|----|----|----|
| `type` | What is returned | Typical use |
| `"quantile"` (default), `"response"` | the conditional tau-quantile; with `tau` supplied, several quantiles of the same fitted distribution | the headline prediction |
| `"mu"` | the conditional tau-quantile, always ignoring `tau` | when the fitted level is wanted inside code that passes `tau` around |
| `"link"` | the linear predictors of every part, each on its own link scale | diagnostics, and plotting on the scale the model is linear in |
| `"parameter"` | the reconstructed five distribution parameters per row | handing the fitted conditional law to `gkwdist` functions |
| `"mean"` | the conditional mean, by 64-node Gauss-Legendre quadrature of the quantile function | comparison with a mean-parametrized fit |
| `"variance"` | the conditional variance, by the same quadrature | assessing conditional dispersion |
| `"density"` | the conditional density evaluated at `at` | likelihood displays, simulated envelopes |
| `"probability"` | the conditional distribution function evaluated at `at` | exceedance probabilities, probability-integral-transform checks |
| `"terms"` | each term's additive contribution to its part's linear predictor | partial-effect plots |

`"mean"` and `"variance"` integrate \\\int_0^1 Q(u)^k \\ du\\ rather
than the density, because \\Q\\ is available in closed form and has no
boundary singularity to work around. They are provided so that a
quantile fit can be put beside a mean-parametrized one; reporting the
conditional mean as though it were the estimand of this model defeats
the point of the parametrization.

## Reading other levels off the fit versus refitting

Supplying `tau` returns quantiles of the **same fitted conditional
distribution** read at other probabilities. This is not the same object
as `gkwqreg(..., tau = tau_new)`, which solves a different estimating
problem and in general returns different coefficients.

- Reading off is a *distributional extrapolation*. It borrows strength
  from the parametric family, costs nothing beyond a call to the
  quantile function, and is monotone in `tau` by construction, so
  predictions read off one fit can never cross.

- Refitting is a *level-specific* estimate. It targets the new level
  directly, and is the honest thing to report when the family may be
  wrong away from the level already fitted.

The two coincide only if the family is correctly specified. The gap
between them is therefore a specification diagnostic in its own right: a
small discrepancy far from the fitted `tau` supports the shape
assumption, and a large one says the parametric tail, rather than the
systematic part of the model, is doing the work. See
[`quantile_process()`](https://evandeilton.github.io/gkwqreg/reference/quantile_process.md)
for the systematic version of this comparison and
[`check_crossing()`](https://evandeilton.github.io/gkwqreg/reference/check_crossing.md)
for what can go wrong once several levels are fitted separately.

## See also

[`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md)
for the model,
[`fitted.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/fitted.gkwqreg.md)
for the in-sample quantiles,
[`marginal_effects()`](https://evandeilton.github.io/gkwqreg/reference/marginal_effects.md)
for effects on the quantile scale,
[`quantile_process()`](https://evandeilton.github.io/gkwqreg/reference/quantile_process.md)
and
[`check_crossing()`](https://evandeilton.github.io/gkwqreg/reference/check_crossing.md)
for several levels at once.

## Examples

``` r
## A Kumaraswamy sample whose conditional MEDIAN follows a logit model.
## Drawn by inverse transform through gkwq_quantile(), so the data-generating
## process is exactly the one the model assumes.
set.seed(2024)
n   <- 300
x1  <- runif(n, -2, 2)
x2  <- rbinom(n, 1, 0.5)
mu  <- plogis(0.3 + 0.9 * x1 - 0.5 * x2)      # true conditional median
bt  <- log1p(-0.5) / log1p(-mu^2)             # beta anchoring Q(0.5) = mu
y   <- gkwq_quantile(runif(n), alpha = 2, beta = bt)
dat <- data.frame(y = y, x1 = x1, x2 = factor(x2))

fit <- gkwqreg(y ~ x1 + x2 | x1, data = dat, tau = 0.5, family = "kw")

## -- the default prediction is a conditional quantile --------------------
nd <- data.frame(x1 = c(-1, 0, 1), x2 = factor(0, levels = c("0", "1")))
predict(fit, nd)
#> [1] 0.3650941 0.5755396 0.7617504
## 0.3651 0.5755 0.7618.  Read this as: among units with x1 = 0 and x2 = 0,
## half are predicted to fall below 0.5755 -- NOT "they average 0.5755".

## -- the arithmetic check ------------------------------------------------
predict(fit, nd, type = "probability", at = predict(fit, nd),
        elementwise = TRUE)
#> [1] 0.5 0.5 0.5
## 0.5 0.5 0.5 exactly: the fitted distribution puts mass tau below its own
## fitted tau-quantile, by construction of the anchor.

## -- other levels of the SAME fit ---------------------------------------
predict(fit, nd, tau = c(0.1, 0.5, 0.9))
#>        tau=0.1   tau=0.5   tau=0.9
#> [1,] 0.1584320 0.3650941 0.5945101
#> [2,] 0.2616780 0.5755396 0.8438056
#> [3,] 0.3746906 0.7617504 0.9668158
## The "tau=0.5" column reproduces the default call above; the other two are
## read off the same conditional law and were never targeted by estimation.

## How much does that extrapolation cost?  Compare with a genuine refit.
fit90 <- gkwqreg(y ~ x1 + x2 | x1, data = dat, tau = 0.9, family = "kw")
cbind(read_off = predict(fit, nd, tau = 0.9)[, 1],
      refit    = predict(fit90, nd))
#>       read_off     refit
#> [1,] 0.5945101 0.6251732
#> [2,] 0.8438056 0.8532645
#> [3,] 0.9668158 0.9529930
##      read_off  refit
## [1,]   0.5945 0.6252
## [2,]   0.8438 0.8533
## [3,]   0.9668 0.9530
## Close, as it should be here: the family is correctly specified. A wide
## gap would indict the shape assumption, not the covariate effects.

## -- the fitted conditional distribution, row by row ---------------------
predict(fit, nd, type = "parameter")   # alpha, beta, gamma, delta, lambda
#>      alpha      beta gamma delta lambda
#> 1 2.198537 5.9985203     1     0      1
#> 2 2.208171 1.9808351     1     0      1
#> 3 2.217847 0.8756712     1     0      1
predict(fit, nd, type = "link")        # eta for every part
#>           mu     alpha
#> 1 -0.5533215 0.7877921
#> 2  0.3044892 0.7921644
#> 3  1.1622998 0.7965368
predict(fit, nd, type = "mean")        # 0.3720 0.5630 0.7123
#> [1] 0.3720496 0.5629987 0.7123465
predict(fit, nd, type = "variance")
#> [1] 0.02691687 0.04633717 0.05042043
## The conditional mean (0.5630) is NOT the conditional median (0.5755):
## the fitted Kumaraswamy law is skewed. This is why fitted() must not be
## read as a mean.

## -- densities and probabilities ----------------------------------------
predict(fit, nd, type = "density", at = c(0.25, 0.50, 0.75))
#> [1] 1.963491 1.490407 1.502078
## A VECTOR of length 3, not a matrix: `at` happens to have exactly one
## element per row, so `elementwise` defaults to TRUE and pairs the ith
## point with the ith row. Say so explicitly when that is not the intent:
predict(fit, nd, type = "density", at = c(0.25, 0.50, 0.75),
        elementwise = FALSE)                                     # 3 x 3
#>           0.25      0.50      0.75
#> [1,] 1.9634907 1.6825462 0.2116031
#> [2,] 0.7817324 1.4904067 1.4739907
#> [3,] 0.3610858 0.8604633 1.5020780
predict(fit, nd, type = "probability", at = 0.5)                 # 3 x 1
#>            0.5
#> [1,] 0.7709818
#> [2,] 0.3831080
#> [3,] 0.1909808

## -- additive term contributions ----------------------------------------
predict(fit, nd, type = "terms")$mu
#>           x1 x21
#> 1 -0.8578107   0
#> 2  0.0000000   0
#> 3  0.8578107   0
```
