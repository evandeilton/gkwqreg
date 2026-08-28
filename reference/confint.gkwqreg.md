# Confidence intervals for a quantile regression fit

Computes confidence intervals for the coefficients of a `"gkwqreg"` fit
by the Wald approximation, by inverting the profile likelihood, or by
the bootstrap.

## Usage

``` r
# S3 method for class 'gkwqreg'
confint(
  object,
  parm,
  level = 0.95,
  method = c("wald", "profile", "boot"),
  R = 200L,
  ...
)
```

## Arguments

- object:

  A `"gkwqreg"` fit, as returned by
  [`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md).

- parm:

  Coefficients to report, given either as a character vector of names
  from `names(coef(object))` – which have the form `"part:term"`, as in
  `"mu:x1"` or `"alpha:(Intercept)"` – or as integer positions into that
  vector. Defaults to every coefficient of every part.

- level:

  Confidence level, a single number in `(0,1)`. The default, `0.95`,
  gives equal-tailed limits at the 2.5th and 97.5th percentiles.

- method:

  How to construct the interval; see the table below. Partial matching
  applies.

- R:

  Number of bootstrap replicates, used only when `method = "boot"`. The
  default of 200 is adequate for a percentile interval at the usual
  levels; raise it for tail levels or for publication.

- ...:

  Passed to
  [`gkwq_boot()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_boot.md)
  when `method = "boot"`, which is where `type` and `seed` are set.
  Ignored by the other methods.

## Value

A numeric matrix with one row per entry of `parm` and two columns, the
lower and upper limits. Row names are the coefficient names; column
names are the percentages of the limits, for example `"2.5 %"` and
`"97.5 %"`. Entries that could not be computed are `NA`.

## Details

The intervals are for coefficients on the **link scale**: under the
default logit link an interval for a `mu` coefficient is an interval for
an effect on the log quantile odds. Exponentiating the limits gives an
interval for the multiplicative effect on the quantile odds
\\\mu\_\tau/(1-\mu\_\tau)\\, which is legitimate because the
transformation is monotone. An interval for the effect on the quantile
itself is a different object, since it involves the link derivative as
well; use
[`marginal_effects()`](https://evandeilton.github.io/gkwqreg/reference/marginal_effects.md)
for that.

## The three methods

|  |  |  |  |
|----|----|----|----|
| `method` | Construction | Cost | Use when |
| `"wald"` (default) | estimate plus or minus `z * se`, from the stored standard errors | free | the sample is large and the log-likelihood is close to quadratic |
| `"profile"` | inversion of the `TMB` profile likelihood, one coefficient at a time | one profile per coefficient | the Wald approximation is suspect and the sample is small |
| `"boot"` | percentile interval from [`gkwq_boot()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_boot.md) | `R` refits | the family may be misspecified, or the profile does not resolve |

Wald intervals use `object$se`, the square roots of the diagonal of the
stored observed-information inverse, and are therefore always symmetric
about the estimate. This matters more here than in an ordinary
generalized linear model: the conditional quantile is not orthogonal to
the remaining parameters in the Cox-Reid sense, so the profile
log-likelihood in a `mu` coefficient inherits curvature from the
nuisance parameters and can be visibly asymmetric in small samples.
Neither `"profile"` nor `"boot"` imposes symmetry.

`"boot"` passes `...` through to
[`gkwq_boot()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_boot.md),
whose default resampling scheme is `type = "pairs"`. That default is the
right one for interval construction: the appeal of parametric quantile
regression is that it buys efficiency under a distributional assumption,
and a parametric bootstrap that re-imposes the same assumption cannot
say anything about whether it holds. Set `type = "parametric"` only when
the family is not in question.

`"profile"` is the most accurate of the three when it succeeds. It is
not always able to interpolate the likelihood cut-off from the profile
it obtains – the profile may simply not extend far enough on one side.
Limits that cannot be obtained come back as `NA` and the function warns,
naming how many coefficients were affected, rather than reporting a
silently wrong number. Fall back on `"wald"` or `"boot"` when that
happens.

## See also

[`vcov.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/vcov.gkwqreg.md)
for the covariance matrices behind the Wald interval,
[`gkwq_boot()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_boot.md)
for the bootstrap,
[`marginal_effects()`](https://evandeilton.github.io/gkwqreg/reference/marginal_effects.md)
for intervals on the scale of the response,
[`summary.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg-extractors.md)
for the accompanying coefficient table.

## Examples

``` r
set.seed(2024)
n   <- 300
x1  <- runif(n, -2, 2)
x2  <- rbinom(n, 1, 0.5)
mu  <- plogis(0.3 + 0.9 * x1 - 0.5 * x2)
bt  <- log1p(-0.5) / log1p(-mu^2)
y   <- gkwq_quantile(runif(n), alpha = 2, beta = bt)
dat <- data.frame(y = y, x1 = x1, x2 = factor(x2))

fit <- gkwqreg(y ~ x1 + x2 | x1, data = dat, tau = 0.5, family = "kw")

## -- Wald intervals for every coefficient of every part ------------------
round(confint(fit), 4)
#>                     2.5 %  97.5 %
#> mu:(Intercept)     0.1390  0.4700
#> mu:x1              0.7660  0.9496
#> mu:x21            -0.7006 -0.3423
#> alpha:(Intercept)  0.6664  0.9179
#> alpha:x1          -0.0945  0.1032
## mu:x1 lies in (0.7660, 0.9496); the true value used to simulate was 0.9.

## -- on the quantile-odds scale ------------------------------------------
round(exp(confint(fit, parm = "mu:x1")), 4)
#>        2.5 % 97.5 %
#> mu:x1 2.1512 2.5847
## (2.1512, 2.5847): a one-unit rise in x1 multiplies the odds of the median
## by between 2.15 and 2.58. The transformation is monotone, so the limits
## carry over; the effect on the median itself is marginal_effects().

## -- a narrower selection and a different level --------------------------
round(confint(fit, parm = c("mu:x1", "mu:x21"), level = 0.9), 4)
#>            5 %    95 %
#> mu:x1   0.7808  0.9349
#> mu:x21 -0.6718 -0.3711

# \donttest{
## -- percentile bootstrap, resampling observations -----------------------
round(confint(fit, parm = c("mu:x1", "mu:x21"), method = "boot",
              R = 100, seed = 7), 4)
#>          2.5 %  97.5 %
#> mu:x1   0.7829  0.9529
#> mu:x21 -0.7060 -0.2651
##          2.5 %  97.5 %
## mu:x1   0.7829  0.9529
## mu:x21 -0.7060 -0.2651
## For mu:x1 these sit almost on top of the Wald limits, as they should when
## the family is correctly specified and n = 300. For the factor mu:x21 the
## upper limit is visibly further out; with R = 100 a percentile limit is
## itself an estimate, so raise R before reading much into a single
## discrepancy of this size.
# }
```
