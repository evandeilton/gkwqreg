# Bootstrap a quantile regression fit

Refits the model on `R` resampled data sets and summarizes the resulting
spread of the coefficients: a bootstrap standard error and a percentile
interval for each one, together with the full replicate matrix and its
covariance. It is the tool to reach for when the model-based Wald
standard errors are in doubt – a small sample, a weakly identified
family, an information matrix that is barely positive definite, or
simply a distributional assumption one would rather not lean on.

## Usage

``` r
gkwq_boot(object, R = 200L, type = c("pairs", "parametric"), seed = NULL, ...)
```

## Arguments

- object:

  A `"gkwqreg"` fit. It must have been fitted with `model = TRUE`, which
  is the default, because the retained model frame is what gets
  resampled; otherwise the function stops with an explanatory message.

- R:

  Number of bootstrap replicates. The default of 200 is adequate for a
  standard error, which is a scale estimate and converges quickly.
  Percentile intervals estimate the tails of the replicate distribution
  and need considerably more, of the order of 1000 or beyond.

- type:

  Resampling scheme, one of `"pairs"` or `"parametric"`. `"pairs"`, the
  default, resamples whole observations – rows of the model frame,
  response and covariates together – with replacement. `"parametric"`
  holds the covariates fixed and redraws only the responses, from the
  fitted conditional distributions, via
  [`simulate.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/simulate.gkwqreg.md).

- seed:

  Optional seed, passed to
  [`base::set.seed()`](https://rdrr.io/r/base/Random.html). Unlike
  [`simulate.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/simulate.gkwqreg.md),
  the ambient random number stream is **not** restored afterwards, so a
  seeded call does change `.Random.seed`.

- ...:

  Unused; present for future extension.

## Value

An object of class `"gkwq_boot"`: a list with components

- `replicates`:

  an `R` by \\p\\ matrix of replicate coefficient vectors, with column
  names from the fit and `NA` rows where a replicate failed.

- `vcov`:

  the \\p\\ by \\p\\ bootstrap covariance matrix \\\widehat V\\,
  computed over the surviving replicates only.

- `R`:

  the number of replicates requested.

- `n_ok`:

  the number that converged and entered the summaries.

- `type`:

  the resampling scheme used.

- `estimate`:

  the coefficients of the original fit.

- `se`:

  bootstrap standard errors, the square root of `diag(vcov)`.

- `percentile`:

  a \\p\\ by 2 matrix of percentile interval endpoints at 0.025 and
  0.975, one row per coefficient.

- `tau`,`family`,`anchor`:

  copied from the fit, for the printed header.

The [`print()`](https://rdrr.io/r/base/print.html) method shows the
original estimates beside the bootstrap standard errors and the
percentile interval.

## Details

Each replicate re-evaluates the original call on the resampled data.
Only `subset` is stripped from that call, having already been applied in
building the model frame. Everything else that is indexed by observation
is carried across and reindexed with the resampled rows: prior
`weights`, an `offset` argument, and
[`offset()`](https://rdrr.io/r/stats/offset.html) terms written into any
part of the formula. A weighted or offset fit is therefore bootstrapped
as the estimator it actually is.

## Why pairs is the default

The whole bargain of parametric quantile regression is that it buys
efficiency, and interpretable coefficients, by assuming a conditional
distribution. A parametric bootstrap re-imposes that very assumption at
every replicate: the responses it draws are by construction exactly what
the fitted family says they should be, since
[`simulate.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/simulate.gkwqreg.md)
samples from the fitted distributions themselves. Its intervals
therefore describe the sampling behaviour of the estimator *in a world
where the family is correct*, and can say nothing whatever about whether
that is the world the data came from. If the distributional assumption
is the thing in question, the parametric bootstrap is constitutionally
unable to answer.

The pairs bootstrap makes no such assumption. It resamples observations
and so requires only that they be exchangeable draws from the joint
distribution of response and covariates. Its spread reflects whatever
generated the data, correctly specified or not, which is why it is the
honest default. The price is some efficiency when the family really is
right, and a heavier tail in small samples, where a resample can happen
to contain too few distinct covariate patterns to identify the model – a
factor variable with a rare level is the usual culprit, and shows up as
a low `n_ok`.

Use `type = "parametric"` when the family is *not* what is being
questioned: to check the accuracy of the Wald approximation for a family
one is willing to assume, to study the estimator's behaviour in a
simulation, or to build envelopes for a diagnostic plot.

## Failed replicates

A replicate whose refit throws an error, or which the optimizer does not
certify as converged, contributes a row of `NA` to `replicates` and is
excluded from the covariance, the standard errors and the percentile
intervals. The count that survived is reported as `n_ok` and printed
alongside `R`. That ratio is itself a diagnostic: a fit that loses a
noticeable share of its replicates is telling you that the likelihood
surface is fragile under mild perturbation of the data, and the
surviving replicates are then a selected sample rather than a random
one, so the resulting intervals are optimistic. The function stops
outright if fewer than two replicates survive.

## What the output means

With \\\hat\theta^{(r)}\\ the coefficient vector from replicate \\r\\
and \\B\\ the number that converged, the reported covariance is the
ordinary sample covariance of the surviving replicates,

\$\$\widehat{V} = \frac{1}{B - 1} \sum\_{r} \left(\hat\theta^{(r)} -
\bar{\theta}^{\*}\right) \left(\hat\theta^{(r)} -
\bar{\theta}^{\*}\right)^{\top},\$\$

with \\\bar\theta^{\*}\\ the mean over surviving replicates; `se` is the
square root of its diagonal. The percentile interval for coefficient
\\j\\ is the pair of empirical quantiles of \\\hat\theta_j^{(r)}\\ at
0.025 and 0.975, which is the basic percentile method and is not
corrected for bias or acceleration.

Note that `estimate` is the coefficient vector of the *original* fit,
not the mean of the replicates. The difference between the two estimates
the bootstrap bias, and a percentile interval that sits noticeably
off-centre from `estimate` is a sign that the sampling distribution is
skewed and that the symmetric Wald interval is the wrong shape.

## See also

[`simulate.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/simulate.gkwqreg.md)
for the draws behind `type = "parametric"`;
[`vcov.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/vcov.gkwqreg.md)
and
[`confint.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/confint.gkwqreg.md)
for the model-based alternatives;
[`gkwq_control()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_control.md)
for the optimizer settings that govern each refit.

## Examples

``` r
set.seed(11)
n  <- 400
x  <- runif(n, -2, 2)
mu <- plogis(0.4 + 1.1 * x)
y  <- gkwdist::rkw(n, alpha = 2, beta = log1p(-0.9) / log1p(-mu^2))
d  <- data.frame(y = y, x = x)

fit <- gkwqreg(y ~ x, data = d, tau = 0.9, family = "kw")

# \donttest{
## The default: resample observations. No distributional assumption is
## re-imposed, so this is the version that can disagree with the model-based
## standard errors when the family is wrong.
bp <- gkwq_boot(fit, R = 200, seed = 1)
bp
#> 
#> Bootstrap for a kw quantile regression (tau = 0.9, anchor beta)
#> pairs bootstrap, 200 of 200 replicates converged
#> 
#>                   Estimate Boot SE   2.5%  97.5%
#> mu:(Intercept)      0.3972  0.0586 0.2763 0.5127
#> mu:x                1.0582  0.0496 0.9694 1.1681
#> alpha:(Intercept)   0.6960  0.0416 0.6192 0.7897
##                   Estimate Boot SE   2.5%  97.5%
## mu:(Intercept)      0.3972  0.0586 0.2763 0.5127
## mu:x                1.0582  0.0496 0.9694 1.1681
## alpha:(Intercept)   0.6960  0.0416 0.6192 0.7897

## The parametric alternative, for comparison only. Its replicates are drawn
## from the fitted Kumaraswamy laws themselves, so it assumes precisely what
## a specification check would want to test.
bq <- gkwq_boot(fit, R = 200, type = "parametric", seed = 1)

## Here the family is in fact correct -- the data were generated from it --
## so all three sets of standard errors agree closely. That agreement is
## itself the evidence: under misspecification it is the pairs column that
## remains trustworthy.
round(cbind(wald       = sqrt(diag(vcov(fit))),
            pairs      = bp$se,
            parametric = bq$se), 4)
#>                     wald  pairs parametric
#> mu:(Intercept)    0.0603 0.0586     0.0567
#> mu:x              0.0505 0.0496     0.0466
#> alpha:(Intercept) 0.0470 0.0416     0.0479
##                     wald  pairs parametric
## mu:(Intercept)    0.0603 0.0586     0.0567
## mu:x              0.0505 0.0496     0.0466
## alpha:(Intercept) 0.0470 0.0416     0.0479

## Every replicate converged, so nothing was silently discarded.
c(requested = bp$R, converged = bp$n_ok)
#> requested converged 
#>       200       200 
## requested converged
##       200       200
# }
```
