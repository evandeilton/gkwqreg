# Compare the seven families at one quantile level

Refits the same model under each family in turn and tabulates the
results side by side, reporting AIC, BIC and check loss together. The
three criteria answer different questions and can genuinely disagree;
the sections below explain when that happens and which number to act on.

## Usage

``` r
compare_families(
  object,
  families = c("kw", "ekw", "kkw", "bkw", "gkw", "mc", "beta"),
  ...
)
```

## Arguments

- object:

  A fitted `"gkwqreg"` model whose call is reused. Everything is held
  fixed except the family: the same formula, data, quantile level,
  links, weights, offsets and control settings are used for every refit.

- families:

  Character vector of families to try, any subset of
  `c("kw", "ekw", "kkw", "bkw", "gkw", "mc", "beta")`. Defaults to all
  seven. Restricting it is often sensible: `"gkw"` is weakly identified
  in every parametrization and warns accordingly, and a family that
  cannot represent the data at all merely costs time.

- ...:

  Currently unused; present for future extension. Arguments given here
  do not reach
  [`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md).
  To vary anything other than the family, change it in `object` and call
  this function again.

## Value

A data frame with one row per entry of `families`, ordered by increasing
`AIC` with failures last, and with row names reset. Its columns are

- `family`:

  the family fitted in that row.

- `anchor`:

  the anchor that family defaulted to; `NA` on failure.

- `df`:

  the number of estimated coefficients, \\p\\.

- `logLik`:

  the maximized log-likelihood, \\\ell(\hat\theta)\\.

- `AIC`:

  \\-2\ell(\hat\theta) + 2p\\.

- `BIC`:

  \\-2\ell(\hat\theta) + p\log n\\.

- `pinball`:

  the in-sample mean check loss \\L\_\tau\\.

- `converged`:

  `TRUE` when the optimizer reported convergence.

## Details

The stored call of `object` is re-evaluated once per family, with
`family` replaced and `anchor` deleted, in the environment from which
`compare_families()` was itself called. Two practical consequences
follow. The data and any other objects named in the call must be visible
from that environment, exactly as they were when `object` was fitted.
And because the anchor is dropped rather than carried over, **each
family is fitted under its own default anchor**: `beta` for the five
Kumaraswamy members, `lambda` for `"mc"`, `gamma` for `"beta"`. That is
not a lapse. Not every anchor is admissible in every family – `"mc"` and
`"beta"` each admit exactly one – so insisting on a common anchor would
simply fail. It does mean the rows can differ in more than the family,
which is why `anchor` is a column of the result and worth reading.

A family that fails to fit is not silently dropped. It contributes a row
of `NA` values with `converged = FALSE`, so an absent family is visible
as an absence rather than as nothing at all. The `converged` column
carries the optimizer's own verdict, and a row reporting `FALSE` should
be excluded from the ranking rather than interpreted: an unconverged
fit's log-likelihood is a lower bound on its maximum, so its AIC is not
comparable with the rest.

## The three criteria

With \\\ell(\hat\theta)\\ the maximized log-likelihood, \\p\\ the number
of estimated coefficients and \\n\\ the number of observations,

\$\$\mathrm{AIC} = -2\ell(\hat\theta) + 2p, \qquad \mathrm{BIC} =
-2\ell(\hat\theta) + p \log n.\$\$

Both are computed from the whole conditional density and differ only in
how hard they penalize dimension: \\\log n\\ exceeds 2 for any \\n \ge
8\\, so BIC prefers smaller models and is consistent for the true model
when it is among the candidates, whereas AIC targets predictive accuracy
of the density and tends to over-select.

The `pinball` column is the mean check loss at the fitted quantiles
\\\hat Q_i(\tau)\\,

\$\$L\_\tau = \frac{\sum\_{i=1}^n w_i \\(y_i - \hat Q_i(\tau)) \left\\
\tau - \mathbf{1}(y_i \< \hat Q_i(\tau)) \right\\} {\sum\_{i=1}^n
w_i},\$\$

with \\w_i\\ the prior weights. The summand is \\\tau\\(y_i - \hat
Q_i)\\ when the observation lies above the fitted quantile and \\(1 -
\tau)(\hat Q_i - y_i)\\ when it lies below, so the loss is asymmetric in
exactly the proportion \\\tau : 1 - \tau\\ that defines the quantile.
Its population minimizer is the true conditional \\\tau\\-quantile and
nothing else; see
[`pinball()`](https://evandeilton.github.io/gkwqreg/reference/pinball.md).

## Why AIC and check loss can disagree

They measure different things, and the disagreement is not a paradox but
the point. AIC and BIC score the entire conditional density: how well
the model describes the response at every point of its support, the bulk
of the distribution included. Check loss scores one feature of that
density, the conditional \\\tau\\-quantile, and is indifferent to
everything else. A family can capture the bulk of the distribution well,
which is what earns it the likelihood, and still place the
\\\tau\\-quantile worse than a rival that fits the bulk less well –
especially at an extreme level, where the quantile is determined by a
region of the support that contributes little to the likelihood.

This is not hypothetical. In the case study shipped with the package,
[`vignette("gkwqreg-case-study")`](https://evandeilton.github.io/gkwqreg/articles/gkwqreg-case-study.md),
AIC ranks the `"beta"` family first while out-of-sample check loss ranks
it last, behind every Kumaraswamy member. The example below reproduces a
milder version of the same reversal on simulated data.

The resolution is to let the question decide. If the object of inference
is the conditional quantile, the check loss is the criterion that
matches the goal, and AIC is answering a question that was not asked.

## Compute the check loss out of sample

One qualification is essential. The `pinball` column returned here is
computed **in sample**, on the same observations that were used to fit,
and an in-sample loss rewards flexibility for its own sake: a larger
family can always drive it down a little by tracking noise. It is
reported because it is free and because it exposes gross failures, not
because it settles anything.

A check loss used to *choose* a family must be evaluated on data the fit
has not seen – a holdout sample or a cross-validation, both shown in the
example below – via
[`pinball()`](https://evandeilton.github.io/gkwqreg/reference/pinball.md)
with its `newdata` argument. The out-of-sample check loss is the number
to act on for a quantile question.

## See also

[`pinball()`](https://evandeilton.github.io/gkwqreg/reference/pinball.md)
for the check loss, in or out of sample;
[`anova.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/anova.gkwqreg.md)
for a formal test along a nesting chain;
[`vuong_test()`](https://evandeilton.github.io/gkwqreg/reference/vuong_test.md)
for a formal test between two non-nested fits;
[`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md)
for the family and anchor arguments themselves.

## Examples

``` r
## Data from an EKW distribution, split into a training and a holdout part.
set.seed(2026)
n  <- 500
x  <- runif(n, -2, 2)
mu <- plogis(0.3 + 1.2 * x)
b  <- log1p(-0.9^(1 / 4)) / log1p(-mu^2)      # beta anchor solve at tau = 0.9
y  <- gkwdist::rekw(n, alpha = 2, beta = b, lambda = 4)
d  <- data.frame(y = y, x = x)
train <- d[1:350, ]
test  <- d[351:500, ]

fit <- gkwqreg(y ~ x, data = train, tau = 0.9, family = "kw")
cmp <- compare_families(fit, families = c("kw", "ekw", "kkw", "beta"))
cmp
#>   family anchor df   logLik       AIC       BIC    pinball converged
#> 1    ekw   beta  4 391.9588 -775.9176 -760.4858 0.01504029      TRUE
#> 2    kkw   beta  5 392.7234 -775.4469 -756.1572 0.01500626      TRUE
#> 3     kw   beta  3 383.0529 -760.1058 -748.5320 0.01507971      TRUE
#> 4   beta  gamma  3 297.2013 -588.4027 -576.8289 0.01994690      TRUE
##   family anchor df   logLik       AIC       BIC    pinball converged
## 1    ekw   beta  4 391.9588 -775.9176 -760.4858 0.01504029      TRUE
## 2    kkw   beta  5 392.7234 -775.4469 -756.1572 0.01500626      TRUE
## 3     kw   beta  3 383.0529 -760.1058 -748.5320 0.01507971      TRUE
## 4   beta  gamma  3 297.2013 -588.4027 -576.8289 0.01994690      TRUE

## Now recompute the check loss on the holdout, which is the criterion that
## matches a quantile question.
holdout <- vapply(cmp$family, function(fm)
  pinball(gkwqreg(y ~ x, data = train, tau = 0.9, family = fm),
          newdata = test),
  numeric(1))
cbind(cmp[, c("family", "AIC", "pinball")], holdout = round(holdout, 5))
#>      family       AIC    pinball holdout
#> ekw     ekw -775.9176 0.01504029 0.01859
#> kkw     kkw -775.4469 0.01500626 0.01871
#> kw       kw -760.1058 0.01507971 0.01854
#> beta   beta -588.4027 0.01994690 0.02362
##      family       AIC    pinball holdout
## ekw     ekw -775.9176 0.01504029 0.01859
## kkw     kkw -775.4469 0.01500626 0.01871
## kw       kw -760.1058 0.01507971 0.01854
## beta   beta -588.4027 0.01994690 0.02362
##
## The orderings differ. AIC puts "kw" third and the in-sample check loss puts
## it third as well, both rewarding the larger families -- which is what extra
## parameters are for. On data the fits have not seen, the smallest family
## predicts the 0.9-quantile best. Only "beta" is last by every criterion,
## and it is the one family here that is genuinely excluded.
```
