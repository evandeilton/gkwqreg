# Pinball (check) loss of a fitted quantile regression

Evaluates, on the fitting data or on a held-out sample, the loss
function that a quantile estimate actually targets. This is the
criterion to use when the question is which family to fit and the answer
will be acted on as a statement about a quantile.

## Usage

``` r
pinball(object, newdata = NULL, y = NULL)
```

## Arguments

- object:

  A `"gkwqreg"` fit, or a `"gkwqregs"` container of fits at several
  quantile levels.

- newdata:

  Optional data frame on which to evaluate the loss. When supplied, the
  fitted quantiles are obtained by
  `predict(object, newdata, type = "quantile")`, so `newdata` must
  contain every covariate used by any part of the model. When `NULL`
  (the default) the loss is computed in sample, from
  `object$fitted.values` and the stored response, and is then identical
  to the value cached in `object$pinball`.

- y:

  Optional numeric response to pair with `newdata`. If omitted, the
  response is taken from the column of `newdata` named on the left-hand
  side of the model formula; if that column is absent, an error asks for
  `y` explicitly. Ignored when `newdata` is `NULL`.

## Value

A single number for a `"gkwqreg"` fit: the weighted mean check loss, on
the scale of the response. For a `"gkwqregs"` container, a numeric
vector with one entry per quantile level, named as the container's
elements are (`"tau=0.1"` and so on).

## Details

Write \\\hat{Q}\_\tau(x_i)\\ for the fitted conditional
\\\tau\\-quantile of observation \\i\\ and \\u_i = y_i -
\hat{Q}\_\tau(x_i)\\ for its residual. The pinball, or check, loss is
the asymmetrically weighted absolute deviation

\$\$\rho\_\tau(u) \\=\\ u \\ \left( \tau - \mathbf{1}\\u \< 0\\
\right),\$\$

equal to \\\tau u\\ when the observation lies above the fitted quantile
and to \\-(1 - \tau) u\\ when it lies below. `pinball()` returns the
weighted sample mean of that loss,

\$\$L\_\tau \\=\\ \frac{\sum_i w_i \\ \rho\_\tau(y_i -
\hat{Q}\_\tau(x_i))}{\sum_i w_i},\$\$

with \\w_i\\ the prior weights supplied to
[`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md).
Weights apply only to the in-sample calculation; when `newdata` is given
every row counts once.

Underprediction is charged \\\tau\\ per unit and overprediction \\1 -
\tau\\ per unit, so at \\\tau = 0.9\\ falling short of an observation
costs nine times as much as overshooting it by the same amount. The loss
carries the units of the response and is therefore numerically small for
a response confined to \\(0,1)\\; only *differences* between competing
fits evaluated on the *same* rows are interpretable.

## Why check loss and not AIC

The population minimiser of \\E\[\rho\_\tau(Y - q)\]\\ over \\q\\ is the
\\\tau\\-quantile of \\Y\\. This is the defining property of the check
function (Koenker and Bassett, 1978) and the reason it is the natural
scoring rule for a quantile forecast (Gneiting, 2011). Check loss
therefore scores exactly the quantity this package estimates, and
nothing else.

AIC scores something materially different: the whole conditional
density. A family can describe the bulk of the distribution well, be
rewarded for it by the likelihood, and still place the quantile of
interest worse than a rival that fits the bulk less well. **The two
criteria can rank families differently, and when they do it is not a
symptom of anything having gone wrong.** In the case study shipped with
this package
([`vignette("gkwqreg-case-study")`](https://evandeilton.github.io/gkwqreg/articles/gkwqreg-case-study.md))
AIC ranks family `"beta"` first while out-of-sample check loss ranks it
last, behind every Kumaraswamy member. When the inferential target is a
quantile, act on the check loss; when it is the conditional distribution
as a whole, AIC is the relevant criterion. Deciding which of the two
applies is a modelling decision, not a computational one.

## Why out of sample

In-sample check loss is not a model-selection criterion. It carries no
penalty for dimension, so a family with more free parameters can track
the observed rows more closely whether or not the extra flexibility
corresponds to anything real, and comparing it with a sub-family in
sample mostly rewards that flexibility. Evaluate on rows the model has
not seen: pass a held-out `newdata`, or average `pinball()` over
cross-validation folds. The in-sample value remains useful as a
description of fit – and it is what
[`compare_families()`](https://evandeilton.github.io/gkwqreg/reference/compare_families.md)
reports, alongside AIC and BIC, which do carry penalties – but it should
not be the basis for choosing between families.

## References

Koenker, R. and Bassett, G. (1978). Regression quantiles.
*Econometrica*.

Gneiting, T. (2011). Making and evaluating point forecasts. *Journal of
the American Statistical Association*.

## See also

[`compare_families()`](https://evandeilton.github.io/gkwqreg/reference/compare_families.md),
which tabulates AIC, BIC and the in-sample check loss for every family
at once;
[`residuals.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/residuals.gkwqreg.md)
with `type = "check"` for the per-observation contributions to
\\L\_\tau\\;
[`vuong_test()`](https://evandeilton.github.io/gkwqreg/reference/vuong_test.md)
for a formal test between two non-nested families.

## Examples

``` r
## A response in (0,1) whose dispersion, not only its location, depends on x.
## Fitting `y ~ x` regresses the quantile alone and leaves the shape constant,
## so every family below is misspecified, to a different degree and in a
## different way -- which is exactly when the choice of criterion matters.
set.seed(1)
n <- 600
x <- runif(n, -2, 2)
mu <- plogis(0.3 + 1.0 * x)
phi <- exp(1 + 0.8 * x)
y <- rbeta(n, phi * mu, phi * (1 - mu))
d <- data.frame(y = y, x = x)
train <- d[1:300, ]
test <- d[301:600, ]

fit <- gkwqreg(y ~ x, data = train, tau = 0.5, family = "kw")

## In sample this is the value the likelihood machinery already cached.
c(pinball(fit), fit$pinball)
#> [1] 0.1116671 0.1116671
#> 0.1116671 0.1116671

## Out of sample, and identical to the definition applied by hand.
e <- test$y - predict(fit, newdata = test, type = "quantile")
c(pinball(fit, newdata = test), mean(e * (0.5 - (e < 0))))
#> [1] 0.1034218 0.1034218
#> 0.1034218 0.1034218

## AIC and out-of-sample check loss need not agree.
tab <- t(sapply(c("kw", "ekw", "bkw", "beta", "mc"), function(f) {
  m <- suppressWarnings(gkwqreg(y ~ x, data = train, tau = 0.5, family = f))
  c(AIC = AIC(m), in_sample = pinball(m), out_of_sample = pinball(m, newdata = test))
}))
round(tab, 4)
#>            AIC in_sample out_of_sample
#> kw   -291.5234    0.1117        0.1034
#> ekw  -480.8364    0.1003        0.0934
#> bkw  -417.9420    0.1010        0.0938
#> beta -484.7745    0.1012        0.0940
#> mc   -482.7745    0.1012        0.0940
#>            AIC in_sample out_of_sample
#> kw   -291.5234    0.1117        0.1034
#> ekw  -480.8364    0.1003        0.0934
#> bkw  -417.9420    0.1010        0.0938
#> beta -484.7745    0.1012        0.0940
#> mc   -482.7745    0.1012        0.0940

rownames(tab)[which.min(tab[, "AIC"])]            # "beta"
#> [1] "beta"
rownames(tab)[which.min(tab[, "out_of_sample"])]  # "ekw"
#> [1] "ekw"

## AIC prefers "beta"; the median it delivers is the fourth best of five.
## For a question about the median, "ekw" is the fit to report.

## One value per level from a container of independent fits.
fits <- gkwqreg(y ~ x, data = train, tau = c(0.1, 0.5, 0.9), family = "kw")
round(pinball(fits, newdata = test), 4)
#> tau=0.1 tau=0.5 tau=0.9 
#>  0.0550  0.1034  0.0389 
#> tau=0.1 tau=0.5 tau=0.9
#>  0.0550  0.1034  0.0389
## The loss shrinks towards the tails because rho_tau is bounded by the
## smaller of tau and 1 - tau; these three numbers are NOT comparable with
## one another, only with a rival fit at the same level.
```
