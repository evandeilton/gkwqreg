# Fitted conditional quantiles

**Returns the fitted conditional \\\tau\\-quantile, not the conditional
mean.** This is the single most consequential difference between a
`"gkwqreg"` object and the mean-parametrized `gkwreg` fits it otherwise
resembles, and the class deliberately does not inherit from that one so
that the confusion cannot happen by silent dispatch.

## Usage

``` r
# S3 method for class 'gkwqreg'
fitted(object, type = c("quantile", "mean", "parameter"), ...)
```

## Arguments

- object:

  A `"gkwqreg"` fit, as returned by
  [`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md).

- type:

  Which fitted quantity to return; see the table below. Partial matching
  applies.

- ...:

  Currently unused; present for consistency with the generic.

## Value

For `type = "quantile"` and `type = "mean"`, a numeric vector of length
`nobs(object)`, aligned with the rows retained after `na.action`. For
`type = "parameter"`, a data frame with `nobs(object)` rows and the five
columns `alpha`, `beta`, `gamma`, `delta`, `lambda` in `gkwdist` order.

## Details

The default value is \$\$\hat{y}\_i \\=\\ \widehat{Q}(\tau \mid x_i)
\\=\\ \hat\mu\_\tau(x_i),\$\$ the fitted \\\tau\\-quantile for
observation \\i\\, taken straight from the quantile part of the linear
predictor through its inverse link. It is not a prediction of \\y_i\\ in
the least-squares sense and should not be averaged, differenced against
\\y_i\\, or fed to a residual routine that assumes a mean. Use
[`residuals.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/residuals.gkwqreg.md)
for residuals defined on the right scale.

The operational meaning is a coverage statement: a fraction \\\tau\\ of
the observations should fall at or below their own fitted value. That
check is one line and worth running, and
[`summary()`](https://rdrr.io/r/base/summary.html) reports it as the
empirical coverage. A conditional mean satisfies no such property, which
is what the examples make concrete.

## Values of `type`

|  |  |  |
|----|----|----|
| `type` | What is returned | Shape |
| `"quantile"` (default) | the fitted conditional tau-quantile | numeric vector of length `nobs(object)` |
| `"mean"` | the conditional mean by 64-node Gauss-Legendre quadrature of the fitted quantile function | numeric vector of length `nobs(object)` |
| `"parameter"` | the reconstructed distribution parameters, one row per observation | data frame with columns `alpha`, `beta`, `gamma`, `delta`, `lambda` |

`type = "mean"` exists so that a quantile fit can be placed beside a
mean-parametrized one, and so that a mean-based diagnostic borrowed from
elsewhere can be reproduced deliberately rather than by accident. It is
never the default, and it is never what the model was estimated to get
right.

`type = "parameter"` returns the fitted conditional law itself. For each
row the anchored parameter has been solved so that the distribution has
\\\hat\mu\_\tau(x_i)\\ as its exact \\\tau\\-quantile, so these values
can be handed directly to the `gkwdist` density, distribution and
random-number functions.

## See also

[`predict.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/predict.gkwqreg.md)
for the same quantities at new covariate values and for further quantile
levels,
[`residuals.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/residuals.gkwqreg.md)
for quantile residuals,
[`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md)
for the model.

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

## Fit the ninth decile, where quantile and mean are far apart.
fit90 <- gkwqreg(y ~ x1 + x2 | x1, data = dat, tau = 0.9, family = "kw")

## -- what fitted() actually promises -------------------------------------
mean(dat$y <= fitted(fit90))                  # 0.8900, target tau = 0.9
#> [1] 0.89
mean(dat$y <= fitted(fit90, type = "mean"))   # 0.5067
#> [1] 0.5066667
## Ninety per cent of the sample lies below its own fitted value, as a 0.9
## quantile must. The conditional mean covers only about half. Reading
## fitted() as a mean would misstate the estimand by four deciles here.

mean(fitted(fit90))                # 0.6991
#> [1] 0.6991453
mean(fitted(fit90, "mean"))        # 0.4736
#> [1] 0.4736314

## -- the fitted conditional distribution ---------------------------------
pv <- fitted(fit90, type = "parameter")
head(pv, 3)
#>      alpha      beta gamma delta lambda
#> 1 2.186951 0.8493235     1     0      1
#> 2 2.176294 7.5514114     1     0      1
#> 3 2.183712 1.1587222     1     0      1
## Each row is a Kumaraswamy law (gamma = 1, delta = 0, lambda = 1) whose
## 0.9-quantile is exactly the corresponding fitted value:
all.equal(gkwq_quantile(0.9, pv$alpha, pv$beta, pv$gamma, pv$delta,
                        pv$lambda),
          fitted(fit90))
#> [1] TRUE
## TRUE. beta varies across rows because it is the anchor: it is recomputed
## from the fitted quantile rather than estimated.
```
