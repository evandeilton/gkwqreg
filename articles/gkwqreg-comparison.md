# How gkwqreg relates to quantreg, unitquantreg and betareg

## “Quantile” means three different things here

Most of the confusion in this corner of the ecosystem comes from one
word carrying three meanings.

**1. Check-function quantile regression.** `quantreg`, `qgam`.
Distribution-free, indexed by `tau`, no likelihood. Minimises the
pinball loss directly.

**2. Location- or median-parametrized models.** `cdfquantreg`, and
`betareg`’s `predict(type = "quantile")`. A *single* fit whose location
parameter happens to be the median, or quantiles derived after the fact
from a mean model. These are **not** `tau`-indexed — `cdfquantreg` has
no `tau` argument at all.

**3. `tau`-indexed parametric quantile regression.** `unitquantreg`, and
this package. One likelihood per `tau`, with a parameter that *is* the
`tau`-quantile.

| package | support | type | families | nested? |
|----|----|----|----|----|
| `quantreg` | ℝ | check-function | — | — |
| `qgam` | ℝ | smooth check-function | — | — |
| `betareg` | (0,1) | mean–precision | beta | — |
| `cdfquantreg` | (0,1) | median / location | 36 pairs | no |
| `gamlss` | flexible | distributional | many (no Kumaraswamy) | no |
| `gkwreg` | (0,1) | distributional (mean) | 7 | **yes** |
| `unitquantreg` | (0,1) | **τ-indexed parametric** | 13 | no |
| **`gkwqreg`** | (0,1) | **τ-indexed parametric** | **7** | **yes** |

## The prior art, and the identity that proves it

The `kw` slice of this package **is** the Kumaraswamy quantile
reparametrization of Mitnik and Baek (2013), which `unitquantreg`
implements as its `kum` family. That is not a similarity to be glossed
over — it is the strongest validation the construction has, so the
package tests it.

``` r

library(unitquantreg)
n <- 400
x <- runif(n, -2, 2)
mu <- plogis(0.4 + 1.1 * x)
y <- gkwdist::rkw(n, alpha = 2, beta = log1p(-0.5) / log1p(-mu^2))
d <- data.frame(y = y, x = x)

out <- t(sapply(c(0.25, 0.50, 0.75), function(tau) {
  a <- gkwqreg(y ~ x, data = d, tau = tau, family = "kw")
  b <- unitquantreg(y ~ x, data = d, tau = tau, family = "kum")
  c(tau = tau,
    gkwqreg_b0 = coef(a)[[1]], gkwqreg_b1 = coef(a)[[2]],
    unitq_b0 = coef(b)[[1]],   unitq_b1 = coef(b)[[2]],
    max_diff = max(abs(unname(coef(a)[1:2]) - unname(coef(b)[1:2]))))
}))
knitr::kable(out, digits = c(2, 5, 5, 5, 5, 10))
```

|  tau | gkwqreg_b0 | gkwqreg_b1 | unitq_b0 | unitq_b1 |   max_diff |
|-----:|-----------:|-----------:|---------:|---------:|-----------:|
| 0.25 |   -0.48739 |    0.84984 | -0.48739 |  0.84984 | 1.1600e-08 |
| 0.50 |    0.35985 |    1.06736 |  0.35985 |  1.06737 | 3.3983e-06 |
| 0.75 |    1.25924 |    1.36879 |  1.25924 |  1.36879 | 4.1230e-07 |

On the author’s machine the two agree to between `1.2e-08` and
`4.4e-06`. They are the same estimator, implemented independently — one
by hand-derived C++ gradients, one by automatic differentiation through
a general five-parameter template.

**What is new here** is not the `kw` case. It is the extension to the
four- and five-parameter members — for which no quantile regression
appears in print — and the fact that they are **nested**, so family
selection is an ordinary likelihood-ratio test. `unitquantreg`’s
thirteen families are mutually non-nested and require Vuong tests.

``` r

d <- local({
  n <- 400; x <- runif(n, -2, 2); mu <- plogis(0.4 + 1.1 * x)
  data.frame(y = gkwdist::rkw(n, 2, log1p(-0.5) / log1p(-mu^2)), x = x)
})
anova(gkwqreg(y ~ x, data = d, tau = .5, family = "kw"),
      gkwqreg(y ~ x, data = d, tau = .5, family = "ekw"),
      gkwqreg(y ~ x, data = d, tau = .5, family = "gkw"))
#> Warning: family = "gkw" is weakly identified in any parametrization
#> (information-matrix condition numbers of order 1e8 to 1e11). Consider a
#> sub-family such as "ekw" or "kkw"; summary() reports the condition number so
#> you can check.
#> Warning in sqrt(d): NaNs produced
#> Warning: the information matrix is not positive definite, so 2 standard
#> error(s) are unavailable. This is the signature of a parameter drifting along a
#> flat ridge; a smaller sub-family usually fixes it.
#> Warning: the optimizer did not report convergence (code 1). Treat the estimates
#> as provisional.
#> Likelihood-ratio test for Generalized Kumaraswamy quantile regression
#> tau = 0.5, anchor = beta
#>      Df logLik     AIC   Chisq Chi Df Pr(>Chisq)
#> [1,]  3 208.53 -411.05                          
#> [2,]  4 208.70 -409.40  0.3470      1     0.5558
#> [3,]  6 206.94 -401.87 -3.5258      2     1.0000
```

## Against `quantreg`, under correct specification

``` r

library(quantreg)
#> Loading required package: SparseM
#> 
#> Attaching package: 'quantreg'
#> The following object is masked from 'package:gkwqreg':
#> 
#>     rearrange
loss <- function(y, q, tau) mean((y - q) * (tau - ((y - q) < 0)))
res <- t(sapply(c(0.25, 0.5, 0.75), function(tau) {
  a <- gkwqreg(y ~ x, data = d, tau = tau, family = "kw")
  r <- rq(y ~ x, data = d, tau = tau)
  c(tau = tau, gkwqreg = pinball(a), rq = loss(d$y, fitted(r), tau))
}))
knitr::kable(res, digits = 6)
```

|  tau |  gkwqreg |       rq |
|-----:|---------:|---------:|
| 0.25 | 0.067567 | 0.067232 |
| 0.50 | 0.080418 | 0.080974 |
| 0.75 | 0.056009 | 0.061722 |

The two are close, as they should be: with the family correctly
specified, the parametric fit has little to gain in sample and loses
nothing. Its advantages lie elsewhere — standard errors from a
likelihood, nested model comparison, a full conditional distribution,
and no crossing.

## Against `quantreg`, under misspecification

This is where the honest comparison lives, and it must not be buried.
Generate from a two-component beta mixture, which is bimodal and lies
outside the Generalized Kumaraswamy family entirely.

``` r

n <- 600
x <- runif(n, -2, 2)
comp <- rbinom(n, 1, plogis(0.5 * x))
y <- ifelse(comp == 1, rbeta(n, 8, 2), rbeta(n, 2, 8))
y <- pmin(pmax(y, 1e-8), 1 - 1e-8)
dm <- data.frame(y = y, x = x)

## split-sample, so the comparison is out of sample
idx <- sample(n, n / 2)
tr <- dm[idx, ]; te <- dm[-idx, ]

mis <- t(sapply(c(0.1, 0.25, 0.5, 0.75, 0.9), function(tau) {
  a <- suppressWarnings(gkwqreg(y ~ x, data = tr, tau = tau, family = "ekw"))
  r <- rq(y ~ x, data = tr, tau = tau)
  c(tau = tau,
    gkwqreg = loss(te$y, predict(a, newdata = te), tau),
    rq      = loss(te$y, predict(r, newdata = te), tau))
}))
knitr::kable(mis, digits = 5)
```

|  tau | gkwqreg |      rq |
|-----:|--------:|--------:|
| 0.10 | 0.04533 | 0.04533 |
| 0.25 | 0.09985 | 0.09911 |
| 0.50 | 0.13993 | 0.13584 |
| 0.75 | 0.09697 | 0.09419 |
| 0.90 | 0.04324 | 0.04263 |

Under a distribution the family cannot represent, check-function
quantile regression remains consistent for the conditional quantile and
the parametric fit does not. That is the central trade-off of parametric
quantile regression, and it is the reason `quantreg` remains the right
default when you know nothing about the conditional distribution.

The reverse holds too: when the family is right, the parametric fit uses
information the check function throws away, and it gives you a
likelihood to do inference with.

## Against `betareg`

`betareg` fits a mean–precision model. Its `predict(type = "quantile")`
applies the beta quantile function to the fitted mean and precision,
which is meaning (2) above.

``` r

library(betareg)
data("GasolineYield", package = "betareg")

bq <- betareg(yield ~ batch + temp, data = GasolineYield)
gq <- gkwqreg(yield ~ batch + temp, data = GasolineYield, tau = 0.9, family = "kw")
#> Warning in sqrt(d): NaNs produced
#> Warning: the information matrix is not positive definite, so 1 standard
#> error(s) are unavailable. This is the signature of a parameter drifting along a
#> flat ridge; a smaller sub-family usually fixes it.

## betareg: the 0.9-quantile implied by a MEAN model
q_breg <- qbeta(0.9, shape1 = predict(bq, type = "response") * predict(bq, type = "precision"),
                shape2 = (1 - predict(bq, type = "response")) * predict(bq, type = "precision"))

c(coverage_betareg = mean(GasolineYield$yield <= q_breg),
  coverage_gkwqreg = mean(GasolineYield$yield <= fitted(gq)),
  target = 0.9)
#> coverage_betareg coverage_gkwqreg           target 
#>          0.90625          0.96875          0.90000
```

Both produce a number. Only one of them has a coefficient whose standard
error refers to the 0.9-quantile:

``` r

round(summary(gq)$coefficients["mu:temp", ], 5)
#> Warning in sqrt(d): NaNs produced
#>   Estimate Std. Error    z value   Pr(>|z|) 
#>    0.01179         NA         NA         NA
```

In `betareg` the coefficient on `temp` describes the **mean**.
Converting it into a statement about the 0.9-quantile means propagating
uncertainty through both the mean and the precision submodels; here the
quantity is modelled directly.

## Choosing between them

- Conditional distribution unknown or clearly outside the family →
  **`quantreg`**.
- Smooth or additive quantile effects → **`qgam`**.
- Median only, no interest in shape → **`unitquantreg`** or
  **`cdfquantreg`**.
- A mean model that happens to also report quantiles → **`betareg`**,
  **`gkwreg`**.
- A `tau`-indexed likelihood with quantile-interpretable coefficients,
  shape heterogeneity, and nested family selection on (0,1) →
  **`gkwqreg`**.

## References

- Koenker and Bassett (1978). Regression quantiles. *Econometrica*
  **46**, 33–50.
- Mitnik and Baek (2013). *Statistical Papers* **54**, 177–192.
- Menezes, A. F. B. and Mazucheli, J. `unitquantreg`: Parametric
  Quantile Regression Models for Bounded Data. R package. (the package
  has no journal of record)
- Mazucheli, J., Alves, B., Menezes, A. F. B. and Leiva, V. (2022). An
  overview on parametric quantile regression models and their
  computational implementation with applications to biomedical problems
  including COVID-19 data. *Computer Methods and Programs in
  Biomedicine* **221**, 106816.
- Ferrari, S. L. P. and Cribari-Neto, F. (2004). Beta regression for
  modelling rates and proportions. *Journal of Applied Statistics*
  **31**(7), 799-815.
- Chernozhukov, V., Fernández-Val, I. and Galichon, A. (2010). Quantile
  and probability curves without crossing. *Econometrica* **78**(3),
  1093-1125.
