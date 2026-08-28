# Covariance matrix of a quantile regression fit

Returns the estimated covariance matrix of the full stacked coefficient
vector – the quantile part and every nuisance part together – either
from the observed information alone or in sandwich form.

## Usage

``` r
# S3 method for class 'gkwqreg'
vcov(object, type = c("expected", "observed", "sandwich"), ...)
```

## Arguments

- object:

  A `"gkwqreg"` fit, as returned by
  [`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md).

- type:

  Which covariance estimator to return; see the table below. Partial
  matching applies. All three require that the fit carried out the
  Hessian computation, that is,
  `control = gkwq_control(hessian = TRUE)`, which is the default;
  otherwise the function errors with instructions to refit.

- ...:

  Currently unused; present for consistency with the generic.

## Value

A symmetric numeric matrix of dimension `p` by `p`, where
`p = object$npar` is the total number of estimated coefficients across
all parts. Both dimnames are `names(coef(object))`, whose entries have
the form `"part:term"` – `"mu:(Intercept)"`, `"mu:x1"`, `"alpha:x1"` and
so on – so blocks can be extracted by
[`grep()`](https://rdrr.io/r/base/grep.html) on the part prefix.

## Details

All three types are built from the same observed information matrix
\$\$\hat{H} \\=\\ \left. \frac{\partial^2 \\ (-\ell)} {\partial\theta \\
\partial\theta^{\top}} \right\|\_{\theta = \hat\theta},\$\$ evaluated by
[`stats::optimHess()`](https://rdrr.io/r/stats/optim.html) applied to
the exact automatic- differentiation gradient supplied by `TMB` and
stored in the fit at estimation time. Two routes that look shorter are
deliberately not taken. The Hessian is never read off `obj$he()`, which
would require second-order differentiation through the incomplete-beta
atomic; and it is never obtained from the naive \\J^{\top} H J\\
transformation of the unreparametrized Hessian, which omits the
curvature term \\\sum_k g_k \\ \nabla^2 \theta_k\\ and is simply wrong
once the parameters vary with covariates.

## The three types

|  |  |  |
|----|----|----|
| `type` | Estimator | Appropriate when |
| `"expected"` (default) | \\\hat{H}^{-1}\\ | the family is taken to be correctly specified |
| `"observed"` | \\\hat{H}^{-1}\\, identical to the above | provided for interface compatibility |
| `"sandwich"` | \\\hat{H}^{-1} \left(\sum_i s_i s_i^{\top}\right) \hat{H}^{-1}\\ | the conditional distribution may be misspecified |

`"expected"` and `"observed"` return **the same matrix**. Nothing here
evaluates a Fisher information by expectation: the two names are kept so
that code written for other regression classes runs unchanged, and the
label is recorded in [`summary()`](https://rdrr.io/r/base/summary.html)
for the record. If the distinction matters to the argument being made,
say "observed information" and mean it.

The sandwich form replaces the middle of the expression by the empirical
outer product of the per-observation scores \\s_i\\ returned by
[`estfun.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/estfun.gkwqreg.md).
It is consistent for the covariance of the maximum likelihood estimator
when the family is wrong but the quantile equation is right, at the cost
of a noisier estimate in small samples. A large gap between the two is
itself informative: under correct specification the information equality
makes them agree up to sampling error, so a systematic discrepancy is
evidence against the family.

## Why the whole matrix matters

The conditional quantile is **not** orthogonal to the remaining
parameters in the Cox-Reid sense. The information matrix has
non-negligible off-diagonal blocks linking \\\beta\_\mu\\ to the
nuisance coefficients, so uncertainty in the shape parameters propagates
into the quantile coefficients and vice versa. Two consequences follow.
Any linear or non-linear function that mixes parts – a contrast, a
marginal effect, a predicted quantile – must be given the full matrix,
not a per-part block. And per-part standard errors reported in isolation
understate what is actually known jointly; use
[`confint.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/confint.gkwqreg.md)
with `method = "profile"` or `"boot"` when the Wald approximation is
doing too much work.

## See also

[`estfun.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/estfun.gkwqreg.md)
and
[`bread.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/bread.gkwqreg.md)
for the ingredients of the sandwich and for interoperability with the
sandwich package,
[`confint.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/confint.gkwqreg.md)
for intervals that do not rely on the Wald approximation,
[`summary.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg-extractors.md)
which uses this to build its coefficient table.

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

V <- vcov(fit)
dim(V)                                   # 5 x 5: three mu, two alpha
#> [1] 5 5
rownames(V)                              # "mu:(Intercept)" ... "alpha:x1"
#> [1] "mu:(Intercept)"    "mu:x1"             "mu:x21"           
#> [4] "alpha:(Intercept)" "alpha:x1"         
identical(V, vcov(fit, type = "observed"))   # TRUE, by construction
#> [1] TRUE

## -- the quantile block is NOT orthogonal to the nuisance block ----------
round(cov2cor(V), 3)
#>                   mu:(Intercept) mu:x1 mu:x21 alpha:(Intercept) alpha:x1
#> mu:(Intercept)             1.000 0.359 -0.711             0.319    0.107
#> mu:x1                      0.359 1.000  0.037             0.191    0.468
#> mu:x21                    -0.711 0.037  1.000            -0.009    0.045
#> alpha:(Intercept)          0.319 0.191 -0.009             1.000    0.426
#> alpha:x1                   0.107 0.468  0.045             0.426    1.000
## cor(mu:(Intercept), alpha:(Intercept)) = 0.319 and
## cor(mu:x1, alpha:x1) = 0.468.  Uncertainty about the shape parameter is
## not separable from uncertainty about the median; this is why anything
## mixing parts must be handed the whole matrix.

## -- model-based versus sandwich standard errors -------------------------
Vs <- vcov(fit, type = "sandwich")
round(cbind(model = sqrt(diag(V)), sandwich = sqrt(diag(Vs)),
            ratio = sqrt(diag(Vs) / diag(V))), 4)
#>                    model sandwich  ratio
#> mu:(Intercept)    0.0844   0.0882 1.0448
#> mu:x1             0.0468   0.0467 0.9963
#> mu:x21            0.0914   0.1018 1.1138
#> alpha:(Intercept) 0.0642   0.0583 0.9086
#> alpha:x1          0.0504   0.0468 0.9274
## Ratios between 0.91 and 1.11 here.  The family is correctly specified in
## this simulation, so the information equality holds and the two agree up
## to sampling error, exactly as it should.

## -- a Wald table built by hand from whichever matrix you trust ----------
est <- coef(fit); se <- sqrt(diag(Vs))
round(cbind(Estimate = est, `Std. Error` = se, z = est / se,
            `Pr(>|z|)` = 2 * pnorm(-abs(est / se))), 4)
#>                   Estimate Std. Error       z Pr(>|z|)
#> mu:(Intercept)      0.3045     0.0882  3.4512   0.0006
#> mu:x1               0.8578     0.0467 18.3826   0.0000
#> mu:x21             -0.5215     0.1018 -5.1223   0.0000
#> alpha:(Intercept)   0.7922     0.0583 13.5857   0.0000
#> alpha:x1            0.0044     0.0468  0.0935   0.9255
```
