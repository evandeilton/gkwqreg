# Per-observation score contributions

Returns the empirical estimating functions of the fit: the gradient of
each observation's log-likelihood contribution with respect to the full
stacked coefficient vector. This is the "meat" of the sandwich
covariance, and supplying it in the conventional shape is what makes
[`sandwich::sandwich()`](https://zeileis.codeberg.page/sandwich/reference/sandwich.html),
[`sandwich::vcovCL()`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html)
and [`lmtest::coeftest()`](https://rdrr.io/pkg/lmtest/man/coeftest.html)
work on `"gkwqreg"` fits with no further glue.

## Usage

``` r
estfun.gkwqreg(x, ...)
```

## Arguments

- x:

  A `"gkwqreg"` fit, as returned by
  [`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md),
  still carrying its `TMB` object in `x$obj`.

- ...:

  Currently unused.

## Value

A numeric matrix with `nobs(x)` rows and `length(coef(x))` columns, the
columns named as `names(coef(x))` (`"part:term"`). Row `i` holds the
score contribution of observation `i`.

## Details

The `i`th row is \$\$s_i(\hat\theta) \\=\\ \left. \frac{\partial \\ \log
f(y_i \mid x_i; \theta)}{\partial \theta} \right\|\_{\theta =
\hat\theta},\$\$ obtained by applying
[`numDeriv::jacobian()`](https://rdrr.io/pkg/numDeriv/man/jacobian.html)
to the vector of per-observation log-likelihood contributions that the
`TMB` object reports. Numerical differentiation of the reported
likelihood is used rather than a second automatic-differentiation pass
because the tape returns the summed objective; the per-observation
decomposition lives in the report, not on the tape.

At the maximum the column sums are zero up to optimizer and
differentiation tolerance, which is a useful convergence check in its
own right: column sums that are not small say the reported optimum is
not a stationary point.

The computation goes through the fit's `TMB` object, held in
`object$obj`. If that component is absent – it is dropped by anything
that strips the fit down for storage – the function stops with a message
asking for a refit in the current session, rather than returning numbers
it cannot stand behind.

## What does not work

[`sandwich::vcovHC()`](https://zeileis.codeberg.page/sandwich/reference/vcovHC.html)
is **not** supported and will error with "cannot match dimension of
model.matrix and estfun". It needs working residuals and a single design
matrix whose columns line up one-for-one with the scores. This is a
multi-part model with one design matrix *per part*, so
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) returns
the `mu` block while `estfun()` spans every coefficient of every part,
and the two cannot be made to line up. There is no sensible thing for a
heteroskedasticity-consistent correction to do here. Use
`vcov(object, type = "sandwich")`, which is the same estimator computed
correctly, or
[`sandwich::vcovCL()`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html)
for clustered data.

## See also

[`bread.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/bread.gkwqreg.md)
for the other half of the sandwich,
[`vcov.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/vcov.gkwqreg.md)
for the assembled estimator.

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

S <- estfun.gkwqreg(fit)
dim(S)                          # 300 x 5: one row per observation
#> [1] 300   5
round(colSums(S), 6)            # all zero to 1e-4: a stationary point
#>    mu:(Intercept)             mu:x1            mu:x21 alpha:(Intercept) 
#>         -0.000142         -0.000045         -0.000103          0.000013 
#>          alpha:x1 
#>          0.000110 

## -- the sandwich, assembled by hand -------------------------------------
B <- bread.gkwqreg(fit)                       # n * H^{-1}
Vs <- B %*% crossprod(S) %*% B / nobs(fit)^2
all.equal(Vs, vcov(fit, type = "sandwich"), check.attributes = FALSE)
#> [1] TRUE
## TRUE.  vcov(type = "sandwich") is exactly H^{-1} (sum s_i s_i') H^{-1},
## and is also what sandwich::sandwich() returns for this class.

## -- the scores carry the clustering information -------------------------
grp <- rep(1:30, each = 10)
M   <- crossprod(rowsum(S, grp))              # clustered meat
Vcl <- (B %*% M %*% B) / nobs(fit)^2
round(sqrt(diag(Vcl)), 4)
#>    mu:(Intercept)             mu:x1            mu:x21 alpha:(Intercept) 
#>            0.0739            0.0415            0.1052            0.0537 
#>          alpha:x1 
#>            0.0443 
## The same construction sandwich::vcovCL() performs; shown here so the
## estimator is legible without taking a dependency on that package.
```
