# Bread matrix for the sandwich covariance

Returns the inverse observed information scaled by the sample size, in
the normalization the sandwich package expects. Together with
[`estfun.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/estfun.gkwqreg.md)
this makes
[`sandwich::sandwich()`](https://zeileis.codeberg.page/sandwich/reference/sandwich.html),
[`sandwich::vcovCL()`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html)
and [`lmtest::coeftest()`](https://rdrr.io/pkg/lmtest/man/coeftest.html)
work on `"gkwqreg"` fits without any further glue.

## Usage

``` r
bread.gkwqreg(x, ...)
```

## Arguments

- x:

  A `"gkwqreg"` fit that carries an observed information matrix, that
  is, one fitted with `gkwq_control(hessian = TRUE)` (the default).
  Otherwise the function errors.

- ...:

  Currently unused.

## Value

A symmetric numeric matrix of dimension `p` by `p`, where
`p = length(coef(x))`, with dimnames inherited from `vcov(x)`.

## Details

The value is \$\$\widehat{\mathrm{bread}} \\=\\ n \\ \hat{H}^{-1} \\=\\
\left(\frac{1}{n}\hat{H}\right)^{-1},\$\$ the inverse of the *average*
observed information, which is the quantity that converges to a fixed
matrix as \\n\\ grows. That normalization is what lets sandwich assemble
\\\widehat{\mathrm{bread}} \\ \widehat{\mathrm{meat}} \\
\widehat{\mathrm{bread}} / n\\ and recover the same estimator as
`vcov(object, type = "sandwich")`.

[`sandwich::vcovHC()`](https://zeileis.codeberg.page/sandwich/reference/vcovHC.html)
is not supported for this class; see
[`estfun.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/estfun.gkwqreg.md)
for why.

## See also

[`estfun.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/estfun.gkwqreg.md),
[`vcov.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/vcov.gkwqreg.md).

## Examples

``` r
set.seed(2024)
n   <- 300
x1  <- runif(n, -2, 2)
mu  <- plogis(0.3 + 0.9 * x1)
bt  <- log1p(-0.5) / log1p(-mu^2)
y   <- gkwq_quantile(runif(n), alpha = 2, beta = bt)
dat <- data.frame(y = y, x1 = x1)

fit <- gkwqreg(y ~ x1, data = dat, tau = 0.5, family = "kw")

B <- bread.gkwqreg(fit)
all.equal(B / nobs(fit), vcov(fit), check.attributes = FALSE)   # TRUE
#> [1] TRUE

## The bread alone gives the model-based standard errors back:
round(sqrt(diag(B / nobs(fit))), 4)
#>    mu:(Intercept)             mu:x1 alpha:(Intercept) 
#>            0.0665            0.0460            0.0642 
round(sqrt(diag(vcov(fit))), 4)
#>    mu:(Intercept)             mu:x1 alpha:(Intercept) 
#>            0.0665            0.0460            0.0642 
```
