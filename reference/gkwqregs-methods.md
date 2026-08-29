# Methods for a set of quantile regression fits

Accessors for the `"gkwqregs"` container returned by
[`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md)
when `tau` has length greater than one. The container holds one entirely
independent fit per quantile level in `object$fits`, and these methods
arrange their output with one column per level.

## Usage

``` r
# S3 method for class 'gkwqregs'
coef(object, ...)

# S3 method for class 'gkwqregs'
fitted(object, ...)

# S3 method for class 'gkwqregs'
logLik(object, ...)

# S3 method for class 'gkwqregs'
summary(object, ...)

# S3 method for class 'gkwqregs'
predict(object, newdata = NULL, ...)

# S3 method for class 'gkwqregs'
residuals(object, ...)
```

## Arguments

- object:

  A `"gkwqregs"` container, as returned by
  [`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md)
  with a vector-valued `tau`.

- ...:

  Passed on to the corresponding method for each individual fit.

- newdata:

  Used by [`predict()`](https://rdrr.io/r/stats/predict.html) only: an
  optional data frame in which to evaluate every fit, with the same
  requirements as
  [`predict.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/predict.gkwqreg.md).
  `NULL`, the default, uses the estimation data.

## Value

- [`coef()`](https://rdrr.io/r/stats/coef.html): a numeric matrix with
  one row per coefficient and one column per quantile level, the columns
  named by level.

- [`fitted()`](https://rdrr.io/r/stats/fitted.values.html): a numeric
  matrix with one row per observation and one column per level, holding
  each fit's conditional quantiles.

- [`summary()`](https://rdrr.io/r/base/summary.html): a list of
  `"summary.gkwqreg"` objects, one per level.

- [`predict()`](https://rdrr.io/r/stats/predict.html): a numeric matrix
  with one row per prediction and one column per level, from each fit in
  turn.

- [`residuals()`](https://rdrr.io/r/stats/residuals.html): a numeric
  matrix with one row per observation and one column per level.

- [`logLik()`](https://rdrr.io/r/stats/logLik.html): never returns; it
  raises an error explaining why.

## Details

Levels are stored sorted and de-duplicated, and each fit is estimated on
its own. Nothing ties the levels together, which is what makes quantile
crossing possible and worth testing for; see
[`check_crossing()`](https://evandeilton.github.io/gkwqreg/reference/check_crossing.md)
and
[`rearrange()`](https://evandeilton.github.io/gkwqreg/reference/rearrange.md).

[`logLik()`](https://rdrr.io/r/stats/logLik.html) deliberately fails.
The container holds one log-likelihood per level, and those are neither
comparable across levels nor additive: summing them would treat the same
observations as independent evidence several times over, and comparing
them would compare answers to different questions. Take the likelihood
of an individual fit instead, as the error message says.

## See also

[`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md),
[`quantile_process()`](https://evandeilton.github.io/gkwqreg/reference/quantile_process.md)
for the quantile process on a grid,
[`check_crossing()`](https://evandeilton.github.io/gkwqreg/reference/check_crossing.md)
for crossing diagnostics,
[`predict.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/predict.gkwqreg.md)
for reading several levels off a *single* fit.

## Examples

``` r
set.seed(2024)
n   <- 300
x1  <- runif(n, -2, 2)
mu  <- plogis(0.3 + 0.9 * x1)
bt  <- log1p(-0.5) / log1p(-mu^2)
y   <- gkwq_quantile(runif(n), alpha = 2, beta = bt)
dat <- data.frame(y = y, x1 = x1)

fits <- gkwqreg(y ~ x1, data = dat, tau = c(0.25, 0.5, 0.75),
                family = "kw")
fits
#> 
#> 3 Generalized Kumaraswamy quantile regressions (family: kw, anchor: beta)
#> 
#> Call:
#> gkwqreg(formula = y ~ x1, data = dat, tau = c(0.25, 0.5, 0.75), 
#>     family = "kw")
#> 
#> Coefficients by tau:
#>                      0.25   0.50   0.75
#> mu:(Intercept)    -0.3503 0.4088 1.2208
#> mu:x1              0.6959 0.8630 1.0947
#> alpha:(Intercept)  0.7935 0.7908 0.7818

round(coef(fits), 4)
#>                      0.25   0.50   0.75
#> mu:(Intercept)    -0.3503 0.4088 1.2208
#> mu:x1              0.6959 0.8630 1.0947
#> alpha:(Intercept)  0.7935 0.7908 0.7818
## One column per level. The mu:x1 row is the covariate effect on the log
## quantile odds at each level; its drift across columns is what a quantile
## regression is for.

Q <- fitted(fits)
dim(Q)
#> [1] 300   3
colnames(Q)                   # the levels, formatted: "0.25" "0.50" "0.75"
#> [1] "0.25" "0.50" "0.75"
colMeans(dat$y <= Q)          # 0.25, 0.46, 0.73: coverage at each level
#> 0.25 0.50 0.75 
#> 0.25 0.46 0.73 

## Each level is a separate fit, so nothing forces the columns to be
## ordered; check rather than assume.
mean(Q[, 1] <= Q[, 2] & Q[, 2] <= Q[, 3])   # 1: no crossing in this sample
#> [1] 1
```
