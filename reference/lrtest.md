# Likelihood-ratio test for nested quantile regression fits

A thin wrapper on
[`anova.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/anova.gkwqreg.md),
offered under the name that users of
[`lmtest::lrtest()`](https://rdrr.io/pkg/lmtest/man/lrtest.html) will
reach for first. Behaviour, guards and return value are identical: the
fits are sorted by dimension, each is tested against the one above it by
the statistic \\LR = 2\\\ell_1(\hat\theta_1) - \ell_0(\hat\theta_0)\\\\
referred to a chi-squared distribution on the difference in dimension,
and comparisons across quantile levels or across anchors are refused.

## Usage

``` r
lrtest(object, ...)

# S3 method for class 'gkwqreg'
lrtest(object, ...)
```

## Arguments

- object:

  A `"gkwqreg"` fit.

- ...:

  One or more further `"gkwqreg"` fits, at the same quantile level, with
  the same anchor and the same number of observations.

## Value

A data frame of class `"anova.gkwqreg"`, exactly as returned by
[`anova.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/anova.gkwqreg.md);
see that help page for the columns and the printed layout.

## Details

`lrtest()` is a generic defined by this package, so that `gkwqreg` need
not depend on `lmtest` merely to offer a familiar name. It is a
*different* generic from
[`lmtest::lrtest()`](https://rdrr.io/pkg/lmtest/man/lrtest.html), and
the two mask one another when both packages are attached: whichever was
attached last wins.

The distinction matters, because the two behave differently on these
objects. If `lmtest` is attached after `gkwqreg`, a bare call to
`lrtest()` reaches `lmtest`'s default method, which knows nothing about
quantile levels or anchors. Handed two fits at different levels it will
happily compute a difference of log-likelihoods and report a significant
result – the very number that
[`anova.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/anova.gkwqreg.md)
refuses to produce, and for the reasons given there. Call
`gkwqreg::lrtest()` or
[`anova.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/anova.gkwqreg.md)
explicitly if there is any doubt about which generic is in scope; both
carry the guards.

## See also

[`anova.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/anova.gkwqreg.md)
for the full account of the test, the nesting of the families and the
two refusals;
[`vuong_test()`](https://evandeilton.github.io/gkwqreg/reference/vuong_test.md)
for the non-nested case.

## Examples

``` r
## Data with a genuine third shape parameter: lambda = 4, so "kw" is too
## small and "ekw" is right.
set.seed(2026)
n  <- 500
x  <- runif(n, -2, 2)
mu <- plogis(0.3 + 1.2 * x)
b  <- log1p(-0.5^(1 / 4)) / log1p(-mu^2)
y  <- gkwdist::rekw(n, alpha = 2, beta = b, lambda = 4)
d  <- data.frame(y = y, x = x)

m_kw  <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
m_ekw <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "ekw")

lrtest(m_kw, m_ekw)
#> Likelihood-ratio test for Generalized Kumaraswamy quantile regression
#> tau = 0.5, anchor = beta
#>      Df logLik     AIC  Chisq Chi Df Pr(>Chisq)    
#> [1,]  3 482.58 -959.16                             
#> [2,]  4 496.90 -985.80 28.641      1  8.713e-08 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##      Df logLik     AIC   Chisq Chi Df Pr(>Chisq)
## [1,]  3 482.58 -959.16
## [2,]  4 496.90 -985.80  28.641      1  8.713e-08

## Identical to anova(), which is what it calls.
identical(lrtest(m_kw, m_ekw), anova(m_kw, m_ekw))
#> [1] TRUE
## [1] TRUE

## The guards travel with it: a cross-level comparison is refused here too.
m_tau9 <- gkwqreg(y ~ x, data = d, tau = 0.9, family = "kw")
try(lrtest(m_kw, m_tau9))
#> Error : models at different quantile levels cannot be compared: each level is a separate likelihood. Levels seen: 0.5, 0.9.
```
