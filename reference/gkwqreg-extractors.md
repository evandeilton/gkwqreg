# Extractor methods for a quantile regression fit

The standard model accessors, adapted to a model whose coefficients are
split across several parts.
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`terms()`](https://rdrr.io/r/stats/terms.html) and
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) therefore
take a `part` argument; the rest behave as their generics lead one to
expect.

## Usage

``` r
# S3 method for class 'gkwqreg'
coef(object, part = NULL, ...)

# S3 method for class 'gkwqreg'
logLik(object, ...)

# S3 method for class 'gkwqreg'
nobs(object, ...)

# S3 method for class 'gkwqreg'
AIC(object, ..., k = 2)

# S3 method for class 'gkwqreg'
BIC(object, ...)

# S3 method for class 'gkwqreg'
family(object, ...)

# S3 method for class 'gkwqreg'
formula(x, ...)

# S3 method for class 'gkwqreg'
terms(x, part = "mu", ...)

# S3 method for class 'gkwqreg'
model.frame(formula, ...)

# S3 method for class 'gkwqreg'
model.matrix(object, part = "mu", ...)

# S3 method for class 'gkwqreg'
getCall(x, ...)

# S3 method for class 'gkwqreg'
update(object, formula., ..., evaluate = TRUE)

# S3 method for class 'gkwqreg'
summary(object, level = 0.95, vcov_type = NULL, ...)
```

## Arguments

- object, x, formula:

  A `"gkwqreg"` fit, as returned by
  [`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md).
  The three names correspond to the argument names of the respective
  generics.

- part:

  Which model part to extract, named as in
  [`gkwq_parts()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_parts.md):
  `"mu"` for the conditional quantile, then the family's remaining
  parameters. Partial matching applies.
  [`coef()`](https://rdrr.io/r/stats/coef.html) accepts a vector of
  parts and returns a named list in that case, and `part = NULL` (its
  default) returns the whole stacked coefficient vector.
  [`terms()`](https://rdrr.io/r/stats/terms.html) and
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) take
  exactly one part and default to `"mu"`.

- ...:

  For [`update()`](https://rdrr.io/r/stats/update.html), further named
  arguments replacing those of the original call. Unused by the other
  methods.

- k:

  The penalty per parameter in
  [`AIC()`](https://rdrr.io/r/stats/AIC.html). The default `k = 2` gives
  the usual Akaike criterion; `k = log(nobs(object))` reproduces
  [`BIC()`](https://rdrr.io/r/stats/AIC.html).

- formula.:

  A formula update specification for
  [`update()`](https://rdrr.io/r/stats/update.html), in the sense of
  [`stats::update.formula()`](https://rdrr.io/r/stats/update.formula.html),
  applied to the multi-part model formula.

- evaluate:

  If `TRUE` (the default)
  [`update()`](https://rdrr.io/r/stats/update.html) evaluates the
  modified call and returns the new fit; if `FALSE` it returns the
  unevaluated call, which is useful for inspection or for deferred
  evaluation.

- level:

  Confidence level recorded by
  [`summary()`](https://rdrr.io/r/base/summary.html) and used when its
  result is printed or post-processed.

- vcov_type:

  Which covariance estimator
  [`summary()`](https://rdrr.io/r/base/summary.html) should use for its
  standard errors, one of `"expected"`, `"observed"` or `"sandwich"`;
  see
  [`vcov.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/vcov.gkwqreg.md).
  `NULL`, the default, uses the type recorded in the fit's control
  object. If the requested type cannot be computed,
  [`summary()`](https://rdrr.io/r/base/summary.html) falls back to the
  stored matrix rather than failing.

## Value

- [`coef()`](https://rdrr.io/r/stats/coef.html): a named numeric vector.
  The whole stacked vector, with names `"part:term"`, when `part` is
  `NULL`; a single part's coefficients, with bare term names, when one
  part is given; a named list of such vectors when several are given.

- [`logLik()`](https://rdrr.io/r/stats/logLik.html): an object of class
  `"logLik"` with attributes `df` (the number of estimated coefficients)
  and `nobs`.

- [`nobs()`](https://rdrr.io/r/stats/nobs.html): a single integer, the
  number of observations used in the fit.

- [`AIC()`](https://rdrr.io/r/stats/AIC.html),
  [`BIC()`](https://rdrr.io/r/stats/AIC.html): a single number.

- [`family()`](https://rdrr.io/r/stats/family.html): an object of class
  `"gkwq_family"`, a list with components `family`, `anchor`, `tau`,
  `parts` and `link`, with a `print` method.

- [`formula()`](https://rdrr.io/r/stats/formula.html): the multi-part
  model formula as supplied.

- [`terms()`](https://rdrr.io/r/stats/terms.html): the `"terms"` object
  of one part.

- [`model.frame()`](https://rdrr.io/r/stats/model.frame.html): the
  stored model frame.

- [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html): the
  design matrix of one part, with an `"assign"` attribute and contrasts
  as recorded at fit time.

- [`getCall()`](https://rdrr.io/r/stats/update.html): the matched call
  that produced the fit.

- [`update()`](https://rdrr.io/r/stats/update.html): a new `"gkwqreg"`
  fit, or the unevaluated call when `evaluate = FALSE`.

- [`summary()`](https://rdrr.io/r/base/summary.html): an object of class
  `"summary.gkwqreg"` carrying the coefficient table, the fit
  statistics, the pinball loss, the pseudo-R1, the empirical coverage
  and the information condition number, with a `print` method.

## Details

A `"gkwqreg"` object carries one linear predictor per modelled
parameter. The first part is always `mu`, the conditional
\\\tau\\-quantile; the remaining parts are the family's free parameters
with the anchored one removed, as reported by
[`gkwq_parts()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_parts.md).
Coefficients are stored in a single stacked vector whose names have the
form `"part:term"`, so `coef(object)` returns everything at once and
`coef(object, part = "mu")` returns just the quantile block with its
bare term names.

[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`BIC()`](https://rdrr.io/r/stats/AIC.html) describe the model at
**one** quantile level. They are comparable across families, anchors and
covariate sets fitted at the same `tau`, and are not comparable across
levels: different levels answer different questions of the same data and
their likelihoods are not on a common scale.

[`family()`](https://rdrr.io/r/stats/family.html) returns a small
description object rather than a
[`stats::family()`](https://rdrr.io/r/stats/family.html) object, because
a Generalized Kumaraswamy quantile model has no single variance function
or canonical link to report. It records the family, the anchor, the
quantile level, the part names and the link used by each part, and
prints them compactly.

Some accessors depend on what the fit was asked to retain.
[`model.frame()`](https://rdrr.io/r/stats/model.frame.html) requires
`model = TRUE` (the default) and errors otherwise;
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) uses the
stored design matrices when the fit was made with `x = TRUE`, and
otherwise rebuilds them from the model frame, so it needs one or the
other. Each error message names the argument to set on the refit.

## What [`update()`](https://rdrr.io/r/stats/update.html) does and does not carry over

[`update()`](https://rdrr.io/r/stats/update.html) modifies the recorded
call and re-evaluates it, so any argument of
[`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md)
can be replaced by name – `tau`, `family`, `anchor`, `data`, `control`
and so on. Changing `tau` refits at the new level, which is a different
model and not the same thing as reading another level off the present
fit; see
[`predict.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/predict.gkwqreg.md)
for that distinction. Changing `family` or `anchor` may change the
number of formula parts the model expects, in which case the formula
must be updated in the same call.

## See also

[`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md)
for the model,
[`gkwq_parts()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_parts.md)
for the part contract,
[`vcov.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/vcov.gkwqreg.md),
[`confint.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/confint.gkwqreg.md),
[`fitted.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/fitted.gkwqreg.md).

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

## -- coefficients, whole and by part -------------------------------------
round(coef(fit), 4)                 # names are "part:term"
#>    mu:(Intercept)             mu:x1            mu:x21 alpha:(Intercept) 
#>            0.3045            0.8578           -0.5215            0.7922 
#>          alpha:x1 
#>            0.0044 
round(coef(fit, part = "mu"), 4)    # just the conditional-median equation
#> (Intercept)          x1         x21 
#>      0.3045      0.8578     -0.5215 
names(coef(fit, part = c("mu", "alpha")))
#> [1] "mu"    "alpha"

## -- what the model is -----------------------------------------------------
family(fit)
#> Generalized Kumaraswamy quantile family
#>   family kw | tau 0.5 | anchor beta
#>   parts  mu | alpha
gkwq_parts("kw")                    # the part contract this fit obeys
#> [1] "mu"    "alpha"
formula(fit)
#> y ~ x1 + x2 | x1
#> <environment: 0x56542f2790e8>

## -- fit statistics, all at this one quantile level ----------------------
c(logLik = as.numeric(logLik(fit)), df = attr(logLik(fit), "df"),
  n = nobs(fit), AIC = AIC(fit), BIC = BIC(fit))
#>    logLik        df         n       AIC       BIC 
#>  153.7594    5.0000  300.0000 -297.5188 -278.9999 
all.equal(AIC(fit, k = log(nobs(fit))), BIC(fit))   # TRUE
#> [1] TRUE

## -- design matrices differ by part --------------------------------------
fitx <- gkwqreg(y ~ x1 + x2 | x1, data = dat, tau = 0.5, family = "kw",
                x = TRUE)
colnames(model.matrix(fitx, part = "mu"))       # intercept, x1, x21
#> [1] "(Intercept)" "x1"          "x21"        
colnames(model.matrix(fitx, part = "alpha"))    # intercept, x1
#> [1] "(Intercept)" "x1"         
## This is why sandwich::vcovHC() cannot work here: there is no single
## design matrix whose columns match the score matrix.

## -- summary, with a covariance estimator of your choosing ---------------
summary(fit, vcov_type = "sandwich")
#> 
#> Generalized Kumaraswamy quantile regression
#> family: kw   tau: 0.5   anchor: beta
#> 
#> Call:
#> gkwqreg(formula = y ~ x1 + x2 | x1, data = dat, tau = 0.5, family = "kw")
#> 
#> Quantile residuals:
#>     Min      1Q  Median      3Q     Max 
#> -2.2285 -0.6799 -0.0364  0.6648  2.9745 
#> 
#> Conditional 0.5-quantile (link logit) -- coefficients are effects on the LOG QUANTILE ODDS log(mu/(1-mu)):
#>                Estimate Std. Error z value Pr(>|z|)    
#> mu:(Intercept)  0.30449    0.08823   3.451 0.000558 ***
#> mu:x1           0.85781    0.04666  18.383  < 2e-16 ***
#> mu:x21         -0.52147    0.10180  -5.122 3.02e-07 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> alpha (link log):
#>                   Estimate Std. Error z value Pr(>|z|)    
#> alpha:(Intercept) 0.792164   0.058309  13.586   <2e-16 ***
#> alpha:x1          0.004372   0.046762   0.094    0.926    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Anchored parameter: beta (computed from the quantile, not estimated)
#> Log-likelihood: 153.8 on 5 Df   AIC: -297.5   BIC: -279
#> Pinball loss: 0.07067   Pseudo-R1: 0.3959
#> Empirical coverage: 0.5233 (target 0.5)
#> Covariance: sandwich   Information condition number: 22.17

## -- update() re-evaluates the call --------------------------------------
getCall(fit)
#> gkwqreg(formula = y ~ x1 + x2 | x1, data = dat, tau = 0.5, family = "kw")
update(fit, . ~ x1 | x1, evaluate = FALSE)      # inspect before running
#> gkwqreg(formula = y ~ x1 | x1, data = dat, tau = 0.5, family = "kw")
fit_simpler <- update(fit, . ~ x1 | x1)
c(full = AIC(fit), simpler = AIC(fit_simpler))
#>      full   simpler 
#> -297.5188 -267.4854 
```
