# Marginal effects on the quantile scale

Translates the coefficients of the conditional quantile, which are
effects on the *log quantile odds*, into effects on the quantile itself,
with delta-method standard errors and confidence intervals.

## Usage

``` r
marginal_effects(
  object,
  variables = NULL,
  at = c("ame", "mem", "observed"),
  newdata = NULL,
  level = 0.95,
  vcov. = NULL,
  ...
)
```

## Arguments

- object:

  A `"gkwqreg"` fit, as returned by
  [`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md).

- variables:

  Character vector naming the covariates to report, using the
  coefficient names of the quantile part (`"x1"`, `"x21"` for the second
  level of a factor `x2`, and so on). Defaults to every non-intercept
  covariate in the quantile part. Names not present are dropped
  silently; if nothing remains the function errors.

- at:

  Where to evaluate the derivative: `"ame"` averages it over the rows,
  `"mem"` evaluates it at the average linear predictor, and `"observed"`
  returns the untouched per-row vector for each covariate. See the table
  below.

- newdata:

  Optional data frame over which to average, defaulting to the
  estimation data. Useful for reporting the effect for a subpopulation,
  or for a counterfactual covariate distribution, without refitting.

- level:

  Confidence level for the reported interval. The interval is Wald on
  the effect scale, `estimate` \\\pm\\ `z * std.error`, and is ignored
  when `at = "observed"`.

- vcov.:

  Optional covariance matrix for the full stacked coefficient vector,
  dimensioned and named as `vcov(object)`. Typical choices are
  `vcov(object, type = "sandwich")` or a clustered matrix such as
  [`sandwich::vcovCL()`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html).
  When `NULL` the model-based `vcov(object)` is used, and if that is
  unavailable the standard errors come back as `NA` rather than the
  function failing.

- ...:

  Currently unused.

## Value

For `at = "ame"` or `"mem"`, an object of class `"gkwq_meff"`: a list
with components `table` (a data frame with one row per covariate and
columns `variable`, `effect`, `std.error`, `lower`, `upper`), `at`,
`tau`, `level`, `link`, and `overlap` (the covariates that also enter a
nuisance part). It has a `print` method; extract `x$table` for further
computation.

For `at = "observed"`, an object of class `"gkwq_meff_observed"`: a list
with components `effects` (an `n` by `length(variables)` numeric matrix
of per-row derivatives, columns named after the covariates), `at` and
`tau`. There is no print method for this class; use `x$effects`
directly.

## Details

Under the default logit link the quantile part of the model is \$\$\log
\frac{\mu\_\tau(x)}{1 - \mu\_\tau(x)} \\=\\ x^{\top}\beta\_\mu ,\$\$ so
\\\exp(\beta_j)\\ is the multiplicative effect of a one-unit increase in
\\x_j\\ on the **quantile odds** \\\mu\_\tau/(1-\mu\_\tau)\\. It is
*not* the odds ratio of an event, and it is *not* an effect on a mean.
Reported on its own it is close to uninterpretable, because a fixed
multiplicative change in the odds displaces the quantile by an amount
that depends entirely on where that quantile already sits.

What is interpretable on the scale of the response is the derivative
\$\$\frac{\partial Q(\tau \mid x)}{\partial x_j} \\=\\ \beta_j \\
\frac{d\mu}{d\eta} \\=\\ \beta_j \\ \mu\_\tau (1 - \mu\_\tau) \quad
\text{(logit)},\$\$ which is what this function reports. For links other
than the logit the factor \\d\mu/d\eta\\ changes accordingly and
everything else is unaltered.

For a factor the derivative is still the quantity reported, evaluated at
the dummy column: it is the *marginal* effect of that indicator on the
quantile, which for a binary regressor is an accurate local
approximation to the discrete contrast rather than the contrast itself.

## Where the effect is evaluated

|  |  |  |
|----|----|----|
| `at` | What is reported | Shape |
| `"ame"` (default) | average marginal effect: the derivative averaged over the rows | one number per covariate, with a standard error |
| `"mem"` | marginal effect at the mean: the derivative evaluated at the average linear predictor | one number per covariate, with a standard error |
| `"observed"` | the whole vector of per-row derivatives | an `n` by `p` matrix, no standard errors |

`"ame"` is the default because it answers the question actually being
asked of a bounded response: how much does the quantile move, on average
over the population at hand? `"mem"` is cheaper to explain but describes
a unit whose covariates are all at their mean, which need not resemble
anybody. `"observed"` is the raw material for a plot of how the effect
varies across the covariate space; its column means reproduce the
`"ame"` estimates exactly.

## How the standard errors are obtained

For `at = "ame"` the reported effect is
\$\$\widehat{\mathrm{AME}}\_j(\beta) \\=\\ \beta_j \cdot
\frac{1}{n}\sum\_{i=1}^{n} g'\\\left(x_i^{\top}\beta\_\mu +
o_i\right),\$\$ and for `at = "mem"` the same expression with the
averaging done inside, \\g'(\bar\eta)\\. **Both** factors are functions
of the coefficient vector: the slope \\\beta_j\\ directly, and the
averaged derivative through \\\eta = X\beta\_\mu\\. The standard error
is therefore the exact delta method taken over the whole stacked
coefficient vector,
\$\$\widehat{\mathrm{Var}}\big(\widehat{\mathrm{AME}}\big) \\=\\ J \\ V
\\ J^{\top}, \qquad J \\=\\ \partial\\\widehat{\mathrm{AME}} /
\partial\beta^{\top},\$\$ with \\J\\ evaluated numerically by
[`numDeriv::jacobian()`](https://rdrr.io/pkg/numDeriv/man/jacobian.html)
and \\V\\ taken from `vcov.` when supplied and from `vcov(object)`
otherwise.

Freezing \\g'(\eta)\\ at its sample average and differentiating only
\\\beta_j\\ – the obvious shortcut, which yields
\\\mathrm{s.e.}(\beta_j) \times \overline{g'}\\ – is not the delta
method. It discards the second term, which partially cancels the first:
raising \\\beta_j\\ spreads \\\eta\\, and a more spread \\\eta\\ has a
smaller \\\overline{g'}\\ once it is away from the mode of the link
derivative. The shortcut consequently *overstates* the standard error,
badly for a strong covariate; the examples show the discrepancy on a fit
where it is a factor of about 1.6.

Whether the interval deserves to be believed is a separate question from
whether it is computed correctly. Pass
`vcov. = vcov(object, type = "sandwich")` if the conditional
distribution may be misspecified, or a clustered covariance matrix if
the data are grouped; the delta method is applied to whatever matrix is
supplied.

## Covariates that also enter a nuisance part

If \\x_j\\ also appears in the equation for one of the remaining
parameters, the effect on the quantile is still exactly \\\beta_j \\
d\mu/d\eta\\, with no extra term. This is a genuine convenience of the
parametrization rather than an approximation: \\\mu\_\tau\\ *is* the
\\\tau\\-quantile, so a nuisance part alters the spread and shape of the
conditional distribution while leaving the quantile being modelled
untouched. The print method says so when it detects the overlap. Note
that this is a statement about the \\\tau\\-quantile only: the same
covariate does move *other* quantiles through the nuisance part.

## See also

[`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md)
for the model and the interpretation of its coefficients,
[`predict.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/predict.gkwqreg.md)
for predictions,
[`vcov.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/vcov.gkwqreg.md)
for the covariance matrices that can be passed to `vcov.`.

## Examples

``` r
## Same simulated design as predict.gkwqreg(): a conditional median that
## follows a logit model, drawn by inverse transform.
set.seed(2024)
n   <- 300
x1  <- runif(n, -2, 2)
x2  <- rbinom(n, 1, 0.5)
mu  <- plogis(0.3 + 0.9 * x1 - 0.5 * x2)
bt  <- log1p(-0.5) / log1p(-mu^2)
y   <- gkwq_quantile(runif(n), alpha = 2, beta = bt)
dat <- data.frame(y = y, x1 = x1, x2 = factor(x2))

fit <- gkwqreg(y ~ x1 + x2 | x1, data = dat, tau = 0.5, family = "kw")

## -- what the coefficient says, and what it does NOT say -----------------
exp(coef(fit, part = "mu"))
#> (Intercept)          x1         x21 
#>    1.355932    2.357993    0.593646 
## x1 -> 2.3580.  A one-unit rise in x1 multiplies the odds of the median,
## m/(1-m), by 2.36.  That is not an event odds ratio and not a mean effect.

me <- marginal_effects(fit)
me
#> 
#> Marginal effects on the conditional 0.5-quantile (averaged over the sample)
#> dQ(tau|x)/dx_j, logit link
#> 
#>  variable  effect std.error   lower    upper
#>        x1  0.1683  0.005758  0.1571  0.17962
#>       x21 -0.1023  0.017803 -0.1372 -0.06744
#> 
#> Note: x1 also appear(s) in a nuisance part. The effect on the
#>       quantile is still exactly the above: mu IS the quantile, so the
#>       nuisance part changes the spread, not the quantile modelled.
## x1  effect  0.16834  (s.e. 0.005758)
## x21 effect -0.10233  (s.e. 0.017803)
## A one-unit rise in x1 lifts the conditional MEDIAN by 0.168 on average,
## on the (0,1) scale of the response.  The factor-2.36 odds statement and
## the 0.168 response-scale statement describe the same coefficient.

## -- the two are linked by the chain rule --------------------------------
b   <- coef(fit)[["mu:x1"]]
muh <- fitted(fit)                       # fitted conditional medians
b * mean(muh * (1 - muh))                # 0.1683372
#> [1] 0.1683372
me$table$effect[me$table$variable == "x1"]   # identical
#> [1] 0.1683372

## -- the effect is not constant across the covariate space ---------------
obs <- marginal_effects(fit, at = "observed")
range(obs$effects[, "x1"])               # 0.0901 to 0.2144
#> [1] 0.09005712 0.21443906
colMeans(obs$effects)                    # reproduces the AME exactly
#>         x1        x21 
#>  0.1683372 -0.1023339 
## The AME is an average of effects that differ by more than twofold across
## the sample, which is precisely why the odds-scale coefficient is constant
## and the response-scale effect is not.

## -- at the covariate means instead --------------------------------------
marginal_effects(fit, at = "mem")
#> 
#> Marginal effects on the conditional 0.5-quantile (at covariate means)
#> dQ(tau|x)/dx_j, logit link
#> 
#>  variable  effect std.error   lower    upper
#>        x1  0.2143   0.01188  0.1910  0.23757
#>       x21 -0.1303   0.02286 -0.1751 -0.08547
#> 
#> Note: x1 also appear(s) in a nuisance part. The effect on the
#>       quantile is still exactly the above: mu IS the quantile, so the
#>       nuisance part changes the spread, not the quantile modelled.
## Larger in absolute value here (0.2143 for x1): the average linear
## predictor sits nearer the middle of the logistic curve, where mu(1-mu)
## attains its maximum, than the typical observation does.

## -- the delta method versus the shortcut --------------------------------
se_shortcut <- sqrt(diag(vcov(fit)))[["mu:x1"]] * mean(muh * (1 - muh))
se_exact    <- me$table$std.error[me$table$variable == "x1"]
c(shortcut = se_shortcut, exact = se_exact, ratio = se_shortcut / se_exact)
#>    shortcut       exact       ratio 
#> 0.009191738 0.005758051 1.596327912 
## shortcut 0.009192, exact 0.005758, ratio 1.596.  Freezing mu(1-mu) at its
## average ignores that it, too, depends on the coefficients, and the two
## dependencies partially cancel.

## -- robust and subpopulation variants -----------------------------------
marginal_effects(fit, vcov. = vcov(fit, type = "sandwich"))
#> 
#> Marginal effects on the conditional 0.5-quantile (averaged over the sample)
#> dQ(tau|x)/dx_j, logit link
#> 
#>  variable  effect std.error   lower    upper
#>        x1  0.1683  0.005836  0.1569  0.17978
#>       x21 -0.1023  0.020095 -0.1417 -0.06295
#> 
#> Note: x1 also appear(s) in a nuisance part. The effect on the
#>       quantile is still exactly the above: mu IS the quantile, so the
#>       nuisance part changes the spread, not the quantile modelled.
marginal_effects(fit, variables = "x1", newdata = subset(dat, x2 == "1"))
#> 
#> Marginal effects on the conditional 0.5-quantile (averaged over the sample)
#> dQ(tau|x)/dx_j, logit link
#> 
#>  variable effect std.error  lower  upper
#>        x1 0.1675  0.006011 0.1557 0.1793
#> 
#> Note: x1 also appear(s) in a nuisance part. The effect on the
#>       quantile is still exactly the above: mu IS the quantile, so the
#>       nuisance part changes the spread, not the quantile modelled.
```
