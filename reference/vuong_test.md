# Vuong test for two non-nested quantile regression fits

Tests whether two models fitted to the same responses, at the same
quantile level, are equally close to the true data-generating
distribution. This is the right instrument for the two comparisons a
likelihood-ratio test cannot make: two anchors of the same family, and
two families that do not nest. Both produce models of equal dimension
fitted to the same data, for which no null hypothesis of restriction
exists and no degrees of freedom remain.

## Usage

``` r
vuong_test(object, object2, correction = TRUE)
```

## Arguments

- object, object2:

  Two `"gkwqreg"` fits at the same quantile level and over the same
  observations. The order matters only through the sign of the
  statistic: `object` is model 1 and `object2` is model 2. The two may
  differ in family, in anchor, in the covariates entering any part, or
  in all three.

- correction:

  Logical. If `TRUE`, the default, apply the BIC-style dimension
  correction described below, which penalizes whichever model estimates
  more coefficients. If `FALSE`, compare the raw per-observation
  log-likelihood contributions. The correction is exactly zero when the
  two fits have the same number of coefficients, which is the case for a
  comparison of two anchors.

## Value

An object of class `"gkwq_vuong"`: a list with components

- `statistic`:

  the value of \\Z\\.

- `p.value`:

  the two-sided p-value \\2\\\Phi(-\|Z\|)\\.

- `n`:

  the number of observations, common to both fits.

- `model1`,`model2`:

  character labels of the form `"family / anchor a"` identifying the two
  fits.

- `tau`:

  the common quantile level.

- `correction`:

  the value of `correction` that was used.

The [`print()`](https://rdrr.io/r/base/print.html) method reports the
two labels, the statistic, the p-value, and the verdict at the five per
cent level.

## Details

Let \\f_1\\ and \\f_2\\ denote the two fitted conditional densities and
\\h\\ the unknown true one. Vuong's null hypothesis is that the two
models are *equally close to the truth* in Kullback-Leibler divergence,

\$\$H_0 : \quad E\left\[ \log \frac{f_1(Y_i \mid x_i; \theta_1^{\*})}
{f_2(Y_i \mid x_i; \theta_2^{\*})} \right\] = 0,\$\$

the expectation taken under \\h\\ and \\\theta_k^{\*}\\ being the
pseudo-true parameter of model \\k\\. The alternatives are one-sided in
either direction: \\H_1\\, that model 1 is closer to the truth, and
\\H_2\\, that model 2 is. Crucially, **neither model is assumed
correct**. The test is about relative distance to the truth, so
rejecting \\H_0\\ says that one model is nearer, never that it is right.

## The statistic

Write \\\ell\_{ki}(\hat\theta_k)\\ for the log-likelihood contribution
of observation \\i\\ under model \\k\\, stored in `object$loglik_i`, and
let

\$\$m_i = \ell\_{1i}(\hat\theta_1) - \ell\_{2i}(\hat\theta_2), \qquad i
= 1, \ldots, n.\$\$

With `correction = TRUE` these differences are shifted by the
per-observation share of the dimension penalty,

\$\$\tilde m_i = m_i - \frac{(p_1 - p_2)\log n}{2n},\$\$

where \\p_k\\ is the number of coefficients estimated by model \\k\\;
with `correction = FALSE`, \\\tilde m_i = m_i\\. The reported statistic
is the studentized mean

\$\$Z = \frac{\sqrt{n}\\ \bar{\tilde m}}{s}, \qquad s^2 = \frac{1}{n -
1} \sum\_{i=1}^{n} (\tilde m_i - \bar{\tilde m})^2,\$\$

with \\s\\ the ordinary sample standard deviation, as computed by
[`stats::sd()`](https://rdrr.io/r/stats/sd.html). Under \\H_0\\, and for
strictly non-nested models, \\Z\\ converges in distribution to a
standard normal, so the two-sided p-value is \\2\\\Phi(-\|Z\|)\\.

The denominator is what distinguishes this from a naive comparison of
log-likelihoods. A difference in total log-likelihood says nothing on
its own about sampling variability; \\s\\ measures how consistently one
model beats the other across observations, so a small total advantage
accumulated steadily over every observation can be decisive, while a
large advantage driven by a handful of points is not.

## Reading the sign

The sign of \\Z\\ points at the favoured model, and the p-value says
whether the data support the preference at all.

- \\Z \> 0\\ with a small p-value: model 1, that is `object`, is
  favoured.

- \\Z \< 0\\ with a small p-value: model 2, that is `object2`, is
  favoured.

- a large p-value: the data do not discriminate between the two. This is
  a genuine third conclusion and not a failure to establish a
  preference: it states that both models sit at the same
  Kullback-Leibler distance from the truth as far as these data can
  tell. For an anchor comparison it is the expected outcome once the
  nuisance parameters are saturated, and the example below shows exactly
  that.

The [`print()`](https://rdrr.io/r/base/print.html) method states the
verdict at the conventional five per cent level.

## The dimension correction

The uncorrected statistic compares maximized log-likelihoods, and
therefore rewards the more flexible model for its flexibility alone. The
correction subtracts half the difference of BIC penalties, spread evenly
over the observations. Summing the corrected differences shows what it
buys:

\$\$\sum\_{i=1}^{n} \tilde m_i \\=\\ \ell_1 - \ell_2 - \frac{(p_1 -
p_2)\log n}{2} \\=\\ \frac{\mathrm{BIC}\_2 - \mathrm{BIC}\_1}{2},\$\$

since \\\mathrm{BIC}\_k = -2\ell_k + p_k \log n\\. The corrected
numerator is thus one half of the BIC difference, oriented so that a
positive value favours model 1, and the test is a formal significance
statement about a BIC comparison rather than an informal reading of one.
When \\p_1 = p_2\\ the correction term vanishes identically and the two
settings of `correction` agree to the last digit.

## When this test does not apply

The standard normal reference requires the two models to be *strictly*
non-nested. For a genuinely nested pair, one family restricting another,
the variance \\s^2\\ tends to zero under the null and the limit of \\Z\\
is not normal; the correct reference there is chi-squared, which is what
[`anova.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/anova.gkwqreg.md)
uses. Use [`anova()`](https://rdrr.io/r/stats/anova.html) for a nested
pair and this function for everything else.

The two fits must also share the quantile level and the number of
observations, both of which are checked. Fits at different levels answer
different questions and their likelihoods are not commensurable; fits on
different samples give differences \\m_i\\ that are not paired.

## References

Vuong, Q. H. (1989). Likelihood ratio tests for model selection and
non-nested hypotheses. *Econometrica* **57**, 307-333.

## See also

[`anova.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/anova.gkwqreg.md)
for the nested case;
[`compare_families()`](https://evandeilton.github.io/gkwqreg/reference/compare_families.md)
for a sweep over all seven families by AIC, BIC and check loss;
[`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md)
for why two anchors are two models.

## Examples

``` r
## Data generated under the "beta" anchor of the Kumaraswamy family: alpha is
## constant at 2 and beta is the quantity that varies with the covariate.
set.seed(11)
n  <- 400
x  <- runif(n, -2, 2)
mu <- plogis(0.4 + 1.1 * x)                  # true conditional 0.9-quantile
y  <- gkwdist::rkw(n, alpha = 2, beta = log1p(-0.9) / log1p(-mu^2))
d  <- data.frame(y = y, x = x)

## Two anchors, same family, same level, same data: three coefficients each,
## and neither nested in the other. anova() refuses them; this does not.
f_beta  <- gkwqreg(y ~ x, data = d, tau = 0.9, family = "kw")
f_alpha <- gkwqreg(y ~ x, data = d, tau = 0.9, family = "kw",
                   anchor = "alpha")
vuong_test(f_beta, f_alpha)
#> 
#> Vuong test for non-nested quantile regression fits
#> tau = 0.9, BIC-corrected
#>   model 1: kw / anchor beta
#>   model 2: kw / anchor alpha
#> 
#>   z = 6.316, p-value = 2.679e-10
#>   model 1 is favoured
##   z = 6.316, p-value = 2.679e-10
##   model 1 is favoured
## z is large and positive, so model 1 -- the beta anchor, the one the data
## were generated under -- is favoured decisively.

## Both fits estimate three coefficients, so the dimension correction is
## identically zero and switching it off changes nothing.
identical(vuong_test(f_beta, f_alpha)$statistic,
          vuong_test(f_beta, f_alpha, correction = FALSE)$statistic)
#> [1] TRUE
## [1] TRUE

## Saturating the nuisance parameter -- regressing alpha, respectively beta,
## on the same covariate as the quantile -- makes the anchor choice almost
## immaterial, which is the recommended practice when theory does not
## dictate an anchor.
s_beta  <- gkwqreg(y ~ x | x, data = d, tau = 0.9, family = "kw")
s_alpha <- gkwqreg(y ~ x | x, data = d, tau = 0.9, family = "kw",
                   anchor = "alpha")
vuong_test(s_beta, s_alpha)
#> 
#> Vuong test for non-nested quantile regression fits
#> tau = 0.9, BIC-corrected
#>   model 1: kw / anchor beta
#>   model 2: kw / anchor alpha
#> 
#>   z = 1.62, p-value = 0.1052
#>   neither model is favoured
##   z = 1.62, p-value = 0.1052
##   neither model is favoured
## The evidence has evaporated. Saturation has turned two different models
## into two parametrizations of very nearly the same one.

## Non-nested families work the same way. "kw" restricts alpha and beta while
## "beta" restricts them away entirely, so neither contains the other.
f_bet <- gkwqreg(y ~ x, data = d, tau = 0.9, family = "beta")
vuong_test(f_beta, f_bet)
#> 
#> Vuong test for non-nested quantile regression fits
#> tau = 0.9, BIC-corrected
#>   model 1: kw / anchor beta
#>   model 2: beta / anchor gamma
#> 
#>   z = 6.898, p-value = 5.277e-12
#>   model 1 is favoured
```
