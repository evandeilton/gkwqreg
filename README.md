# gkwqreg

<!-- badges: start -->
[![R-CMD-check](https://github.com/evandeilton/gkwqreg/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/evandeilton/gkwqreg/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R >= 4.1](https://img.shields.io/badge/R-%3E%3D%204.1-blue.svg)](https://cran.r-project.org/)
<!-- badges: end -->

**Parametric quantile regression for the Generalized Kumaraswamy family**

For a response in the open interval $(0,1)$ and a quantile level fixed in advance,
`gkwqreg` reparametrizes the Generalized Kumaraswamy distribution so that one of its
parameters *is* the conditional quantile at that level, and models that parameter
directly. Every coefficient then has a quantile interpretation and a standard error that
refers to one.

This is the sibling of [`gkwreg`](https://github.com/evandeilton/gkwreg), which fits the
same seven families with mean semantics.

## Why this is not `predict(type = "quantile")`

`gkwreg` and `betareg` can already report quantiles: fit a mean-oriented likelihood, then
apply the quantile function to the estimated parameters. That answers *"given the fitted
distribution, what is its 90th percentile?"*

It does not answer *"how do covariates shift the 90th percentile?"*, because no parameter
in that model **is** the 90th percentile. `gkwqreg` closes exactly that gap.

## Installation

```r
# install.packages("remotes")
remotes::install_github("evandeilton/gkwqreg")
```

The TMB template is compiled at install time, so no C++ toolchain is needed at run time.

## Five minutes to your first model

Six steps. Copy each block in turn; everything is self-contained.

### 1. Simulate a bounded response

The response must lie strictly inside `(0,1)` — a share, a rate, an index. Here it is
"share of households with piped-water access", and we build it so the true **0.9-quantile**
is `plogis(0.20 + 1.60*income - 0.50*urban)`.

```r
library(gkwqreg)

set.seed(2026)
n      <- 600
income <- runif(n, 0, 1)
region <- factor(rbinom(n, 1, 0.5), labels = c("rural", "urban"))
mu90   <- plogis(0.20 + 1.60 * income - 0.50 * (region == "urban"))
access <- gkwdist::rkw(n, alpha = 2.5, beta = log1p(-0.90) / log1p(-mu90^2.5))

dat <- data.frame(access, income, region)
```

### 2. Fit

One line. `tau = 0.9` says *which* quantile you are asking about; it is chosen by you and
never estimated.

```r
fit <- gkwqreg(access ~ income + region, data = dat, tau = 0.9, family = "kw")
```

### 3. Read the summary

```r
summary(fit)
```

```
Generalized Kumaraswamy quantile regression
family: kw   tau: 0.9   anchor: beta

Conditional 0.9-quantile (link logit) -- coefficients are effects on the LOG QUANTILE ODDS log(mu/(1-mu)):
               Estimate Std. Error z value Pr(>|z|)
mu:(Intercept)  0.14521    0.07983   1.819   0.0689 .
mu:income       1.73121    0.14868  11.644  < 2e-16 ***
mu:regionurban -0.49843    0.07893  -6.315 2.71e-10 ***

alpha (link log):
                  Estimate Std. Error z value Pr(>|z|)
alpha:(Intercept)  0.90635    0.03894   23.27   <2e-16 ***

Anchored parameter: beta (computed from the quantile, not estimated)
Log-likelihood: 230.6 on 4 Df   AIC: -453.2   BIC: -435.7
Pinball loss: 0.02995   Pseudo-R1: 0.1603
Empirical coverage: 0.895 (target 0.9)
Covariance: expected   Information condition number: 30.99
```

The true values were `0.20`, `1.60`, `-0.50`; every estimate is within about one standard
error. Line by line:

| Line | What it tells you |
|---|---|
| `mu:` block | The quantile model. **These are effects on the log quantile odds**, not on a mean and not on the odds of an event |
| `alpha:` block | A nuisance parameter. It controls the spread of the conditional distribution, not the quantile you are modelling |
| `Anchored parameter` | Which parameter was eliminated to make room for the quantile. A **modeling choice** — see below |
| `Empirical coverage` | Share of observations below the fitted quantile. Should sit at `tau`. **0.895 vs 0.900 — the model is doing its job** |
| `Information condition number` | Conditioning. Small (here 31) is healthy; above `1e8` means the standard errors are unreliable |

### 4. Turn a coefficient into a sentence

`mu:income = 1.731` is on the log-quantile-odds scale, which nobody can interpret directly.
Convert it to the response scale:

```r
marginal_effects(fit)
```

```
Marginal effects on the conditional 0.9-quantile (averaged over the sample)
dQ(tau|x)/dx_j, logit link

    variable  effect std.error   lower    upper
      income  0.3640   0.02764  0.3099  0.41819
 regionurban -0.1048   0.01609 -0.1363 -0.07327
```

Now you can write it down:

> Moving income from 0 to 1 raises the **90th percentile** of access by about **0.36**
> (95% CI 0.31 to 0.42). Urban regions sit about **0.10 lower** at the 90th percentile
> than rural ones (95% CI 0.14 to 0.07 lower).

These describe *conditional quantiles under the fitted model*. They are not causal effects.

### 5. Check you should believe it

Three numbers, all in the summary above, plus one plot:

```r
mean(dat$access <= fitted(fit))        # 0.8950  -- should be ~ tau
fit$cond_number                        # 31.0    -- small is good
mean(residuals(fit, type = "tau-sign"))# -0.0050 -- should be ~ 0

plot(fit, which = c(3, 4, 5))          # Q-Q, observed-vs-fitted, calibration
```

Panel 5 is the one to look at. It bins the residuals `1{y <= Q} - tau` by fitted value;
under a correct model every bin sits at zero. A bin whose interval excludes zero says the
fit is miscalibrated *there*, which no mean-regression diagnostic can tell you.

### 6. Use it

```r
nd <- data.frame(income = c(0.25, 0.75),
                 region = factor("rural", levels = c("rural", "urban")))

predict(fit, newdata = nd)                                     # 0.6406  0.8090
1 - predict(fit, newdata = nd, type = "probability", at = 0.5)  # 0.3227  0.6008
predict(fit, newdata = nd, type = "quantile", tau = 0.5)        # 0.4165  0.5584
```

Read in order: the fitted **90th percentile** at each profile; the probability that access
**exceeds 0.5** there; and the **median** for the same profiles, read off the *same fitted
distribution* without refitting.

That third line is the thing check-function quantile regression cannot do. A `gkwqreg` fit
carries a whole conditional distribution per observation, so any quantile, any exceedance
probability and any moment is available from one fit — and the quantiles it reports can
never cross, because they come from one proper distribution function.

### Where to go next

```r
vignette("gkwqreg")            # the fuller tour
vignette("gkwqreg-anchor")     # what `anchor` means and why it is a modeling choice
vignette("gkwqreg-design")     # the vocabulary and the design decisions
```

## The seven families

They form a nesting **lattice**, not a single chain, so a likelihood-ratio test applies
along a chain of containments. `kw` sits inside all of `ekw`, `kkw`, `bkw` and `gkw`;
`ekw` inside `kkw` and `gkw`; `beta` inside `mc`, `bkw` and `gkw`. But `ekw` and `bkw`
contain neither each other, and `anova()` returns `NA` with a warning for any such pair
rather than reporting a test that does not apply.

| family | free parameters | constraints | default anchor |
|---|---|---|---|
| `kw`   | α, β          | γ=1, δ=0, λ=1 | β |
| `ekw`  | α, β, λ       | γ=1, δ=0      | β |
| `kkw`  | α, β, δ, λ    | γ=1           | β |
| `bkw`  | α, β, γ, δ    | λ=1           | β |
| `gkw`  | α, β, γ, δ, λ | —             | β |
| `mc`   | γ, δ, λ       | α=β=1         | λ |
| `beta` | γ, δ          | α=β=λ=1       | γ |

```r
compare_families(fit, families = c("kw", "ekw", "beta"))
#>   family anchor df   logLik       AIC       BIC    pinball converged
#> 1     kw   beta  4 230.6209 -453.2418 -435.6540 0.02995121      TRUE
#> 2    ekw   beta  5 230.8573 -451.7146 -429.7299 0.02995948      TRUE
#> 3   beta  gamma  4 210.5077 -413.0155 -395.4278 0.03196525      TRUE
```

Compare on `pinball` when choosing out of sample: the pinball loss is what a quantile estimate
actually targets.

## The anchor is a modeling choice

`anchor` names the parameter eliminated in favour of the conditional quantile. It looks
like a reparametrization, and when the nuisance parameters are regressed on at least the
covariates used for the quantile, it is one — the likelihood is identical either way.

It stops being one as soon as the quantile varies with covariates while a nuisance
parameter is held constant. Anchoring on β then asserts *"α constant, βᵢ = f(μᵢ, α)"*;
anchoring on α asserts the reverse. In the design study behind this package the two
differed by **131 in log-likelihood** and by **41% in the coefficient of interest**.

Two anchors give non-nested models of equal dimension. Compare them with `vuong_test()`,
never with `anova()` — which refuses:

```r
fa <- gkwqreg(y ~ x, data = d, tau = 0.9, family = "kw", anchor = "alpha")
anova(fit, fa)
#> Error: models with different anchors cannot be compared by a likelihood-ratio
#> test: they are non-nested models of equal dimension.
vuong_test(fit, fa)
```

The practical advice: regress the nuisance parameters on the same covariates as the
quantile (`y ~ x | x`) and the fit becomes largely anchor-insensitive.

## Across quantile levels

```r
fits <- gkwqreg(y ~ x, data = d, tau = seq(0.1, 0.9, by = 0.1), family = "kw")
plot(quantile_process(fits))     # coefficient paths, in the idiom of plot.rqs
check_crossing(fits)             # crossing across independent fits
rearrange(fits)                  # Chernozhukov et al. (2010) monotonization
```

A single fit cannot cross, and the printed output says why:

```
  rows with a crossing: 0 of 300 (0.00%)
  Zero crossings is guaranteed here, not luck: these are quantiles of
  a single proper distribution function.
```

That guarantee is what a parametric fit buys over independently fitted levels — and the
price is paid under misspecification, where check-function quantile regression stays
consistent and this does not.

## Prior art, stated plainly

The `kw` slice of this package **is** the Kumaraswamy quantile reparametrization of
**Mitnik and Baek (2013)**, which the [`unitquantreg`](https://cran.r-project.org/package=unitquantreg)
package implements as its `kum` family. The general formula here contains that result
exactly as its two-parameter slice, and the test suite checks the identity.

What is new is the extension to the four- and five-parameter members, and the fact that
they are **nested**: `unitquantreg`'s thirteen families are mutually non-nested and need
Vuong tests, while here `kw ⊂ ekw ⊂ gkw` are genuine nestings compared by ordinary LR
tests.

## Where this is not the right tool

- Under **misspecification** of the family, `quantreg::rq` stays consistent for the
  conditional quantile and this does not. That is the central trade-off of parametric
  quantile regression.
- For a **single quantile** from a wide catalogue of two-parameter families, with no need
  for nested family selection or regressions on the shape, `unitquantreg` is simpler and
  established — it is `tau`-indexed exactly as this package is. `cdfquantreg` covers
  location and median models, which are not `tau`-indexed at all.
- For **smooth or additive** quantile effects, `qgam` is the right tool.
- The five-parameter `gkw` is weakly identified in any parametrization; the package warns
  and `summary()` reports the information matrix's condition number.

## Vignettes

```r
vignette("gkwqreg")               # getting started
vignette("gkwqreg-anchor")        # the reparametrization and identifiability
vignette("gkwqreg-comparison")    # versus quantreg, unitquantreg, betareg
vignette("gkwqreg-case-study")    # upper-quantile body-fat modelling
vignette("gkwqreg-design")        # the vocabulary and the four design decisions
```

The last one is worth knowing about. It defines the terms this package commits to
(*anchor*, *saturated*, *part*, *quantile level*) and records the four decisions a reader
is most likely to question — including why the package ships its own incomplete-beta
derivatives instead of using TMB's `qbeta`, whose shape derivatives are silently **zero**
at whole-number shapes. It demonstrates that failure, and this package's fix, live.

## References

Every DOI below was resolved against the registration agency's own record.

**The distribution**

- Kumaraswamy, P. (1980). A generalized probability density function for
  double-bounded random processes. *Journal of Hydrology* **46**(1–2), 79–88.
  [doi:10.1016/0022-1694(80)90036-0](https://doi.org/10.1016/0022-1694(80)90036-0)
- Carrasco, J. M. F., Ferrari, S. L. P. and Cordeiro, G. M. (2010). A New
  Generalized Kumaraswamy Distribution. *arXiv preprint* arXiv:1004.0911.
  [doi:10.48550/arXiv.1004.0911](https://doi.org/10.48550/arXiv.1004.0911)
  (no peer-reviewed version exists)

**The quantile reparametrization**

- Mitnik, P. A. and Baek, S. (2013). The Kumaraswamy distribution:
  median-dispersion re-parameterizations for regression modeling and
  simulation-based estimation. *Statistical Papers* **54**(1), 177–192.
  [doi:10.1007/s00362-011-0417-y](https://doi.org/10.1007/s00362-011-0417-y)
- Bayes, C. L., Bazán, J. L. and de Castro, M. (2017). A quantile parametric
  mixed regression model for bounded response variables. *Statistics and Its
  Interface* **10**(3), 483–493.
  [doi:10.4310/SII.2017.v10.n3.a11](https://doi.org/10.4310/SII.2017.v10.n3.a11)

**Quantile regression**

- Koenker, R. and Bassett, G. (1978). Regression quantiles. *Econometrica*
  **46**(1), 33–50. [doi:10.2307/1913643](https://doi.org/10.2307/1913643)
- Koenker, R. and Machado, J. A. F. (1999). Goodness of fit and related
  inference processes for quantile regression. *JASA* **94**(448), 1296–1310.
  [doi:10.1080/01621459.1999.10473882](https://doi.org/10.1080/01621459.1999.10473882)
- Chernozhukov, V., Fernández-Val, I. and Galichon, A. (2009). Improving point
  and interval estimators of monotone functions by rearrangement. *Biometrika*
  **96**(3), 559–575. [doi:10.1093/biomet/asp030](https://doi.org/10.1093/biomet/asp030)
- Chernozhukov, V., Fernández-Val, I. and Galichon, A. (2010). Quantile and
  probability curves without crossing. *Econometrica* **78**(3), 1093–1125.
  [doi:10.3982/ECTA7880](https://doi.org/10.3982/ECTA7880)

**Related software and the comparison class**

- Menezes, A. F. B. and Mazucheli, J. `unitquantreg`: Parametric Quantile
  Regression Models for Bounded Data. R package.
  [doi:10.32614/CRAN.package.unitquantreg](https://doi.org/10.32614/CRAN.package.unitquantreg)
- Mazucheli, J., Alves, B., Menezes, A. F. B. and Leiva, V. (2022). An overview
  on parametric quantile regression models and their computational
  implementation with applications to biomedical problems including COVID-19
  data. *Computer Methods and Programs in Biomedicine* **221**, 106816.
  [doi:10.1016/j.cmpb.2022.106816](https://doi.org/10.1016/j.cmpb.2022.106816)
- Ferrari, S. L. P. and Cribari-Neto, F. (2004). Beta regression for modelling
  rates and proportions. *Journal of Applied Statistics* **31**(7), 799–815.
  [doi:10.1080/0266476042000214501](https://doi.org/10.1080/0266476042000214501)
- Lopes, J. E. and Bonat, W. H. (2026). `gkwreg`: An R package for generalized
  Kumaraswamy regression models for bounded data. *JOSS* **11**(117), 8991.
  [doi:10.21105/joss.08991](https://doi.org/10.21105/joss.08991)

**Inference and numerics**

- Dunn, P. K. and Smyth, G. K. (1996). Randomized quantile residuals. *JCGS*
  **5**(3), 236–244. [doi:10.1080/10618600.1996.10474708](https://doi.org/10.1080/10618600.1996.10474708)
- Cox, D. R. and Reid, N. (1987). Parameter orthogonality and approximate
  conditional inference. *JRSS-B* **49**(1), 1–18.
  [doi:10.1111/j.2517-6161.1987.tb01422.x](https://doi.org/10.1111/j.2517-6161.1987.tb01422.x)
- Vuong, Q. H. (1989). Likelihood ratio tests for model selection and
  non-nested hypotheses. *Econometrica* **57**(2), 307–333.
  [doi:10.2307/1912557](https://doi.org/10.2307/1912557)
- Boik, R. J. and Robison-Cox, J. F. (1998). Derivatives of the incomplete beta
  function. *Journal of Statistical Software* **3**(1), 1–20.
  [doi:10.18637/jss.v003.i01](https://doi.org/10.18637/jss.v003.i01)
- Zeileis, A. (2006). Object-oriented computation of sandwich estimators.
  *Journal of Statistical Software* **16**(9), 1–16.
  [doi:10.18637/jss.v016.i09](https://doi.org/10.18637/jss.v016.i09)
- Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H. and Bell, B. M. (2016).
  TMB: Automatic differentiation and Laplace approximation. *Journal of
  Statistical Software* **70**(5), 1–21.
  [doi:10.18637/jss.v070.i05](https://doi.org/10.18637/jss.v070.i05)

## License

MIT
