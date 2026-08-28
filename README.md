# gkwqreg

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

## Getting started

```r
library(gkwqreg)

set.seed(1)
n  <- 300
x  <- runif(n, -2, 2)
mu <- plogis(0.4 + 1.1 * x)
y  <- gkwdist::rkw(n, alpha = 2, beta = log1p(-0.5) / log1p(-mu^2))
d  <- data.frame(y = y, x = x)

fit <- gkwqreg(y ~ x, data = d, tau = 0.9, family = "kw")
summary(fit)
```

```
Generalized Kumaraswamy quantile regression
family: kw   tau: 0.9   anchor: beta

Quantile residuals:
    Min      1Q  Median      3Q     Max
-2.9359 -0.5983  0.0174  0.6741  2.8844

Conditional 0.9-quantile (link logit) -- coefficients are effects on the LOG QUANTILE ODDS log(mu/(1-mu)):
               Estimate Std. Error z value Pr(>|z|)
mu:(Intercept)   2.3382     0.1434   16.30   <2e-16 ***
mu:x             1.8560     0.1008   18.42   <2e-16 ***

alpha (link log):
                  Estimate Std. Error z value Pr(>|z|)
alpha:(Intercept)  0.68066    0.06448   10.56   <2e-16 ***

Anchored parameter: beta (computed from the quantile, not estimated)
Log-likelihood: 149.2 on 3 Df   AIC: -292.4   BIC: -281.3
Pinball loss: 0.0277   Pseudo-R1: 0.358
Empirical coverage: 0.8867 (target 0.9)
Covariance: expected   Information condition number: 14.44
```

Note the last two lines: the fit puts 88.7% of observations below its fitted quantile
against a target of 90%, and the information matrix is well conditioned.

### `fitted()` returns quantiles, not means

This is the one place a `gkwreg` user will get burned, so the class deliberately does not
inherit from `gkwreg` — a missing method errors rather than silently applying mean
semantics.

```r
head(fitted(fit))                  # conditional 0.9-quantiles
head(fitted(fit, type = "mean"))   # E[Y|x], the escape hatch
```

### Reading the coefficients

Under the default logit link, `exp(1.856) = 6.4` is the multiplicative effect of a
one-unit increase in `x` on the **odds of the conditional 0.9-quantile** — not on the odds
of an event, and not on a mean. For the effect on the quantile itself:

```r
marginal_effects(fit)
#>  variable effect std.error  lower  upper
#>         x   0.21    0.0114 0.1876 0.2323
```

## The seven families

They form a genuine nesting, so family selection is an ordinary likelihood-ratio test.

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
#>   family anchor df    logLik       AIC    pinball converged
#> 1     kw   beta  3 149.21283 -292.4257 0.02769714      TRUE
#> 2    ekw   beta  4 149.36135 -290.7227 0.02769820      TRUE
#> 3   beta  gamma  3  65.15314 -124.3063 0.03894886      TRUE
```

Compare on `pinball` when choosing out of sample: check loss is what a quantile estimate
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
- For a **median-only** model with no interest in shape, `unitquantreg` and `cdfquantreg`
  are simpler and established.
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

- Carrasco, Ferrari and Cordeiro (2010). A generalized Kumaraswamy distribution.
- Mitnik and Baek (2013). The Kumaraswamy distribution: median-dispersion
  re-parameterizations. *Statistical Papers* **54**, 177–192.
- Mazucheli, Alves, Menezes and Leiva (2022). `unitquantreg`. *CMPB* **221**, 106816.
- Dunn and Smyth (1996). Randomized quantile residuals.
- Chernozhukov, Fernández-Val and Galichon (2010). Quantile and probability curves
  without crossing.
- Boik and Robinson-Cox (1998). Derivatives of the incomplete beta function. *JSS* **3**(1).
- Lopes and Bonat (2026). `gkwreg`. *JOSS*. doi:10.21105/joss.08991

## License

MIT
