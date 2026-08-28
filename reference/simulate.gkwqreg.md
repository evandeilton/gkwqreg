# Simulate responses from a fitted quantile regression

Draws new responses from each observation's fitted conditional
distribution. Because the anchored parameter is reconstructed from the
fitted conditional quantile before anything is drawn, every simulated
sample carries the fitted conditional `tau`-quantile by construction,
exactly and not approximately. That property is what makes this function
the engine behind the parametric bootstrap of
[`gkwq_boot()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_boot.md)
and the simulated envelopes of
[`plot.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/plot.gkwqreg.md).

## Usage

``` r
# S3 method for class 'gkwqreg'
simulate(object, nsim = 1L, seed = NULL, ...)
```

## Arguments

- object:

  A `"gkwqreg"` fit.

- nsim:

  Number of replicate samples to draw, returned as columns. Defaults to
  1.

- seed:

  Optional seed for reproducibility, handled as in
  [`stats::simulate()`](https://rdrr.io/r/stats/simulate.html): the
  current value of `.Random.seed` is saved, `set.seed(seed)` is called,
  and the original state is restored when the function exits. A seeded
  call is therefore reproducible without disturbing the surrounding
  random number stream.

- ...:

  Unused; present for consistency with the generic.

## Value

A data frame of `nobs(object)` rows and `nsim` columns, named `sim_1`
through `sim_nsim`, the rows in the order of the model frame used for
fitting. Each column is one independent replicate sample. Row names are
the default integer sequence, not those of the original data.

## Details

For observation \\i\\ the fit stores the reconstructed parameter vector
\\(\hat\alpha_i, \hat\beta_i, \hat\gamma_i, \hat\delta_i,
\hat\lambda_i)\\ in `object$parameter_vectors`. It is obtained by
evaluating each modelled parameter at its own linear predictor, then
solving the anchor equation for the eliminated one so that

\$\$Q(\tau; \hat\alpha_i, \hat\beta_i, \hat\gamma_i, \hat\delta_i,
\hat\lambda_i) \\=\\ \hat\mu_i,\$\$

\\\hat\mu_i\\ being the fitted conditional `tau`-quantile returned by
`fitted(object)`. Responses are then drawn independently across
observations from the corresponding Generalized Kumaraswamy laws,

\$\$y_i^{\*} \\\sim\\ \mathrm{GKw}(\hat\alpha_i, \hat\beta_i,
\hat\gamma_i, \hat\delta_i, \hat\lambda_i),\$\$

using
[`gkwdist::rgkw()`](https://evandeilton.github.io/gkwdist/reference/rgkw.html).
Since \\\hat\mu_i\\ *is* the `tau`-quantile of that law, it follows
immediately that

\$\$\Pr\left(y_i^{\*} \le \hat\mu_i\right) = \tau, \qquad i = 1, \ldots,
n,\$\$

whatever the family, the anchor, the covariates or the level. The
example below verifies this numerically, both marginally over the whole
sample and observation by observation. It is the simulation counterpart
of the fact that a single `gkwqreg` fit cannot produce crossing
quantiles: the fitted quantiles are those of one proper distribution
function, and this function samples from that distribution.

## What the draws do and do not capture

The parameters are held at their estimates and treated as known. The
draws therefore reproduce the sampling variation of the *response*
around a fixed fitted model; they say nothing about the estimation
uncertainty in \\\hat\theta\\ itself. Propagating that uncertainty into
the coefficients is what
[`gkwq_boot()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_boot.md)
does, and it is why a parametric bootstrap refits the model on each
simulated sample instead of stopping here.

A corollary worth keeping in view: a sample drawn here is, by
construction, precisely what the fitted family says the data should look
like. Comparing such a sample with the observed data can reveal a poor
*fit* within the family, which is the point of a simulated envelope, but
it cannot arbitrate between families, since every family would pass its
own test.

## See also

[`gkwq_boot()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_boot.md)
for the bootstrap built on this;
[`plot.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/plot.gkwqreg.md)
for simulated envelopes;
[`residuals.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/residuals.gkwqreg.md)
for the residuals they envelope;
[`stats::simulate()`](https://rdrr.io/r/stats/simulate.html) for the
generic and its seed convention.

## Examples

``` r
set.seed(11)
n  <- 400
x  <- runif(n, -2, 2)
mu <- plogis(0.4 + 1.1 * x)
y  <- gkwdist::rkw(n, alpha = 2, beta = log1p(-0.9) / log1p(-mu^2))
d  <- data.frame(y = y, x = x)

fit <- gkwqreg(y ~ x, data = d, tau = 0.9, family = "kw")

## One replicate, of the same length as the data.
s1 <- simulate(fit, nsim = 1, seed = 1)
dim(s1)
#> [1] 400   1
## [1] 400   1

## The defining property: the fitted values are the 0.9-quantiles of the
## distributions being sampled, so about 90 per cent of every simulated
## sample falls below them.
s <- as.matrix(simulate(fit, nsim = 200, seed = 1))
mean(s < fitted(fit))
#> [1] 0.8989625
## [1] 0.8989625

## Observation by observation, not merely on average. Each row of the matrix
## is one observation across the 200 replicates.
summary(rowMeans(s < fitted(fit)))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   0.835   0.885   0.900   0.899   0.915   0.950 
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
##   0.835   0.885   0.900   0.899   0.915   0.950

## Seeding is reproducible and leaves the ambient stream untouched.
identical(simulate(fit, nsim = 2, seed = 1),
          simulate(fit, nsim = 2, seed = 1))
#> [1] TRUE
## [1] TRUE
```
