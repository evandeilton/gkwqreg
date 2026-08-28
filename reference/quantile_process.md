# The estimated quantile process

Fits the model independently at each level of a grid of quantile levels
and collects the resulting coefficient paths \\\hat{\beta}\_j(\tau)\\
together with pointwise confidence bands. The companion [plot
method](https://evandeilton.github.io/gkwqreg/reference/plot.gkwq_process.md)
draws them in the visual idiom of `quantreg`'s `plot.summary.rqs`,
because that picture is the literature's standard summary of a quantile
process and readers should recognise it on sight.

## Usage

``` r
quantile_process(object, taus = seq(0.05, 0.95, by = 0.05), level = 0.95, ...)
```

## Arguments

- object:

  Either a `"gkwqreg"` fit, whose stored call is re-evaluated with `tau`
  replaced by `taus`, or a `"gkwqregs"` container of fits that have
  already been computed and are simply harvested.

- taus:

  Increasing grid of quantile levels in \\(0,1)\\, used only when
  `object` is a single fit. **Silently ignored when `object` is a
  `"gkwqregs"` container**, whose levels are already fixed by its fits.

- level:

  Confidence level for the pointwise bands, `0.95` by default.

- ...:

  Additional arguments merged into the re-evaluated call, hence passed
  to
  [`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md);
  `family`, `anchor` or `control` may be overridden here. Ignored when
  `object` is a container.

## Value

An object of class `"gkwq_process"`, a list with components

- `taus`:

  Numeric vector of the levels, increasing.

- `coef`:

  Coefficient matrix, one row per model coefficient and one column per
  level. Rows are named `"part:term"` (for example `"mu:x"`,
  `"alpha:(Intercept)"`); columns are named by level.

- `se`:

  Standard errors, same shape and dimnames as `coef`.

- `lower`, `upper`:

  Pointwise Wald band limits, `coef` plus or minus
  `qnorm(1 - (1 - level) / 2)` standard errors.

- `parts`:

  Character vector aligned with the rows of `coef`, naming the formula
  part each coefficient belongs to; this is what the plot method's
  `parts` argument selects on.

- `level`:

  The confidence level used for the bands.

- `family`, `anchor`:

  Family and anchor, common to every fit.

- `fits`:

  The list of underlying `"gkwqreg"` fits, so that any per-level
  quantity (log-likelihood, check loss, residuals) remains available.

- `call`:

  The matched call.

A `print` method reports the family, anchor, level range and the rounded
coefficient matrix.

## Details

Each level is a separate maximum-likelihood problem with its own
likelihood and its own parameter vector, exactly as
[`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md)
produces when given a vector `tau`. The process object collects those
separate answers side by side; it does not estimate them jointly, and it
imposes no smoothness or monotonicity across \\\tau\\.

The bands are pointwise Wald intervals computed level by level,
\\\hat{\beta}\_j(\tau) \pm z\_{1 - (1 - \texttt{level})/2}\\
\mathrm{se}\\\hat{\beta}\_j(\tau)\\\\, from each fit's own observed
information. Three consequences follow and none of them should be
glossed over when the picture is presented. The bands are pointwise in
\\\tau\\ and not simultaneous, so a band that excludes zero at one level
out of the nineteen in the default grid is not evidence at the nominal
level for the process as a whole. They ignore the dependence between
levels, which is substantial because every fit uses the same rows. And
they are symmetric on the coefficient scale, which is the scale of the
linear predictor, not of the quantile; use
[`marginal_effects()`](https://evandeilton.github.io/gkwqreg/reference/marginal_effects.md)
for statements about the response.

When `object` is a single fit the call is rebuilt and evaluated in the
caller's frame, so the data referred to by that call must still be
visible from where `quantile_process()` is invoked. Passing a
`"gkwqregs"` container avoids the issue entirely and refits nothing.

## A non-flat coefficient path is not evidence of tail heterogeneity

This is the single most consequential misreading of the picture this
function produces, and under the present parametrization it is a trap
rather than a subtlety.

In check-function quantile regression the coefficient path is flat
exactly when the covariate shifts the whole conditional distribution
without changing its shape, so a sloping path is read, correctly, as
tail heterogeneity. **That reading does not carry over here.** The
quantity plotted is the slope of \\g(\mu\_\tau)\\, the link-transformed
conditional \\\tau\\-quantile, and for a fixed conditional distribution
the map from a covariate to \\g\\ of its \\\tau\\-quantile is a
different function at each \\\tau\\. The path is therefore *generically*
\\\tau\\-varying even when the data-generating process has no tail
heterogeneity whatsoever. It is a property of the link and the
reparametrization, not a finding about the data.

The effect is large, not a second-order curiosity. Simulating from a
correctly specified `kw` model with a **constant** \\\alpha\\, one
covariate, a single index driving the entire conditional distribution
and a true median slope of `1.1` (\\n = 20000\\), the fitted `mu:x`
coefficient runs

|        |       |       |       |       |       |
|--------|-------|-------|-------|-------|-------|
| `tau`  | 0.10  | 0.25  | 0.50  | 0.75  | 0.90  |
| `mu:x` | 0.748 | 0.883 | 1.100 | 1.397 | 1.744 |

The median is recovered exactly. The other four levels are not the same
number and were never supposed to be: a path rising by a factor of 2.3
from the first decile to the ninth is what *no* tail heterogeneity looks
like on this scale. A reader shown only the picture would report that
the covariate "matters more in the upper tail", and would be describing
the parametrization.

What to do instead. The null hypothesis to compare against is not a
horizontal line but the path implied by a homogeneous model, and there
are two honest ways to obtain it. Simulate from the fitted homogeneous
model with
[`simulate.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/simulate.gkwqreg.md),
refit the process on the simulated data, and use that path as the
reference; the `reference` argument of
[`plot.gkwq_process()`](https://evandeilton.github.io/gkwqreg/reference/plot.gkwq_process.md)
takes a value per coefficient for exactly this purpose. Or, better, test
the question directly inside the model: give the shape parameter its own
regression (`y ~ x | x` rather than `y ~ x`) and compare the two nested
fits at a single level with
[`lrtest()`](https://evandeilton.github.io/gkwqreg/reference/lrtest.md).
Tail heterogeneity in this parametrization means the shape parameters
depend on the covariates, and that is a hypothesis with a
likelihood-ratio test attached, not something to be read off a slope.

The shape coefficients' own paths are the useful diagnostic in the plot:
under a homogeneous data-generating process they are flat, because a
constant shape parameter is estimated as the same constant at every
level.

## See also

[`plot.gkwq_process()`](https://evandeilton.github.io/gkwqreg/reference/plot.gkwq_process.md)
for the picture;
[`check_crossing()`](https://evandeilton.github.io/gkwqreg/reference/check_crossing.md)
and
[`rearrange()`](https://evandeilton.github.io/gkwqreg/reference/rearrange.md)
for the monotonicity of independently fitted levels;
[`marginal_effects()`](https://evandeilton.github.io/gkwqreg/reference/marginal_effects.md)
for effects on the scale of the response;
[`lrtest()`](https://evandeilton.github.io/gkwqreg/reference/lrtest.md)
for testing shape heterogeneity properly.

## Examples

``` r
## A correctly specified kw model with a CONSTANT shape parameter: by
## construction there is no tail heterogeneity of any kind here.
set.seed(6)
n <- 1000
x <- runif(n, -2, 2)
mu <- plogis(0.4 + 1.1 * x)                    # true median, logit slope 1.1
y <- gkwdist::rkw(n, alpha = 2, beta = log1p(-0.5) / log1p(-mu^2))
d <- data.frame(y = y, x = x)

fit <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
qp <- quantile_process(fit, taus = c(0.1, 0.25, 0.5, 0.75, 0.9))
round(qp$coef, 3)
#>                     0.10   0.25  0.50  0.75  0.90
#> mu:(Intercept)    -1.180 -0.409 0.436 1.342 2.272
#> mu:x               0.750  0.891 1.116 1.425 1.789
#> alpha:(Intercept)  0.664  0.680 0.688 0.685 0.675
#>                     0.10   0.25  0.50  0.75  0.90
#> mu:(Intercept)    -1.180 -0.409 0.436 1.342 2.272
#> mu:x               0.750  0.891 1.116 1.425 1.789
#> alpha:(Intercept)  0.664  0.680 0.688 0.685 0.675

## Read the two informative rows against each other.
##   mu:x rises from 0.75 to 1.79 and recovers 1.1 at the median. That rise
##   is NOT tail heterogeneity: the simulation has none. It is what a single
##   fixed conditional distribution looks like when read at five levels
##   through a logit link.
##   alpha:(Intercept) is flat to within 0.025, correctly reporting that the
##   shape parameter is the same constant at every level. THAT is the row
##   that carries the heterogeneity question.

## Test the question directly instead of eyeballing the slope of mu:x.
m0 <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")       # alpha constant
m1 <- gkwqreg(y ~ x | x, data = d, tau = 0.5, family = "kw")   # alpha ~ x
lrtest(m0, m1)
#> Likelihood-ratio test for Generalized Kumaraswamy quantile regression
#> tau = 0.5, anchor = beta
#>      Df logLik     AIC  Chisq Chi Df Pr(>Chisq)
#> [1,]  3 521.43 -1036.9                         
#> [2,]  4 521.43 -1034.9 0.0018      1     0.9664
#>      Df logLik     AIC  Chisq Chi Df Pr(>Chisq)
#> [1,]  3 521.43 -1036.9
#> [2,]  4 521.43 -1034.9 0.0018      1     0.9664
## No evidence at all that the shape depends on x -- which is the right
## answer, and the one the sloping mu:x path would have contradicted.

## Harvesting an existing container refits nothing and ignores `taus`.
fits <- gkwqreg(y ~ x, data = d, tau = c(0.25, 0.5, 0.75), family = "kw")
quantile_process(fits, taus = seq(0.1, 0.9, by = 0.1))$taus
#> [1] 0.25 0.50 0.75
#> 0.25 0.50 0.75
```
