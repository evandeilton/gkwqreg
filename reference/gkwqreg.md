# Parametric quantile regression for the Generalized Kumaraswamy family

Fits a regression model for the conditional \\\tau\\-quantile of a
response confined to the open unit interval \\(0,1)\\, at a level
\\\tau\\ that the analyst fixes in advance. The Generalized Kumaraswamy
(GKw) distribution is reparametrized so that one of its five parameters
*is* the conditional \\\tau\\-quantile \\\mu\_\tau\\; that parameter is
then modelled directly through a link function, while the parameters
left free may carry regressions of their own. Estimation is by maximum
likelihood with exact gradients from automatic differentiation.

Seven nested families are available (`kw`, `ekw`, `kkw`, `bkw`, `gkw`,
`mc`, `beta`), so the shape of the conditional distribution is chosen by
model selection rather than assumed. Unlike check-function quantile
regression, a single fit carries a complete conditional distribution:
every other quantile, the density, the mean and the variance are all
available from it.

## Usage

``` r
gkwqreg(
  formula,
  data,
  tau = 0.5,
  family = c("kw", "ekw", "kkw", "bkw", "gkw", "mc", "beta"),
  anchor = NULL,
  link = NULL,
  link_scale = NULL,
  subset = NULL,
  weights = NULL,
  offset = NULL,
  na.action = stats::na.omit,
  contrasts = NULL,
  control = gkwq_control(),
  model = TRUE,
  x = FALSE,
  y = TRUE,
  ...
)
```

## Arguments

- formula:

  A multi-part formula, parts separated by `|`. **Part one is always the
  conditional quantile.** The remaining parts follow the family's
  parameter order with the anchored parameter removed;
  [`gkwq_parts()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_parts.md)
  reports the exact contract. Omitted trailing parts default to `~ 1`.

- data:

  An optional data frame, list or environment holding the variables in
  `formula`. If missing, they are looked up in the environment of
  `formula`.

- tau:

  Quantile level in `(0,1)`. A vector returns a `"gkwqregs"` container
  of independent fits, one per level, fitted from the level nearest the
  median outward with warm starts. Mandatory in the sense that it is
  never estimated: the profile likelihood in `tau` is exactly flat, so
  `tau` indexes the question, not the model.

- family:

  One of `"kw"`, `"ekw"`, `"kkw"`, `"bkw"`, `"gkw"`, `"mc"`, `"beta"`.
  See the table above for the constraints each imposes.

- anchor:

  The parameter eliminated in favour of the conditional quantile, or
  `NULL` for the family default. See the section above; this is a
  modeling argument, not an internal detail.
  [`gkwq_anchors()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_anchors.md)
  lists what a family allows.

- link:

  Link functions, as a single name or a named vector indexed by part.
  The quantile part defaults to `"logit"` and must map to `(0,1)`, so it
  is one of `"logit"`, `"probit"`, `"cauchy"`, `"cloglog"`. The
  remaining parts default to `"log"` and may also be `"identity"`,
  `"sqrt"`, `"inverse"` or `"inverse-square"`. A single *unnamed* link
  applies to the nuisance parts only, never silently to the quantile; an
  unnamed vector as long as
  [`gkwq_parts()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_parts.md)
  applies part by part.

- link_scale:

  Multiplicative scale factors for the inverse links, as a single value
  or a named vector indexed by part. Only the bounded links use them,
  rescaling the image from `(0,1)` to `(0, scale)`; the conditional
  quantile is a probability and is never rescaled, so its scale is
  forced to 1.

- subset:

  An optional expression selecting a subset of rows, evaluated in `data`
  first and then in the calling frame.

- weights:

  Optional prior weights \\w_i\\, entering the likelihood as multipliers
  of the per-observation log-density; the reported check loss is
  weighted to match. May name a column of `data`. Negative weights are
  an error.

- offset:

  An optional numeric offset added to the linear predictor of the
  **conditional quantile**, which is the part an offset is normally
  meant for. Offsets for another part are written as an
  [`offset()`](https://rdrr.io/r/stats/offset.html) term inside that
  part of `formula`, and the two add.

- na.action:

  How to treat missing values in the model frame; defaults to
  [`stats::na.omit()`](https://rdrr.io/r/stats/na.fail.html). Responses
  at exactly 0 or 1 lie outside the support and are an error, never
  silently clamped.

- contrasts:

  An optional list of contrasts for factor covariates, as in
  [`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html).

- control:

  An object from
  [`gkwq_control()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_control.md)
  holding the optimizer, starting-value and covariance settings.

- model, x, y:

  Whether to keep the model frame, the per-part design matrices and the
  response in the fitted object. `y = TRUE` (the default) is required by
  [`residuals.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/residuals.gkwqreg.md),
  [`plot.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/plot.gkwqreg.md)
  and the goodness-of-fit lines of
  [`summary()`](https://rdrr.io/r/base/summary.html); `x = TRUE` avoids
  rebuilding design matrices in
  [`marginal_effects()`](https://evandeilton.github.io/gkwqreg/reference/marginal_effects.md).

- ...:

  Passed to
  [`gkwq_control()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_control.md).
  If any are supplied, the resulting control object **replaces**
  `control` entirely, so use one route or the other rather than both.

## Value

For a single `tau`, an object of class `"gkwqreg"`: a list whose
components include

- `coefficients`:

  the stacked coefficient vector, named `"part:term"`.

- `coef_list`:

  the same coefficients split by part; `coef(fit, part =)` reads from
  it.

- `se`, `vcov`, `hessian`, `cond_number`:

  standard errors, the inverse observed information, the observed
  information itself, and its exact condition number. `cond_number`
  above `1e8` is the signature of a weakly identified fit.

- `fitted.values`:

  the fitted conditional \\\tau\\-quantiles \\\mu\_{\tau i}\\ –
  **quantiles, not means**.

- `linear.predictors`:

  a named list of linear predictors, one per part, offsets included.

- `parameter_vectors`:

  an `n` by 5 data frame of the reconstructed `alpha`, `beta`, `gamma`,
  `delta`, `lambda` per observation, the anchored column included.

- `loglik`, `loglik_i`, `npar`, `nobs`, `aic`, `bic`:

  the maximized log-likelihood (re-evaluated at the reported
  coefficients, never taken from the optimizer's own record), its
  per-observation contributions, the number of estimated coefficients,
  the sample size, and the two information criteria.

- `pinball`:

  the in-sample check loss \\n^{-1}\sum_i (y_i-\mu\_{\tau
  i})\bigl(\tau-\mathbf{1}(y_i\<\mu\_{\tau i})\bigr)\\, the criterion a
  quantile estimate actually targets.

- `family`, `tau`, `anchor`, `parts`, `spec`:

  the family name, the quantile level, the anchored parameter, the part
  names in order, and the full internal family specification.

- `link`, `link_scale`:

  the link name and scale used for each part.

- `call`, `formula`, `terms`, `levels`, `contrasts`:

  the matched call, the multi-part formula, per-part terms objects, the
  factor levels and the contrasts recorded at fit time, all used by
  `predict(newdata =)`.

- `weights`, `offsets`:

  the prior weights and the per-part offsets.

- `convergence`, `message`, `iterations`:

  the optimizer's exit code, message and iteration count. A non-zero
  code triggers a warning at fit time and a note in
  [`summary()`](https://rdrr.io/r/base/summary.html).

- `control`, `start`, `obj`:

  the control object, the starting values actually used, and the TMB
  object. `obj` is what makes
  [`estfun.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/estfun.gkwqreg.md)
  and profile-likelihood intervals possible; it holds external pointers,
  so it is valid only in the session that built it.

- `model`, `x`, `y`:

  present only when the matching argument was `TRUE`: the model frame,
  the list of design matrices, and the (clamped) response.

For a vector `tau`, an object of class `"gkwqregs"` with components
`fits` (a named list of `"gkwqreg"` objects, one per level), `taus`,
`family`, `anchor`, `parts`, `call` and `nobs`. Each level is a separate
likelihood, so the container has no single
[`logLik()`](https://rdrr.io/r/stats/logLik.html); take it from an
individual fit.

## The distribution and its quantile function

For \\0\<y\<1\\, with \\\alpha,\beta,\gamma,\lambda\>0\\ and \\\delta\ge
0\\, the GKw density (Carrasco, Ferrari and Cordeiro, 2010) is

\$\$f(y;\alpha,\beta,\gamma,\delta,\lambda)=
\frac{\lambda\alpha\beta}{B(\gamma,\delta+1)}\\
y^{\alpha-1}\\(1-y^{\alpha})^{\beta-1}\\
\bigl\[1-(1-y^{\alpha})^{\beta}\bigr\]^{\gamma\lambda-1}\\
\bigl\\1-\bigl\[1-(1-y^{\alpha})^{\beta}\bigr\]^{\lambda}\bigr\\^{\delta}.\$\$

Writing \\v=1-(1-y^{\alpha})^{\beta}\\, the distribution function is the
regularized incomplete beta function
\\F(y)=I\_{v^{\lambda}}(\gamma,\delta+1)\\, which inverts in closed
form. With \\z\_\tau=I^{-1}\_\tau(\gamma,\delta+1)\\, that is
`qbeta(tau, gamma, delta + 1)`,

\$\$Q(\tau)=\Bigl\\1-\bigl\[1-z\_\tau^{1/\lambda}\bigr\]^{1/\beta}\Bigr\\^{1/\alpha}.\$\$

Two special cases keep the incomplete beta out of the computation
entirely: \\z\_\tau=\tau\\ when \\\gamma=1,\delta=0\\ (families `kw`,
`ekw`), and \\z\_\tau=1-(1-\tau)^{1/(\delta+1)}\\ when \\\gamma=1\\
(family `kkw`).

## The quantile reparametrization (the anchor)

Fix \\\tau\\ and set \\\mu=Q(\tau)\\. Rearranging the quantile function
gives the master identity

\$\$\mu^{\alpha}=1-\bigl(1-z\_\tau^{1/\lambda}\bigr)^{1/\beta},\$\$

from which any one of \\\alpha\\, \\\beta\\ or \\\lambda\\ can be solved
for in closed form:

\$\$\alpha=\frac{\log\bigl(1-(1-z\_\tau^{1/\lambda})^{1/\beta}\bigr)}{\log\mu},
\qquad
\beta=\frac{\log\bigl(1-z\_\tau^{1/\lambda}\bigr)}{\log(1-\mu^{\alpha})},
\qquad \lambda=\frac{\log
z\_\tau}{\log\bigl(1-(1-\mu^{\alpha})^{\beta}\bigr)}.\$\$

The parameter eliminated this way is the **anchor**. Each solve is a
ratio of two logarithms of numbers in \\(0,1)\\, hence strictly positive
and finite for every \\\mu,\tau\in(0,1)\\ and every admissible value of
the remaining parameters. No box constraint is ever needed, which is
what makes an unconstrained optimizer safe here.

Setting \\\gamma=1\\, \\\delta=0\\, \\\lambda=1\\ gives \\z\_\tau=\tau\\
and reduces the \\\beta\\-solve to
\\\beta=\log(1-\tau)/\log(1-\mu^{\alpha})\\, which is exactly the
Kumaraswamy quantile (median-dispersion) reparametrization of Mitnik and
Baek (2013). The general construction therefore contains the established
two-parameter case as a slice.

\\\gamma\\ and \\\delta\\ admit no closed-form solve: they enter only
through \\z\_\tau\\. Family `beta` is anchored on \\\gamma\\ by a
bracketed root of \\I\_\mu(\gamma,\delta+1)=\tau\\, which always exists
because the left-hand side sweeps all of \\(0,1)\\ as \\\gamma\\ ranges
over \\(0,\infty)\\. The \\\delta\\-solve is inadmissible in every
family: it would require \\\tau\>\mu^{\gamma}\\, which nothing
guarantees.

## The regression and its estimating equations

For observation \\i\\, with \\\tau\\ fixed,

\$\$g\_\mu(\mu\_{\tau i}) = \mathbf{x}\_i^{\top}\boldsymbol\beta\_\mu,
\qquad g_k(\theta\_{ki}) =
\mathbf{x}\_{ki}^{\top}\boldsymbol\beta_k,\$\$

where \\g\_\mu\\ maps \\(0,1)\to\mathbb{R}\\ (logit by default) and each
\\g_k\\ maps \\(0,\infty)\to\mathbb{R}\\ (log by default), one equation
per non-anchored parameter \\k\\. The anchored parameter is *computed*
from \\(\mu\_{\tau i},\tau,\theta_i)\\ by the closed form above; it is
never a free coefficient. The log-likelihood is the ordinary GKw
log-likelihood evaluated at the reconstructed parameter vector,

\$\$\ell(\boldsymbol\beta)=\sum\_{i=1}^{n} w_i \\ \log
f\bigl(y_i;\alpha_i,\beta_i,\gamma_i,\delta_i,\lambda_i\bigr),\$\$

and the score \\\partial\ell/\partial\boldsymbol\beta\\ is obtained by
automatic differentiation through the anchor solve, so the chain rule
through the reparametrization is exact rather than approximated.
Maximization uses
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html) or
[`stats::optim()`](https://rdrr.io/r/stats/optim.html); see
[`gkwq_control()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_control.md).

Standard errors come from the observed information,
[`stats::optimHess()`](https://rdrr.io/r/stats/optim.html) applied to
the automatic-differentiation gradient. The naive \\J^\top H J\\
transformation of the unreparametrized Hessian is **wrong** in a
regression because it omits the curvature term \\\sum_k g_k
\nabla^2\theta_k\\, and is never used. See
[`vcov.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/vcov.gkwqreg.md).

`tau` is fixed, never estimated. The profile likelihood in \\\tau\\ is
exactly flat – for any other level there is a parameter value
reproducing the same distribution – so \\\tau\\ indexes the *question
asked of the data*, not the model fitted to it.

## The multi-part formula

The right-hand side is split into parts by `|`, as in the Formula
package. The contract is fixed and is reported by
[`gkwq_parts()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_parts.md):

1.  **Part one is always the conditional \\\tau\\-quantile**, `mu`.

2.  The remaining parts follow the family's own parameter order
    (\\\alpha\\, \\\beta\\, \\\gamma\\, \\\delta\\, \\\lambda\\) **with
    the anchored parameter removed**, because that one is computed
    rather than estimated.

3.  Omitted trailing parts default to `~ 1`, so `y ~ x` and `y ~ x | 1`
    are the same model.

4.  Supplying more parts than the family has is an error, not a silent
    truncation; the message names the parts the model expects.

Because the anchor changes which parameter disappears, it also changes
what the later parts mean. Always confirm with
`gkwq_parts(family, anchor)`:

    gkwq_parts("kw")            # "mu" "alpha"
    gkwq_parts("kw", "alpha")   # "mu" "beta"
    gkwq_parts("gkw")           # "mu" "alpha" "gamma" "delta" "lambda"

So for `family = "gkw"`, `y ~ x1 | x2 | 1 | x3` puts `x1` in the
quantile, `x2` in \\\alpha\\, an intercept in \\\gamma\\, `x3` in
\\\delta\\, and an intercept in \\\lambda\\. Offsets for a specific part
are written as an [`offset()`](https://rdrr.io/r/stats/offset.html) term
inside that part.

## The seven families

Every family is a constrained GKw, so they nest: `kw` inside `ekw`
inside `kkw` inside `gkw`, and `bkw` and `beta` inside `gkw` as well.
Nested pairs can be compared by a likelihood-ratio test
([`anova.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/anova.gkwqreg.md));
non-nested pairs by AIC, BIC or
[`vuong_test()`](https://evandeilton.github.io/gkwqreg/reference/vuong_test.md).

|  |  |  |  |  |
|----|----|----|----|----|
| family | free parameters | constraints | default anchor | other anchors |
| `kw` | \\\alpha,\beta\\ | \\\gamma=1,\\ \delta=0,\\ \lambda=1\\ | \\\beta\\ | \\\alpha\\ |
| `ekw` | \\\alpha,\beta,\lambda\\ | \\\gamma=1,\\ \delta=0\\ | \\\beta\\ | \\\alpha\\, \\\lambda\\ |
| `kkw` | \\\alpha,\beta,\delta,\lambda\\ | \\\gamma=1\\ | \\\beta\\ | \\\alpha\\, \\\lambda\\ |
| `bkw` | \\\alpha,\beta,\gamma,\delta\\ | \\\lambda=1\\ | \\\beta\\ | \\\alpha\\ |
| `gkw` | \\\alpha,\beta,\gamma,\delta,\lambda\\ | none | \\\beta\\ | \\\alpha\\, \\\lambda\\ |
| `mc` | \\\gamma,\delta,\lambda\\ | \\\alpha=\beta=1\\ | \\\lambda\\ | none |
| `beta` | \\\gamma,\delta\\ | \\\alpha=\beta=\lambda=1\\ | \\\gamma\\ | none |

`mc` fixes \\\alpha=\beta=1\\, so neither is available and the
\\\lambda\\-solve is forced; it collapses to the exact form
\\\lambda=\log z\_\tau/\log\mu\\. `beta` fixes
\\\alpha=\beta=\lambda=1\\, leaving \\\gamma\\ as the only admissible
anchor and the only one solved numerically.
[`gkwq_anchors()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_anchors.md)
lists the admissible anchors for a family, the default first.

Two identifiability facts are enforced rather than left to the user.
When \\\delta=0\\ with \\\gamma\\ free, \\I_u(\gamma,1)=u^{\gamma}\\ and
only the product \\\gamma\lambda\\ is identified, so a
\\\lambda\\-anchor there is refused with an error. And `family = "gkw"`
is weakly identified in *any* parametrization, with information-matrix
condition numbers of order \\10^{8}\\ to \\10^{11}\\; it warns, and
[`summary()`](https://rdrr.io/r/base/summary.html) prints the condition
number so the warning can be checked.

## The anchor is a modeling choice

`anchor` names the parameter eliminated in favour of the conditional
quantile. It looks like a reparametrization and, when every remaining
parameter is regressed on at least the covariates used for `mu`, it is
one: the likelihood is identical whichever anchor is chosen. It stops
being one as soon as `mu` varies with covariates while a nuisance
parameter is held constant. Anchoring on `beta` then asserts "`alpha`
constant, `beta_i = f(mu_i, alpha)`"; anchoring on `alpha` asserts the
reverse, and the two trace different paths through parameter space. In
the design study the two differed by 131 in log-likelihood and by 41% in
the coefficient of interest.

Two anchors on the same data give non-nested models of equal dimension:
compare them by AIC, BIC or Vuong's test, **not** by a likelihood-ratio
test.
[`anova.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/anova.gkwqreg.md)
refuses the comparison rather than printing a number that would mean
nothing. Regressing the nuisance parameters on the same covariates as
`mu` makes the fit largely anchor-insensitive, and is recommended
whenever theory does not dictate an anchor.

The two closed-form solves also fail in **opposite tails**: the
\\\beta\\-solve overflows as \\\mu\to0\\ and the \\\alpha\\-solve as
\\\mu\to1\\. Neither dominates, so where the response mass lies is a
legitimate reason to prefer one; see the `eps_mu` entry of
[`gkwq_control()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_control.md)
for the arithmetic.

## Interpreting the coefficients

Under the default logit link for the quantile part, \\\exp(\beta\_{\mu
j})\\ is the multiplicative effect of a one-unit increase in \\x_j\\ on
the **odds of the conditional \\\tau\\-quantile**,
\\\mu\_\tau/(1-\mu\_\tau)\\. It is *not* an effect on the odds of an
event, and *not* an effect on a mean.
[`summary()`](https://rdrr.io/r/base/summary.html) labels the block
accordingly.

On the scale of the quantile itself the effect is

\$\$\frac{\partial Q\_\tau(x)}{\partial x_j} = \beta\_{\mu
j}\\\frac{d\mu}{d\eta} = \beta\_{\mu j}\\\mu\_\tau(1-\mu\_\tau)
\quad\text{(logit link)},\$\$

reported with delta-method standard errors by
[`marginal_effects()`](https://evandeilton.github.io/gkwqreg/reference/marginal_effects.md).
A covariate that also appears in a nuisance part does not change this
expression: `mu` *is* the quantile, so a nuisance part alters the spread
of the conditional distribution, not the quantile being modelled. That
separation is a genuine convenience of this parametrization.

Nuisance coefficients are on the log scale by default, so
[`exp()`](https://rdrr.io/r/base/Log.html) of one is a multiplicative
effect on the corresponding shape parameter.

## What [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns

**[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns the
fitted conditional `tau`-quantile, not the conditional mean.** This is
the single most important difference from `gkwreg`, which fits the same
seven families with mean semantics. A `type = "mean"` escape hatch
exists but is never the default. For the same reason the class
deliberately does not inherit from `"gkwreg"`: a missing method must
error rather than silently dispatch mean semantics on a quantile object.

## References

Carrasco, J. M. F., Ferrari, S. L. P. and Cordeiro, G. M. (2010). A
generalized Kumaraswamy distribution.

Chernozhukov, V., Fernandez-Val, I. and Galichon, A. (2010). Quantile
and probability curves without crossing. *Econometrica* **78**,
1093-1125.

Cox, D. R. and Reid, N. (1987). Parameter orthogonality and approximate
conditional inference. *Journal of the Royal Statistical Society B*
**49**, 1-39.

Koenker, R. and Bassett, G. (1978). Regression quantiles. *Econometrica*
**46**, 33-50.

Mitnik, P. A. and Baek, S. (2013). The Kumaraswamy distribution:
median-dispersion re-parameterizations. *Statistical Papers* **54**,
177-192.
[doi:10.1007/s00362-011-0417-y](https://doi.org/10.1007/s00362-011-0417-y)

## See also

[`gkwq_parts()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_parts.md)
and
[`gkwq_anchors()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_anchors.md)
for the formula and anchor contracts;
[`gkwq_control()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_control.md)
for optimizer, starting-value and covariance settings;
[`predict.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/predict.gkwqreg.md),
[`fitted.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/fitted.gkwqreg.md)
and
[`residuals.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/residuals.gkwqreg.md)
for what comes out of a fit;
[`marginal_effects()`](https://evandeilton.github.io/gkwqreg/reference/marginal_effects.md)
for effects on the quantile scale;
[`vcov.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/vcov.gkwqreg.md)
and
[`confint.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/confint.gkwqreg.md)
for inference;
[`anova.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/anova.gkwqreg.md),
[`vuong_test()`](https://evandeilton.github.io/gkwqreg/reference/vuong_test.md)
and
[`compare_families()`](https://evandeilton.github.io/gkwqreg/reference/compare_families.md)
for model choice;
[`quantile_process()`](https://evandeilton.github.io/gkwqreg/reference/quantile_process.md),
[`check_crossing()`](https://evandeilton.github.io/gkwqreg/reference/check_crossing.md),
[`rearrange()`](https://evandeilton.github.io/gkwqreg/reference/rearrange.md)
and
[`pinball()`](https://evandeilton.github.io/gkwqreg/reference/pinball.md)
for quantile-specific tooling;
[`plot.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/plot.gkwqreg.md)
for diagnostics;
[`gkwq_quantile()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_quantile.md)
for the closed-form quantile function on its own.

Other model fitting:
[`gkwq_control()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_control.md)

## Examples

``` r
## -------------------------------------------------------------------------
## 1. A Kumaraswamy median regression, end to end
## -------------------------------------------------------------------------
## Simulate a response whose conditional MEDIAN is logit-linear in x1 and x2.
## The generator is the beta anchor itself: a Kumaraswamy with
## beta = log(1 - tau) / log(1 - mu^alpha) has its tau-quantile exactly at mu.
set.seed(2024)
n  <- 300
x1 <- runif(n, -2, 2)
x2 <- rbinom(n, 1, 0.5)
mu <- plogis(0.4 + 1.1 * x1 - 0.6 * x2)      # the true conditional median
y  <- gkwdist::rkw(n, alpha = 2, beta = log(1 - 0.5) / log(1 - mu^2))
dat <- data.frame(y = y, x1 = x1, x2 = x2)

fit <- gkwqreg(y ~ x1 + x2, data = dat, tau = 0.5, family = "kw")
summary(fit)
#> 
#> Generalized Kumaraswamy quantile regression
#> family: kw   tau: 0.5   anchor: beta
#> 
#> Call:
#> gkwqreg(formula = y ~ x1 + x2, data = dat, tau = 0.5, family = "kw")
#> 
#> Quantile residuals:
#>     Min      1Q  Median      3Q     Max 
#> -2.2296 -0.6736 -0.0368  0.6642  2.9883 
#> 
#> Conditional 0.5-quantile (link logit) -- coefficients are effects on the LOG QUANTILE ODDS log(mu/(1-mu)):
#>                Estimate Std. Error z value Pr(>|z|)    
#> mu:(Intercept)  0.39579    0.08649   4.576 4.73e-06 ***
#> mu:x1           1.05114    0.04346  24.189  < 2e-16 ***
#> mu:x2          -0.61862    0.09053  -6.834 8.28e-12 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> alpha (link log):
#>                   Estimate Std. Error z value Pr(>|z|)    
#> alpha:(Intercept)  0.79052    0.05808   13.61   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Anchored parameter: beta (computed from the quantile, not estimated)
#> Log-likelihood: 191.9 on 4 Df   AIC: -375.9   BIC: -361.1
#> Pinball loss: 0.06714   Pseudo-R1: 0.4726
#> Empirical coverage: 0.5267 (target 0.5)
#> Covariance: expected   Information condition number: 20.93
## The "mu" block should recover (0.4, 1.1, -0.6) up to sampling error, and
## the "alpha" block log(2) = 0.693 -- alpha is modelled on the log scale.
## The footer reports "Anchored parameter: beta", the reminder that beta was
## computed from the fitted median rather than estimated, plus the pinball
## loss and the empirical coverage, which should sit near tau.

## -------------------------------------------------------------------------
## 2. Reading the coefficients: quantile odds, not event odds
## -------------------------------------------------------------------------
exp(coef(fit, part = "mu"))
#> (Intercept)          x1          x2 
#>   1.4855519   2.8609146   0.5386851 
## exp(b_x1) should be close to exp(1.1) = 3.0: a one-unit rise in x1
## multiplies the ODDS OF THE CONDITIONAL MEDIAN, mu/(1 - mu), by about 3.
## It is not the odds of an event, and not an effect on a mean.

## On the quantile scale itself the effect is b_j * mu * (1 - mu), i.e. in
## probability units, averaged over the sample:
marginal_effects(fit)
#> 
#> Marginal effects on the conditional 0.5-quantile (averaged over the sample)
#> dQ(tau|x)/dx_j, logit link
#> 
#>  variable  effect std.error   lower    upper
#>        x1  0.1879  0.003905  0.1802  0.19552
#>        x2 -0.1106  0.016020 -0.1420 -0.07917

## -------------------------------------------------------------------------
## 3. Fitted values are quantiles
## -------------------------------------------------------------------------
head(fitted(fit))                     # conditional medians
#> [1] 0.8596617 0.2736838 0.7602655 0.6480307 0.4004488 0.7760286
head(fitted(fit, type = "mean"))      # conditional means, by quadrature
#> [1] 0.7860434 0.2822266 0.7109591 0.6230794 0.4057070 0.7229222
mean(dat$y <= fitted(fit))            # empirical coverage; target is tau
#> [1] 0.5266667

## -------------------------------------------------------------------------
## 4. Prediction: one fit, a whole conditional distribution
## -------------------------------------------------------------------------
nd <- data.frame(x1 = c(-1, 0, 1), x2 = 0)
predict(fit, newdata = nd)                        # conditional medians
#> [1] 0.3417838 0.5976749 0.8095252
## Other quantiles read off the SAME fitted distribution -- these cannot
## cross, because they are quantiles of one proper distribution function:
predict(fit, newdata = nd, tau = c(0.1, 0.5, 0.9))
#>        tau=0.1   tau=0.5   tau=0.9
#> [1,] 0.1481773 0.3417838 0.5604450
#> [2,] 0.2732389 0.5976749 0.8639080
#> [3,] 0.4091304 0.8095252 0.9827776

# \donttest{
## -------------------------------------------------------------------------
## 5. The formula grammar
## -------------------------------------------------------------------------
gkwq_parts("kw")             # "mu" "alpha" -- part 1 is always the quantile
#> [1] "mu"    "alpha"
gkwq_parts("kw", "alpha")    # "mu" "beta"  -- the anchor is what disappears
#> [1] "mu"   "beta"
gkwq_parts("gkw")            # five parts, in family order
#> [1] "mu"     "alpha"  "gamma"  "delta"  "lambda"

## Give the shape parameter a regression of its own. Omitted trailing parts
## default to ~ 1, so `y ~ x1 + x2` is `y ~ x1 + x2 | 1`.
fit_sat <- gkwqreg(y ~ x1 + x2 | x1 + x2, data = dat, tau = 0.5,
                   family = "kw")
anova(fit, fit_sat)          # nested, same anchor: an ordinary LR test
#> Likelihood-ratio test for Generalized Kumaraswamy quantile regression
#> tau = 0.5, anchor = beta
#>      Df logLik     AIC  Chisq Chi Df Pr(>Chisq)
#> [1,]  4 191.94 -375.88                         
#> [2,]  6 192.52 -373.03 1.1561      2      0.561

## -------------------------------------------------------------------------
## 6. The anchor is a modeling choice
## -------------------------------------------------------------------------
fb <- gkwqreg(y ~ x1, data = dat, tau = 0.5, family = "kw", anchor = "beta")
fa <- gkwqreg(y ~ x1, data = dat, tau = 0.5, family = "kw", anchor = "alpha")
c(beta = as.numeric(logLik(fb)), alpha = as.numeric(logLik(fa)))
#>      beta     alpha 
#> 169.37780  76.00107 
## Same number of parameters, different likelihood: these are NON-NESTED
## models, so compare them by AIC/BIC or a Vuong test, never by an LR test.
vuong_test(fb, fa)
#> 
#> Vuong test for non-nested quantile regression fits
#> tau = 0.5, BIC-corrected
#>   model 1: kw / anchor beta
#>   model 2: kw / anchor alpha
#> 
#>   z = 6.602, p-value = 4.065e-11
#>   model 1 is favoured

## Saturating the nuisance part makes the choice nearly immaterial, which is
## the recommended default when theory does not dictate an anchor:
fbs <- gkwqreg(y ~ x1 | x1, data = dat, tau = 0.5, family = "kw",
               anchor = "beta")
fas <- gkwqreg(y ~ x1 | x1, data = dat, tau = 0.5, family = "kw",
               anchor = "alpha")
c(beta = coef(fbs)[["mu:x1"]], alpha = coef(fas)[["mu:x1"]])
#>     beta    alpha 
#> 1.079045 1.066362 

## -------------------------------------------------------------------------
## 7. A grid of quantile levels
## -------------------------------------------------------------------------
fits <- gkwqreg(y ~ x1, data = dat, tau = c(0.1, 0.25, 0.5, 0.75, 0.9),
                family = "kw")
round(coef(fits), 3)   # one column per level: the estimated quantile process
#>                     0.10   0.25  0.50  0.75  0.90
#> mu:(Intercept)    -1.454 -0.736 0.028 0.817 1.594
#> mu:x1              0.771  0.880 1.053 1.287 1.558
#> alpha:(Intercept)  0.660  0.665 0.664 0.655 0.641
## The slope on x1 growing with tau is heteroskedasticity that a mean
## regression would average away.
## Independently fitted levels are not ordered by construction; here none
## cross, but nothing guaranteed that, which is exactly what this reports:
check_crossing(fits)
#> 
#> Quantile crossing check
#> levels: 0.10, 0.25, 0.50, 0.75, 0.90
#> source: independently fitted models, one per level
#> 
#>   rows with a crossing: 0 of 300 (0.00%)
#>   worst violation     : 0.000e+00

## -------------------------------------------------------------------------
## 8. Choosing a family
## -------------------------------------------------------------------------
fit_ekw <- gkwqreg(y ~ x1 + x2, data = dat, tau = 0.5, family = "ekw")
anova(fit, fit_ekw)    # kw is nested in ekw: one degree of freedom
#> Likelihood-ratio test for Generalized Kumaraswamy quantile regression
#> tau = 0.5, anchor = beta
#>      Df logLik     AIC  Chisq Chi Df Pr(>Chisq)
#> [1,]  4 191.94 -375.88                         
#> [2,]  5 192.06 -374.13 0.2522      1     0.6155
## Check loss is the criterion a quantile fit targets, so read it beside
## AIC -- and read the `converged` column before either of them, since a
## richer family that failed to converge can still post the best logLik:
compare_families(fit, families = c("kw", "ekw", "kkw", "beta"))
#>   family anchor df   logLik       AIC       BIC    pinball converged
#> 1     kw   beta  4 191.9378 -375.8756 -361.0604 0.06713950      TRUE
#> 2    ekw   beta  5 192.0639 -374.1278 -355.6089 0.06712997      TRUE
#> 3    kkw   beta  6 192.3088 -372.6177 -350.3950 0.06709513     FALSE
#> 4   beta  gamma  4  78.2969 -148.5938 -133.7787 0.07321204      TRUE
# }
```
