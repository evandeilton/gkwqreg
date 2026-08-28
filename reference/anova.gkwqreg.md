# Likelihood-ratio comparison of nested quantile regression fits

Compares two or more fitted quantile regressions that stand in a nesting
relation, by the classical likelihood-ratio test, and tabulates the
result in the idiom of
[`stats::anova()`](https://rdrr.io/r/stats/anova.html). Because the
seven Generalized Kumaraswamy families are genuine parameter
restrictions of one another, choosing among them is here an ordinary
hypothesis test. It does not require the Vuong machinery that a
collection of mutually non-nested families would force upon it, and
which packages built on such a collection have no way to avoid.

## Usage

``` r
# S3 method for class 'gkwqreg'
anova(object, ...)
```

## Arguments

- object:

  A `"gkwqreg"` fit.

- ...:

  One or more further `"gkwqreg"` fits. Every fit must use the same
  quantile level `tau`, the same anchor, and the same number of
  observations; the section *What this method refuses* explains why each
  requirement is enforced rather than merely advised. The order of the
  arguments does not matter: the fits are sorted internally by number of
  estimated coefficients.

## Value

A data frame of class `"anova.gkwqreg"`, inheriting from `"anova"` and
`"data.frame"`, with one row per fit sorted by increasing `Df` and a
`"heading"` attribute recording the common `tau` and anchor. Its columns
are

- `family`:

  the family of the fit in that row.

- `Df`:

  \\p\\, the number of estimated coefficients.

- `logLik`:

  \\\ell(\hat\theta)\\, the maximized log-likelihood.

- `AIC`:

  \\-2\ell(\hat\theta) + 2p\\.

- `Chisq`:

  the likelihood-ratio statistic against the row above; `NA` in the
  first row.

- `Chi Df`:

  the difference in `Df` against the row above; `NA` in the first row.

- `Pr(>Chisq)`:

  the upper-tail chi-squared probability; `NA` in the first row and
  wherever `Chi Df` is not positive.

The [`print()`](https://rdrr.io/r/base/print.html) method formats the
table through
[`stats::printCoefmat()`](https://rdrr.io/r/stats/printCoefmat.html) and
prepends the heading.

## Details

Write \\\ell(\hat\theta)\\ for the maximized log-likelihood of a fit and
\\p\\ for its number of estimated coefficients, reported in the `Df`
column. The anchored parameter is not counted in \\p\\: it is computed
from the conditional quantile by the anchor solve, not estimated.

For a model \\M_0\\ nested inside a model \\M_1\\, with \\p_0 \< p_1\\,
the statistic reported in the `Chisq` column is

\$\$LR = 2 \left\\ \ell_1(\hat\theta_1) - \ell_0(\hat\theta_0)
\right\\,\$\$

whose null distribution, under \\M_0\\ and the usual regularity
conditions, is asymptotically chi-squared on \\p_1 - p_0\\ degrees of
freedom. The columns `Chi Df` and `Pr(>Chisq)` hold \\p_1 - p_0\\ and
the upper-tail probability \\\Pr(\chi^2\_{p_1 - p_0} \> LR)\\. The `AIC`
column reports \\-2\ell(\hat\theta) + 2p\\ for each row; unlike the test
itself, it is comparable down the whole table and requires no nesting.

With more than two fits the table is **sequential**: after sorting by
dimension, each row is tested against the row immediately above it,
never against the first row. Nothing verifies that consecutive rows
really are nested; that remains the analyst's judgement, and the map in
the next section is the guide.

## The families are genuinely nested

Every family is the five-parameter `gkw` distribution with some of its
parameters held at neutral values, so each is a restriction of the ones
above it in the lattice.

|  |  |  |
|----|----|----|
| family | free parameters | held fixed |
| `"gkw"` | `alpha`, `beta`, `gamma`, `delta`, `lambda` | none |
| `"bkw"` | `alpha`, `beta`, `gamma`, `delta` | \\\lambda = 1\\ |
| `"kkw"` | `alpha`, `beta`, `delta`, `lambda` | \\\gamma = 1\\ |
| `"ekw"` | `alpha`, `beta`, `lambda` | \\\gamma = 1\\, \\\delta = 0\\ |
| `"kw"` | `alpha`, `beta` | \\\gamma = 1\\, \\\delta = 0\\, \\\lambda = 1\\ |
| `"mc"` | `gamma`, `delta`, `lambda` | \\\alpha = \beta = 1\\ |
| `"beta"` | `gamma`, `delta` | \\\alpha = \beta = \lambda = 1\\ |

The restrictions compose into chains, and it is along a chain that this
method is valid:

1.  `"kw"` is `"ekw"` with \\\lambda = 1\\; `"ekw"` is `"kkw"` with
    \\\delta = 0\\; `"kkw"` is `"gkw"` with \\\gamma = 1\\.

2.  `"kw"` is also `"bkw"` with \\\gamma = 1\\ and \\\delta = 0\\, and
    `"bkw"` is `"gkw"` with \\\lambda = 1\\.

3.  `"beta"` is `"mc"` with \\\lambda = 1\\, and `"mc"` is `"gkw"` with
    \\\alpha = \beta = 1\\; `"beta"` is likewise `"bkw"` with \\\alpha =
    \beta = 1\\.

The reparametrization does not disturb any of this. The anchor solve
eliminates the same parameter from both members of a pair, so a
restriction on the original family is a restriction on the anchored one
and \\LR\\ is the same quantity it always was.

Two consequences are worth stating. First, a larger family cannot attain
a smaller maximum log-likelihood than a family nested inside it, so the
`logLik` column is non-decreasing down the table as a matter of
mathematics, not of luck; the package's test suite checks that
invariant. A `logLik` that *falls* as `Df` rises is therefore not
evidence about the data but evidence that one of the fits has not
reached its maximum, and the offending fit should be refitted or
discarded. Second, not every pair of families is nested: `"kw"` and
`"mc"`, for instance, restrict disjoint sets of parameters and neither
contains the other. Such a pair belongs to
[`vuong_test()`](https://evandeilton.github.io/gkwqreg/reference/vuong_test.md)
or
[`compare_families()`](https://evandeilton.github.io/gkwqreg/reference/compare_families.md),
not here.

A practical corollary of the anchor guard is that `"mc"` admits only the
`lambda` anchor and `"beta"` only the `gamma` anchor, whereas every
Kumaraswamy member defaults to `beta`. Those two families therefore
cannot reach this method in company with a Kumaraswamy member at all,
even where the nesting exists on paper.

## A caution at the boundary

The restriction \\\delta = 0\\ – the step from `"kkw"` down to `"ekw"`,
or from `"bkw"` down to `"kw"` – lies on the boundary of the parameter
space: under the default log link \\\delta\\ is strictly positive and
reaches zero only in the limit. For such a restriction the null
distribution of \\LR\\ is not chi-squared on one degree of freedom but a
mixture with an atom at zero (Chernoff 1954; Self and Liang 1987), so
the p-value reported here is conservative and errs towards keeping the
smaller model. Its signature in practice is an \\LR\\ of essentially
zero with a p-value close to one, as in the third row of the example
below. Restrictions setting \\\gamma\\ or \\\lambda\\ to 1 are interior
and need no such caveat.

## What this method refuses, and why

Two comparisons raise an error instead of returning a number. Both
refusals are deliberate, and in both cases the arithmetic would have
produced a perfectly presentable statistic that means nothing.

**Fits at different quantile levels.** The level `tau` indexes the
question being asked, not the model: it never enters the parameter
vector, and the profile likelihood in `tau` is exactly flat. Each level
therefore defines its own likelihood, and two likelihoods for two
different questions are not comparable at all – neither nested, nor
non-nested, but incommensurable. A difference of log-likelihoods between
a `tau = 0.5` fit and a `tau = 0.9` fit is a difference between the fit
of one model to one target and another model to another target. Note
that a generic likelihood-ratio routine written for other model classes
has no way to know this and will happily report a significant result;
see
[`lrtest()`](https://evandeilton.github.io/gkwqreg/reference/lrtest.md).

**Fits with different anchors.** Two anchors on the same data give
*non-nested models of equal dimension*. Neither is a restriction of the
other, so no null hypothesis relates them; and \\p_1 - p_0 = 0\\ leaves
no degrees of freedom against which to refer a statistic. The difference
in log-likelihood between them can nonetheless be substantial, because
the anchor is a modeling choice and not a relabelling whenever the
nuisance parameters are unsaturated (see
[`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md)).
The correct instrument is
[`vuong_test()`](https://evandeilton.github.io/gkwqreg/reference/vuong_test.md),
which is built for models of equal dimension that do not nest; AIC and
BIC are also legitimate.

A third case is warned about rather than refused. If two fits of *equal*
dimension reach the test by some other route – same family, same anchor,
same level, but different covariates of the same count – a warning is
issued and the corresponding `Chisq` and `Pr(>Chisq)` entries are `NA`.

## See also

[`lrtest()`](https://evandeilton.github.io/gkwqreg/reference/lrtest.md)
for the same test under the name `lmtest` users expect;
[`vuong_test()`](https://evandeilton.github.io/gkwqreg/reference/vuong_test.md)
for non-nested comparisons, including two anchors;
[`compare_families()`](https://evandeilton.github.io/gkwqreg/reference/compare_families.md)
to sweep every family at once;
[`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md)
for what the anchor is and why it is a modeling argument.

## Examples

``` r
## A response generated from an EKW distribution: the true lambda is 4, so
## the two-parameter Kumaraswamy family is genuinely wrong, and its larger
## relatives are genuinely right.
set.seed(2026)
n  <- 500
x  <- runif(n, -2, 2)
mu <- plogis(0.3 + 1.2 * x)          # the true conditional median
b  <- log1p(-0.5^(1 / 4)) / log1p(-mu^2)   # the beta anchor solve at tau = 0.5
y  <- gkwdist::rekw(n, alpha = 2, beta = b, lambda = 4)
d  <- data.frame(y = y, x = x)

m_kw  <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
m_ekw <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "ekw")
m_kkw <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kkw")

## The chain kw inside ekw inside kkw, tested one step at a time.
anova(m_kw, m_ekw, m_kkw)
#> Likelihood-ratio test for Generalized Kumaraswamy quantile regression
#> tau = 0.5, anchor = beta
#>      Df logLik     AIC   Chisq Chi Df Pr(>Chisq)    
#> [1,]  3 482.58 -959.16                              
#> [2,]  4 496.90 -985.80 28.6406      1  8.713e-08 ***
#> [3,]  5 496.91 -983.81  0.0093      1      0.923    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##      Df logLik     AIC   Chisq Chi Df Pr(>Chisq)
## [1,]  3 482.58 -959.16
## [2,]  4 496.90 -985.80 28.6406      1  8.713e-08
## [3,]  5 496.91 -983.81  0.0093      1      0.923
##
## Row 2 tests lambda = 1, that is kw inside ekw, and rejects it decisively:
## the extra shape parameter is real, exactly as the simulation was built.
## Row 3 tests delta = 0, that is ekw inside kkw, and finds nothing -- note
## the statistic of essentially zero, the boundary signature discussed above.
## Note also that logLik never falls as Df rises, as the nesting requires.

## The refusals. Neither of these returns a number.
m_tau9 <- gkwqreg(y ~ x, data = d, tau = 0.9, family = "kw")
try(anova(m_kw, m_tau9))       # different levels: separate likelihoods
#> Error : models at different quantile levels cannot be compared: each level is a separate likelihood. Levels seen: 0.5, 0.9.

m_alpha <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw",
                   anchor = "alpha")
try(anova(m_kw, m_alpha))      # different anchors: non-nested, equal dimension
#> Error : models with different anchors cannot be compared by a likelihood-ratio test: they are non-nested models of equal dimension. Compare them by AIC, BIC or a Vuong test instead. Anchors seen: beta, alpha.

## The comparison anova() refused is exactly the one vuong_test() makes.
vuong_test(m_kw, m_alpha)
#> 
#> Vuong test for non-nested quantile regression fits
#> tau = 0.5, BIC-corrected
#>   model 1: kw / anchor beta
#>   model 2: kw / anchor alpha
#> 
#>   z = 7.642, p-value = 2.137e-14
#>   model 1 is favoured
##   z = 7.642, p-value = 2.137e-14 -- model 1, the beta anchor, is favoured.
```
