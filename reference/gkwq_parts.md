# The multi-part formula contract for a family and anchor

Reports, in order, the names of the right-hand formula parts that
[`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md)
will read for a given family and anchor. The contract is exported as a
function rather than merely described in prose, so that user code, error
messages and this manual cannot disagree about it.

## Usage

``` r
gkwq_parts(family = "kw", anchor = NULL)
```

## Arguments

- family:

  Character, one of `"kw"`, `"ekw"`, `"kkw"`, `"bkw"`, `"gkw"`, `"mc"`,
  `"beta"`. Partial matching is applied. An unrecognised name is an
  error.

- anchor:

  Character or `NULL`. The parameter eliminated in favour of the
  conditional quantile. `NULL`, the default, selects the family default,
  which is the first element of
  [`gkwq_anchors()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_anchors.md).
  An anchor the family does not admit raises an error naming those it
  does.

## Value

A character vector naming the formula parts in order. Its length is one
plus the number of free parameters remaining after the anchor is
removed, and its first element is always `"mu"`.

## Details

A `gkwqreg` model formula carries one right-hand part per modelled
parameter, the parts separated by `|`, as in `y ~ x1 + x2 | z | 1`. Two
rules fix the meaning of every part.

1.  **Part one is always `mu`**, the conditional \\\tau\\-quantile of
    the response, \\\mu\_\tau(x) = Q(\tau \mid x)\\. This holds for
    every family and every anchor without exception. It is the parameter
    that the reparametrization creates, and the reason the model is a
    quantile regression rather than a mean regression.

2.  **The remaining parts follow the family's own parameter order, with
    the anchored parameter removed.** The anchor is the parameter
    eliminated in favour of \\\mu\\; having been solved for, it is no
    longer free, so it receives no formula part. See
    [`gkwq_anchors()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_anchors.md).

The Kumaraswamy family, for instance, has parameters \\(\alpha,
\beta)\\. Anchoring on \\\beta\\, its default, leaves `mu | alpha`;
anchoring on \\\alpha\\ leaves `mu | beta`. These are different models,
not different labels for one model, unless every remaining part is
regressed on the same covariates as `mu`.

## Omitted and excess parts

Trailing parts that are omitted default to `~ 1`, an intercept only, so
`y ~ x` and `y ~ x | 1` describe the same Kumaraswamy model. Supplying
**more** parts than the model has is an error, and the message names the
parts the model expected. This is deliberate: silently discarding the
surplus would conceal a genuine mistake, most often a formula written
for a different anchor or a different family.

## The family parameter sets

Each family holds some of the five Generalized Kumaraswamy parameters
\\(\alpha, \beta, \gamma, \delta, \lambda)\\ at structural constants and
leaves the rest free. The free set, listed in the family's own order, is
what rule two refers to.

|  |  |  |
|----|----|----|
| family | free parameters | held fixed |
| `"kw"` | `alpha`, `beta` | \\\gamma = 1\\, \\\delta = 0\\, \\\lambda = 1\\ |
| `"ekw"` | `alpha`, `beta`, `lambda` | \\\gamma = 1\\, \\\delta = 0\\ |
| `"kkw"` | `alpha`, `beta`, `delta`, `lambda` | \\\gamma = 1\\ |
| `"bkw"` | `alpha`, `beta`, `gamma`, `delta` | \\\lambda = 1\\ |
| `"gkw"` | `alpha`, `beta`, `gamma`, `delta`, `lambda` | none |
| `"mc"` | `gamma`, `delta`, `lambda` | \\\alpha = 1\\, \\\beta = 1\\ |
| `"beta"` | `gamma`, `delta` | \\\alpha = 1\\, \\\beta = 1\\, \\\lambda = 1\\ |

Removing the anchor from the free set and prefixing `mu` gives the
parts. Under the default anchor the number of parts equals the number of
free parameters, since one parameter is removed and `mu` is added.

## See also

[`gkwq_anchors()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_anchors.md)
for which anchors a family admits and why;
[`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md)
for the model these parts describe.

## Examples

``` r
## Rule one: the first part is always the conditional quantile.
gkwq_parts("kw")
#> [1] "mu"    "alpha"
## [1] "mu"    "alpha"

## Rule two: the anchored parameter is dropped from the family's own order.
gkwq_parts("kw", anchor = "alpha")
#> [1] "mu"   "beta"
## [1] "mu"   "beta"

## The five-parameter family under each anchor it admits.
for (a in gkwq_anchors("gkw")) {
  cat(sprintf("%-7s %s\n", a, paste(gkwq_parts("gkw", a), collapse = " | ")))
}
#> beta    mu | alpha | gamma | delta | lambda
#> alpha   mu | beta | gamma | delta | lambda
#> lambda  mu | alpha | beta | gamma | delta
## beta    mu | alpha | gamma | delta | lambda
## alpha   mu | beta | gamma | delta | lambda
## lambda  mu | alpha | beta | gamma | delta

## The whole contract at default anchors, at a glance.
for (f in c("kw", "ekw", "kkw", "bkw", "gkw", "mc", "beta")) {
  cat(sprintf("%-5s %s\n", f, paste(gkwq_parts(f), collapse = " | ")))
}
#> kw    mu | alpha
#> ekw   mu | alpha | lambda
#> kkw   mu | alpha | delta | lambda
#> bkw   mu | alpha | gamma | delta
#> gkw   mu | alpha | gamma | delta | lambda
#> mc    mu | gamma | delta
#> beta  mu | delta
## kw    mu | alpha
## ekw   mu | alpha | lambda
## kkw   mu | alpha | delta | lambda
## bkw   mu | alpha | gamma | delta
## gkw   mu | alpha | gamma | delta | lambda
## mc    mu | gamma | delta
## beta  mu | delta

## Trailing parts may be omitted and default to ~ 1.
set.seed(1)
dat <- data.frame(x = runif(50), y = runif(50, 0.2, 0.8))
f1 <- gkwqreg(y ~ x,     data = dat, tau = 0.5, family = "kw")
f2 <- gkwqreg(y ~ x | 1, data = dat, tau = 0.5, family = "kw")
all.equal(coef(f1), coef(f2))
#> [1] TRUE
## [1] TRUE

## Supplying more parts than the model has is an error, not a silent drop,
## and the message names the parts that were expected.
err <- try(gkwqreg(y ~ x | x | x, data = dat, tau = 0.5, family = "kw"),
           silent = TRUE)
cat(conditionMessage(attr(err, "condition")), "\n")
#> formula has 3 right-hand parts but family ‘kw’ with anchor ‘beta’ takes at most 2: mu | alpha.
#>   Parts are separated by `|`, in the order given by gkwq_parts("kw", "beta"). 
## formula has 3 right-hand parts but family 'kw' with anchor 'beta' takes
##   at most 2: mu | alpha.
##   Parts are separated by `|`, in the order given by gkwq_parts("kw", "beta").
```
