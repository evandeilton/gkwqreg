# Admissible anchors for a family

The *anchor* is the parameter eliminated in favour of the conditional
quantile when a family is reparametrized at a level \\\tau\\ fixed in
advance. `gkwq_anchors()` reports which of a family's parameters may
play that role, the first being the default that
[`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md)
uses when `anchor = NULL`.

## Usage

``` r
gkwq_anchors(family = "kw")
```

## Arguments

- family:

  Character, one of `"kw"`, `"ekw"`, `"kkw"`, `"bkw"`, `"gkw"`, `"mc"`,
  `"beta"`. Partial matching is applied. An unrecognised name is an
  error.

## Value

A character vector of the anchor names admissible for `family`, ordered
with the family default first. Supplying any other value as the `anchor`
argument of
[`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md)
or
[`gkwq_parts()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_parts.md)
raises an error that names these.

## Details

Fix a level \\\tau \in (0,1)\\ and let \\\mu\\ denote the target
conditional quantile. Anchoring means solving

\$\$Q(\tau; \alpha, \beta, \gamma, \delta, \lambda) = \mu\$\$

for one parameter and substituting the solution back into the
likelihood. What remains is a distribution indexed by \\\mu\\ together
with the parameters left free, in which \\\mu\\ *is* the conditional
\\\tau\\-quantile and can therefore be modelled directly through a link
function. The anchored parameter is no longer free and takes no formula
part; see
[`gkwq_parts()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_parts.md).

A parameter is admissible as an anchor only when that equation has a
unique solution for it at every \\(\mu, \tau)\\ in \\(0,1)^2\\ and for
every admissible value of the parameters left free. Three of the five
parameters admit a closed-form solve. Writing \\z\_\tau\\ for
`qbeta(tau, gamma, delta + 1)`,

\$\$\beta = \frac{\log(1 - z\_\tau^{1/\lambda})}{\log(1 - \mu^\alpha)},
\qquad \alpha = \frac{\log\\1 - (1 -
z\_\tau^{1/\lambda})^{1/\beta}\\}{\log \mu}, \qquad \lambda = \frac{\log
z\_\tau}{\log\\1 - (1 - \mu^\alpha)^\beta\\}.\$\$

Each is a ratio of two logarithms of quantities lying in \\(0,1)\\ and
is consequently strictly positive and finite for every \\\mu\\ and
\\\tau\\ in \\(0,1)\\. No box constraint on the anchored parameter is
ever needed, which is what makes an unconstrained optimizer safe for
this model.

## Which families admit which anchors

|  |  |  |
|----|----|----|
| family | admissible anchors (default first) | how the anchor is obtained |
| `"kw"` | `"beta"`, `"alpha"` | closed form |
| `"ekw"` | `"beta"`, `"alpha"`, `"lambda"` | closed form |
| `"kkw"` | `"beta"`, `"alpha"`, `"lambda"` | closed form |
| `"bkw"` | `"beta"`, `"alpha"` | closed form |
| `"gkw"` | `"beta"`, `"alpha"`, `"lambda"` | closed form |
| `"mc"` | `"lambda"` | closed form, \\\lambda = \log z\_\tau / \log \mu\\ |
| `"beta"` | `"gamma"` | one-dimensional root find |

A parameter that a family holds fixed as a structural constant cannot be
an anchor, since it is not free to be solved for, and that alone
accounts for most of the table. The `"mc"` family fixes \\\alpha = \beta
= 1\\, so neither is available and the \\\lambda\\-solve is forced; with
those two substitutions the general expression above collapses to the
elementary form \\\lambda = \log z\_\tau / \log \mu\\. The `"beta"`
family fixes \\\alpha = \beta = \lambda = 1\\, leaving only \\\gamma\\
and \\\delta\\ as candidates.

## Why gamma is admissible and delta never is

In the `"beta"` family both free parameters enter the quantile solely
through \\z\_\tau\\, so \\Q(\tau) = \mu\\ reduces to

\$\$I\_\mu(\gamma, \delta + 1) = \tau,\$\$

where \\I\\ is the regularized incomplete beta function, that is
`pbeta`. There is no closed-form inverse, but admissibility is a
question of existence and uniqueness rather than of closed form.

Read as an equation in \\\gamma\\, the left-hand side is continuous and
strictly decreasing, tending to \\1\\ as \\\gamma \to 0\\ and to \\0\\
as \\\gamma \to \infty\\. It therefore sweeps the whole of \\(0,1)\\ and
a root exists, uniquely, for every \\\tau \in (0,1)\\ and every \\\mu\\.
A bracketed one-dimensional root find is all that is required, and
\\\gamma\\ is admissible.

Read instead as an equation in \\\delta\\, the left-hand side is
increasing, but \\\delta\\ is bounded below by zero and
\\I\_\mu(\gamma, 1) = \mu^\gamma\\. Over the admissible range \\\delta
\in \[0, \infty)\\ the left-hand side spans only \\\[\mu^\gamma, 1)\\. A
\\\delta\\-solve therefore exists only when \\\tau \> \mu^\gamma\\, a
condition that nothing in the data or the model guarantees: a low level
combined with a large fitted quantile simply breaks it. This is why
\\\gamma\\ is the only anchor the `"beta"` family admits, and why
\\\delta\\ is an anchor in no family at all – wherever it is free it is
bounded below by zero and enters only through \\z\_\tau\\, so the same
range restriction recurs.

## The anchor is a modelling choice

When every free parameter is regressed on the same covariates as
\\\mu\\, the choice of anchor leaves the likelihood unchanged and is a
genuine reparametrization. It ceases to be innocuous as soon as \\\mu\\
varies with covariates while a free parameter is held constant, because
the anchor then decides which parameter absorbs that variation. Two
anchors applied to the same data give non-nested models of equal
dimension: compare them by AIC, BIC or
[`vuong_test()`](https://evandeilton.github.io/gkwqreg/reference/vuong_test.md),
never by a likelihood-ratio test. See
[`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md)
for the full discussion.

## See also

[`gkwq_parts()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_parts.md)
for the formula parts that follow from a choice of anchor;
[`gkwq_quantile()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_quantile.md)
for the quantile function being solved;
[`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md)
for fitting.

## Examples

``` r
## Every family, default anchor first.
for (f in c("kw", "ekw", "kkw", "bkw", "gkw", "mc", "beta")) {
  cat(sprintf("%-5s %s\n", f, paste(gkwq_anchors(f), collapse = ", ")))
}
#> kw    beta, alpha
#> ekw   beta, alpha, lambda
#> kkw   beta, alpha, lambda
#> bkw   beta, alpha
#> gkw   beta, alpha, lambda
#> mc    lambda
#> beta  gamma
## kw    beta, alpha
## ekw   beta, alpha, lambda
## kkw   beta, alpha, lambda
## bkw   beta, alpha
## gkw   beta, alpha, lambda
## mc    lambda            <- alpha and beta are fixed at 1
## beta  gamma             <- alpha, beta and lambda are fixed at 1

## An inadmissible anchor is an error that names the admissible ones, so the
## registry never has to be memorised.
err <- try(gkwq_parts("mc", anchor = "beta"), silent = TRUE)
cat(conditionMessage(attr(err, "condition")), "\n")
#> anchor ‘beta’ is not available for family ‘mc’; choose one of ‘lambda’. 
## anchor 'beta' is not available for family 'mc'; choose one of 'lambda'.

## The forced lambda-solve for "mc", in closed form. With alpha = beta = 1
## the quantile is z_tau^(1/lambda), so lambda = log z_tau / log mu.
tau <- 0.30; mu <- 0.60; g <- 2; d <- 1.5
lam <- log(qbeta(tau, g, d + 1)) / log(mu)
gkwq_quantile(tau, alpha = 1, beta = 1, gamma = g, delta = d, lambda = lam)
#> [1] 0.6
## [1] 0.6
## The tau-quantile is exactly the target mu: that is what anchoring buys.

## Why gamma is admissible for the "beta" family. Solving
## I_mu(gamma, delta + 1) = tau in gamma always has a root, because the left
## side sweeps all of (0,1) as gamma ranges over (0, Inf):
round(pbeta(0.4, c(0.01, 0.1, 0.5, 1, 2, 5, 20), 2), 4)
#> [1] 0.9968 0.9672 0.8222 0.6400 0.3520 0.0410 0.0000
## [1] 0.9968 0.9672 0.8222 0.6400 0.3520 0.0410 0.0000

## Why delta is not. Over the admissible range delta >= 0 the same expression
## spans only [mu^gamma, 1), so any tau below mu^gamma is unreachable:
round(pbeta(0.4, 1.5, c(0, 1, 5, 50) + 1), 4)
#> [1] 0.2530 0.4807 0.9049 1.0000
## [1] 0.2530 0.4807 0.9049 1.0000
round(0.4^1.5, 4)
#> [1] 0.253
## [1] 0.253    <- the floor, attained at delta = 0
```
