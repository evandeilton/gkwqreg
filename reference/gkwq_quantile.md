# The Generalized Kumaraswamy quantile function, in closed form

Evaluates \\Q(\tau)\\, the \\\tau\\-quantile of the five-parameter
Generalized Kumaraswamy distribution on the unit interval. The
distribution function inverts analytically, so no root finding is
involved. That closed form is what makes the fixed-level
reparametrization behind
[`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md)
possible at all: because \\Q(\tau)\\ can be written down, one parameter
can be solved for so that the conditional quantile becomes a parameter
itself.

## Usage

``` r
gkwq_quantile(tau, alpha = 1, beta = 1, gamma = 1, delta = 0, lambda = 1)
```

## Arguments

- tau:

  Numeric vector of quantile levels, strictly inside \\(0,1)\\. Values
  outside \\\[0,1\]\\ produce `NaN`.

- alpha, beta:

  Strictly positive shape parameters of the Kumaraswamy kernel.
  \\\alpha\\ acts on the left tail through \\y^\alpha\\ and \\\beta\\ on
  the right tail through \\(1 - y^\alpha)^\beta\\.

- gamma, delta:

  Shape parameters of the beta generator, which enter only through
  \\z\_\tau\\: \\\gamma\\ strictly positive, \\\delta\\ non-negative.
  Their defaults, \\\gamma = 1\\ and \\\delta = 0\\, remove the
  generator.

- lambda:

  Strictly positive exponentiation parameter, acting on \\v^\lambda\\.
  Its default \\\lambda = 1\\ removes it.

## Value

A numeric vector of quantiles in \\(0,1)\\, as long as the longest
argument; all arguments are recycled to that length under the usual R
rules. Inadmissible or missing inputs propagate as `NaN` or `NA`.

## Details

The Generalized Kumaraswamy distribution function at \\y \in (0,1)\\ is

\$\$F(y) = I\_{v^\lambda}(\gamma, \delta + 1), \qquad v = 1 - (1 -
y^\alpha)^\beta,\$\$

where \\I\\ denotes the regularized incomplete beta function, that is
`pbeta`. Inverting this expression one layer at a time, and writing
\\z\_\tau = I^{-1}\_\tau(\gamma, \delta + 1)\\, which in R is
`qbeta(tau, gamma, delta + 1)`, gives

\$\$Q(\tau) = \left\[\\ 1 - \left\\ 1 - z\_\tau^{1/\lambda}
\right\\^{1/\beta} \\\right\]^{1/\alpha}.\$\$

Two special cases of \\z\_\tau\\ are worth knowing, because the
sub-families that enjoy them never touch the incomplete beta function at
all. When \\\gamma = 1\\ and \\\delta = 0\\ the generator is the
identity and \\z\_\tau = \tau\\; this covers the `"kw"` and `"ekw"`
families. When \\\gamma = 1\\ with \\\delta\\ free, \\z\_\tau = 1 - (1 -
\tau)^{1/(\delta + 1)}\\; this covers `"kkw"`.

## Numerical behaviour

The expression above is a cascade of two nested complements, each of the
form \\1 - w\\ with \\w\\ close to one in the regions that matter. It is
therefore evaluated entirely in the log domain, using a two-branch
stable evaluation of \\\log(1 - e^{x})\\ for \\x \< 0\\ (Maechler 2012),
which switches between `log(-expm1(x))` and `log1p(-exp(x))` according
to the magnitude of \\x\\. Naive evaluation loses precision precisely
where a bounded-response model needs it most: when \\\beta\\ is large
the inner complement approaches one, a regime the default `beta` anchor
visits routinely.

Checked against
[`gkwdist::qgkw`](https://evandeilton.github.io/gkwdist/reference/qgkw.html)
over a 1701-point parameter grid – seven levels of \\\tau\\ spanning
0.01 to 0.99, with each of the five parameters taking three values – the
largest absolute discrepancy is `2.7e-14`.

## See also

[`gkwq_anchors()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_anchors.md)
for the solves that invert this identity in each parameter in turn;
[`gkwq_parts()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_parts.md)
for the family parameter sets;
[`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md)
for the regression model built on it.

## Examples

``` r
## The median of a Kumaraswamy(2, 3) variate. With gamma = 1, delta = 0 and
## lambda = 1 the generator drops out and Q collapses to the Kumaraswamy
## quantile [1 - (1 - tau)^(1/beta)]^(1/alpha).
gkwq_quantile(0.5, alpha = 2, beta = 3)
#> [1] 0.454202
## [1] 0.454202
all.equal(gkwq_quantile(0.4, alpha = 2, beta = 3),
          (1 - (1 - 0.4)^(1 / 3))^(1 / 2))
#> [1] TRUE
## [1] TRUE

## Q is increasing in tau, as any quantile function must be.
round(gkwq_quantile(c(0.05, 0.25, 0.5, 0.75, 0.95), alpha = 2, beta = 3), 4)
#> [1] 0.1302 0.3024 0.4542 0.6083 0.7947
## [1] 0.1302 0.3024 0.4542 0.6083 0.7947

## Arguments recycle: one level, three shapes.
round(gkwq_quantile(0.5, alpha = c(1, 2, 4), beta = 3), 4)
#> [1] 0.2063 0.4542 0.6739
## [1] 0.2063 0.4542 0.6739

## Q genuinely inverts F, to machine precision.
y <- c(0.05, 0.30, 0.62, 0.95)
p <- gkwdist::pgkw(y, 2, 3, 1.5, 0.5, 0.8)
max(abs(gkwq_quantile(p, 2, 3, 1.5, 0.5, 0.8) - y))
#> [1] 2.831069e-14
## [1] 2.831069e-14

## Agreement with the reference implementation across the parameter space:
## 1701 combinations of level and the five parameters.
g <- expand.grid(t = c(.01, .1, .25, .5, .75, .9, .99), a = c(.5, 1, 2.5),
                 b = c(.7, 1.5, 3), gm = c(.8, 1, 2), d = c(0, .5, 2),
                 L = c(.6, 1, 2))
max(abs(gkwq_quantile(g$t, g$a, g$b, g$gm, g$d, g$L) -
        gkwdist::qgkw(g$t, g$a, g$b, g$gm, g$d, g$L)))
#> [1] 2.738261e-14
## [1] 2.738261e-14

## The identity that anchoring inverts. Fix a level and a target quantile,
## then choose beta so that Q(tau) hits the target exactly. This is the
## "beta" anchor of the "kw" family, the default used by gkwqreg().
tau <- 0.5; mu <- 0.7; alpha <- 2
beta <- log1p(-tau) / log1p(-mu^alpha)
beta
#> [1] 1.029409
## [1] 1.029409
gkwq_quantile(tau, alpha = alpha, beta = beta)
#> [1] 0.7
## [1] 0.7   -- the median is the target mu, exactly
```
