# Changelog

## gkwqreg 0.1.0

First release. Fixed-level parametric quantile regression for the seven
Generalized Kumaraswamy families.

### Features

- [`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md)
  fits the conditional `tau`-quantile directly, for `kw`, `ekw`, `kkw`,
  `bkw`, `gkw`, `mc` and `beta`. A vector `tau` returns a `gkwqregs`
  container of independent fits.
- `anchor` selects which parameter is eliminated in favour of the
  quantile. It is a documented modeling argument, not an internal
  detail: two anchors give non-nested models of equal dimension unless
  the nuisance parameters are saturated.
- Full S3 surface.
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) returns
  conditional quantiles, not means; the class does not inherit from
  `gkwreg`, so a missing method errors rather than applying mean
  semantics.
- Quantile-specific diagnostics: `"tau-sign"` and `"check"` residuals, a
  calibration panel in
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), and
  empirical coverage in
  [`summary()`](https://rdrr.io/r/base/summary.html).
- [`marginal_effects()`](https://evandeilton.github.io/gkwqreg/reference/marginal_effects.md),
  [`pinball()`](https://evandeilton.github.io/gkwqreg/reference/pinball.md),
  [`check_crossing()`](https://evandeilton.github.io/gkwqreg/reference/check_crossing.md),
  [`rearrange()`](https://evandeilton.github.io/gkwqreg/reference/rearrange.md),
  [`quantile_process()`](https://evandeilton.github.io/gkwqreg/reference/quantile_process.md),
  [`vuong_test()`](https://evandeilton.github.io/gkwqreg/reference/vuong_test.md),
  [`compare_families()`](https://evandeilton.github.io/gkwqreg/reference/compare_families.md),
  [`gkwq_boot()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_boot.md),
  [`simulate()`](https://rdrr.io/r/stats/simulate.html),
  [`lrtest()`](https://evandeilton.github.io/gkwqreg/reference/lrtest.md).
- Wald, likelihood-ratio, profile, bootstrap and sandwich inference.
  `estfun()` and `bread()` are registered when `sandwich` is installed,
  so
  [`sandwich::sandwich()`](https://zeileis.codeberg.page/sandwich/reference/sandwich.html),
  [`sandwich::vcovCL()`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html)
  and
  [`lmtest::coeftest()`](https://rdrr.io/pkg/lmtest/man/coeftest.html)
  work directly.
  [`sandwich::vcovHC()`](https://zeileis.codeberg.page/sandwich/reference/vcovHC.html)
  is not supported: it needs working residuals and a single design
  matrix aligned with the scores, which a multi-part model does not
  have. Use `vcov(type = "sandwich")` instead.

### Numerical notes

- Ships its own incomplete-beta atomic (`inst/include/gkwq_atomic.hpp`).
  TMB’s `qbeta` reports shape derivatives of exactly zero at
  whole-number shapes, because `tiny_ad/beta/toms708.cpp:341` assigns a
  literal to an AD variable; an optimizer reading that gradient leaves
  `delta` at its starting value and reports convergence. See
  [`vignette("gkwqreg-design")`](https://evandeilton.github.io/gkwqreg/articles/gkwqreg-design.md),
  which demonstrates both the failure and the fix.
- One TMB template for all seven families, compiled at install time. No
  C++ toolchain is needed at run time and no template is compiled during
  `R CMD check`.
- The log-likelihood is a log-domain cascade, never routed through
  [`gkwdist::dgkw`](https://evandeilton.github.io/gkwdist/reference/dgkw.html),
  which returns `-Inf` for responses within about `1e-10` of 1 where the
  true log-density is finite.
- No branch on `asDouble()` of anything parameter-dependent. Such a
  branch is frozen into the tape at the starting values, so the gradient
  the optimizer later reads belongs to a different function. Test `V9`
  checks the gradient far from the taping point.
- Observed information from `optimHess` on the AD gradient, never
  `obj$he()` and never the naive `J'HJ` sandwich of the unreparametrized
  Hessian.
- Weights and offsets are applied inside the likelihood.
