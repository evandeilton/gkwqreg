# gkwqreg 0.1.0

First release. Fixed-level parametric quantile regression for the seven Generalized
Kumaraswamy families.

## Features

* `gkwqreg()` fits the conditional `tau`-quantile directly, for `kw`, `ekw`, `kkw`,
  `bkw`, `gkw`, `mc` and `beta`. A vector `tau` returns a `gkwqregs` container of
  independent fits.
* `anchor` selects which parameter is eliminated in favour of the quantile. It is a
  documented modeling argument, not an internal detail: two anchors give non-nested
  models of equal dimension unless the nuisance parameters are saturated.
* Full S3 surface. `fitted()` returns conditional quantiles, not means; the class does
  not inherit from `gkwreg`, so a missing method errors rather than applying mean
  semantics.
* Quantile-specific diagnostics: `"tau-sign"` and `"check"` residuals, a calibration
  panel in `plot()`, and empirical coverage in `summary()`.
* `marginal_effects()`, `pinball()`, `check_crossing()`, `rearrange()`,
  `quantile_process()`, `vuong_test()`, `compare_families()`, `gkwq_boot()`,
  `simulate()`, `lrtest()`.
* Wald, likelihood-ratio, profile, bootstrap and sandwich inference. `estfun()` and `bread()` are
  registered when `sandwich` is installed, so `sandwich::vcovHC()` and
  `lmtest::coeftest()` work directly.

## Numerical notes

* Ships its own incomplete-beta atomic (`inst/include/gkwq_atomic.hpp`). TMB's `qbeta`
  reports shape derivatives of exactly zero at whole-number shapes, because
  `tiny_ad/beta/toms708.cpp:341` assigns a literal to an AD variable; an optimizer reading
  that gradient leaves `delta` at its starting value and reports convergence. See
  `vignette("gkwqreg-design")`, which demonstrates both the failure and the fix.
* One TMB template for all seven families, compiled at install time. No C++ toolchain is
  needed at run time and no template is compiled during `R CMD check`.
* The log-likelihood is a log-domain cascade, never routed through `gkwdist::dgkw`, which
  returns `-Inf` for responses within about `1e-10` of 1 where the true log-density is
  finite.
* No branch on `asDouble()` of anything parameter-dependent. Such a branch is frozen into
  the tape at the starting values, so the gradient the optimizer later reads belongs to a
  different function. Test `V9` checks the gradient far from the taping point.
* Observed information from `optimHess` on the AD gradient, never `obj$he()` and never the
  naive `J'HJ` sandwich of the unreparametrized Hessian.
* Weights and offsets are applied inside the likelihood.
