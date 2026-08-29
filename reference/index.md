# Package index

## Fitting the model

The entry point, and the pieces that configure a fit.

- [`gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg.md)
  : Parametric quantile regression for the Generalized Kumaraswamy
  family
- [`gkwq_control()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_control.md)
  : Tuning parameters for a Generalized Kumaraswamy quantile regression
- [`gkwq_parts()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_parts.md)
  : The multi-part formula contract for a family and anchor
- [`gkwq_anchors()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_anchors.md)
  : Admissible anchors for a family
- [`gkwq_quantile()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_quantile.md)
  : The Generalized Kumaraswamy quantile function, in closed form

## Extractors and summaries

The standard S3 surface, for a single fit and for a set of fits over
several quantile levels.

- [`coef(`*`<gkwqreg>`*`)`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg-extractors.md)
  [`logLik(`*`<gkwqreg>`*`)`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg-extractors.md)
  [`nobs(`*`<gkwqreg>`*`)`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg-extractors.md)
  [`AIC(`*`<gkwqreg>`*`)`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg-extractors.md)
  [`BIC(`*`<gkwqreg>`*`)`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg-extractors.md)
  [`family(`*`<gkwqreg>`*`)`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg-extractors.md)
  [`formula(`*`<gkwqreg>`*`)`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg-extractors.md)
  [`terms(`*`<gkwqreg>`*`)`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg-extractors.md)
  [`model.frame(`*`<gkwqreg>`*`)`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg-extractors.md)
  [`model.matrix(`*`<gkwqreg>`*`)`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg-extractors.md)
  [`getCall(`*`<gkwqreg>`*`)`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg-extractors.md)
  [`update(`*`<gkwqreg>`*`)`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg-extractors.md)
  [`summary(`*`<gkwqreg>`*`)`](https://evandeilton.github.io/gkwqreg/reference/gkwqreg-extractors.md)
  : Extractor methods for a quantile regression fit
- [`coef(`*`<gkwqregs>`*`)`](https://evandeilton.github.io/gkwqreg/reference/gkwqregs-methods.md)
  [`fitted(`*`<gkwqregs>`*`)`](https://evandeilton.github.io/gkwqreg/reference/gkwqregs-methods.md)
  [`logLik(`*`<gkwqregs>`*`)`](https://evandeilton.github.io/gkwqreg/reference/gkwqregs-methods.md)
  [`summary(`*`<gkwqregs>`*`)`](https://evandeilton.github.io/gkwqreg/reference/gkwqregs-methods.md)
  [`predict(`*`<gkwqregs>`*`)`](https://evandeilton.github.io/gkwqreg/reference/gkwqregs-methods.md)
  [`residuals(`*`<gkwqregs>`*`)`](https://evandeilton.github.io/gkwqreg/reference/gkwqregs-methods.md)
  : Methods for a set of quantile regression fits
- [`fitted(`*`<gkwqreg>`*`)`](https://evandeilton.github.io/gkwqreg/reference/fitted.gkwqreg.md)
  : Fitted conditional quantiles

## Prediction and marginal effects

- [`predict(`*`<gkwqreg>`*`)`](https://evandeilton.github.io/gkwqreg/reference/predict.gkwqreg.md)
  : Predictions from a fixed-level quantile regression fit
- [`marginal_effects()`](https://evandeilton.github.io/gkwqreg/reference/marginal_effects.md)
  : Marginal effects on the quantile scale

## Inference

Wald, likelihood-ratio, profile, bootstrap and sandwich routes, plus the
tests for comparing fits.

- [`vcov(`*`<gkwqreg>`*`)`](https://evandeilton.github.io/gkwqreg/reference/vcov.gkwqreg.md)
  : Covariance matrix of a quantile regression fit
- [`confint(`*`<gkwqreg>`*`)`](https://evandeilton.github.io/gkwqreg/reference/confint.gkwqreg.md)
  : Confidence intervals for a quantile regression fit
- [`anova(`*`<gkwqreg>`*`)`](https://evandeilton.github.io/gkwqreg/reference/anova.gkwqreg.md)
  : Likelihood-ratio comparison of nested quantile regression fits
- [`lrtest()`](https://evandeilton.github.io/gkwqreg/reference/lrtest.md)
  : Likelihood-ratio test for nested quantile regression fits
- [`vuong_test()`](https://evandeilton.github.io/gkwqreg/reference/vuong_test.md)
  : Vuong test for two non-nested quantile regression fits
- [`compare_families()`](https://evandeilton.github.io/gkwqreg/reference/compare_families.md)
  : Compare the seven families at one quantile level
- [`gkwq_boot()`](https://evandeilton.github.io/gkwqreg/reference/gkwq_boot.md)
  : Bootstrap a quantile regression fit
- [`estfun.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/estfun.gkwqreg.md)
  : Per-observation score contributions
- [`bread.gkwqreg()`](https://evandeilton.github.io/gkwqreg/reference/bread.gkwqreg.md)
  : Bread matrix for the sandwich covariance

## Residuals and diagnostics

- [`residuals(`*`<gkwqreg>`*`)`](https://evandeilton.github.io/gkwqreg/reference/residuals.gkwqreg.md)
  : Residuals from a fitted quantile regression
- [`pinball()`](https://evandeilton.github.io/gkwqreg/reference/pinball.md)
  : Pinball (check) loss of a fitted quantile regression

## The quantile process

Fitting over a grid of quantile levels, detecting crossing and repairing
it.

- [`quantile_process()`](https://evandeilton.github.io/gkwqreg/reference/quantile_process.md)
  : The estimated quantile process
- [`check_crossing()`](https://evandeilton.github.io/gkwqreg/reference/check_crossing.md)
  [`qcrossing()`](https://evandeilton.github.io/gkwqreg/reference/check_crossing.md)
  : Detect quantile crossing
- [`rearrange()`](https://evandeilton.github.io/gkwqreg/reference/rearrange.md)
  : Monotone rearrangement of crossing quantile curves

## Plots

- [`plot(`*`<gkwqreg>`*`)`](https://evandeilton.github.io/gkwqreg/reference/plot.gkwqreg.md)
  : Diagnostic plots for a quantile regression fit
- [`plot(`*`<gkwq_process>`*`)`](https://evandeilton.github.io/gkwqreg/reference/plot.gkwq_process.md)
  [`plot(`*`<gkwqregs>`*`)`](https://evandeilton.github.io/gkwqreg/reference/plot.gkwq_process.md)
  : Plot a quantile process

## Simulation

- [`simulate(`*`<gkwqreg>`*`)`](https://evandeilton.github.io/gkwqreg/reference/simulate.gkwqreg.md)
  : Simulate responses from a fitted quantile regression
