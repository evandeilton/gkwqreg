## The registry is the single source of truth the C++ template reads. These
## tests pin the contract that both sides depend on.

test_that("V11: family nesting holds -- kw is ekw with lambda at 1", {
  d <- sim_kw(n = 300)
  f_kw <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "kw")
  f_ekw <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "ekw")
  ## ekw nests kw, so it cannot fit worse, and on kw data it should not fit
  ## much better either.
  expect_gte(f_ekw$loglik, f_kw$loglik - 1e-6)
  expect_lt(f_ekw$loglik - f_kw$loglik, 5)
  a <- anova(f_kw, f_ekw)
  expect_s3_class(a, "anova.gkwqreg")
  expect_equal(a$`Chi Df`[2], 1)
})

test_that("gkwq_parts() puts the quantile first and drops the anchor", {
  expect_equal(gkwq_parts("kw"), c("mu", "alpha"))
  expect_equal(gkwq_parts("kw", "alpha"), c("mu", "beta"))
  expect_equal(gkwq_parts("ekw"), c("mu", "alpha", "lambda"))
  expect_equal(gkwq_parts("kkw"), c("mu", "alpha", "delta", "lambda"))
  expect_equal(gkwq_parts("bkw"), c("mu", "alpha", "gamma", "delta"))
  expect_equal(gkwq_parts("gkw"), c("mu", "alpha", "gamma", "delta", "lambda"))
  expect_equal(gkwq_parts("mc"), c("mu", "gamma", "delta"))
  expect_equal(gkwq_parts("beta"), c("mu", "delta"))
  for (fam in ALL_FAMILIES) expect_identical(gkwq_parts(fam)[1], "mu")
})

test_that("the registry and the template agree on the family description", {
  info <- gq(".gkwq_family_info")
  for (fam in ALL_FAMILIES) {
    for (anc in gkwq_anchors(fam)) {
      s <- info(fam, anc)
      ## par_id names exactly the non-quantile parts, in order.
      tail_parts <- s$parts[-1L]
      expect_equal(sum(s$par_id > 0), length(tail_parts),
                   info = sprintf("%s/%s", fam, anc))
      expect_length(s$par_id, 4L)
      expect_length(s$fixed_val, 5L)
      ## z_mode must match what the family fixes.
      expect_equal(s$z_mode,
                   if (identical(s$fixed$gamma, 1) && identical(s$fixed$delta, 0)) 1L
                   else if (identical(s$fixed$gamma, 1)) 2L else 3L,
                   info = sprintf("%s/%s z_mode", fam, anc))
    }
  }
})

test_that("supplying more formula parts than the model has is an error", {
  d <- sim_kw(n = 100)
  expect_error(gkwqreg(y ~ x | x | x, data = d, family = "kw"),
               "right-hand parts")
})

test_that("an unavailable anchor is refused by name", {
  expect_error(gkwq_parts("kw", "lambda"), "not available")
  expect_error(gkwq_parts("mc", "beta"), "not available")
  expect_error(gkwq_parts("beta", "delta"), "not available")
})

test_that("padding a formula with `~ 1` parts keeps its environment", {
  ## `y ~ x` is padded to `y ~ x | 1` for a two-part family. If the rebuilt
  ## formula lost its environment, a variable living in the caller rather than
  ## in `data` would resolve in one form and not the other.
  d <- sim_kw(n = 120)
  outside <- d$x                      # deliberately NOT a column of `dd`
  dd <- data.frame(y = d$y)
  f1 <- gkwqreg(y ~ outside, data = dd, tau = 0.5, family = "kw")
  f2 <- gkwqreg(y ~ outside | 1, data = dd, tau = 0.5, family = "kw")
  expect_equal(coef(f1), coef(f2), tolerance = 1e-8)
  expect_equal(f1$loglik, f2$loglik, tolerance = 1e-8)
})

test_that("omitted trailing parts default to an intercept", {
  d <- sim_kw(n = 150)
  a <- gkwqreg(y ~ x, data = d, tau = 0.5, family = "ekw")
  b <- gkwqreg(y ~ x | 1 | 1, data = d, tau = 0.5, family = "ekw")
  expect_equal(coef(a), coef(b), tolerance = 1e-6)
})

test_that("the whole nesting lattice holds: a bigger family never fits worse", {
  ## kw subset ekw subset kkw subset gkw, and kw subset bkw subset gkw. Each
  ## containment forces logLik to be non-decreasing, so a violation means the
  ## optimizer failed rather than that the model is worse.
  ##
  ## This is a real regression test: an intercept-only pre-fit used to run off
  ## along kkw's flat ridge to alpha ~ 1e-12 and lambda ~ 1e13, and the full fit
  ## inherited that start and landed 30 log-likelihood units below ekw.
  d <- sim_kw(n = 300, seed = 7)
  ll <- vapply(c("kw", "ekw", "kkw", "bkw", "gkw"), function(fam) {
    suppressWarnings(gkwqreg(y ~ x, data = d, tau = 0.5, family = fam))$loglik
  }, numeric(1))

  tol <- 1e-4                      # optimizer slack, not modelling slack
  expect_gte(ll[["ekw"]], ll[["kw"]] - tol)
  expect_gte(ll[["kkw"]], ll[["ekw"]] - tol)
  expect_gte(ll[["gkw"]], ll[["kkw"]] - tol)
  expect_gte(ll[["bkw"]], ll[["kw"]] - tol)
  expect_gte(ll[["gkw"]], ll[["bkw"]] - tol)
})

test_that("a bad pre-fit is discarded rather than inherited", {
  ## The starting value is chosen by the objective it actually produces, so a
  ## pre-fit that lands on a ridge cannot make the final fit worse than the
  ## plain start would have.
  d <- sim_kw(n = 250, seed = 7)
  auto <- suppressWarnings(gkwqreg(y ~ x, data = d, tau = 0.5, family = "kkw"))
  ols  <- suppressWarnings(gkwqreg(y ~ x, data = d, tau = 0.5, family = "kkw",
                                   control = gkwq_control(start_method = "ols")))
  expect_gte(auto$loglik, ols$loglik - 1e-4)
})

test_that("anova() refuses to report a test the lattice does not support", {
  ## The seven families form a lattice, not a chain: ekw and bkw both contain
  ## kw and neither contains the other, and 9 of the 21 pairs are unordered
  ## like that. anova() used to run the likelihood-ratio test on them anyway,
  ## with no warning and an ordinary-looking p-value.
  set.seed(41); n <- 200; x <- stats::runif(n)
  d <- data.frame(y = stats::qbeta(stats::runif(n), 2 + x, 3), x = x)
  fit <- function(f) suppressWarnings(gkwqreg(y ~ x, d, tau = 0.5, family = f))

  f_kw <- fit("kw"); f_ekw <- fit("ekw"); f_bkw <- fit("bkw")

  ## A genuine chain still reports a test.
  ok <- anova(f_kw, f_ekw)
  expect_true(is.finite(ok$`Pr(>Chisq)`[2L]))

  ## ekw and bkw are not ordered by containment.
  expect_warning(bad <- anova(f_ekw, f_bkw), "not nested")
  expect_true(is.na(bad$`Pr(>Chisq)`[2L]))
  expect_true(is.finite(bad$Chisq[2L]))   # the statistic is still shown

  ## The relation is read off the registry, so it cannot drift from the specs.
  nested <- gq(".gkwq_family_nested")
  expect_true(nested("kw", "ekw"));  expect_true(nested("kw", "gkw"))
  expect_true(nested("beta", "mc")); expect_true(nested("mc", "gkw"))
  expect_false(nested("ekw", "bkw")); expect_false(nested("bkw", "ekw"))
  expect_false(nested("ekw", "kw"))   # containment is strict and directed
})

test_that("a likelihood that falls as Df rises is not reported as p = 1", {
  ## pchisq() of a negative statistic returns 1, so a larger model that failed
  ## to converge printed as a comfortable "keep the smaller one". The fits are
  ## constructed by hand rather than hunted for, so the case is the same on
  ## every platform.
  set.seed(41); n <- 150; x <- stats::runif(n)
  d <- data.frame(y = stats::qbeta(stats::runif(n), 2 + x, 3), x = x)
  a <- suppressWarnings(gkwqreg(y ~ x, d, tau = 0.5, family = "kw"))
  b <- suppressWarnings(gkwqreg(y ~ x, d, tau = 0.5, family = "ekw"))
  b$loglik <- a$loglik - 5           # the larger model, fitting worse

  expect_warning(out <- anova(a, b), "falls as the dimension rises")
  expect_lt(out$Chisq[2L], 0)
  expect_true(is.na(out$`Pr(>Chisq)`[2L]))
})

test_that("a link is checked against the range of the parameter it carries", {
  ## .GKWQ_POS_LINKS was defined for the nuisance parts and never consulted, so
  ## a unit link on a shape passed in silence: link = c(alpha = "logit")
  ## confines alpha to (0,1) and the fit converges onto the boundary at 1,
  ## reporting nothing. Both ends of the rule are checked here.
  set.seed(9); n <- 80; x <- stats::runif(n)
  d <- data.frame(y = stats::qbeta(stats::runif(n), 2 + x, 3), x = x)

  for (lk in c("logit", "probit", "cauchy", "cloglog")) {
    expect_error(gkwqreg(y ~ x, d, tau = 0.5, family = "kw",
                         link = stats::setNames(lk, "alpha")),
                 "must map to \\(0, Inf\\)", info = lk)
  }
  for (lk in c("log", "sqrt", "identity")) {
    expect_no_error(suppressWarnings(
      gkwqreg(y ~ x, d, tau = 0.5, family = "kw",
              link = stats::setNames(lk, "alpha"))))
  }

  ## The quantile keeps its own rule, in the other direction.
  expect_error(gkwqreg(y ~ x, d, tau = 0.5, family = "kw", link = c(mu = "log")),
               "must map to \\(0,1\\)")
})
