## Shared fixtures. The design study's own pilot: Kumaraswamy, beta anchor,
## logit(mu) = 0.4 + 1.1 x, alpha = 2, so the truth is known exactly.
sim_kw <- function(n = 300, tau = 0.5, b0 = 0.4, b1 = 1.1, alpha = 2,
                   seed = 42) {
  set.seed(seed)
  x <- stats::runif(n, -2, 2)
  mu <- stats::plogis(b0 + b1 * x)
  y <- gkwdist::rkw(n, alpha = alpha, beta = log1p(-tau) / log1p(-mu^alpha))
  data.frame(y = y, x = x)
}

## Reach the package internals the tests are allowed to check against.
gq <- function(nm) get(nm, envir = asNamespace("gkwqreg"))

ALL_FAMILIES <- c("kw", "ekw", "kkw", "bkw", "gkw", "mc", "beta")
