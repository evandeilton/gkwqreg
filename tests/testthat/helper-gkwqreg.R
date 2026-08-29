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

## A bare TMB object for one family, at a chosen anchor and level. Shared by
## the gradient tests and the anchor tests, so it lives here rather than in
## whichever test file happens to load first.
make_obj <- function(fam, d, tau = 0.5, par = NULL, anchor = NULL) {
  info <- gq(".gkwq_family_info"); md_f <- gq(".gkwq_model_data")
  val <- gq(".gkwq_validate_y"); lks <- gq(".gkwq_links")
  scl <- gq(".gkwq_link_scales"); tdat <- gq(".gkwq_tmb_data")
  tpar <- gq(".gkwq_tmb_params")
  spec <- info(fam, anchor)
  md <- md_f(y ~ x, d, spec$parts, NULL, stats::na.omit, NULL, NULL, NULL,
             fam, spec$anchor)
  yv <- val(md$y, 1e-10)
  lk <- lks(NULL, spec$parts); sc <- scl(NULL, lk, spec$parts)
  st <- lapply(spec$parts, function(p) numeric(ncol(md$X[[p]])))
  names(st) <- spec$parts
  TMB::MakeADFun(tdat(yv, md, spec, tau, lk, sc, gkwq_control()),
                 tpar(st, md, spec), DLL = "gkwqreg", silent = TRUE)
}
