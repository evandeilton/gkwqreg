## ---------------------------------------------------------------------------
## The family registry.
##
## This is the single source of truth for what a family is.  src/gkwqreg.cpp
## switches on nothing: the family reaches the tape entirely as data (par_id,
## fixed_val, z_mode, anchor_code, delta_is_zero), so adding a family costs an
## entry here rather than a branch in C++.
## ---------------------------------------------------------------------------

## Parameter identifiers shared with the template.
.GKWQ_PAR_ID <- c(alpha = 1L, beta = 2L, gamma = 3L, delta = 4L, lambda = 5L)

## Anchor identifiers shared with the template.  4 (delta) is deliberately
## absent: see .gkwq_family_spec() for why the delta-solve is inadmissible.
.GKWQ_ANCHOR_ID <- c(beta = 1L, alpha = 2L, lambda = 3L, gamma = 5L)

.GKWQ_FAMILIES <- c("kw", "ekw", "kkw", "bkw", "gkw", "mc", "beta")

.gkwq_family_spec <- function(family) {
  switch(family,
    kw = list(
      code = 1L, full = c("alpha", "beta"),
      fixed = list(gamma = 1, delta = 0, lambda = 1),
      anchors = c("beta", "alpha"), default_anchor = "beta"
    ),
    ekw = list(
      code = 2L, full = c("alpha", "beta", "lambda"),
      fixed = list(gamma = 1, delta = 0),
      anchors = c("beta", "alpha", "lambda"), default_anchor = "beta"
    ),
    kkw = list(
      code = 3L, full = c("alpha", "beta", "delta", "lambda"),
      fixed = list(gamma = 1),
      anchors = c("beta", "alpha", "lambda"), default_anchor = "beta"
    ),
    bkw = list(
      code = 4L, full = c("alpha", "beta", "gamma", "delta"),
      fixed = list(lambda = 1),
      anchors = c("beta", "alpha"), default_anchor = "beta"
    ),
    gkw = list(
      code = 5L, full = c("alpha", "beta", "gamma", "delta", "lambda"),
      fixed = list(),
      anchors = c("beta", "alpha", "lambda"), default_anchor = "beta"
    ),
    mc = list(
      ## alpha = beta = 1 leaves neither available, so the lambda-solve is
      ## forced.  It collapses to the exact form lambda = log z / log mu.
      code = 6L, full = c("gamma", "delta", "lambda"),
      fixed = list(alpha = 1, beta = 1),
      anchors = "lambda", default_anchor = "lambda"
    ),
    beta = list(
      ## The one family with no closed-form solve: gamma and delta enter only
      ## through z_tau = qbeta(tau, gamma, delta+1).  gamma is the only
      ## admissible anchor.  I_mu(gamma, delta+1) sweeps all of (0,1) as gamma
      ## ranges over (0, Inf), so a root always exists; the delta-solve instead
      ## requires tau > mu^gamma, which nothing guarantees.
      code = 7L, full = c("gamma", "delta"),
      fixed = list(alpha = 1, beta = 1, lambda = 1),
      anchors = "gamma", default_anchor = "gamma"
    ),
    stop("unknown family: ", sQuote(family), call. = FALSE)
  )
}

.gkwq_family_info <- function(family = "kw", anchor = NULL) {
  family <- match.arg(family, .GKWQ_FAMILIES)
  spec <- .gkwq_family_spec(family)
  spec$family <- family

  if (is.null(anchor)) {
    anchor <- spec$default_anchor
  } else {
    anchor <- as.character(anchor)[1L]
    if (!anchor %in% spec$anchors) {
      stop(sprintf(
        "anchor %s is not available for family %s; choose one of %s.",
        sQuote(anchor), sQuote(family),
        paste(sQuote(spec$anchors), collapse = ", ")
      ), call. = FALSE)
    }
  }
  spec$anchor <- anchor
  spec$anchor_code <- .GKWQ_ANCHOR_ID[[anchor]]

  ## THE part-order contract: the conditional quantile first, then the
  ## family's own order with the anchored parameter removed.
  spec$parts <- c("mu", setdiff(spec$full, anchor))
  spec$positions <- stats::setNames(seq_along(spec$parts), spec$parts)

  ## --- the data-driven family description handed to the tape ---------------
  ## par_id[j] says which parameter part position j+1 carries; 0 = unused.
  par_id <- integer(4L)
  tail_parts <- spec$parts[-1L]
  if (length(tail_parts)) {
    par_id[seq_along(tail_parts)] <- .GKWQ_PAR_ID[tail_parts]
  }
  spec$par_id <- par_id

  ## fixed_val holds the family's structural constants in (a,b,g,d,L) order.
  ## Slots for modelled or anchored parameters are overwritten in the tape.
  fixed_val <- c(alpha = 1, beta = 1, gamma = 1, delta = 0, lambda = 1)
  for (nm in names(spec$fixed)) fixed_val[[nm]] <- spec$fixed[[nm]]
  spec$fixed_val <- unname(fixed_val)

  gamma_fixed_one <- isTRUE(spec$fixed$gamma == 1)
  delta_fixed_zero <- isTRUE(spec$fixed$delta == 0)

  ## z_mode picks how log z_tau is obtained.  The elementary cases are not an
  ## optimisation: they keep the incomplete beta out of the tape entirely,
  ## which is what makes kw, ekw and kkw carry zero differentiation risk.
  spec$z_mode <- if (gamma_fixed_one && delta_fixed_zero) {
    1L
  } else if (gamma_fixed_one) {
    2L
  } else {
    3L
  }
  spec$delta_is_zero <- as.integer(delta_fixed_zero)
  spec$needs_qbeta <- spec$z_mode == 3L || anchor == "gamma"
  spec
}

#' The multi-part formula contract for a family and anchor
#'
#' Reports, in order, the names of the right-hand formula parts that
#' [gkwqreg()] will read for a given family and anchor. The contract is
#' exported as a function rather than merely described in prose, so that user
#' code, error messages and this manual cannot disagree about it.
#'
#' @details
#' A `gkwqreg` model formula carries one right-hand part per modelled
#' parameter, the parts separated by `|`, as in `y ~ x1 + x2 | z | 1`. Two
#' rules fix the meaning of every part.
#'
#' 1. **Part one is always `mu`**, the conditional \eqn{\tau}-quantile of the
#'    response, \eqn{\mu_\tau(x) = Q(\tau \mid x)}. This holds for every family
#'    and every anchor without exception. It is the parameter that the
#'    reparametrization creates, and the reason the model is a quantile
#'    regression rather than a mean regression.
#' 2. **The remaining parts follow the family's own parameter order, with the
#'    anchored parameter removed.** The anchor is the parameter eliminated in
#'    favour of \eqn{\mu}; having been solved for, it is no longer free, so it
#'    receives no formula part. See [gkwq_anchors()].
#'
#' The Kumaraswamy family, for instance, has parameters
#' \eqn{(\alpha, \beta)}. Anchoring on \eqn{\beta}, its default, leaves
#' `mu | alpha`; anchoring on \eqn{\alpha} leaves `mu | beta`. These are
#' different models, not different labels for one model, unless every remaining
#' part is regressed on the same covariates as `mu`.
#'
#' @section Omitted and excess parts:
#' Trailing parts that are omitted default to `~ 1`, an intercept only, so
#' `y ~ x` and `y ~ x | 1` describe the same Kumaraswamy model. Supplying
#' **more** parts than the model has is an error, and the message names the
#' parts the model expected. This is deliberate: silently discarding the
#' surplus would conceal a genuine mistake, most often a formula written for a
#' different anchor or a different family.
#'
#' @section The family parameter sets:
#' Each family holds some of the five Generalized Kumaraswamy parameters
#' \eqn{(\alpha, \beta, \gamma, \delta, \lambda)} at structural constants and
#' leaves the rest free. The free set, listed in the family's own order, is
#' what rule two refers to.
#'
#' | family | free parameters | held fixed |
#' |:--|:--|:--|
#' | `"kw"` | `alpha`, `beta` | \eqn{\gamma = 1}, \eqn{\delta = 0}, \eqn{\lambda = 1} |
#' | `"ekw"` | `alpha`, `beta`, `lambda` | \eqn{\gamma = 1}, \eqn{\delta = 0} |
#' | `"kkw"` | `alpha`, `beta`, `delta`, `lambda` | \eqn{\gamma = 1} |
#' | `"bkw"` | `alpha`, `beta`, `gamma`, `delta` | \eqn{\lambda = 1} |
#' | `"gkw"` | `alpha`, `beta`, `gamma`, `delta`, `lambda` | none |
#' | `"mc"` | `gamma`, `delta`, `lambda` | \eqn{\alpha = 1}, \eqn{\beta = 1} |
#' | `"beta"` | `gamma`, `delta` | \eqn{\alpha = 1}, \eqn{\beta = 1}, \eqn{\lambda = 1} |
#'
#' Removing the anchor from the free set and prefixing `mu` gives the parts.
#' Under the default anchor the number of parts equals the number of free
#' parameters, since one parameter is removed and `mu` is added.
#'
#' @param family Character, one of `"kw"`, `"ekw"`, `"kkw"`, `"bkw"`, `"gkw"`,
#'   `"mc"`, `"beta"`. Partial matching is applied. An unrecognised name is an
#'   error.
#' @param anchor Character or `NULL`. The parameter eliminated in favour of the
#'   conditional quantile. `NULL`, the default, selects the family default,
#'   which is the first element of [gkwq_anchors()]. An anchor the family does
#'   not admit raises an error naming those it does.
#'
#' @return A character vector naming the formula parts in order. Its length is
#'   one plus the number of free parameters remaining after the anchor is
#'   removed, and its first element is always `"mu"`.
#'
#' @seealso [gkwq_anchors()] for which anchors a family admits and why;
#'   [gkwqreg()] for the model these parts describe.
#'
#' @examples
#' ## Rule one: the first part is always the conditional quantile.
#' gkwq_parts("kw")
#' ## [1] "mu"    "alpha"
#'
#' ## Rule two: the anchored parameter is dropped from the family's own order.
#' gkwq_parts("kw", anchor = "alpha")
#' ## [1] "mu"   "beta"
#'
#' ## The five-parameter family under each anchor it admits.
#' for (a in gkwq_anchors("gkw")) {
#'   cat(sprintf("%-7s %s\n", a, paste(gkwq_parts("gkw", a), collapse = " | ")))
#' }
#' ## beta    mu | alpha | gamma | delta | lambda
#' ## alpha   mu | beta | gamma | delta | lambda
#' ## lambda  mu | alpha | beta | gamma | delta
#'
#' ## The whole contract at default anchors, at a glance.
#' for (f in c("kw", "ekw", "kkw", "bkw", "gkw", "mc", "beta")) {
#'   cat(sprintf("%-5s %s\n", f, paste(gkwq_parts(f), collapse = " | ")))
#' }
#' ## kw    mu | alpha
#' ## ekw   mu | alpha | lambda
#' ## kkw   mu | alpha | delta | lambda
#' ## bkw   mu | alpha | gamma | delta
#' ## gkw   mu | alpha | gamma | delta | lambda
#' ## mc    mu | gamma | delta
#' ## beta  mu | delta
#'
#' ## Trailing parts may be omitted and default to ~ 1.
#' set.seed(1)
#' dat <- data.frame(x = runif(50), y = runif(50, 0.2, 0.8))
#' f1 <- gkwqreg(y ~ x,     data = dat, tau = 0.5, family = "kw")
#' f2 <- gkwqreg(y ~ x | 1, data = dat, tau = 0.5, family = "kw")
#' all.equal(coef(f1), coef(f2))
#' ## [1] TRUE
#'
#' ## Supplying more parts than the model has is an error, not a silent drop,
#' ## and the message names the parts that were expected.
#' err <- try(gkwqreg(y ~ x | x | x, data = dat, tau = 0.5, family = "kw"),
#'            silent = TRUE)
#' cat(conditionMessage(attr(err, "condition")), "\n")
#' ## formula has 3 right-hand parts but family 'kw' with anchor 'beta' takes
#' ##   at most 2: mu | alpha.
#' ##   Parts are separated by `|`, in the order given by gkwq_parts("kw", "beta").
#' @export
gkwq_parts <- function(family = "kw", anchor = NULL) {
  .gkwq_family_info(family, anchor)$parts
}

#' Admissible anchors for a family
#'
#' The *anchor* is the parameter eliminated in favour of the conditional
#' quantile when a family is reparametrized at a level \eqn{\tau} fixed in
#' advance. `gkwq_anchors()` reports which of a family's parameters may play
#' that role, the first being the default that [gkwqreg()] uses when
#' `anchor = NULL`.
#'
#' @details
#' Fix a level \eqn{\tau \in (0,1)} and let \eqn{\mu} denote the target
#' conditional quantile. Anchoring means solving
#'
#' \deqn{Q(\tau; \alpha, \beta, \gamma, \delta, \lambda) = \mu}
#'
#' for one parameter and substituting the solution back into the likelihood.
#' What remains is a distribution indexed by \eqn{\mu} together with the
#' parameters left free, in which \eqn{\mu} *is* the conditional
#' \eqn{\tau}-quantile and can therefore be modelled directly through a link
#' function. The anchored parameter is no longer free and takes no formula
#' part; see [gkwq_parts()].
#'
#' A parameter is admissible as an anchor only when that equation has a unique
#' solution for it at every \eqn{(\mu, \tau)} in \eqn{(0,1)^2} and for every
#' admissible value of the parameters left free. Three of the five parameters
#' admit a closed-form solve. Writing \eqn{z_\tau} for
#' `qbeta(tau, gamma, delta + 1)`,
#'
#' \deqn{\beta = \frac{\log(1 - z_\tau^{1/\lambda})}{\log(1 - \mu^\alpha)},
#'   \qquad
#'   \alpha = \frac{\log\{1 - (1 - z_\tau^{1/\lambda})^{1/\beta}\}}{\log \mu},
#'   \qquad
#'   \lambda = \frac{\log z_\tau}{\log\{1 - (1 - \mu^\alpha)^\beta\}}.}
#'
#' Each is a ratio of two logarithms of quantities lying in \eqn{(0,1)} and is
#' consequently strictly positive and finite for every \eqn{\mu} and \eqn{\tau}
#' in \eqn{(0,1)}. No box constraint on the anchored parameter is ever needed,
#' which is what makes an unconstrained optimizer safe for this model.
#'
#' @section Which families admit which anchors:
#'
#' | family | admissible anchors (default first) | how the anchor is obtained |
#' |:--|:--|:--|
#' | `"kw"` | `"beta"`, `"alpha"` | closed form |
#' | `"ekw"` | `"beta"`, `"alpha"`, `"lambda"` | closed form |
#' | `"kkw"` | `"beta"`, `"alpha"`, `"lambda"` | closed form |
#' | `"bkw"` | `"beta"`, `"alpha"` | closed form |
#' | `"gkw"` | `"beta"`, `"alpha"`, `"lambda"` | closed form |
#' | `"mc"` | `"lambda"` | closed form, \eqn{\lambda = \log z_\tau / \log \mu} |
#' | `"beta"` | `"gamma"` | one-dimensional root find |
#'
#' A parameter that a family holds fixed as a structural constant cannot be an
#' anchor, since it is not free to be solved for, and that alone accounts for
#' most of the table. The `"mc"` family fixes \eqn{\alpha = \beta = 1}, so
#' neither is available and the \eqn{\lambda}-solve is forced; with those two
#' substitutions the general expression above collapses to the elementary form
#' \eqn{\lambda = \log z_\tau / \log \mu}. The `"beta"` family fixes
#' \eqn{\alpha = \beta = \lambda = 1}, leaving only \eqn{\gamma} and
#' \eqn{\delta} as candidates.
#'
#' @section Why gamma is admissible and delta never is:
#' In the `"beta"` family both free parameters enter the quantile solely
#' through \eqn{z_\tau}, so \eqn{Q(\tau) = \mu} reduces to
#'
#' \deqn{I_\mu(\gamma, \delta + 1) = \tau,}
#'
#' where \eqn{I} is the regularized incomplete beta function, that is `pbeta`.
#' There is no closed-form inverse, but admissibility is a question of existence
#' and uniqueness rather than of closed form.
#'
#' Read as an equation in \eqn{\gamma}, the left-hand side is continuous and
#' strictly decreasing, tending to \eqn{1} as \eqn{\gamma \to 0} and to \eqn{0}
#' as \eqn{\gamma \to \infty}. It therefore sweeps the whole of \eqn{(0,1)} and
#' a root exists, uniquely, for every \eqn{\tau \in (0,1)} and every \eqn{\mu}.
#' A bracketed one-dimensional root find is all that is required, and
#' \eqn{\gamma} is admissible.
#'
#' Read instead as an equation in \eqn{\delta}, the left-hand side is
#' increasing, but \eqn{\delta} is bounded below by zero and
#' \eqn{I_\mu(\gamma, 1) = \mu^\gamma}. Over the admissible range
#' \eqn{\delta \in [0, \infty)} the left-hand side spans only
#' \eqn{[\mu^\gamma, 1)}. A \eqn{\delta}-solve therefore exists only when
#' \eqn{\tau > \mu^\gamma}, a condition that nothing in the data or the model
#' guarantees: a low level combined with a large fitted quantile simply breaks
#' it. This is why \eqn{\gamma} is the only anchor the `"beta"` family admits,
#' and why \eqn{\delta} is an anchor in no family at all -- wherever it is free
#' it is bounded below by zero and enters only through \eqn{z_\tau}, so the same
#' range restriction recurs.
#'
#' @section The anchor is a modelling choice:
#' When every free parameter is regressed on the same covariates as \eqn{\mu},
#' the choice of anchor leaves the likelihood unchanged and is a genuine
#' reparametrization. It ceases to be innocuous as soon as \eqn{\mu} varies with
#' covariates while a free parameter is held constant, because the anchor then
#' decides which parameter absorbs that variation. Two anchors applied to the
#' same data give non-nested models of equal dimension: compare them by AIC, BIC
#' or [vuong_test()], never by a likelihood-ratio test. See [gkwqreg()] for the
#' full discussion.
#'
#' @inheritParams gkwq_parts
#'
#' @return A character vector of the anchor names admissible for `family`,
#'   ordered with the family default first. Supplying any other value as the
#'   `anchor` argument of [gkwqreg()] or [gkwq_parts()] raises an error that
#'   names these.
#'
#' @seealso [gkwq_parts()] for the formula parts that follow from a choice of
#'   anchor; [gkwq_quantile()] for the quantile function being solved;
#'   [gkwqreg()] for fitting.
#'
#' @examples
#' ## Every family, default anchor first.
#' for (f in c("kw", "ekw", "kkw", "bkw", "gkw", "mc", "beta")) {
#'   cat(sprintf("%-5s %s\n", f, paste(gkwq_anchors(f), collapse = ", ")))
#' }
#' ## kw    beta, alpha
#' ## ekw   beta, alpha, lambda
#' ## kkw   beta, alpha, lambda
#' ## bkw   beta, alpha
#' ## gkw   beta, alpha, lambda
#' ## mc    lambda            <- alpha and beta are fixed at 1
#' ## beta  gamma             <- alpha, beta and lambda are fixed at 1
#'
#' ## An inadmissible anchor is an error that names the admissible ones, so the
#' ## registry never has to be memorised.
#' err <- try(gkwq_parts("mc", anchor = "beta"), silent = TRUE)
#' cat(conditionMessage(attr(err, "condition")), "\n")
#' ## anchor 'beta' is not available for family 'mc'; choose one of 'lambda'.
#'
#' ## The forced lambda-solve for "mc", in closed form. With alpha = beta = 1
#' ## the quantile is z_tau^(1/lambda), so lambda = log z_tau / log mu.
#' tau <- 0.30; mu <- 0.60; g <- 2; d <- 1.5
#' lam <- log(qbeta(tau, g, d + 1)) / log(mu)
#' gkwq_quantile(tau, alpha = 1, beta = 1, gamma = g, delta = d, lambda = lam)
#' ## [1] 0.6
#' ## The tau-quantile is exactly the target mu: that is what anchoring buys.
#'
#' ## Why gamma is admissible for the "beta" family. Solving
#' ## I_mu(gamma, delta + 1) = tau in gamma always has a root, because the left
#' ## side sweeps all of (0,1) as gamma ranges over (0, Inf):
#' round(pbeta(0.4, c(0.01, 0.1, 0.5, 1, 2, 5, 20), 2), 4)
#' ## [1] 0.9968 0.9672 0.8222 0.6400 0.3520 0.0410 0.0000
#'
#' ## Why delta is not. Over the admissible range delta >= 0 the same expression
#' ## spans only [mu^gamma, 1), so any tau below mu^gamma is unreachable:
#' round(pbeta(0.4, 1.5, c(0, 1, 5, 50) + 1), 4)
#' ## [1] 0.2530 0.4807 0.9049 1.0000
#' round(0.4^1.5, 4)
#' ## [1] 0.253    <- the floor, attained at delta = 0
#' @export
gkwq_anchors <- function(family = "kw") {
  spec <- .gkwq_family_spec(match.arg(family, .GKWQ_FAMILIES))
  c(spec$default_anchor, setdiff(spec$anchors, spec$default_anchor))
}

## ---------------------------------------------------------------------------
## Identifiability guards (SPEC 1.6).
## ---------------------------------------------------------------------------

.gkwq_check_identifiability <- function(spec, warn = TRUE) {
  ## (b) delta = 0 confounds gamma and lambda.  With delta = 0 the incomplete
  ## beta collapses to I_u(gamma, 1) = u^gamma, so F(y) = v^(gamma*lambda) and
  ## only the PRODUCT is identified.  Anchoring on lambda then leaves the
  ## quantile coefficients themselves unidentified, which is an error rather
  ## than a warning.
  gamma_free <- "gamma" %in% spec$full
  delta_zero <- isTRUE(spec$fixed$delta == 0)
  if (spec$anchor == "lambda" && delta_zero && gamma_free) {
    stop(
      "anchor = \"lambda\" is not identified for family ", sQuote(spec$family),
      ": delta is fixed at 0 while gamma is free, so only the product ",
      "gamma * lambda is identified and the quantile coefficients are not. ",
      "Use anchor = \"beta\" or anchor = \"alpha\".",
      call. = FALSE
    )
  }

  ## (c) the full gkw model is weakly identified in ANY parametrization,
  ## with information-matrix condition numbers of order 1e8 to 1e11.
  if (warn && spec$family == "gkw") {
    warning(
      "family = \"gkw\" is weakly identified in any parametrization ",
      "(information-matrix condition numbers of order 1e8 to 1e11). ",
      "Consider a sub-family such as \"ekw\" or \"kkw\"; summary() reports the ",
      "condition number so you can check.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
