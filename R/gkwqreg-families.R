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

#' Formula parts for a family and anchor
#'
#' The multi-part formula contract, exported so that users, error messages and
#' documentation cannot disagree about it. Part one is always `mu`, the
#' conditional quantile; the remaining parts follow the family's own parameter
#' order with the anchored parameter removed.
#'
#' @param family One of `"kw"`, `"ekw"`, `"kkw"`, `"bkw"`, `"gkw"`, `"mc"`,
#'   `"beta"`.
#' @param anchor The parameter eliminated in favour of the conditional quantile,
#'   or `NULL` for the family default. See [gkwqreg()].
#'
#' @return A character vector naming the formula parts in order.
#'
#' @examples
#' gkwq_parts("kw")            # "mu" "alpha"
#' gkwq_parts("kw", "alpha")   # "mu" "beta"
#' gkwq_parts("gkw")           # "mu" "alpha" "gamma" "delta" "lambda"
#' @export
gkwq_parts <- function(family = "kw", anchor = NULL) {
  .gkwq_family_info(family, anchor)$parts
}

#' Available anchors for a family
#'
#' @inheritParams gkwq_parts
#' @return A character vector of admissible anchors, the first being the default.
#' @examples
#' gkwq_anchors("gkw")
#' gkwq_anchors("mc")
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
