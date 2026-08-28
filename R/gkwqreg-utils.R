## ---------------------------------------------------------------------------
## Links, formulas, model frames.
##
## Reimplemented rather than borrowed from gkwreg: every equivalent there is a
## dot-prefixed internal, so reuse would mean `gkwreg:::` (a check NOTE CRAN
## rejects for new submissions) or growing that package's API to serve a
## sibling.
## ---------------------------------------------------------------------------

## Link codes match gkwreg's .convert_links_to_int, so link names stay portable
## between the two packages.
.GKWQ_LINKS <- c(log = 1L, logit = 2L, probit = 3L, cauchy = 4L, cloglog = 5L,
                 identity = 6L, sqrt = 7L, inverse = 8L, `inverse-square` = 9L)

## Links mapping to (0,1): admissible for the conditional quantile.
.GKWQ_UNIT_LINKS <- c("logit", "probit", "cauchy", "cloglog")
## Links mapping to (0, Inf): admissible for the remaining parameters.
.GKWQ_POS_LINKS <- c("log", "identity", "sqrt", "inverse", "inverse-square")

.gkwq_links <- function(link, parts) {
  out <- stats::setNames(
    c("logit", rep("log", length(parts) - 1L)), parts
  )
  if (!is.null(link)) {
    if (is.null(names(link))) {
      if (length(link) == 1L) {
        ## A single unnamed link applies to the nuisance parts only: silently
        ## overriding the quantile's link is never what a user means.
        out[parts[-1L]] <- link
      } else if (length(link) == length(parts)) {
        out[] <- link
      } else {
        stop("`link` must be named, or of length 1 or ", length(parts),
             " (one per part: ", paste(parts, collapse = ", "), ").",
             call. = FALSE)
      }
    } else {
      unknown <- setdiff(names(link), parts)
      if (length(unknown)) {
        stop("`link` names unknown parts: ", paste(sQuote(unknown), collapse = ", "),
             ". This model has parts ", paste(sQuote(parts), collapse = ", "), ".",
             call. = FALSE)
      }
      out[names(link)] <- unlist(link)
    }
  }
  bad <- setdiff(out, names(.GKWQ_LINKS))
  if (length(bad)) {
    stop("unknown link(s): ", paste(sQuote(bad), collapse = ", "),
         ". Available: ", paste(names(.GKWQ_LINKS), collapse = ", "), ".",
         call. = FALSE)
  }
  ## The conditional quantile lives in (0,1) and its link must respect that.
  if (!out[["mu"]] %in% .GKWQ_UNIT_LINKS) {
    stop("link for the conditional quantile must map to (0,1); got ",
         sQuote(out[["mu"]]), ". Use one of ",
         paste(sQuote(.GKWQ_UNIT_LINKS), collapse = ", "), ".", call. = FALSE)
  }
  out
}

.gkwq_link_scales <- function(link_scale, links, parts) {
  out <- stats::setNames(rep(1, length(parts)), parts)
  if (!is.null(link_scale)) {
    if (is.null(names(link_scale))) {
      out[] <- rep_len(link_scale, length(parts))
    } else {
      out[names(link_scale)] <- unlist(link_scale)
    }
  }
  out[["mu"]] <- 1  # the quantile is a probability; never rescaled
  out
}

## Inverse link and its derivative, mirroring inv_link() in src/gkwqreg.cpp.
.gkwq_linkinv <- function(eta, link, scale = 1) {
  switch(link,
    log = exp(eta),
    logit = scale * stats::plogis(eta),
    probit = scale * stats::pnorm(eta),
    cauchy = scale * stats::pcauchy(eta),
    cloglog = scale * (1 - exp(-exp(eta))),
    identity = eta,
    sqrt = eta^2,
    inverse = 1 / eta,
    `inverse-square` = 1 / sqrt(eta),
    stop("unknown link: ", sQuote(link), call. = FALSE)
  )
}

.gkwq_linkfun <- function(mu, link, scale = 1) {
  switch(link,
    log = log(mu),
    logit = stats::qlogis(mu / scale),
    probit = stats::qnorm(mu / scale),
    cauchy = stats::qcauchy(mu / scale),
    cloglog = log(-log1p(-mu / scale)),
    identity = mu,
    sqrt = sqrt(mu),
    inverse = 1 / mu,
    `inverse-square` = 1 / mu^2,
    stop("unknown link: ", sQuote(link), call. = FALSE)
  )
}

## d mu / d eta -- the chain-rule factor in marginal effects.
.gkwq_mu_eta <- function(eta, link, scale = 1) {
  switch(link,
    log = exp(eta),
    logit = scale * stats::dlogis(eta),
    probit = scale * stats::dnorm(eta),
    cauchy = scale * stats::dcauchy(eta),
    cloglog = scale * exp(eta - exp(eta)),
    identity = rep_len(1, length(eta)),
    sqrt = 2 * eta,
    inverse = -1 / eta^2,
    `inverse-square` = -0.5 * eta^(-1.5),
    stop("unknown link: ", sQuote(link), call. = FALSE)
  )
}

## ---------------------------------------------------------------------------
## Formula handling.
## ---------------------------------------------------------------------------

## Supplying MORE parts than the model has is an error naming the expected
## parts. gkwreg silently ignores the extras, which hides a real mistake.
.gkwq_check_parts <- function(fobj, parts, family, anchor) {
  n_rhs <- length(attr(fobj, "rhs"))
  if (n_rhs > length(parts)) {
    stop(sprintf(
      "formula has %d right-hand parts but family %s with anchor %s takes at most %d: %s.\n  Parts are separated by `|`, in the order given by gkwq_parts(\"%s\", \"%s\").",
      n_rhs, sQuote(family), sQuote(anchor), length(parts),
      paste(parts, collapse = " | "), family, anchor
    ), call. = FALSE)
  }
  n_rhs
}

## Build the model frame once from the union of all parts, then one design
## matrix per part. Parts beyond those supplied default to ~ 1.
.gkwq_model_data <- function(formula, data, parts, subset, na.action, weights,
                             offset, contrasts, family, anchor) {
  fobj <- Formula::as.Formula(formula)
  n_rhs <- .gkwq_check_parts(fobj, parts, family, anchor)

  ## Pad missing parts with ~ 1 so the frame covers every part we will fit.
  ## The rebuilt formula must keep the ORIGINAL environment: variables that live
  ## in the caller rather than in `data` are resolved through it, and dropping it
  ## would make `y ~ x` behave differently from `y ~ x | 1`.
  if (n_rhs < length(parts)) {
    env <- environment(fobj)
    dep <- function(e) paste(deparse(e, width.cutoff = 500L), collapse = " ")
    rhs <- vapply(seq_len(n_rhs),
                  function(k) dep(formula(fobj, lhs = 0, rhs = k)[[2L]]),
                  character(1))
    rhs <- c(rhs, rep("1", length(parts) - n_rhs))
    fobj <- Formula::as.Formula(stats::as.formula(
      paste(dep(fobj[[2L]]), "~", paste(rhs, collapse = " | ")), env = env))
  }

  ## Only pass weights/offset when they exist: passing NULL still makes
  ## model.frame look the symbol up, and in the formula's environment `weights`
  ## finds stats::weights.
  args <- list(formula = fobj, data = data, na.action = na.action,
               drop.unused.levels = TRUE)
  if (!is.null(subset)) args$subset <- subset
  if (!is.null(weights)) args$weights <- weights
  ## `offset` is deliberately NOT passed to model.frame. It would arrive back
  ## through model.offset() *in addition to* the per-part offset() terms already
  ## collected below, and the conditional quantile would receive it twice.
  mf <- do.call(stats::model.frame, args)

  y <- stats::model.response(mf, "numeric")
  n <- length(y)

  X <- vector("list", length(parts))
  terms_list <- vector("list", length(parts))
  offs <- vector("list", length(parts))
  names(X) <- names(terms_list) <- names(offs) <- parts

  for (k in seq_along(parts)) {
    mtk <- stats::terms(fobj, data = data, rhs = k, lhs = 0)
    X[[k]] <- stats::model.matrix(mtk, mf, contrasts.arg = contrasts)
    terms_list[[k]] <- mtk
    ## Read this part's offset BY NAME out of the model frame that already
    ## exists. Rebuilding a frame here -- model.frame(mtk, mf) -- re-evaluates
    ## the terms' predvars against mf, and mf holds a column literally named
    ## "log(z)" rather than a column "z", so ANY transformed term
    ## (log(), I(), poly(), splines) died with "object 'z' not found".
    offs[[k]] <- .gkwq_part_offset(mtk, mf, n)
  }

  w <- stats::model.weights(mf)
  if (is.null(w)) w <- rep(1, n)
  if (any(w < 0)) stop("negative weights are not allowed.", call. = FALSE)

  ## An `offset` ARGUMENT, as opposed to an offset() term inside the formula,
  ## attaches to the conditional quantile: it is the part an offset is normally
  ## meant for. offset() terms were already collected per part above, so this
  ## adds the argument only, and adds it once.
  if (!is.null(offset)) {
    ov <- as.numeric(offset)
    if (length(ov) != n) {
      stop("`offset` has length ", length(ov), " but there are ", n,
           " observations.", call. = FALSE)
    }
    offs[[1L]] <- offs[[1L]] + ov
  }

  list(y = y, X = X, weights = as.numeric(w), offsets = offs,
       terms = terms_list, mf = mf, formula = fobj,
       levels = .getXlevels(stats::terms(fobj, data = data), mf),
       contrasts = lapply(X, function(m) attr(m, "contrasts")))
}

.getXlevels <- function(Terms, m) stats::.getXlevels(Terms, m)

.gkwq_validate_y <- function(y, eps_y) {
  if (!is.numeric(y)) stop("the response must be numeric.", call. = FALSE)
  if (!length(y)) {
    stop("no observations remain after subsetting and na.action.", call. = FALSE)
  }
  if (anyNA(y)) stop("the response contains NA after na.action.", call. = FALSE)
  if (any(y <= 0 | y >= 1)) {
    n_bad <- sum(y <= 0 | y >= 1)
    stop(sprintf(
      "the response must lie strictly in (0,1); %d value%s at or outside the boundary. Exact 0s and 1s are outside the family's support and cannot be clamped away honestly.",
      n_bad, if (n_bad == 1L) "" else "s"), call. = FALSE)
  }
  pmin(pmax(y, eps_y), 1 - eps_y)
}

## Coefficient names are "part:term", and the order here IS the order the tape
## stacks its ADREPORTs, so sdreport's covariance needs no reordering.
.gkwq_coef_names <- function(X, parts) {
  unlist(lapply(parts, function(p) {
    cn <- colnames(X[[p]])
    if (is.null(cn) || !length(cn)) character(0) else paste0(p, ":", cn)
  }), use.names = FALSE)
}

## Offsets declared as offset() terms inside one formula part. attr(, "offset")
## indexes the part's own variable list, so the expressions are matched by name
## against the columns of the already-built model frame; nothing is re-evaluated.
.gkwq_part_offset <- function(mtk, mf, n) {
  idx <- attr(mtk, "offset")
  if (is.null(idx) || !length(idx)) return(rep(0, n))
  vars <- as.character(attr(mtk, "variables"))[-1L]
  cols <- match(vars[idx], names(mf))
  cols <- cols[!is.na(cols)]
  if (!length(cols)) return(rep(0, n))
  rowSums(vapply(mf[cols], as.numeric, numeric(n)))
}

`%||%` <- function(x, y) if (is.null(x)) y else x
