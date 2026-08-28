## Register estfun() and bread() for the sandwich package's generics only when
## sandwich is actually installed. It sits in Suggests, so we cannot import its
## generics unconditionally, but making sandwich::vcovHC(), sandwich::vcovCL()
## and lmtest::coeftest() work out of the box is worth the small dance.
.onLoad <- function(libname, pkgname) {
  if (requireNamespace("sandwich", quietly = TRUE)) {
    registerS3method("estfun", "gkwqreg", estfun.gkwqreg,
                     envir = asNamespace("sandwich"))
    registerS3method("bread", "gkwqreg", bread.gkwqreg,
                     envir = asNamespace("sandwich"))
  }
  ## lmtest exports its own lrtest() generic, which masks ours whenever lmtest
  ## is attached after this package. Its default method would then run a naive
  ## likelihood-ratio test on two gkwqreg fits -- including across different
  ## quantile levels or anchors, where the test is meaningless and
  ## anova.gkwqreg() deliberately refuses. Registering our method against their
  ## generic makes the guard hold whichever lrtest() the user reaches.
  if (requireNamespace("lmtest", quietly = TRUE)) {
    registerS3method("lrtest", "gkwqreg", lrtest.gkwqreg,
                     envir = asNamespace("lmtest"))
  }
  invisible()
}
