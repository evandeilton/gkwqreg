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
  invisible()
}
