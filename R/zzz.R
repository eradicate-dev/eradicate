.onUnload <- function (libpath) {
  library.dynam.unload("eradicate", libpath)
}
