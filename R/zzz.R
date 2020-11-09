
.onUnload <- function (libpath) {
  library.dynam.unload("eradicate", libpath)
  library.dynam.unload("unmarked", libpath)
}
