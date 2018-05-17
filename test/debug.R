fixwinpath <- function() {
    PATH <- Sys.getenv("PATH")
    PATH <- paste0(R.home(), "/bin/x64;", PATH)
    PATH <- paste0("c:/Rtools/mingw_64/bin;", PATH)
    Sys.setenv(PATH=PATH)
}
fixwinpath()
setwd("C:/Users/VANBE/Desktop/post-doc/DATA/CCAM/test/")
TMB :: gdbsource ('CCAM.r' , interactive = TRUE )
