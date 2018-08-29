##' onload function
##' @rdname onLoad
##' @details stuff that need's to be done when loading the package
.onLoad <- function(libname=find.package("CCAM"), pkgname="CCAM") {
    my_theme <- theme_bw()+theme(strip.background = element_blank())
    theme_set(my_theme)
    update_geom_defaults("line",   list(size = 2))
    assign("ccamcol", c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499"), globalenv())
}

##' check availability of certain objects
##' @param x object (character)
##' @details x can be either MP (management procedure) or IE (implementation error)
##' @export
avail <- function(x){
    ClassFilter <- function(y) inherits(get(y), x)
    Objs <- Filter( ClassFilter, ls("package:CCAM"))
    return(Objs)
}

##' copy a list (MPs or OMs)
##' @param x list to copy
##' @param n number of times to copy. If vector than the number of copies equals the sum of n (with the i'th element corresponding to name i)
##' @param name name basis for the list copies
##' @details list names are name plus a number from 1 to n, OM/MP label is automatically adapted
##' @export
copy <- function(x,n,name){
    nt <- sum(n)
    l <- rep(list(x),nt)
    mynames <- paste0(rep(name,n),unlist(sapply(n,function(x) seq(1,x))))

    l <- lapply(1:nt,function(y){myl <- x
                                ix <- grep('label', names(myl), value=TRUE)
                                myl[[ix]]<-mynames[y]
                                return(myl)})
    names(l) <- mynames
    list2env(l, envir = .GlobalEnv)
}

##' stick function
##' @param x value of jy/yref
##' @param a parameter that scales the catch limit (hight of the curve)
##' @param b parameter that scales the catch limit (slope after j0)
##' @param c parameter determining the initial slope
##' @param j0 parameter determining the bending point
##' @details basis for an index based HCR (as in redfish)
##' @export
stick <- function(x,a,b,c,j0){
    pen <- numeric(length(x))
    pen[x<j0] <- c*(x[x<j0]-j0)^2
    C <- a+b*(x-j0)-pen
    C[C<0] <- 0
    return(C)
}

##' geometric mean
##' @param x vector
##' @param na.rm a logical value indicating whether NA values should be stripped before the computation proceeds.
##' @details Geometric mean of a vector. Indicates the central tendency or typical value of a set of numbers by using the product of their values (as opposed to the arithmetic mean which uses their sum). The geometric mean is defined as the nth root of the product of n numbers.
##' @export
gmean = function(x, na.rm=TRUE){
    if(all(is.na(x))) return(NA)
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x[!is.na(x)]))
}

##' save pngs in one line
##' @param plot plot
##' @param wd working directory in which to plase image
##' @param name file name
##' @param dim dimensions of png in cm
##' @param ... arguments passed to png
##' @importFrom gridExtra grid.arrange
##' @export
savepng <- function(plot,wd, name,dim,...){
    png(file=paste0(wd,name,".png"),units="cm",res=300,width=dim[1],height=dim[2],...)
    if("ggplot" %in% class(plot)) grid.arrange(plot) else plot
    dev.off()
}
