##' IEnorm1
##' @param TAC TAC based on MP
##' @param i numeric between 0 and time series length-1 indicating the timestep
##' @rdname IEnorm1
##' @details adds catch to TAC based on normal distribution
##' @export
IEnorm1 <- function(TAC,i){
    IEmean=rep(6000,100)
    IEsd=IEmean/3
    IE=rnorm(1,IEmean[i+1],IEmean[i+1])
    ret=TAC+IE
    return(ret)
}
class(IEnorm1) <- append(class(IEnorm1),"IE")

##' IEnorm1y2
##' @param TAC TAC based on MP
##' @param i numeric between 0 and time series length-1 indicating the timestep
##' @rdname IEnorm1
##' @details adds catch to TAC based on normal distribution
##' @export
IEnorm1y2 <- function(TAC,i){
    IEmean=rep(6000,100)
    IEsd=IEmean/3
    IE=rnorm(1,IEmean[i+1],IEmean[i+1])
    ret=TAC+IE
    return(ret)
}
class(IEnorm1y2) <- append(class(IEnorm1y2),"IE")

##' IEnorm2
##' @param TAC TAC based on MP
##' @param i numeric between 0 and time series length-1 indicating the timestep
##' @rdname IEnorm2
##' @details adds catch to TAC based on normal distribution
##' @export
IEnorm2 <- function(TAC,i){
    IEmean=c(rep(4000,3),rep(3000,100))
    IEsd=IEmean/3
    IE=rnorm(1,IEmean[i+1],IEmean[i+1])
    ret=TAC+IE
    return(ret)
}
class(IEnorm2) <- append(class(IEnorm2),"IE")

##' IEnorm3
##' @param TAC TAC based on MP
##' @param i numeric between 0 and time series length-1 indicating the timestep
##' @rdname IEnorm3
##' @details adds catch to TAC based on normal distribution
##' @export
IEnorm3 <- function(TAC,i){
    IEmean=c(Reduce(function(v, x) .8*v , x=numeric(3),  init=6000, accumulate=TRUE)[-1],rep(3000,100))
    IEsd=IEmean/3
    IE=rnorm(1,IEmean[i+1],IEmean[i+1])
    ret=TAC+IE
    return(ret)
}
class(IEnorm3) <- append(class(IEnorm3),"IE")

##' IEnorm4
##' @param TAC TAC based on MP
##' @param i numeric between 0 and time series length-1 indicating the timestep
##' @rdname IEnorm4
##' @details adds catch to TAC based on normal distribution
##' @export
IEnorm4 <- function(TAC,i){
    IEmean=rep(6000*1.2,100)
    IEsd=IEmean/3
    IE=rnorm(1,IEmean[i+1],IEmean[i+1])
    ret=TAC+IE
    return(ret)
}
class(IEnorm4) <- append(class(IEnorm4),"IE")

##' IEnorm5
##' @param TAC TAC based on MP
##' @param i numeric between 0 and time series length-1 indicating the timestep
##' @rdname IEnorm4
##' @details adds catch to TAC based on normal distribution
##' @export
IEnorm5 <- function(TAC,i){
    IEmean=rep(6000*0.8,100)
    IEsd=IEmean/3
    IE=rnorm(1,IEmean[i+1],IEmean[i+1])
    ret=TAC+IE
    return(ret)
}
class(IEnorm5) <- append(class(IEnorm5),"IE")

##' IEgamma1
##' @param TAC TAC based on MP
##' @param i numeric between 0 and time series length-1 indicating the timestep
##' @rdname IEgamma1
##' @details adds catch to TAC based on normal distribution
##' @export
IEgamma1 <- function(TAC,i){
    IEpar1=0.8
    IEpar2=15
    IE=rgamma(1,IEpar1,IEpar2)
    ret=TAC*(1-IE)
    return(ret) #hist(rgamma(1000000,0.8,15),500)
}
class(IEgamma1) <- append(class(IEgamma1),"IE")

