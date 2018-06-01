##' IEnorm1
##' @param TAC TAC based on MP
##' @param i numeric between 0 and time series length-1 indicating the timestep
##' @rdname IEnorm1
##' @details adds catch to TAC based on normal distribution: constant 6000 (sd 2000)
##' @export
IEnorm1 <- function(TAC,i){
    IEmean=rep(6000,100)
    IEsd=IEmean/4
    IE=max(0,rnorm(1,IEmean[i],IEmean[i]))
    ret=TAC+IE
    return(ret)
}
class(IEnorm1) <- append(class(IEnorm1),"IE")

##' IEnorm2
##' @param TAC TAC based on MP
##' @param i numeric between 0 and time series length-1 indicating the timestep
##' @rdname IEnorm2
##' @details adds catch to TAC based on normal distribution: constant 7200 (sd 2400)
##' @export
IEnorm2 <- function(TAC,i){
    IEmean=rep(7200,100)
    IEsd=IEmean/4
    IE=max(0,rnorm(1,IEmean[i],IEmean[i]))
    ret=TAC+IE
    return(ret)
}
class(IEnorm2) <- append(class(IEnorm2),"IE")

##' IEnorm3
##' @param TAC TAC based on MP
##' @param i numeric between 0 and time series length-1 indicating the timestep
##' @rdname IEnorm3
##' @details adds catch to TAC based on normal distribution: gradualy reduction first 3 years to 3000
##' @export
IEnorm3 <- function(TAC,i){
    IEmean=c(Reduce(function(v, x) .8*v , x=numeric(3),  init=6000, accumulate=TRUE)[-1],rep(3000,97))
    IEsd=IEmean/4
    IE=max(0,rnorm(1,IEmean[i],IEmean[i]))
    ret=TAC+IE
    return(ret)
}
class(IEnorm3) <- append(class(IEnorm3),"IE")

##' IEnorm4
##' @param TAC TAC based on MP
##' @param i numeric between 0 and time series length-1 indicating the timestep
##' @rdname IEnorm4
##' @details adds catch to TAC based on normal distribution: gradual reductino first 6 years to 1000
##' @export
IEnorm4 <- function(TAC,i){
    IEmean=c(Reduce(function(v, x) .75*v , x=numeric(6),  init=6000, accumulate=TRUE)[-1],rep(1000,94))
    IEsd=IEmean/4
    IE=max(0,rnorm(1,IEmean[i],IEmean[i]))
    ret=TAC+IE
    return(ret)
}
class(IEnorm4) <- append(class(IEnorm4),"IE")

##' IEnorm5
##' @param TAC TAC based on MP
##' @param i numeric between 0 and time series length-1 indicating the timestep
##' @rdname IEnorm4
##' @details adds catch to TAC based on normal distribution: constant 4800 (20% decrease)
##' @export
IEnorm5 <- function(TAC,i){
    IEmean=rep(6000*0.8,100)
    IEsd=IEmean/4
    IEmax(0,rnorm(1,IEmean[i],IEmean[i]))
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
    IEpar1=0.6
    IEpar2=15
    IE=rgamma(1,IEpar1,IEpar2)
    ret=TAC*(1-IE)
    return(ret) #hist(1-rgamma(1000000,0.6,15),500, freq = FALSE, ylab="%",xlab='Proportion of TAC attained',main='')
}
class(IEgamma1) <- append(class(IEgamma1),"IE")

##' IEnormgamma1
##' @param TAC TAC based on MP
##' @param i numeric between 0 and time series length-1 indicating the timestep
##' @rdname IEnormgamma1
##' @details adds catch to TAC based on normal distribution: constant 6000 (sd 2000) + gamma
##' @export
IEnormgamma1 <- function(TAC,i){
    IEmean=rep(6000,100)
    IEsd=IEmean/4
    IE=max(0,rnorm(1,IEmean[i],IEmean[i]))
    IEpar1=0.6
    IEpar2=15
    IEg=rgamma(1,IEpar1,IEpar2)
    ret=(TAC+IE)*(1-IEg)
    return(ret)
}
class(IEnormgamma1) <- append(class(IEnormgamma1),"IE")

##' IEnormgamma2
##' @param TAC TAC based on MP
##' @param i numeric between 0 and time series length-1 indicating the timestep
##' @rdname IEnormgamma2
##' @details adds catch to TAC based on normal distribution: constant 7200 (sd 2400) + gamma
##' @export
IEnormgamma2 <- function(TAC,i){
    IEmean=rep(7200,100)
    IEsd=IEmean/4
    IE=rnorm(1,IEmean[i],IEmean[i])
    IEpar1=0.6
    IEpar2=15
    IEg=rgamma(1,IEpar1,IEpar2)
    ret=(TAC+IE)*(1-IEg)
    return(ret)
}
class(IEnormgamma2) <- append(class(IEnormgamma2),"IE")

##' IEnormgamma3
##' @param TAC TAC based on MP
##' @param i numeric between 0 and time series length-1 indicating the timestep
##' @rdname IEnormgamma3
##' @details adds catch to TAC based on normal distribution: gradualy reduction first 3 years to 3000 + gamma
##' @export
IEnormgamma3 <- function(TAC,i){
    IEmean=c(Reduce(function(v, x) .8*v , x=numeric(3),  init=6000, accumulate=TRUE)[-1],rep(3000,97))
    IEsd=IEmean/4
    IE=max(0,rnorm(1,IEmean[i],IEmean[i]))
    IEpar1=0.6
    IEpar2=15
    IEg=rgamma(1,IEpar1,IEpar2)
    ret=(TAC+IE)*(1-IEg)
    return(ret)
}
class(IEnormgamma3) <- append(class(IEnormgamma3),"IE")

##' IEnormgamma4
##' @param TAC TAC based on MP
##' @param i numeric between 0 and time series length-1 indicating the timestep
##' @rdname IEnormgamma4
##' @details adds catch to TAC based on normal distribution: gradual reductino first 6 years to 1000 + gamma
##' @export
IEnormgamma4 <- function(TAC,i){
    IEmean=c(Reduce(function(v, x) .75*v , x=numeric(6),  init=6000, accumulate=TRUE)[-1],rep(1000,94))
    IEsd=IEmean/4
    IE=max(0,rnorm(1,IEmean[i],IEmean[i]))
    IEpar1=0.6
    IEpar2=15
    IEg=rgamma(1,IEpar1,IEpar2)
    ret=(TAC+IE)*(1-IEg)
    return(ret)
}
class(IEnormgamma4) <- append(class(IEnormgamma4),"IE")

##' IEnormgamma5
##' @param TAC TAC based on MP
##' @param i numeric between 0 and time series length-1 indicating the timestep
##' @rdname IEnormgamma4
##' @details adds catch to TAC based on normal distribution: constant 4800 (20% decrease) + gamma
##' @export
IEnormgamma5 <- function(TAC,i){
    IEmean=rep(6000*0.8,100)
    IEsd=IEmean/4
    IE=max(0,rnorm(1,IEmean[i],IEmean[i]))
    IEpar1=0.6
    IEpar2=15
    IEg=rgamma(1,IEpar1,IEpar2)
    ret=(TAC+IE)*(1-IEg)
    return(ret)
}
class(IEnormgamma5) <- append(class(IEnormgamma5),"IE")

##' IEnothing
##' @param TAC TAC based on MP
##' @param i numeric between 0 and time series length-1 indicating the timestep
##' @rdname IEnothing
##' @details keeps TAC as it is
##' @export
IEnothing <- function(TAC,i){
    return(TAC)
}
class(IEnothing) <- append(class(IEnothing),"IE")
