###########################################################################################################
## Functions for different implementation errors                                                         ##
## Function names should start with 'IE'                                                                 ##
## There are two types (should be in naming):                                                            ##
## - 'indep': IE independent of SSB that year                                                            ##
## - 'dep': dependent of SSB that year                                                                   ##
## all IE's should provide matrices of size nosim x nofutureyears                                        ##
## for simplicity there is a IE named IEnothing that assumes IE=0                                        ##
###########################################################################################################

############################## helper function ############################################################

##' cumsum.bounded
##' @param x vector of values to cumsum
##' @param lower lower bound
##' @param upper upper bound
##' @rdname cumsum.bounded
##' @details cumsum function but with bounds. Used to create a bounded random walk.
##' @export
cumsum.bounded <- function(x, lower = 0, upper = 500) {
    bsum <- function(x, y) min(upper, max(lower, x+y))
    if (length(x) > 1) Reduce(bsum, x, acc = TRUE) else x
}

############################## independent IEs ############################################################

##' IEindep6000
##' @param x number of simulations
##' @param y number of timestep
##' @param seed set.seed
##' @rdname IEindep6000
##' @details generates matrix nosim x nyears with implementation errors that will be added to the TAC
##' @export
IEindep6000 <- function(x,y,seed=NULL){
    IEmean=rep(6000,y)
    IEsd=IEmean/4
    if(!is.null(seed)) set.seed(seed)
    ret=mapply(function(mu,sigma){rnorm(mu,sigma,n=x)},mu=IEmean,sigma=IEsd)
    return(ret)
}
class(IEindep6000) <- append(class(IEindep6000),"IE")

##' IEindep7200
##' @param x number of simulations
##' @param y number of timestep
##' @param seed set.seed
##' @rdname IEindep7200
##' @details generates matrix nosim x nyears with implementation errors that will be added to the TAC
##' @export
IEindep7200 <- function(x,y,seed=NULL){
    IEmean=rep(7200,y)
    IEsd=IEmean/4
    if(!is.null(seed)) set.seed(seed)
    ret=mapply(function(mu,sigma){rnorm(mu,sigma,n=x)},mu=IEmean,sigma=IEsd)
    return(ret)
}
class(IEindep7200) <- append(class(IEindep7200),"IE")

##' IEindepdecr
##' @param x number of simulations
##' @param y number of timestep
##' @param seed set.seed
##' @rdname IEindepdecr
##' @details generates matrix nosim x nyears with implementation errors that will be added to the TAC00
##' @export
IEindepdecr <- function(x,y,seed=NULL){
    IEmean=c(Reduce(function(v, x) .8*v , x=numeric(3),  init=6000, accumulate=TRUE)[-1],rep(3000,y-3))[1:y]
    IEsd=IEmean/4
    ret=mapply(function(mu,sigma){rnorm(mu,sigma,n=x)},mu=IEmean,sigma=IEsd)
    return(ret)
}
class(IEindepdecr) <- append(class(IEindepdecr),"IE")

##' IEindepdecrsteep
##' @param x number of simulations
##' @param y number of timestep
##' @param seed set.seed
##' @rdname IEindepdecrsteep
##' @details generates matrix nosim x nyears with implementation errors that will be added to the TAC
##' @export
IEindepdecrsteep<- function(x,y,seed=NULL){
    IEmean=c(Reduce(function(v, x) .75*v , x=numeric(6),  init=6000, accumulate=TRUE)[-1],rep(1000,y-6))[1:y]
    IEsd=IEmean/4
    if(!is.null(seed)) set.seed(seed)
    ret=mapply(function(mu,sigma){rnorm(mu,sigma,n=x)},mu=IEmean,sigma=IEsd)
    return(ret)
}
class(IEindepdecrsteep) <- append(class(IEindepdecrsteep),"IE")

##' IEindep4800
##' @param x number of simulations
##' @param y number of timestep
##' @param seed set.seed
##' @rdname IEindep4800
##' @details generates matrix nosim x nyears with implementation errors that will be added to the TAC
##' @export
IEindep4800 <- function(x,y,seed=NULL){
    IEmean=rep(6000*0.8,y)
    IEsd=IEmean/4
    if(!is.null(seed)) set.seed(seed)
    ret=mapply(function(mu,sigma){rnorm(mu,sigma,n=x)},mu=IEmean,sigma=IEsd)
    return(ret)
}
class(IEindep4800) <- append(class(IEindep4800),"IE")

############################## dependent IEs ############################################################

##' IEdep2550
##' @param x number of simulations
##' @param y number of timestep
##' @param seed set.seed
##' @param ssb vector of ssb values
##' @param i future year index
##' @rdname IEdep2550
##' @details generates matrix of size nosim x i
##' @export
IEdep2550<- function(x,y,seed=NULL,ssb,i){
    # USA.base <- 7316.08 # tail(ctUSA,1)
    # USA.sigma <- 7815.438 #sd(diff(ctUSA))
    # USA.AC <- 0.826 #acf(ctUSA)$acf[2,1,1]
    # if(!is.null(seed)) set.seed(seed)
    # IErw <- t(mapply(function(x){cumsum.bounded(c(rnorm(1,0.375,0.08),rnorm(y-1,0,0.08)),0.25,0.50)},x=1:x)) #matplot(t(IErw),type='l')
    # IEdev <- mapply(function(mu,sigma){rnorm(mu,sigma,n=x)},mu=rep(0,y),sigma=rep(USA.sigma,y))
    #
    # if(!exists('USA.catch')) USA.catch <<- as.matrix(rep(USA.base,x))
    #
    # cnew <- USA.AC*USA.catch[,i]+(USA.catch[,i]+IEdev[,i])*sqrt(1-USA.AC^2)
    # cnew[cnew<0] <-0
    # print(cnew)
    # cnew <- sapply(cnew, function(y) min(max(y,ssb*0.15),ssb*0.35)) # bound values to 15 % - 35% true biomass
    # print(cnew)
    # USA.catch <<- cbind(USA.catch,cnew)
    # ret <- USA.catch[,2:ncol(USA.catch),drop=FALSE]*IErw[,1:i,drop=FALSE]
    # if(i==y) {attr(ret, "ctUSA") <- USA.catch
    # rm(USA.catch)
    # }
    # return(ret)

    if(!is.null(seed)) set.seed(seed)
    IErw <- t(mapply(function(x){cumsum.bounded(c(runif(1,0.25,0.5),rnorm(y-1,0,0.08)),0.25,0.50)},x=1:x)) #matplot(t(IErw),type='l')

    s <- c(1396, 1361, 938, 1320, 1644, 1998, 2724, 3891, 3929,
      4364, 4049, 2406, 2006, 1336, 1042, 1974, 2712, 1377, 1605, 1990,
      2683, 6150.98, 4520.68, 6806.87, 8273.29, 9345.23, 13860.37,
      16341.85, 15573.79, 16503.07, 33955.3, 30625.25, 14198.05, 9216.26,
      16141.38, 10109.36, 19297.08, 17613.46, 15436.33, 14302.86, 7481.7,
      14890.45, 28115.28, 35399.09, 57647.13, 44320.53, 58465.28, 26338.47,
      23337.45, 23443.41, 10911.59, 1612.87, 6132.86, 5343.28, 7099.19,
      7013.45, 7316.08)

    if(!exists('smat')) smat <<- matrix(s,x,length(s),byrow = T)

    cnew <-  apply(smat,1,function(x){
           pr <-  predict(ar(x,order.max=1))
           ne <-rnorm(1,pr$pred,pr$se)
    })
    cnew <- mapply(function(x,y){min(max(x,y*0.05),y*0.3)},cnew,ssb) # bound values to 0.5 % - 30% true biomass
    cnew <- mapply(function(x,y){min(max(x,y*0),y*2)},cnew,smat[,ncol(smat)]*2) # also bound them so can at most double
    cnew[cnew>20000] <- 20000

    smat <<- cbind(smat,cnew)
    ret <- smat[,(length(s)+1):ncol(smat),drop=FALSE]*IErw[,1:i,drop=FALSE]
    if(i==y) {attr(ret, "ctUSA") <- smat;rm(smat,envir=.GlobalEnv)}
    return(ret)
}
class(IEdep2550) <- append(class(IEdep2550),"IE")

##' IEdep5075
##' @param x number of simulations
##' @param y number of timestep
##' @param seed set.seed
##' @param ssb vector of ssb values
##' @param i future year index
##' @rdname IEdep5075
##' @details generates matrix of size nosim x i
##' @export
IEdep5075<- function(x,y,seed=NULL,ssb,i){
    if(!is.null(seed)) set.seed(seed)
    IErw <- t(mapply(function(x){cumsum.bounded(c(runif(1,0.50,0.75),rnorm(y-1,0,0.08)),0.50,0.75)},x=1:x)) #matplot(t(IErw),type='l')

    s <- c(1396, 1361, 938, 1320, 1644, 1998, 2724, 3891, 3929,
           4364, 4049, 2406, 2006, 1336, 1042, 1974, 2712, 1377, 1605, 1990,
           2683, 6150.98, 4520.68, 6806.87, 8273.29, 9345.23, 13860.37,
           16341.85, 15573.79, 16503.07, 33955.3, 30625.25, 14198.05, 9216.26,
           16141.38, 10109.36, 19297.08, 17613.46, 15436.33, 14302.86, 7481.7,
           14890.45, 28115.28, 35399.09, 57647.13, 44320.53, 58465.28, 26338.47,
           23337.45, 23443.41, 10911.59, 1612.87, 6132.86, 5343.28, 7099.19,
           7013.45, 7316.08)

    if(!exists('smat')) smat <<- matrix(s,x,length(s),byrow = T)

    cnew <-  apply(smat,1,function(x){
        pr <-  predict(ar(x,order.max=1))
        ne <-rnorm(1,pr$pred,pr$se)
    })
    cnew <- mapply(function(x,y){min(max(x,y*0.05),y*0.3)},cnew,ssb) # bound values to 0.5 % - 30% true biomass
    cnew <- mapply(function(x,y){min(max(x,y*0),y*2)},cnew,smat[,ncol(smat)]*2) # also bound them so can at most double
    cnew[cnew>20000] <- 20000

    smat <<- cbind(smat,cnew)
    ret <- smat[,(length(s)+1):ncol(smat),drop=FALSE]*IErw[,1:i,drop=FALSE]
    if(i==y) {attr(ret, "ctUSA") <- smat;rm(smat,envir=.GlobalEnv)}
    return(ret)
}
class(IEdep5075) <- append(class(IEdep5075),"IE")

############################## NO IE ############################################################

##' IEnothing
##' @param x number of simulations
##' @param y number of timestep
##' @param seed set.seed
##' @rdname IEnothing
##' @details keeps TAC as it is
##' @export
IEnothing <- function(x,y,seed=NULL){
    return(matrix(0,nrow=x,ncol=y))
}
class(IEnothing) <- append(class(IEnothing),"IE")
