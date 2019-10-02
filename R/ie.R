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
    IEsd=IEmean/8
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
    IEsd=IEmean/8
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
    IEmean=c(Reduce(function(v, x) .8*v , x=numeric(3),  init=5500, accumulate=TRUE)[-1],rep(3000,y-3))[1:y]
    IEsd=IEmean/8
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
    IEsd=IEmean/8
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
    IEsd=IEmean/8
    if(!is.null(seed)) set.seed(seed)
    ret=mapply(function(mu,sigma){rnorm(mu,sigma,n=x)},mu=IEmean,sigma=IEsd)
    return(ret)
}
class(IEindep4800) <- append(class(IEindep4800),"IE")


##' IEindep2019
##' @param x number of simulations
##' @param y number of timestep
##' @param seed set.seed
##' @rdname IEindep2019
##' @details generates matrix nosim x nyears with implementation errors that will be added to the TAC00
##' @export
IEindep2019 <- function(x,y,seed=NULL){
    IEmean=c(Reduce(function(v, x) .83*v , x=numeric(3),  init=5500, accumulate=TRUE)[-1],rep(3000,y-3))[1:y]
    IEsd=IEmean/8
    if(!is.null(seed)) set.seed(seed)
    ret=mapply(function(mu,sigma){rnorm(mu,sigma,n=x)},mu=IEmean,sigma=IEsd)
    return(ret)
}
class(IEindep2019) <- append(class(IEindep2019),"IE")

##' IEindep500
##' @param x number of simulations
##' @param y number of timestep
##' @param seed set.seed
##' @rdname IEindep500
##' @details generates matrix nosim x nyears with implementation errors that will be added to the TAC00
##' @export
IEindep500 <- function(x,y,seed=NULL){
    IEmean=c(Reduce(function(v, x) .76*v , x=numeric(5),  init=2000, accumulate=TRUE)[-1],rep(500,y-5))[1:y]
    IEsd=IEmean/8
    if(!is.null(seed)) set.seed(seed)
    ret=mapply(function(mu,sigma){rnorm(mu,sigma,n=x)},mu=IEmean,sigma=IEsd)
    return(ret)
}
class(IEindep500) <- append(class(IEindep500),"IE")

##' IEindep1000
##' @param x number of simulations
##' @param y number of timestep
##' @param seed set.seed
##' @rdname IEindep1000
##' @details generates matrix nosim x nyears with implementation errors that will be added to the TAC00
##' @export
IEindep1000 <- function(x,y,seed=NULL){
    IEmean=c(Reduce(function(v, x) .81*v , x=numeric(5),  init=3000, accumulate=TRUE)[-1],rep(1000,y-5))[1:y]
    IEsd=IEmean/8
    if(!is.null(seed)) set.seed(seed)
    ret=mapply(function(mu,sigma){rnorm(mu,sigma,n=x)},mu=IEmean,sigma=IEsd)
    return(ret)
}
class(IEindep1000) <- append(class(IEindep1000),"IE")

##' IEindep2000
##' @param x number of simulations
##' @param y number of timestep
##' @param seed set.seed
##' @rdname IEindep2000
##' @details generates matrix nosim x nyears with implementation errors that will be added to the TAC00
##' @export
IEindep2000 <- function(x,y,seed=NULL){
    IEmean=c(Reduce(function(v, x) .87*v , x=numeric(5),  init=4100, accumulate=TRUE)[-1],rep(2000,y-2))[1:y]
    IEsd=IEmean/8
    if(!is.null(seed)) set.seed(seed)
    ret=mapply(function(mu,sigma){rnorm(mu,sigma,n=x)},mu=IEmean,sigma=IEsd)
    return(ret)
}
class(IEindep2000) <- append(class(IEindep2000),"IE")

##' IEindep3000
##' @param x number of simulations
##' @param y number of timestep
##' @param seed set.seed
##' @rdname IEindep3000
##' @details generates matrix nosim x nyears with implementation errors that will be added to the TAC00
##' @export
IEindep3000 <- function(x,y,seed=NULL){
    IEmean=c(Reduce(function(v, x) .83*v , x=numeric(3),  init=5500, accumulate=TRUE)[-1],rep(3000,y-3))[1:y]
    IEsd=IEmean/8
    if(!is.null(seed)) set.seed(seed)
    ret=mapply(function(mu,sigma){rnorm(mu,sigma,n=x)},mu=IEmean,sigma=IEsd)
    return(ret)
}
class(IEindep3000) <- append(class(IEindep3000),"IE")

##' IEindep4000
##' @param x number of simulations
##' @param y number of timestep
##' @param seed set.seed
##' @rdname IEindep4000
##' @details generates matrix nosim x nyears with implementation errors that will be added to the TAC00
##' @export
IEindep4000 <- function(x,y,seed=NULL){
    IEmean=c(Reduce(function(v, x) .93*v , x=numeric(3),  init=5000, accumulate=TRUE)[-1],rep(4000,y-3))[1:y]
    IEsd=IEmean/8
    if(!is.null(seed)) set.seed(seed)
    ret=mapply(function(mu,sigma){rnorm(mu,sigma,n=x)},mu=IEmean,sigma=IEsd)
    return(ret)
}
class(IEindep4000) <- append(class(IEindep4000),"IE")

############################## dependent IEs ############################################################

##' IEdep0025
##' @param x number of simulations
##' @param y number of timestep
##' @param seed set.seed
##' @param ssb vector of ssb values
##' @param TAC TAC values
##' @param i future year index
##' @rdname IEdep0025
##' @details generates matrix of size nosim x i
##' @export
IEdep0025<- function(x,y,seed=NULL,ssb,TAC,i,past=NULL){
    if(!is.null(seed)) set.seed(seed)
    IErw <- t(mapply(function(x){cumsum.bounded(c(runif(1,0,0.25),rnorm(y-1,0,0.08)),0,0.20)},x=1:x)) #matplot(t(IErw),type='l')
    if ("ctUSA" %in% ls(envir = .GlobalEnv)) {
        s <- get("ctUSA", envir = .GlobalEnv)[,1]
    } else {
        print('To use IEdep a ctUSA objects need to be available in the global environment')
    }
    if(!exists('smat')){
        smat <<- matrix(s,x,length(s),byrow = T)
    }

    cnew <-  apply(smat,1,function(x){
        pr <- predict(ar(x,order.max=1))
        ne <- rnorm(1,pr$pred,pr$se)
    })
    cnew <- mapply(function(x,y){min(max(x,y*0.05),y*0.3)},cnew,ssb) # bound values to 0.5 % - 30% true biomass
    cnew <- mapply(function(x,y){min(max(x,y*0),y*2)},cnew,smat[,ncol(smat)]*2) # also bound them so can at most double
    cnew[cnew>20000] <- 20000

    smat <<- cbind(smat,cnew)
    ret <- smat[,(length(s)+1):ncol(smat),drop=FALSE]*IErw[,1:i,drop=FALSE]
    if(i==y) {attr(ret, "ctUSA") <- smat;rm(smat,envir=.GlobalEnv)}
    return(ret)
}
class(IEdep0025) <- append(class(IEdep0025),"IE")

##' IEdep2550
##' @param x number of simulations
##' @param y number of timestep
##' @param seed set.seed
##' @param ssb vector of ssb values
##' @param TAC TAC values
##' @param i future year index
##' @rdname IEdep2550
##' @details generates matrix of size nosim x i
##' @export
IEdep2550<- function(x,y,seed=NULL,ssb,TAC,i,past=NULL){

    if(!is.null(seed)) set.seed(seed)
    IErw <- t(mapply(function(x){cumsum.bounded(c(runif(1,0.25,0.5),rnorm(y-1,0,0.08)),0.25,0.50)},x=1:x)) #matplot(t(IErw),type='l')
    if ("ctUSA" %in% ls(envir = .GlobalEnv)) {
        s <- get("ctUSA", envir = .GlobalEnv)[,1]
    } else {
        print('To use IEdep a ctUSA objects need to be available in the global environment')
    }
    if(!exists('smat')){

        smat <<- matrix(s,x,length(s),byrow = T)
    }

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
##' @param TAC TAC values
##' @param i future year index
##' @rdname IEdep5075
##' @details generates matrix of size nosim x i
##' @export
IEdep5075<- function(x,y,seed=NULL,ssb,TAC,i,past=NULL){
    if(!is.null(seed)) set.seed(seed)
    IErw <- t(mapply(function(x){cumsum.bounded(c(runif(1,0.50,0.75),rnorm(y-1,0,0.08)),0.50,0.75)},x=1:x)) #matplot(t(IErw),type='l')
    if ("ctUSA" %in% ls(envir = .GlobalEnv)) {
        s <- get("ctUSA", envir = .GlobalEnv)[,1]
    } else {
        print('To use IEdep a ctUSA objects need to be available in the global environment')
    }
    if(!exists('smat')){
        smat <<- matrix(s,x,length(s),byrow = T)
    }
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

##' IEdep75100
##' @param x number of simulations
##' @param y number of timestep
##' @param seed set.seed
##' @param ssb vector of ssb values
##' @param TAC TAC values
##' @param i future year index
##' @rdname IEdep75100
##' @details generates matrix of size nosim x i
##' @export
IEdep75100<- function(x,y,seed=NULL,ssb,TAC,i,past=NULL){
    if(!is.null(seed)) set.seed(seed)
    IErw <- t(mapply(function(x){cumsum.bounded(c(runif(1,0.75,1),rnorm(y-1,0,0.08)),0.75,1)},x=1:x)) #matplot(t(IErw),type='l')
    if ("ctUSA" %in% ls(envir = .GlobalEnv)) {
        s <- get("ctUSA", envir = .GlobalEnv)[,1]
    } else {
        print('To use IEdep a ctUSA objects need to be available in the global environment')
    }
    if(!exists('smat')){

        smat <<- matrix(s,x,length(s),byrow = T)
    }
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
class(IEdep75100) <- append(class(IEdep75100),"IE")

##' IEdepcopy
##' @param x number of simulations
##' @param y number of timestep
##' @param seed set.seed
##' @param ssb vector of ssb values
##' @param TAC TAC values
##' @param i future year index
##' @rdname IEdepcopy
##' @details generates matrix of size nosim x i. US copies Canadian TAC and catches 25-50% of northern contingent.
##' @export
IEdepcopy<- function(x,y,seed=NULL,ssb,TAC,i,past=NULL){
    if(i==1){
        if(!is.null(seed)) set.seed(seed)
        IErw <- t(mapply(function(x){cumsum.bounded(c(runif(1,0.25,0.5),rnorm(y-1,0,0.08)),0.25,0.50)},x=1:x)) #matplot(t(IErw),type='l')
    }else{
        IErw <- attr(past,'IErw')
    }
    past[,i] <-  TAC*IErw[,i]
    attr(past,'IErw') <- IErw
    return(past)
}
class(IEdepcopy) <- append(class(IEdepcopy),"IE")


############################## NO IE or constant IE ############################################################

##' IEconstant
##' @param x number of simulations
##' @param y number of timestep
##' @param value value of IE (same unit as catch)
##' @rdname IEconstant
##' @details keeps TAC as it is
##' @export
IEconstant <- function(x,y,value=0){
    return(matrix(value,nrow=x,ncol=y))
}
class(IEconstant) <- append(class(IEconstant),"IE")
