##' MPeggsurvey
##' @param data list with data objects
##' @param TAC.base TAC of previous year
##' @rdname MPeggsurvey
##' @details simply changes TAC in function of change in I (last year/year before). Change can not be bigger smaller than 0.5x or 2x.
##' @export
MPeggsurvey<-function(data, TAC.base){ #nosim
    change <- unlist(lapply(data,function(x){ Ix <- tail(x$logobs[which(x$aux[,2]==3),1],2)
                             change <- Ix[2]/Ix[1]
                             change <- min(max(change,0.5),2)
                             return(change)}))
    TAC <- TAC.base*change
    return(TAC)
}
class(MPeggsurvey) <- append(class(MPeggsurvey),"MP")

##' MPeggsurveytrail3interim
##' @param data list with data objects
##' @param TAC.base TAC of previous year
##' @param i year in the future
##' @rdname MPeggsurveytrail3interim
##' @details simply changes TAC in function of change in I (last year/ geometric mean 3 years before). Change can not be bigger smaller than 0.5x or 2x.
##' @export
MPeggsurveytrail3interim<-function(data,TAC.base,i){
    if(!i %% 2 == 1){
        TAC <- TAC.base
    }else{
    change <- unlist(lapply(data,function(x){ Ix <- tail(x$logobs[which(x$aux[,2]==3),1],4)
                                            change <- Ix[4]/gmean(Ix[1:3])
                                            change <- min(max(change,0.5),2)
                                            return(change)}))
    TAC <- TAC.base*change
    }
    return(TAC)
}
class(MPeggsurveytrail3interim) <- append(class(MPeggsurveytrail3interim),"MP")

##' MPeggsurveytrail3
##' @param data list with data objects
##' @param TAC.base TAC of previous year
##' @param i year in the future
##' @rdname MPeggsurveytrail3
##' @details simply changes TAC in function of change in I (last year/ geometric mean 3 years before). Change can not be bigger smaller than 0.5x or 2x.
##' @export
MPeggsurveytrail3<-function(data,TAC.base){
        change <- unlist(lapply(data,function(x){ Ix <- tail(x$logobs[which(x$aux[,2]==3),1],4)
                                                change <- Ix[4]/gmean(Ix[1:3])
                                                change <- min(max(change,0.5),2)
                                                return(change)}))
        TAC <- TAC.base*change
    return(TAC)
}
class(MPeggsurveytrail3) <- append(class(MPeggsurveytrail3),"MP")

##' MPeggsurveytarget
##' @param data list with data objects
##' @param TAC.base TAC of previous year
##' @param i year in the future
##' @rdname MPeggsurveytarget
##' @details simply changes TAC in function of change in I (last year/ geometric mean 3 years before). Change can not be bigger smaller than 0.5x or 2x.
##' @export
MPeggsurveytarget<-function(data,TAC.base){
    Ihist <- cbind(I=exp(data[[1]]$logobs[,1]),data[[1]]$aux)
    Ihist <- Ihist[which(Ihist[,3]==3),c(1,2)]

    Itarget <- gmean(Ihist[which(Ihist[,2] %in% 1979:1994),1])
    I0 <- gmean(Ihist[which(Ihist[,2] %in% c(1996:2000,2004:2008)),1])
    TACtarget <- 30000
    TACmin <- 8000
    w = TACmin/TACtarget
    Irecent <- exp(unlist(lapply(data,function(x){gmean(tail(x$logobs[which(x$aux[,2]==3),1],3))})))

    ItoTac <- function(x){
        if(x<I0){
            TAC <- w*TACtarget*(x/I0)^2
        }else{
            TAC <- TACtarget*(w+(1-w)*((x-I0)/(Itarget-I0)))
        }
        TAC <- ifelse(TAC<TACmin,0,TAC)
        return(TAC)
    }

    TAC <- sapply(Irecent, ItoTac)
    if(length(TAC)==1) TAC <- rep(TAC,length(TAC.base))
    return(TAC)

    # par(mfrow=c(2,1),mar=c(4,4,1,2))
    # plot(Ihist[,2],Ihist[,1],type='o',pch=16,ylab='SSB (t)',xlab='Year', yaxs='i', xaxs='i')
    # points(Ihist[which(Ihist[,2] %in% 1979:1994),2],Ihist[which(Ihist[,2] %in% 1979:1994),1],col='green',pch=16)
    # points(Ihist[which(Ihist[,2] %in% c(1996:2000,2004:2008)),2],Ihist[which(Ihist[,2] %in% c(1996:2000,2004:2008)),1],col='red',pch=16)
    # points(Ihist[which(Ihist[,2] %in% c(2014:2016)),2],Ihist[which(Ihist[,2] %in% c(2014:2016)),1],col='grey',pch=16)
    #
    # x <- seq(0,Itarget*1.2,length.out = 500)
    # y <- sapply(x, ItoTac)
    # plot(x,y,type='l',ylab='TAC',xlab='Index',lwd=3, yaxs='i', xaxs='i')
    # abline(v=I0,col='red',lty=2,lwd=2)
    # abline(v=Itarget,col='green',lty=2,lwd=2)
    # abline(h=TACtarget,col='grey',lty=2,lwd=2)
    # abline(h=w*TACtarget,col='grey',lty=2,lwd=2)
    # text(x=(Itarget-I0)/2,y=TACtarget,labels='TAC target',adj=c(0,0),col='grey')
    # text(x=(Itarget-I0)/2,y=TACmin,labels='TAC minimum',adj=c(0,0),col='grey')
    # text(x=Itarget,y=0,labels='I target',adj=c(-0.2,0),col='green')
    # text(x=I0,y==0,labels='I0',adj=c(-0.4,0),col='red')
    # abline(v=Irecent[1],col='darkgrey',lty=3,lwd=2)
}
class(MPeggsurveytarget) <- append(class(MPeggsurveytarget),"MP")

##' MPeggsurveytargetinterim
##' @param data list with data objects
##' @param TAC.base TAC of previous year
##' @param i year in the future
##' @rdname MPeggsurveytargetinterim
##' @details simply changes TAC in function of change in I (last year/ geometric mean 3 years before). Change can not be bigger smaller than 0.5x or 2x.
##' @export
MPeggsurveytargetinterim<-function(data,TAC.base,i){
    if(!i %% 2 == 1){
        TAC <- TAC.base
    }else{
        Ihist <- cbind(I=exp(data[[1]]$logobs[,1]),data[[1]]$aux)
        Ihist <- Ihist[which(Ihist[,3]==3),c(1,2)]
        Itarget <- gmean(Ihist[which(Ihist[,2] %in% 1979:1994),1])
        I0 <- gmean(Ihist[which(Ihist[,2] %in% c(1996:2000,2004:2008)),1])
        TACtarget <- 30000
        TACmin <- 8000
        w = TACmin/TACtarget
        Irecent <- exp(unlist(lapply(data,function(x){gmean(tail(x$logobs[which(x$aux[,2]==3),1],3))})))

        ItoTac <- function(x){
            if(x<I0){
                TAC <- w*TACtarget*(x/I0)^2
            }else{
                TAC <- TACtarget*(w+(1-w)*((x-I0)/(Itarget-I0)))
            }
            TAC <- ifelse(TAC<TACmin,0,TAC)
            return(TAC)
        }

        TAC <- sapply(Irecent, ItoTac)
        if(length(TAC)==1) TAC <- rep(TAC,length(TAC.base))
    }
    return(TAC)
}
class(MPeggsurveytargetinterim) <- append(class(MPeggsurveytargetinterim),"MP")

##' MPf40stages
##' @param fit output from ccam.fit
##' @param data list with data objects
##' @param parameters list of parameters objects
##' @param TAC.base TAC of previous year
##' @param i year in the future
##' @details determine TAC by fiting censored model, determining F40% and health zones and applying F based on the latter.
##' @export
##' @import parallel
MPf40base <- function(fit, data,parameters,TAC.base, i){
    if(!i %% 2 == 1){
        TAC <- TAC.base
    }else{
        ncores <- detectCores()
        cl <- makeCluster(ncores) #set up nodes
        clusterExport(cl, varlist=c("fit","parameters"), envir=environment())
        obsfit <- parLapply(cl,data,function(x){
                    cen <- 'CE' %in% fit$conf$obsLikelihoodFlag
                    obsfit <- CCAM::ccam.fit(x,fit$conf,parameters,silent=FALSE,paracheck=FALSE,phase=ifelse(cen,1,NULL))
                    if(cen){
                        cenpar <- obsfit$pl
                        obsfit <- CCAM::ccam.fit(x,fit$conf,cenpar,silent=FALSE,paracheck=FALSE,phase=2)
                    }
                    return(obsfit)
                    })
        stopCluster(cl)
        class(obsfit) <- 'ccamset'

        # determine reference point and health zones
        catches <- catchtable(obsfit)
        ref <- ypr(obsfit,what=c('f40ssb','f40'))
        ssbref <- unlist(ref[which(names(ref)=='f40ssb')])
        fref <- unlist(ref[which(names(ref)=='f40')])
        Bupper <- 0.8*ssbref
        Blim <- 0.4*ssbref
        ssbp <- ssbtable(obsfit)
        ssbp <- ssbp[which(rownames(ssbp)==max(rownames(ssbp))),1]

        # determine TAC based on the latter
        TAC=do.call('c',lapply(1:length(data), function(x){
            # no fishing
            if(ssbp[x]<Blim[x]){TAC=1}else{
                # maximal fishing
                if(ssbp[x]>Bupper[x]){
                    fut <- forecast(obsfit[[x]],fval=fref[x],nosim = 300,ave.years = max(obsfit[[x]]$data$years)+(-9:0),rec.years = 1969:max(obsfit[[x]]$data$years), rec.meth = 4)
                    Clim <- attr(fut,"tab")[,'catch:median']
                }else{
                    # in between fishing
                    fref <- (fref[x]*(ssbp[x]-Blim[x]))/(Bupper[x]-Blim[x]) # catious zone
                    fut <- forecast(obsfit[[x]],fval=fref[x],nosim = 300,ave.years = max(obsfit[[x]]$data$years)+(-9:0),rec.years = 1969:max(obsfit[[x]]$data$years), rec.meth = 4)
                    Clim <- attr(fut,"tab")[,'catch:median']
                }
                declared <- exp(obsfit[[x]]$data$logobs[which(!is.na(obsfit[[x]]$data$logobs[,2])),1])
                predicted <- exp(obsfit[[x]]$rep$predObs[which(!is.na(obsfit[[x]]$data$logobs[,2]))])
                U <- predicted-declared
                U <- max(c(0,mean(tail(U,3)))) # undeclared catch 3 last years estimated by model, never lower than 0
                TAC <-  Clim-U
            }
            return(TAC)
            }))

    }
   return(TAC)
}
class(MPf40base) <- append(class(MPf40base),"MP")


##' MPspm
##' @param fit output from ccam.fit
##' @param sim simulations of states
##' @param TAC.base TAC of previous state
##' @rdname MPspm
##' @details determine TAC by fiting surplus production model from SPiCT package
##' @export
##' @importFrom  spict fit.spict get.par manage prop.F
MPspm <- function(data,TAC.base,i){
    if(!i %% 2 == 1){
        TAC <- TAC.base
    }else{
        ncores <- detectCores()
        cl <- makeCluster(ncores) #set up nodes
        clusterExport(cl, varlist=c("TAC.base"), envir=environment())
        #obsfit <- lapply(data, function(x){
        obsfit <- parLapply(cl,data,function(x){
            o <- x$logobs
            aux <- x$aux
            dl <- list(obsC = unname(exp(o[which(aux[,2]==1),1])),
                       timeC = unname(aux[which(aux[,2]==1),1]),
                       obsI = unname(exp(o[which(aux[,2]==3),1])),
                       timeI = unname(aux[which(aux[,2]==3),1]))
            res <- spict::fit.spict(dl)
            Bmsy <- spict::get.par('Bmsys', res)[2]
            Blast <- spict::get.par('logBl', res, exp=TRUE)[2]
            Blim <- Bmsy*0.4
            Bupper <- Bmsy*0.8
            if(Blast<Blim){TAC=1}else{
                # maximal fishing
                if(Blast>Bupper){
                    res <- spict::manage(res, scenarios = 3) #fish at msy
                    TAC <- spict::get.par("Cp", res)[2]
                }else{
                    # in between fishing
                    Fmsy <- spict::get.par("logFmsy", res, exp = TRUE)[2]
                    Fref <- (Fmsy*(Blast-Blim))/(Bupper-Blim) # catious zone
                    Flast <- spict::get.par("logF", res, exp = TRUE)[res$inp$indpred[1], 2]
                    Fchange <- Fref/Flast
                    res <- spict::prop.F(Fchange, res$inp, res, which(res$inp$time >= res$inp$manstart))
                    TAC <- spict::get.par("Cp", res)[2]
                }
            }
            return(TAC)
           })
        stopCluster(cl)
        TAC <- do.call('c',obsfit)
        if(length(TAC)==1) TAC <- rep(TAC, length(TAC.base))
    }
    return(TAC)
}
class(MPspm) <- append(class(MPspm),"MP")
