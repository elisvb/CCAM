##' rmvnorm helper function to draw multivariate normal samples
##' @param n the number of samples.
##' @param mu the mean vector.
##' @param Sigma a positive-definite symmetric matrix specifying the covariance matrix.
##' @details Generates samples via the Cholesky decomposition, which is less platform dependent than eigenvalue decomposition.
##' @return If n = 1 a vector of the same length as mu, otherwise an n by length(mu) matrix with one sample in each row.
##' @export
rmvnorm <- function(n = 1, mu, Sigma){
  p <- length(mu)
  if(!all(dim(Sigma) == c(p, p))){
    stop("incompatible arguments")
  }
  idx <- diag(Sigma) > .Machine$double.xmin
  L <- matrix(0,p,p)
  if(any(idx)){
    L[idx,idx] <- chol(Sigma[idx,idx])
  }
  X <- matrix(rnorm(p * n), n)
  X <- drop(mu) + t(X%*%L)
  if(n == 1){
    drop(X)
  }else{
    t(X)
  }
}


##' forecast function to do shortterm
##' @param fit an assessment object of type ccam, as returned from the function ccam.fit
##' @param fscale a vector of f-scales. See details.
##' @param catchval a vector of target catches. See details.
##' @param fval a vector of target f values. See details.
##' @param MP a vector of management strategie names. see avail('MP')
##' @param nosim number of simulations default is 1000
##' @param year.base starting year default last year in assessment. Currently it is only supported to use last assessment year or the year before
##' @param ave.years vector of years to average for weights, maturity, M and such
##' @param rec.years vector of years to use to resample recruitment from
##' @param rec.meth methode to project recruitment
##' @param label optional label to appear in short table
##' @param overwriteSelYears if a vector of years is specified, then the average selectivity of those years is used (not recommended)
##' @param deterministic option to turn all process noise off (not recommended, as it will likely cause bias)
##' @param CZ critical zone (to calculate probability ssb above)
##' @param HZ Healty zone (to calculate probability ssb above)
##' @param IE implementation error (vector). See avail('IE')
##' @param capLower lower limit of TAC (before IE)
##' @param capUpper upper limit of TAC (before IE)
##' @param UL.years vector of years to use to resample upper/lower catch bound ranges from
##' @param TAC.base TAC last year
##' @details There are four ways to specify a scenario; fscale, catchval, fval and MP.
##' @return an object of type ccamforecast
##' @importFrom stats median uniroot quantile
##' @export
forecast <- function(fit, fscale=NULL, catchval=NULL, fval=NULL, MP=NULL,nosim=1000, year.base=max(fit$data$years),
                     ave.years=max(fit$data$years)+(-4:0), rec.years=max(fit$data$years)+(-9:0),rec.meth=5, MPlabel=NULL,
                     OMlabel=NULL,overwriteSelYears=NULL, deterministic=FALSE,CZ=0,HZ=0,IE=NULL,
                     capLower=NULL,capUpper=NULL,UL.years=max(fit$data$years)+(-4:0),TAC.base=NULL){
  parameters <- as.list(match.call())

  args <- c('fval','fscale','catchval','MP')

  if(missing(fscale)&missing(fval)&missing(catchval)&missing(MP)) stop("No scenario is specified")

  defined <- ls()
  passed <- names(as.list(match.call())[-1])

  if (any(!defined %in% passed)) {
      d <- setdiff(defined, passed)
      d <- intersect(d,args)
      for(i in d){
         assign(i,rep(NA,nrow(cbind(fscale, catchval, fval,MP))))
      }
  }

  if(!all(rowSums(!is.na(cbind(fscale, catchval, fval,MP)))==1)){
      stop("For each forecast year exactly one of fscale, catchval, fval or MP must be specified (all others must be set to NA)")
  }

  ny <-length(fscale)

  MPs=avail(x='MP')
  if(all(!is.na(MP)) & any(!MP %in% MPs) & !is.numeric(MP)) stop("Undefined MP (see avail('MP'))")

  IEs=avail(x='IE')
  if(is.null(IE)) IE <- rep(NA,length(fval))
  if(any(!IE %in% c(IEs,NA))) stop("Undefined IE (see avail('IE'))")

  if(!is.null(overwriteSelYears)){
    fromto <- fit$conf$fbarRange-(fit$conf$minAge-1)
    Ftab <- faytable(fit)
    fixedsel <- colMeans(Ftab[as.integer(rownames(Ftab))%in%overwriteSelYears,,drop=FALSE])
    fixedsel <- fixedsel/mean(fixedsel[fromto[1]:fromto[2]])
  }

  getF <<- function(x){
    idx <- fit$conf$keyLogFsta[1,]+1
    nsize <- length(idx)
    ret <- exp(x[nsize+idx])
    ret[idx==0] <- 0
    if(!is.null(overwriteSelYears)){
      fromto <- fit$conf$fbarRange-(fit$conf$minAge-1)
      thisfbar<-mean(ret[fromto[1]:fromto[2]])
      ret<-fixedsel*thisfbar
    }
    ret
  }

  fbar <<-function(x){
    fromto <- fit$conf$fbarRange-(fit$conf$minAge-1)
    mean(getF(x)[fromto[1]:fromto[2]])
  }

  getN <<- function(x){
    idx <- fit$conf$keyLogFsta[1,]+1
    nsize <- length(idx)
    ret <- exp(x[1:nsize])
    ret
  }

  getState <<- function(N,F){
    #k <- fit$conf$keyLogFsta[1,]
    #F <- F[k>=0]
    #k <- unique(k[k>=0]+1)
    x <- log(c(N,F))
    x
  }

  getProcessVar <- function(fit){
    cof <- coef(fit)
    sdN <- exp(cof[names(cof)=="logSdLogN"][fit$conf$keyVarLogN+1])
    sdN[1]<-0
    nN <- length(sdN)
    sdF <- exp(cof[names(cof)=="logSdLogFsta"][fit$conf$keyVarF+1])
    nF <- length(sdF)
    corr <- diag(nF)
    cov <- matrix(0,nrow=nN+nF,ncol=nN+nF)
    cov[1:nN,1:nN] <- diag(sdN^2)
    cov[nN+1:nF,nN+1:nF] <- (sdF%*%t(sdF))*corr
    cov
  }

  getRPS <- function(fit){
      X <- summary(fit)
      n<-nrow(X)
      lag <- fit$conf$minAge
      idxR <- (lag+1):n
      idxS <- 1:(n-lag)
      R<-X[idxR,1]
      S<-X[idxS,4]
      RPS <- R/S
      return(RPS)
  }

  getRec <- function(sim, recpool, rec.meth, deterministic, i, simpara){
      rho <- 0
      sigma <- 0

      Rautocorr <- function(r,mu,sigma,rho){
          eta <- r-mu+sigma^2/2
          delta <- rnorm(length(r),0,sigma)
          etanew <- rho*eta+delta*sqrt(1-rho^2)
          R <- exp(mu+etanew-sigma^2/2)
          return(R)
      }
      BH <- function(x, a, b){a+log(x)-log(1.0+exp(b)*x)}
      RI <- function(x, a, b){a+log(x)-exp(b)*x}  # not implemented yet

      switch(rec.meth,
             "one" <- { # see Johnson 2016  (BH, lag 1 autocorrelation)
                 ssb.before <- apply(sim, 1, ssb, nm=nm, sw=sw, mo=mo, pm=pm, pf=pf)
                 rho <- acf(recpool,plot=FALSE)$acf[2,1,1]
                 mu <- BH(ssb.before,a=simpara[,'rec_loga'],b=simpara[,'rec_logb'])
                 sigma <- exp(simpara[,'logSdLogN'])
                 if(deterministic) sigma <- 0
                 recs <- Rautocorr(sim[,1],mu,sigma,rho)
             },
             "two" <- { # see Johnson 2016 (mean, lag 1 autocorrelation)
                 rho <- acf(recpool,plot=FALSE)$acf[2,1,1]
                 mu <- mean(log(recpool))
                 sigma <- sd(log(recpool))
                 if(deterministic) sigma <- 0
                 recs <- Rautocorr(sim[,1],mu,sigma,rho)
             },
             "three" <- { # see Johnson 2016 (BH, no AC)
                 ssb.before <- apply(sim, 1, ssb, nm=nm, sw=sw, mo=mo, pm=pm, pf=pf)
                 mu <- BH(ssb.before,a=simpara[,'rec_loga'],b=simpara[,'rec_logb'])
                 sigma <- exp(simpara[,'logSdLogN'])
                 if(deterministic) sigma <- 0
                 recs <- Rautocorr(sim[,1],mu,sigma,rho)
             },
             "four" <- { # see Johnson 2016 (mean, no AC)
                 mu <- mean(log(recpool))
                 sigma <- sd(log(recpool))
                 if(deterministic) sigma <- 0
                 recs <- Rautocorr(sim[,1],mu,sigma,rho)
             },
             "five" <- { # sample recs (original SAM method, no autocorrelation)
                 if(deterministic){
                     recs <- rep(mean(recpool),nrow(sim))
                 }else{
                     recs <- sample(recpool, nrow(sim), replace=TRUE)
                  }
             },
             "six" <- { # sample recs, trailing
                 ref <- tail(recpool,i+2)
                 if(deterministic){
                     recs <- rep(mean(ref),nrow(sim))
                 }else{
                     recs <- sample(ref, nrow(sim), replace=TRUE)
                 }
             },
             "seven" <- { # sample RPS (suggestion Daniel, little autocorrelation)
                 ssb.before <- apply(sim, 1, ssb, nm=nm, sw=sw, mo=mo, pm=pm, pf=pf)
                 if(deterministic){
                     recs <- rep(mean(recpool),nrow(sim))*ssb.before
                 }else{
                     recs <- unname(sample(recpool,nrow(sim), replace=TRUE)*ssb.before)
                 }
             },
             "eight" <- { # sample RPS, trailing (suggestion Daniel, little autocorrelation))
                 ssb.before <- apply(sim, 1, ssb, nm=nm, sw=sw, mo=mo, pm=pm, pf=pf)
                 ref <- tail(recpool,i+2)
                 if(deterministic){
                     recs <- rep(mean(ref),nrow(sim))*ssb.before
                 }else{
                     recs <- unname(sample(ref,nrow(sim), replace=TRUE)*ssb.before)
                 }
             },
             {
                 stop('rec.meth should be 1 to 8') #block sampling would be another one
             }
      )
      attr(recs,'rho') <- rho
      attr(recs,'sigma') <- sigma
      return(recs)
  }

  step <- function(x, nm, scale, rec, inyear=FALSE){
    F <- getF(x)
    N <- getN(x)
    if(!inyear){
      Z <- F+nm
      n <- length(N)
      N <- c(rec,N[-n]*exp(-Z[-n])+c(rep(0,n-2),N[n]*exp(-Z[n]))) # recruitment gets samples without replacement!!!
    }
    F <- F*scale
    xx <- getState(N,F)
    return(xx)
  }

  scaleF <- function(x, scale){
    F <- getF(x)*scale
    N <- getN(x)
    xx <- getState(N,F)
    return(xx)
  }

  maxF <- function(x, Flim=2){
      F <- getF(x)
      F[F==0] <- 0.00000001
      if(any(F>Flim)){
          F <- F/max(F)
          F <- F*Flim
      }
      N <- getN(x)
      xx <- getState(N,F)
      return(xx)
  }

  sel <-function(x){
    return(getF(x)/fbar(x))
  }

  caa <<- function(x, nm){
      F <- getF(x)
      Z <- F+nm
      N <- getN(x)
      C <- F/Z*(1-exp(-Z))*N
      return(C)
  }

  catch <<- function(x, nm, cw){
    C <- caa(x,nm)
    return(sum(cw*C))
  }

  ssb <<- function(x, nm, sw, mo, pm, pf){
    F <- getF(x)
    N <- getN(x)*exp(-pm*nm-pf*F)
    return(sum(N*mo*sw))
  }

  if(!all(fit$data$fleetTypes %in% c(3,6))) warning('MPs based on fleettypes different from 3 and 6 need more code')

  # Get final state
  fy <- max(fit$data$years)
  if(year.base==fy){
    est <- fit$sdrep$estY
    cov <- fit$sdrep$covY
  }
  if(year.base==(fy-1)){
    est <- fit$sdrep$estYm1
    cov <- fit$sdrep$covYm1
  }
  if(year.base<(max(fit$data$years)-1)){
    stop("State not saved, so cannot proceed from this year")
  }

  idxN <- (fit$conf$minAge-fit$conf$minAge+1):(fit$conf$maxAge-fit$conf$minAge+1)
  idxF <-  fit$conf$keyLogFsta[1,]+fit$conf$maxAge-fit$conf$minAge+2
  idx <- c(idxN,idxF)
  rm <- !duplicated(idx)

  sim <- rmvnorm(nosim, mu=est[rm], Sigma=cov[rm,rm])
  sim <- sim[,idx]

  if(is.null(overwriteSelYears)){
    if(!all.equal(est,getState(getN(est),getF(est))))stop("Sorry something is wrong here (check code for getN, getF, and getState)")
  }

  if(is.null(overwriteSelYears)){
      if(!all.equal(est,getState(getN(est),getF(est))))stop("Sorry something is wrong here (check code for getN, getF, and getState)")
  }

  # get parameter values (only used when MP)
  if(any(!is.na(MP))){
    simpara <- rmvnorm(nosim, mu=fit$sdrep$par.fixed, Sigma=fit$sdrep$cov.fixed)
    colnames(simpara) <- names(fit$sdrep$par.fixed)
  }

  # recruitment
  if(rec.meth %in% c(7,8)){
     recpool <- getRPS(fit)
  }else{
     rectab <- rectable(fit)
     recpool <- rectab[rownames(rectab)%in%rec.years,1]
  }

  # Upper lower limits
  UL <- exp(fit$data$logobs[!is.na(fit$data$logobs[,2]),])
  pred <- exp(fit$rep$predObs[which(!is.na(fit$data$logobs[,2]))])
  ULtab <- UL/pred
  ULpool <- ULtab[rownames(ULtab)%in%UL.years,]


  # last year stats to compare subsequent years with
  ssbf <- apply(sim, 1, ssb, nm=nm, sw=sw, mo=mo, pm=pm, pf=pf)
  refs <- t(apply(simpara,1,function(x) unlist(ypr(fit,deterministic=FALSE,simpara=x)[c('f40','f40ssb','f40U')])))

  #helper functions
  getThis <- function(x,y){
      x <- fit$data[[x]]
      x <- x[which(rownames(x)==y),]
      x
  }

  Csim <- function(x,sim){  #Nielsen scales the median catch value of all simulations. I scale all individual simulations
      simtmp<<-sim
      fun<-function(s,b){
          simtmp[b,]<<-scaleF(sim[b,], scale=s)
          simcat<-catch(simtmp[b,], nm=nm, cw=cw)
          m <- ifelse(length(x)==1,x,x[b])
          return(m-simcat)
      }
      ff <- sapply(1:nrow(sim), function(b) {
          check <- try(uniroot(fun, c(0,100),b)$root, silent=TRUE)
          if ('try-error' %in% class(check)) simtmp[b,] <- maxF(sim[b,]) #avoid imposibly high F's
          })
      sim <- simtmp
      sim <- t(apply(sim,1,maxF))
      return(sim)
  }

  doAve <- function(x,y)colMeans(x[rownames(x)%in%ave.years,,drop=FALSE])
  ave.sw <- doAve(fit$data$stockMeanWeight)
  ave.cw <- doAve(fit$data$catchMeanWeight)
  ave.mo <- doAve(fit$data$propMat)
  ave.nm <- doAve(fit$data$natMor)
  ave.lf <- doAve(fit$data$landFrac)
  ave.lw <- doAve(fit$data$landMeanWeight)
  ave.pm <- doAve(fit$data$propM)
  ave.pf <- doAve(fit$data$propF)
  ave.en <- doAve(fit$data$env)
  getThisOrAve <- function(x,y, ave){
    if(y %in% rownames(x)){
      ret <- x[which(rownames(x)==y),]
    }else{
      ret <- ave
    }
    ret
  }
  procVar<-getProcessVar(fit)
  simlist<-list()
  recs <- rep(NA, nosim)
  TAC <- rep(TAC.base,nosim)

  for(i in 0:ny){
    y<-year.base+i                       # starts with final year

    sw<-getThisOrAve(fit$data$stockMeanWeight, y, ave.sw) #if final year: finalyear value, otherwise average
    cw<-getThisOrAve(fit$data$catchMeanWeight, y, ave.cw)
    mo<-getThisOrAve(fit$data$propMat, y, ave.mo)
    nm<-getThisOrAve(fit$data$natMor, y, ave.nm)
    lf<-getThisOrAve(fit$data$landFrac, y, ave.lf)
    lw<-getThisOrAve(fit$data$landMeanWeight, y, ave.lw)
    pm<-getThisOrAve(fit$data$propM, y, ave.pm)
    pf<-getThisOrAve(fit$data$propF, y, ave.pf)
    dw<-getThisOrAve(fit$data$disMeanWeight, y, ave.pf)
    pr<-getThisOrAve(fit$data$propMat, y, ave.pf)
    en<-getThisOrAve(fit$data$env, y, ave.en)

    # abundance at the beginning of the year
    if(i!=0) recs <- getRec(sim, recpool, rec.meth, deterministic, i, simpara)
    sim <- t(sapply(1:nrow(sim), function(s) step(sim[s,], nm=nm, rec=recs[s], scale=1, inyear=(i==0)))) # simulations of next year
    if(i!=0){
      if(deterministic) procVar<-procVar*0
      sim <- sim + rmvnorm(nosim, mu=rep(0,nrow(procVar)), Sigma=procVar) # simulations of next year with process variance
    }

    # scale fishing pressure during the year
    if(i!=0){
        if(!is.na(fscale[i])){
            sim<-t(apply(sim, 1, scaleF, scale=fscale[i]))
        }

        if(!is.na(fval[i])){
            curfbar<-median(apply(sim, 1, fbar))
            adj<-fval[i]/curfbar
            sim<-t(apply(sim, 1, scaleF, scale=adj))
        }

        if(!is.na(catchval[i])){
            sim <- Csim(catchval[i])
        }

        if(!is.na(MP[i])){
            ## MPs
            if(is.numeric(MP[i])){
                TAC=MP[i]
            }else{
                TACfun <- match.fun(MP[i])
                args <- names(as.list(TACfun))
                args <- args[args != ""]
                L <- list(data=dats,parameters=paras, TAC.base=TAC, fit=fit,i=i)
                TAC <-do.call(TACfun, L[args])
            }

            ## IE  # still some work to do for the delayed.
            if(!is.na(IE[i])){
                 IEfun=match.fun(IE[i])
                 TAC=IEfun(TAC,i)
            }

            ## caps
            if(!is.null(capLower)) TAC[TAC<capLower]=capLower
            if(!is.null(capUpper)) TAC[TAC>capUpper]=capUpper

            TAC <- round(TAC,0)
            sim <-Csim(TAC,sim) # N depends on F of last year, so only new F needs to be created as N is already correct after step function
        }
    }
    # generate observations during the year if management procedure used (no matter which one)
    if(any(!is.na(MP))){

        #new parameters
        if(i==0) paras <- fit$pl
        if(i>0) paras <- steppar(paras)

        #new data
        if(i==0){
            dats <- list(fit$data)
        }else{
            obs <- newobs(fit, sim, simpara, ULpool,deterministic,nm)
            aux <- cbind(year=y,unique(fit$data$aux[,2:3]))
            if(i==1){
                dats <- lapply(1:nrow(sim),function(x) stepdat(dats[[1]],ob=obs[,,x],aux=aux,
                                                                   pr,sw,cw,nm,lf,dw,lw,pm,pf,en))
            }
            if(i>1){
                dats <- lapply(1:nrow(sim),function(x) stepdat(dats[[x]],ob=obs[,,x],aux=aux,
                                                                   pr,sw,cw,nm,lf,dw,lw,pm,pf,en))
            }
        }
    }

    # get derived quantities
    fbarsim <- apply(sim, 1, fbar)
    catchsim <- apply(sim, 1, catch, nm=nm, cw=cw)
    ssbsim <- apply(sim, 1, ssb, nm=nm, sw=sw, mo=mo, pm=pm, pf=pf)
    recsim <- exp(sim[,1])
    exploitsim <- catchsim/ssbsim
    ssbmsyratiosim <- ssbsim/refs[,'f40ssb']
    fmsyratiosim <- fbarsim/refs[,'f40']
    Umsyratiosim <- exploitsim/refs[,'f40U']
    simlist[[i+1]] <- list(sim=sim, fbar=fbarsim, catch=catchsim, ssb=ssbsim, rec=recsim, TAC=TAC, exploit=exploitsim, ssbmsyratio=ssbmsyratiosim,
                           fmsyratio=fmsyratiosim,Umsyratio=Umsyratiosim,year=y)
  }

  attr(simlist, "fit")<-fit
  class(simlist) <- "ccamforecast"

  collect <- function(x){
    quan <- quantile(x, c(.50,.025,.975))
    c(median=quan[1], low=quan[2], high=quan[3])
  }
  collectprob <- function(x){
      p1 <- signif(sum(simlist[[x]]$ssb>CZ)/nosim,2)
      p2 <- signif(sum(simlist[[x]]$ssb>HZ)/nosim,2)
      p3 <- signif(sum(simlist[[x]]$ssb>ssbf)/nosim,2)
      p4 <- signif(sum(simlist[[x]]$ssb>ssbf*2)/nosim,2)
      c(probCZ=p1, probHZ=p2,probGrowth=p3, probDouble=p4)
  }
  fbar <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$fbar))),3)
  rec <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$rec))))
  ssb <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$ssb))))
  catch <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$catch))))
  exploit <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$exploit))),3)
  ssbmsyratio <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$ssbmsyratio))),3)
  fmsyratio <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$fmsyratio))),3)
  Umsyratio <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$Umsyratio))),3)

  catchcumul <- round(t(apply(apply(do.call('rbind',lapply(simlist, function(x) x$catch)),2,cumsum),1,quantile,c(.50,.025,.975))))
  names(catchcumul) <- names(catch)
  TAC <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$TAC))))
  probs <- do.call(rbind, lapply(1:(ny+1), collectprob))

  tab <- cbind(fbar, rec,ssb,catch,exploit, ssbmsyratio,fmsyratio,Umsyratio,catchcumul,TAC,probs)
  rownames(tab) <- unlist(lapply(simlist, function(xx)xx$year))
  nam <- c("median","low","high")
  colnames(tab) <- c(paste0(rep(c("fbar:","rec:","ssb:","catch:","exploit:","ssbmsyratio:","fmsyratio:","Umsyratio:",'catchcumul:',"TAC:"), each=length(nam)), nam),colnames(probs))
  attr(simlist, "tab")<-tab
  shorttab<-t(tab[,grep("median",colnames(tab))])
  label <- paste(OMlabel,MPlabel,".")
  rownames(shorttab)<-sub(":median","",paste0(label,if(length(label)!=0)":",rownames(shorttab)))

  attr(simlist, "shorttab")<-shorttab
  attr(simlist, "OMlabel") <- OMlabel
  attr(simlist, "MPlabel") <- MPlabel
  attr(simlist, "parameters") <- parameters
  return(simlist)
}
