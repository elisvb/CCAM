##' rmvnorm helper function to draw multivariate normal samples
##' @param n the number of samples.
##' @param mu the mean vector.
##' @param Sigma a positive-definite symmetric matrix specifying the covariance matrix.
##' @param seed optional seed
##' @details Generates samples via the Cholesky decomposition, which is less platform dependent than eigenvalue decomposition.
##' @return If n = 1 a vector of the same length as mu, otherwise an n by length(mu) matrix with one sample in each row.
##' @export
rmvnorm <- function(n = 1, mu, Sigma, seed){
  p <- length(mu)
  if(!all(dim(Sigma) == c(p, p))){
    stop("incompatible arguments")
  }
  idx <- diag(Sigma) > .Machine$double.xmin
  L <- matrix(0,p,p)
  if(any(idx)){
    L[idx,idx] <- chol(Sigma[idx,idx])
  }
  if (!missing(seed)) set.seed(seed)
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
##' @param rec.meth methode to project recruitment. Autocorrelation can be specified as attribute ('AC').
##' @param rec.scale scale recruitment (default is 1)
##' @param label optional label to appear in short table
##' @param overwriteSelYears if a vector of years is specified, then the average selectivity of those years is used (not recommended)
##' @param deterministic option to turn all process noise off (not recommended, as it will likely cause bias)
##' @param CZ critical zone (to calculate probability ssb above)
##' @param HZ Healty zone (to calculate probability ssb above)
##' @param IE implementation error type. See avail('IE')
##' @param capLower lower limit of TAC (before IE)
##' @param capUpper upper limit of TAC (before IE)
##' @param UL.years vector of years to use to resample upper/lower catch bound ranges from
##' @param TAC.base TAC last year
##' @param bio.scale names list with multipliers used to scale future bio values (nm, en, pr,sw,cw,lf,dw,lw,pm,pf)
##' @param deadzone ssb value below which there is no recruitment anymore and the stock goes extinct. Avoids 'dead' stock popping back to live because of recruitment peak.
##' @param Flim Upper limit of F value. defaults to 3.
##' @param fleet fleet for which the future observations are extracted
##' @details There are four ways to specify a scenario; fscale, catchval, fval and MP. There are also multiple recruitment options;
#' \enumerate{
#'   \item parametric BH estimated within the model. If attribute 'AC' not NULL, use autocorrelation estimated outside the model.
#'   \item mean value of rec.pool with sd. If attribute 'AC' not NULL, use autocorrelation estimated outside the model.
#'   \item sample value from rec.pool
#'   \item sample value backwards (for first year one of the two last estimted values)
#'   \item idem as 3 but sample RPS
#'   \item idem as 4 but sample RPS
#' }
#' Set.seed is used with the number of simulations to facilitate comparison between different runs (using e.g. different catchvalues/MPs) for;
#' \enumerate{
#'   \item last year state
#'   \item observation residuals
#'   \item future data (M, proportion mature, weight-at-age, etc.)
#' }
##' @return an object of type ccamforecast
##' @importFrom stats median uniroot quantile
##' @export
forecast <- function(fit, fscale=NULL, catchval=NULL, fval=NULL, MP=NULL,nosim=1000, year.base=max(fit$data$years),
                     ave.years=max(fit$data$years)+(-4:0), rec.years=max(fit$data$years)+(-9:0),rec.meth=4, rec.scale=1, MPlabel=NULL,
                     OMlabel=NULL,overwriteSelYears=NULL, deterministic=FALSE,CZ=0,HZ=0,IE=0,
                     capLower=NULL,capUpper=NULL,UL.years=max(fit$data$years)+(-4:0),TAC.base=0,bio.scale=NULL,verbose=FALSE,deadzone=0, Flim=3,fleet=3){

  #.. check/clean/prepare input.................................................................................................
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

  if(!all(rowSums(!is.na(cbind(fscale, catchval, fval, MP)))==1)){
      stop("For each forecast year exactly one of fscale, catchval, fval or MP must be specified (all others must be set to NULL)")
  }

  ny <-length(fscale)

  MPs <- avail('MP')
  if(all(!is.na(MP)) & any(!MP %in% MPs) & !is.numeric(MP)) stop("Undefined MP (see avail('MP'))")

  #dummy <- lapply(as.list(IE),function(x) if(!x %in% avail('IE') & !is.numeric(x)) stop("Undefined IE (see avail('IE'))"))

  if(!all(fit$data$fleetTypes %in% c(3,6))) warning('MPs based on fleettypes different from 3 and 6 need more code')

  if(!is.null(overwriteSelYears)){
    fromto <- fit$conf$fbarRange-(fit$conf$minAge-1)
    Ftab <- faytable(fit)
    fixedsel <- colMeans(Ftab[as.integer(rownames(Ftab))%in%overwriteSelYears,,drop=FALSE])
    fixedsel <- fixedsel/mean(fixedsel[fromto[1]:fromto[2]])
  }

  # for when model based MP: track the number that are fitted and how many converged
  nfit <- 0
  nfail <- 0
  myenv <- environment()

  #.. functions.................................................................................................
  getF <<- function(x){
    idx <- fit$conf$keySel[1,]+1
    nsize <- length(idx)
    ret <- exp(x[,nsize+idx,drop=FALSE])
    ret[idx==0] <- 0
    if(!is.null(overwriteSelYears)){
      fromto <- fit$conf$fbarRange-(fit$conf$minAge-1)
      thisfbar<-rowMeans(ret[,fromto[1]:fromto[2],drop=FALSE])
      ret<-outer(thisfbar,fixedsel)
    }
    ret
  }

  fbar <<-function(x){
    fromto <- fit$conf$fbarRange-(fit$conf$minAge-1)
    rowMeans(getF(x)[,fromto[1]:fromto[2],drop=FALSE])
  }

  getN <<- function(x){
    idx <- fit$conf$keySel[1,]+1
    nsize <- length(idx)
    ret <- exp(x[,1:nsize,drop=FALSE])
    ret
  }

  getState <<- function(N,F){
    x <- log(cbind(N,F))
    x
  }

  getProcessVar <- function(fit,simpara){
    N <- length(fit$conf$keyVarLogN)*2
    cov <- array(0,dim=c(N,N,nosim))
    for(i in 1:nosim){
        sdN <- exp(simpara[i,which(colnames(simpara)=='logSdLogN')][fit$conf$keyVarLogN+1])
        sdN[1]<-0
        nN <- length(sdN)
        sdF <- exp(simpara[i,which(colnames(simpara)=='logSdLogFsta')][fit$conf$keyVarF+1])
        nF <- length(sdF)
        corr <- diag(nF)
        cov[,,i] <- matrix(0,nrow=nN+nF,ncol=nN+nF)
        cov[1:nN,1:nN,i] <- diag(sdN^2)
        cov[nN+1:nF,nN+1:nF,i] <- (sdF%*%t(sdF))*corr
    }
    return(cov)
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

  getRec <- function(sim, recpool, rec.meth, deterministic, i, simpara, deadzone){

      Rautocorr <- function(r,mu,sigma,rho,seed=NULL){
          eta <- r-mu+sigma^2/2
          if(!is.null(seed)) set.seed(length(r)+seed)
          delta <- rnorm(length(r),0,sigma)
          etanew <- rho*eta+delta*sqrt(1-rho^2)
          R <- exp(mu+etanew-sigma^2/2)
          return(R)
      }
      BH <- function(x, a, b){a+log(x)-log(1.0+exp(b)*x)}
      RI <- function(x, a, b){a+log(x)-exp(b)*x}  # not implemented yet

      ssb.before <- ssb(sim, nm=nm, sw=sw, mo=mo, pm=pm, pf=pf)

      switch(rec.meth,
             "one" <- { # see Johnson 2016  (BH, lag 1 autocorrelation)
                 rho <- attr(rec.meth,'AC')
                 if(is.null(rho)) rho <- acf(recpool,plot=FALSE)$acf[2,1,1]
                 mu <- BH(ssb.before,a=simpara[,'rec_loga'],b=simpara[,'rec_logb'])
                 sigma <- exp(simpara[,'logSdLogN'])
                 if(deterministic) sigma <- 0
                 recs <- Rautocorr(sim[,1],mu,sigma,rho,i)
             },
             "two" <- { # see Johnson 2016 (mean, lag 1 autocorrelation)
                 rho <- attr(rec.meth,'AC')
                 if(is.null(rho)) rho <- acf(recpool,plot=FALSE)$acf[2,1,1]
                 mu <- mean(log(recpool))
                 sd.option <- attr(rec.meth,'sd.option')
                 if(is.null(sd.option)){
                    sigma <- sd(log(recpool))
                 }
                 if(sd.option=='ci'){
                     l <- log(rectab[rownames(rectab)%in%rec.years,])
                     sigma <- max((l[,1]-l[,2])/1.96)
                 }
                 if(deterministic) sigma <- 0
                 recs <- Rautocorr(sim[,1],mu,sigma,rho,i)
             },
             "three" <- { # sample recs (original SAM method, no autocorrelation)
                 if(deterministic){
                     recs <- rep(mean(recpool),nrow(sim))
                 }else{
                     recs <- sample(recpool, nrow(sim), replace=TRUE)
                  }
             },
             "four" <- { # sample recs, trailing
                 ref <- tail(recpool,i+2)
                 if(deterministic){
                     recs <- rep(mean(ref),nrow(sim))
                 }else{
                     recs <- sample(ref, nrow(sim), replace=TRUE)
                 }
             },
             "five" <- { # sample RPS (suggestion Daniel, little autocorrelation, horrible results)
                 if(deterministic){
                     recs <- rep(mean(recpool),nrow(sim))*ssb.before
                 }else{
                     recs <- unname(sample(recpool,nrow(sim), replace=TRUE)*ssb.before)
                 }
             },
             "six" <- { # sample RPS, trailing (suggestion Daniel, little autocorrelation, horrible results)
                 ref <- tail(recpool,i+2)
                 if(deterministic){
                     recs <- rep(mean(ref),nrow(sim))*ssb.before
                 }else{
                     recs <- unname(sample(ref,nrow(sim), replace=TRUE)*ssb.before)
                 }
             },
             {
                 stop('rec.meth should be 1 to 6') #block sampling would be another one
             }
      )
      deadsim[which(ssb.before<deadzone)] <<- TRUE
      recs[deadsim]=1
      return(recs)
  }

  step <- function(x, nm, scale, rec){
    F <- getF(x)
    N <- getN(x)
    Z <- F+nm
    n <- ncol(N)
    N <- cbind(rec,N[,-n]*exp(-Z[,-n])+cbind(matrix(0,ncol=n-2,nrow=nrow(x)),N[,n]*exp(-Z[,n])))
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

  maxF <- function(x, Flim=3){
      F <- getF(x)
      F[F==0] <- 0.00000001

      limit <- function(x,Flim){
          if(any(x>Flim)){
              x <- x/max(x)
              x <- x*Flim
          }
          return(x)
      }
      F <- t(apply(F,1,limit,Flim))
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
    return(rowSums(cw*C))
  }

  ssb <<- function(x, nm, sw, mo, pm, pf){
    F <- getF(x)
    N <- getN(x)*exp(-pm*nm-pf*F)
    return(rowSums(N*mo*sw))
  }

  ssb0 <<- function(x, sw, mo){
      F <- getF(x)
      N <- getN(x)
      return(rowSums(N*mo*sw))
  }

  tep <<- function(x, nm, mo, pm, pf, pfem, fec){
      F <- getF(x)
      N <- getN(x)*exp(-pm*nm-pf*F)
      TEP <- rowSums(N*mo*pfem*fec)
      return(TEP)
  }

  fsb <<- function(x, cw, nm){
      F <- getF(x)
      Z <- F+nm
      fsb <- F/rowSums(F)*getN(x)*exp(-Z/2)*cw
      return(rowSums(fsb))
  }

  getThis <- function(x,y){
      x <- fit$data[[x]]
      x <- x[which(rownames(x)==y),]
      return(x)
  }

  Csim <- function(x,sim,Flim=3){  #Nielsen scales the median catch value of all simulations. I scale all individual simulations
      simtmp<<-sim
      fun<-function(s,b){
          simtmp[b,]<<-scaleF(sim[b,,drop=FALSE], scale=s)
          simcat<-catch(simtmp[b,,drop=FALSE], nm=nm[b,,drop=FALSE], cw=cw[b,,drop=FALSE])
          m <- ifelse(length(x)==1,x,x[b])
          return(m-simcat)
      }
      ff <- sapply(1:nrow(sim), function(b) {
          ch <- try(uniroot(fun, c(0,100),b)$root, silent=TRUE)
          if ('try-error' %in% class(ch)) simtmp[b,] <- maxF(sim[b,,drop=FALSE],Flim) #avoid imposibly high F's
      })
      sim <- simtmp
      sim <- maxF(sim,Flim)
      return(sim)
  }

  doNew <- function(x,ave.years,nosim,deterministic){
      ra <- apply(x,2,range)
      x <- x[rownames(x)%in%ave.years,,drop=FALSE]
      n <- ncol(x)
      xl <- tail(x,1)
      if(ncol(x)>=nrow(x)){warning('An error might occur because the #age classes > ave.years')}
      if(!deterministic){ #last value plus variance (zero variance if cov is positive definite because last ages are too similar)
          for(i in 0:ncol(x)){
              ix <-head(1:n,10-i)
              works <- try(rmvnorm(ny*nosim,rep(0,n-i),cov(x[,ix,drop=FALSE]),seed=nosim), silent=TRUE)
              if(!'try-error' %in% class(works)){
                  delta <- rmvnorm(ny*nosim,rep(0,n-i),cov(x[,ix,drop=FALSE]),seed=nosim)
                  if(i>0) delta <- cbind(delta,matrix(0,ncol=i,nrow=nrow(delta)))
                  break
              }
          }
          xnew <- sweep(delta,2,xl,'+')
          xnew <- sweep(xnew,2,ra[1,],pmax) #keep values between observed, mainly because of propmature (don't get above 1)
          xnew <- sweep(xnew,2,ra[2,],pmin)
      }else{ #average last years
          xnew <- matrix(colMeans(x),nrow=ny*nosim,ncol=n,byrow=TRUE)
      }
      arr <- array(xnew, dim=c(nosim,ny,n),dimnames = list(nosim=1:nosim,
                                                           year=max(as.numeric(dimnames(x)[[1]]))+1:ny,
                                                           age=dimnames(x)[[2]]))
      return(arr)
  }

  getThisOrNew <- function(x, year, new){
      if(y %in% rownames(x)){
          ret <- x[which(rownames(x)==y),]
          ret <- matrix(ret,nrow=dim(new)[1],ncol=dim(new)[3],byrow=TRUE)
      }else{
          ret <- new[,as.character(year),]
      }
      ret
  }

  pair<-function(x,y){
      0.5*(x+y)*(x+y+1)+x
  }

  #.. create pools for recruitment and observation upper and lower limits.................................................................................................
  # recruitment
  if(rec.meth %in% c(5,6)){
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

  #.. get final state and parameter values .................................................................................................
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

  names(est) <- gsub('par.','',names(est))
  rmsel <- which(names(est)=='logitSel') #high colinearity with lastF so cov matrix not positive definite.
  if(length(rmsel)>0) {est <- est[-rmsel];cov <- cov[-rmsel,-rmsel]}

  idxN <- (fit$conf$minAge-fit$conf$minAge+1):(fit$conf$maxAge-fit$conf$minAge+1)
  idxF <-  fit$conf$keySel[fit$data$noYears,]+fit$conf$maxAge-fit$conf$minAge+2
  idxP <- seq(from=max(idxF)+1,by=1,length=length(est)-length(c(idxN,idxF)))
  idx <- c(idxN,idxF,idxP)
  rm <- !duplicated(idx)

  sim <- rmvnorm(nosim, mu=est[rm], Sigma=cov[rm,rm], seed=nosim) #does not work if hessian not positive definite (check sdreport)

  simpara <- sim[,idxP]
  colnames(simpara) <- gsub('par.','',names(est[rm])[idxP])

  sim <- sim[,c(idxN,idxF)]

  #.. get last year reference points .................................................................................................
  #refs <- t(apply(simpara,1,function(x) unlist(ypr(fit)[c('f40','f40ssb','f40U')]))) # to try stochastically
  refs <- matrix(unlist(ypr(fit)[c('f40','f40ssb','f40U')]),ncol=3,nrow=nosim,byrow=TRUE,dimnames=list(sim=1:nosim,ref=c('f40','f40ssb','f40U')))

  #.. get process variance................................................................................................
  # precalculated so for every forecast the same, if identical number of simulations are used
  procVar <- getProcessVar(fit, simpara)
  if(deterministic) procVar[] <- 0
  procsim <- sapply(1:ny,function(y) t(sapply(1:nosim,function(x) rmvnorm(1,mu=rep(0,nrow(procVar)), Sigma=procVar[,,x],seed=nosim+pair(y,x)))),simplify='array')

  #.. get observation residuals.................................................................................................
  # precalculated so for every forecast the same, if identical number of simulations are used
  if(any(!is.na(MP))){
      sds <- cbind(0.0001,exp(simpara[,which(colnames(simpara)=='logSdLogObs')]))
      aux <- unique(fit$data$aux[,2:3])
      aux[aux<0] <- 1
      idx <- apply(aux,1,function(x){fit$conf$keyVarObs[x[1],x[2]]})+2
      sds <- sds[,idx]
      if(deterministic) sds[] <- 0
      nobs <- nrow(aux)
      ressim <- array(dim=c(nobs,2,nosim,ny))
      for(y in 1:ny){
          for(i in 1:nosim){
              set.seed(pair(y,i))
              ressim[,1,i,y] <- rnorm(nobs,rep(0,nobs),sds[i,])
              ce <- which(fit$conf$obsLikelihoodFlag=='CE')
              if(length(ce)>0){
                  idxce <- which(aux[,1]==ce)
                  if(deterministic){
                      set.seed(pair(y,i))
                      ressim[idxce,1,i,y] <- log(mean(ULpool[,1],1))
                      ressim[idxce,2,i,y] <- log(mean(ULpool[,2],1))
                  }else{
                      set.seed(pair(y,i))
                      ressim[idxce,1,i,y] <- log(sample(ULpool[,1],1))
                      ressim[idxce,2,i,y] <- log(sample(ULpool[,2],1))
                  }
              }
          }
      }
  }

  #.. get implementation errors (if SSB independent, otherwise 0).................................................................................................
  IElist <- lapply(as.list(IE),function(x){
        if(grepl('IEindep',x)){                           #IEs independent of SSB that can be calculated directly
            IEfun <- match.fun(as.character(x))
            IEmatrix <- IEfun(nosim,ny,seed=nosim)
            return(IEmatrix)
        }
        if(is.numeric(x)){                              #constant deterministic IEs
            return(IEconstant(nosim,ny,x))
        }
      return(IEconstant(nosim,ny))                      #0
    })
  names(IElist) <- IE


  #.. get future data.................................................................................................
   new.sw <- doNew(fit$data$stockMeanWeight,ave.years,nosim,deterministic)
   new.sw0 <- doNew(fit$data$stockStartWeight,ave.years,nosim,deterministic)
   new.cw <- doNew(fit$data$catchMeanWeight,ave.years,nosim,deterministic)
   new.mo <- doNew(fit$data$propMat,ave.years,nosim,deterministic)
   new.nm <- doNew(fit$data$natMor,ave.years,nosim,deterministic)
   new.lf <- doNew(fit$data$landFrac,ave.years,nosim,deterministic)
   new.lw <- doNew(fit$data$landMeanWeight,ave.years,nosim,deterministic)
   new.pm <- doNew(fit$data$propM,ave.years,nosim,deterministic)
   new.pf <- doNew(fit$data$propF,ave.years,nosim,deterministic)
   new.dw <- doNew(fit$data$disMeanWeight,ave.years,nosim,deterministic)
   new.pfem <- doNew(fit$data$propFemale,ave.years,nosim,deterministic)
   new.fec <- doNew(fit$data$fec,ave.years,nosim,deterministic)
   new.en <- doNew(fit$data$env,ave.years,nosim,deterministic)

   #.. Start forecast.................................................................................................
   if(verbose)print('start predicting')
  deadsim <<- rep(FALSE,nosim)  #keep track of simulations that went extinct
  TAC <- rep(TAC.base,nosim)
  Ctrue <- TAC
  simlist <- list()

  for(i in 0:ny){
    counteryear <<- i
    y<-year.base+i
    if(verbose) print(y)

    # get the data values (for the last year the actual data and otherwise the pre-calculated predictions)
    sw<-getThisOrNew(fit$data$stockMeanWeight, y, new.sw)
    sw0<-getThisOrNew(fit$data$stockStartWeight, y, new.sw0)
    cw<-getThisOrNew(fit$data$catchMeanWeight, y, new.cw)
    mo<-getThisOrNew(fit$data$propMat, y, new.mo)
    nm<-getThisOrNew(fit$data$natMor, y,new.nm)
    lf<-getThisOrNew(fit$data$landFrac, y, new.lf)
    lw<-getThisOrNew(fit$data$landMeanWeight, y, new.lw)
    pm<-getThisOrNew(fit$data$propM, y, new.pm)
    pf<-getThisOrNew(fit$data$propF, y, new.pf)
    dw<-getThisOrNew(fit$data$disMeanWeight, y, new.dw)
    pfem<-getThisOrNew(fit$data$propFemale, y, new.pfem)
    fec<-getThisOrNew(fit$data$fec, y, new.fec)
    en<-getThisOrNew(fit$data$env, y, new.en)

    # if not the final year get new data, state and observations
    if(i!=0) {
        ## scale any of the above values (weight, propmat, env, nm, etc.). variability is kept.
        if(!is.null(bio.scale)){
            envir <- environment()
            dummy <- lapply(1:length(bio.scale),function(x){
                bio <- get(names(bio.scale)[x])*bio.scale[[x]]
                assign(names(bio.scale)[x],bio,envir = envir)
            })
        }

        ## new abundance and F
        recs <- getRec(sim, recpool, rec.meth, deterministic, i, simpara, deadzone)
        recs <- recs*rec.scale
        sim <- step(sim, nm=nm, rec=recs, scale=1)
        sim <- sim + procsim[,,i]

        ## scale fishing pressure during the year (could be skipped if stock is considered dead...)
        # F based
        if(!is.na(fscale[i])){
            sim <- scaleF(sim,scale=fscale[i])
        }

        if(!is.na(fval[i])){
            adj <- fval[i]/fbar(sim)
            sim <- scaleF(sim,scale=adj)
        }
        # C based
        if(!is.na(catchval[i])){
                TAC <- rep(catchval[i],nosim)
        }
        if(!is.na(MP[i])){
                TACfun <- match.fun(MP[i])
                args <- names(as.list(TACfun))
                args <- args[args != ""]
                L <- list(data=datasim, parameters=parasim, conf=confsim, TAC.base=TAC, fit=fit,i=i)
                TAC <-do.call(TACfun, L[args])
                TAC <- round(TAC,0)
        }

        if(!is.null(capLower)) TAC[TAC<capLower] <- capLower
        if(!is.null(capUpper)) TAC[TAC>capUpper] <- capUpper

        if(any(grepl('IEdep',unlist(IE)))){
            IEfun <- match.fun(as.character(IE[grep('IEdep',IE)]))
            IElist[[grep('IEdep',unlist(IE))]] <- IEfun(nosim,ny,seed=nosim,ssb=ssb0(sim, sw=sw, mo=mo),TAC=TAC,i=i,past=IElist[[grep('IEdep',unlist(IE))]]) #update IE
        }

        Ctrue <- TAC + colSums(do.call('rbind',lapply(IElist,function(x) x[,i]))) # add sum of all IEs
        Ctrue <- round(Ctrue,0)

        if(!is.na(catchval[i]) | !is.na(MP[i])){
            sim <- Csim(Ctrue,sim,Flim)
        }

        #sim2[order(as.numeric(rownames(sim2))),] #if MPs only on alive ones
    }



    # generate data/parameters during the year if management procedure used next year (no matter which one)
    if(any(!is.na(MP))){
        if(i==0){
            parasim <- fit$pl
            confsim <- fit$conf
            datasim <-  list(fit$data)
        }else{
            parasim <- steppar(parasim)
            confsim <- stepconf(confsim)
            obssim <- newobs(fit, sim, simpara, ULpool, deterministic, nm=nm, sw=sw, mo=mo, pm=pm, pf=pf, cw=cw, pfem=pfem, fec=fec)
            obssim <- obssim+ressim[,,,i]
            datasim <- lapply(1:nrow(sim),function(x){stepdat(datasim[[ifelse(i==1,1,x)]],
                                                              obs=obssim[,,x],
                                                              aux=cbind(year=y,unique(fit$data$aux[,2:3])),
                                                              mo[x,],sw[x,],sw0[x,],cw[x,],nm[x,],lf[x,],dw[x,],
                                                              lw[x,],pm[x,],pf[x,],pfem[x,],fec[x,],en[x])})

        }
    }

    # get derived quantities
    if(verbose) print(nrow(sim))
    fbarsim <- fbar(sim)
    catchsim <- catch(sim, nm=nm, cw=cw)
    ssbsim <- ssb(sim, nm=nm, sw=sw, mo=mo, pm=pm, pf=pf)
    ssb0sim <- ssb0(sim, sw=sw0, mo=mo)
    recsim <- exp(sim[,1])
    ssbmsyratiosim <- ssbsim/refs[,'f40ssb']
    fmsyratiosim <- fbarsim/refs[,'f40']
    CZ <- refs[,'f40ssb']*0.4
    HZ <- refs[,'f40ssb']*0.8
    IEperc <- (Ctrue-TAC)/Ctrue*100
    age <- rowSums(sweep(getN(sim)  ,2,1:10,'*'))/rowSums(getN(sim))
    if(any(!is.na(MP))){
        index <- unlist(lapply(datasim,function(x){exp(x$logobs[x$idx1[fleet,ncol(x$idx1)]+1,1])})) #values of last year
    }else{
        index <-NA
    }
    simlist[[i+1]] <- list(sim=sim, fbar=fbarsim, catch=catchsim, ssb=ssbsim, rec=recsim, age=age, TAC=TAC, ssbmsyratio=ssbmsyratiosim,
                           fmsyratio=fmsyratiosim,CZ=CZ,HZ=HZ,year=y,extinctions=length(deadsim[deadsim]),ssb0=ssb0sim,
                           IEperc=IEperc,index=index)

    if(verbose)print(tail(data.frame(size=sort( sapply(ls(),function(x){object.size(get(x))}),decreasing=TRUE)*1e-9,unit='GB'),10))
  }
  if(verbose)print('summarise...')

  attr(simlist, "fit")<-fit
  class(simlist) <- "ccamforecast"
  out <<- simlist

  collect <- function(x,...){
    quan <- quantile(x, c(.50,.025,.975),...)
    c(median=quan[1], low=quan[2], high=quan[3])
  }
  collectprob <- function(x){
      probCZ <- signif(sum(simlist[[x]]$ssb>simlist[[x]]$CZ)/nosim,2)
      probHZ <- signif(sum(simlist[[x]]$ssb>simlist[[x]]$HZ)/nosim,2)
      c(probCZ=probCZ, probHZ=probHZ)
  }
  fbar <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$fbar))),3)
  rec <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$rec))))
  ssb <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$ssb))))
  ssb0 <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$ssb0))))
  age <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$age))))
  catch <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$catch))))
  ssbmsyratio <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$ssbmsyratio))),3)
  fmsyratio <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$fmsyratio))),3)
  Umsyratio <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$Umsyratio))),3)
  extinctions <- round(do.call(rbind, lapply(simlist, function(xx) xx$extinctions/nosim*100)))
  IEperc <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$IEperc,na.rm=TRUE))))

  catchcumul <- round(t(apply(apply(do.call('rbind',lapply(simlist, function(xx) xx$catch)),2,cumsum),1,quantile,c(.50,.025,.975))))
  names(catchcumul) <- names(catch)

  TACrel <- do.call('rbind',lapply(simlist, function(x) x$TAC))
  TACrel <- TACrel[-1,]/TACrel[-nrow(TACrel),]
  TACrel <- rbind(NA,TACrel)
  TACrel <- t(apply(TACrel,1, quantile, c(.50,.025,.975),na.rm=T))

  TAC <- round(do.call(rbind, lapply(simlist, function(xx)collect(xx$TAC))))
  probs <- do.call(rbind, lapply(1:(ny+1), collectprob))

  tab <- cbind(fbar, rec,ssb,catch,age,ssb0, ssbmsyratio,fmsyratio,catchcumul,TAC,TACrel,IEperc,probs,extinctions)
  rownames(tab) <- unlist(lapply(simlist, function(xx)xx$year))
  nam <- c("median","low","high")
  colnames(tab) <- c(paste0(rep(c("fbar:","rec:","ssb:","catch:",'age:',"ssb0:","ssbmsyratio:","fmsyratio:",'catchcumul:',"TAC:","TACrel:",'IEperc'), each=length(nam)), nam),colnames(probs),'percExtinct')
  attr(simlist, "tab")<-tab
  shorttab<-t(tab[,grep("median",colnames(tab))])
  label <- paste(OMlabel,MPlabel,".")
  rownames(shorttab)<-sub(":median","",paste0(label,if(length(label)!=0)":",rownames(shorttab)))

  attr(simlist, "shorttab") <- shorttab
  attr(simlist, "OMlabel") <- OMlabel
  attr(simlist, "MPlabel") <- MPlabel
  attr(simlist, "parameters") <- parameters
  attr(simlist, "IE") <- lapply(IElist,function(x) {ies <- rbind(apply(x,2,quantile, c(.50,.025,.975,0.75,0.25)),x)})
  attr(simlist, "MPmessage") <- paste0('%model fails:',nfail/nfit*100)
  return(simlist)
}
