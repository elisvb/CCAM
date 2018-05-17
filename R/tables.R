##' Table helper
##' @param fit returned object from ccam.fit
##' @param what quoted name of what to extract
##' @param x rownames of table
##' @param trans function to be applied
##' @param ... extra arguments not currently used
##' @details ...
tableit <-function (fit, what, x=fit$data$years, trans=function(x)x,...){
    UseMethod("tableit")
}
##' @rdname tableit
##' @method tableit ccam
##' @export
tableit.ccam <- function (fit, what, x=fit$data$years, trans=function(x)x,...){
   idx<-names(fit$sdrep$value)==what
   y<-fit$sdrep$value[idx]
   ci<-y+fit$sdrep$sd[idx]%o%c(-2,2)
   ret<-trans(cbind(y,ci))
   rownames(ret)<-x
   colnames(ret)<-c("Estimate","Low","High")
   return(ret)
}

##' @rdname tableit
##' @method tableit ccamset
##' @export
tableit.ccamset <- function (fit, what, x=NULL, trans=function(x)x,...){
    if(is.null(x)) x <- unique(unlist(lapply(1:length(fit),function(x) fit[[1]]$data$years)))
    ret <- tableset(fit, tableit, what, x, trans,...)
    return(ret)
}

##' tableset
##' @param fit list of returned objects from ccam.fit
##' @param fun table function to be applied to each list element
##' @export
tableset <- function(fit, fun, ...){
    na <- 1:length(fit)
    tabs <- lapply(na,function(x) {
        tab <- fun(fit[[x]],...)
        tab <- cbind(tab,fit=x)
        return(tab)})
    ret <- do.call('rbind',tabs)
    return(ret)
}

##' SSB table
##' @param  fit ...
##' @param ... extra arguments not currently used
##' @details ...
##' @export
ssbtable<-function(fit,...){
    UseMethod("ssbtable")
}
##' @rdname ssbtable
##' @method ssbtable default
##' @export
ssbtable.default <- function(fit,...){
   ret<-tableit(fit, "logssb", trans=exp,...)
   return(ret)
}

##' TSB table
##' @param  fit ...
##' @param ... extra arguments not currently used
##' @details ...
##' @export
tsbtable<-function(fit,...){
    UseMethod("tsbtable")
}
##' @rdname tsbtable
##' @method tsbtable default
##' @export
tsbtable.default <- function(fit,...){
   ret<-tableit(fit, "logtsb", trans=exp,...)
   return(ret)
}

##' Fbar table
##' @param  fit ...
##' @param ... extra arguments not currently used
##' @details ...
##' @export
fbartable<-function(fit,...){
    UseMethod("fbartable")
}
##' @rdname fbartable
##' @method fbartable default
##' @export
fbartable.default <- function(fit,...){
   ret<-tableit(fit, "logfbar", trans=exp)
   return(ret)
}

##' Recruit table
##' @param  fit ...
##' @param ... extra arguments not currently used
##' @details ...
##' @export
rectable<-function(fit){
    UseMethod("rectable")
}
##' @rdname rectable
##' @method rectable default
##' @export
rectable.default <- function(fit,...){
   ret<-tableit(fit, "logR", trans=exp)
   return(ret)
}

##' Catch table
##' @param  fit ...
##' @param obs.show logical add a column with catch sum of product rowsums(C*W)
##' @param ... extra arguments not currently used
##' @details ...
##' @export
catchtable<-function(fit, obs.show=FALSE,...){
    UseMethod("catchtable")
}
##' @rdname catchtable
##' @method catchtable ccam
##' @export
catchtable.ccam <- function(fit, obs.show=FALSE,...){
   CW <- fit$data$catchMeanWeight
   xx <- as.integer(rownames(CW))
   ret <- tableit(fit, x=xx, "logCatch", trans=exp)
   if(obs.show){
     aux <- fit$data$aux
     logobs <- fit$data$logobs
     .goget <- function(y,a){
       ret <- exp(logobs[aux[,"fleet"]==1 & aux[,"year"]==y & aux[,"age"]==a])
       ifelse(length(ret)==0,0,ret)
      }
      sop<-rowSums(outer(rownames(CW), colnames(CW), Vectorize(.goget))*CW, na.rm=TRUE)
      ret<-cbind(ret,sop.catch=sop)
   }
   return(ret)
}
##' @rdname catchtable
##' @method catchtable ccamset
##' @export
catchtable.ccamset <- function(fit, obs.show=FALSE,...){
    return(tableset(fit, fun=catchtable, ...))
}

##' N table
##' @param  fit ...
##' @param ... extra arguments not currently used
##' @details ...
##' @export
ntable <- function(fit,...){
    UseMethod("ntable")
}
##' @rdname ntable
##' @method ntable ccam
##' @export
ntable.ccam <- function(fit,...){
   ret <- exp(t(fit$pl$logN))
   colnames(ret) <- fit$conf$minAge:fit$conf$maxAge
   rownames(ret) <- fit$data$years
   return(ret)
}
##' @rdname ntable
##' @method ntable ccamset
##' @export
ntable.ccamset <- function(fit,...){
    return(tableset(fit, fun=ntable, ...))
}

##' F-at-age table
##' @param  fit ...
##' @param ... extra arguments not currently used
##' @details ...
##' @export
faytable <- function(fit,...){
    UseMethod("faytable")
}
##' @rdname faytable
##' @method faytable ccam
##' @export
faytable.ccam <- function(fit,...){
   sel = c(invlogit(fit$pl$logitSel),1)[fit$conf$keyLogFsta[1,]+1]
   ret = outer(exp(fit$pl$logFy),sel)
   # idx <- fit$conf$keyLogFsta[1,]+2
   # ret <- cbind(NA,exp(t(fit$pl$logF)))[,idx]
   # ret[,idx==0] <- 0
   colnames(ret) <- fit$conf$minAge:fit$conf$maxAge
   rownames(ret) <- fit$data$years
   return(ret)
}
##' @rdname faytable
##' @method faytable ccam
##' @export
faytable.ccamset <- function(fit,...){
    return(tableset(fit, fun=faytable, ...))
}

##' parameter table
##' @param  fit ...
##' @param ... extra arguments not currently used
##' @details ...
##' @export
partable <- function(fit,...){
    UseMethod("partable")
}
##' @rdname partable
##' @method partable ccam
##' @export
partable.ccam <- function(fit,...){
  param <- coef(fit)
  nam <- names(param)
  dup <- duplicated(nam)
  namadd <- rep(0, length(nam))
  for (i in 2:length(dup)) {
    if(dup[i])namadd[i] <- namadd[i - 1] + 1
  }
  nam <- paste(nam, namadd, sep = "_")
  ret<-cbind(param, attr(param,"sd"))
  ex<-exp(ret[,1])
  lo<-exp(ret[,1]-2*ret[,2])
  hi<-exp(ret[,1]+2*ret[,2])
  ret<-cbind(ret,ex,lo,hi)
  colnames(ret)<-c("par", "sd(par)", "exp(par)", "Low", "High")
  rownames(ret)<-nam
  return(ret)
}
##' @rdname partable
##' @method partable ccamset
##' @export
partable.ccamset <- function(fit,...){
    return(tableset(fit, fun=partable, ...))
}

##' model table
##' @param fits A ccam fit as returned from the ccam.fit function, or a collection c(fit1, fit2, ...) of such fits
##' @param ... extra arguments not currently used
##' @details ...
##' @importFrom stats AIC pchisq
##' @export
modeltable <- function(fits,...){
    UseMethod("modeltable")
}
##' @rdname modeltable
##' @method modeltable ccam
##' @export
modeltable.ccam <- function(fits,...){
    modeltable(c(fits))
}
##' @rdname modeltable
##' @method modeltable ccamset
##' @export
modeltable.ccamset <- function(fits,...){
    if(!is.null(attr(fits,"fit"))){
      fits[[length(fits)+1]] <- attr(fits,"fit")
      fits <- fits[c(length(fits),1:(length(fits)-1))]
    }
    fits <- fits[!sapply(fits, is.null)]
    if(is.null(names(fits))){
        nam <- paste("M", 1:length(fits), sep="")
    }else{
        nam <-ifelse(names(fits)=="",paste("M", 1:length(fits), sep=""), names(fits))
    }
    logL <- sapply(fits, logLik)
    npar <- sapply(fits, function(f)attr(logLik(f),"df"))
    aic <- sapply(fits, AIC)
    res <- cbind("log(L)"=logL, "#par"=npar, "AIC"=aic)
    rownames(res) <- nam
    o <- 1:length(fits)
    if(length(fits)==2){
        o <- order(npar, decreasing=TRUE)
        if(npar[o[1]]>npar[o[2]]){
            df <- npar[o[1]]>npar[o[2]]
            D <- 2*(logL[o[1]]-logL[o[2]])
            P <- 1-pchisq(D,df)
            cnam <- paste0("Pval( ",nam[o[1]]," -> ",nam[o[2]], " )")
            res <- cbind(res, c(NA, P)[o])
            colnames(res)[ncol(res)] <- cnam
        }
    }
    return(res[o,,drop=FALSE])
}
##' @rdname modeltable
##' @method modeltable forecastset
##' @export
modeltable.forecastset <- function(fits,...){
    ret <- lapply(1:length(fit),function(x) attr(fit[[x]],"tab"))
    l <- unlist(lapply(ret, nrow))
    ret <- do.call('rbind',ret)
    nam <- names(fit)
    if(is.null(nam)) nam <- 1:length(fit)
    ret <- cbind(ret,run=rep(nam,l))
    return(ret)
}

##' extract predictions for a given variable
##' @param x A forecast object
##' @param ... extra arguments
##' @details ...
##' @export
extract <- function(x,...){
 UseMethod("extract")
}
##' @rdname extract
##' @method extract ccamforecast
##' @export
extract.ccamforecast <- function(x,what){
    tab <- attr(x,'tab')
    ret <- tab[,what]
    return(ret)
}
##' @rdname extract
##' @method extract forecastset
##' @export
extract.forecastset <- function(x,what){
    ret <- do.call('rbind',lapply(1:length(x),function(y) {d <- as.data.frame(extract(x[[y]],what=what))
                                                           d$id <- y
                                                           d$year <- as.numeric(rownames(d))
                                                           return(d)}))
    ret <- as.data.frame(ret)
    return(ret)
}

##' Yield per recruit calculation
##' @param fit the object returned from ccam.fit
##' @param Flimit Upper limit for Fbar
##' @param Fdelta increments on the Fbar axis
##' @param aveYears Number of years back to use when calculating averages (selection, weights, ...)
##' @param ageLimit Oldest age used (should be high)
##' @param ... extra arguments not currently used
##' @export
ypr<-function(fit, Flimit=2, Fdelta=0.01, aveYears=min(15,length(fit$data$years)), ageLimit=100,...){
    UseMethod("ypr")
}
##' @rdname ypr
##' @method ypr ccam
##' @export
ypr.ccam <- function(fit, Flimit=2, Fdelta=0.01, aveYears=min(15,length(fit$data$years)), ageLimit=100,rec.years=fit$data$years,deterministic=TRUE,simpara=NULL,...){
  last.year.used=max(fit$data$years)
  idxno<-which(fit$data$years==last.year.used)

  extend<-function(x,len=100){
    ret<-numeric(len)
    ret[1:length(x)]<-x
    ret[-c(1:length(x))]<-x[length(x)]
    ret
  }

  if(deterministic){
      ave.sl <- fit$pl$logitSel
  }else{
      if(is.null(simpara)) {
          simpara <- rmvnorm(1, mu=fit$sdrep$par.fixed, Sigma=fit$sdrep$cov.fixed)
          names(simpara) <- names(fit$sdrep$par.fixed)
      }
      ave.sl <- simpara[which(names(simpara)=='logitSel')]
  }

  ave.sl<-c(invlogit(ave.sl),1)[fit$conf$keyLogFsta[1,]+1]

  ave.sw<-colMeans(fit$data$stockMeanWeight[(idxno-aveYears+1):idxno,,drop=FALSE])
  ave.cw<-colMeans(fit$data$catchMeanWeight[(idxno-aveYears+1):(idxno-1),,drop=FALSE])
  ave.pm<-colMeans(fit$data$propMat[(idxno-aveYears+1):idxno,,drop=FALSE])
  ave.nm<-colMeans(fit$data$natMor[(idxno-aveYears+1):idxno,,drop=FALSE])
  ave.lf<-colMeans(fit$data$landFrac[(idxno-aveYears+1):(idxno-1),,drop=FALSE])
  ave.cw.land<-colMeans(fit$data$landMeanWeight[(idxno-aveYears+1):(idxno-1),,drop=FALSE])

  N<-numeric(ageLimit)
  N[1]<-1.0
  M<-extend(ave.nm)
  sw<-extend(ave.sw)
  cw<-extend(ave.cw.land)
  pm<-extend(ave.pm)
  lf<-extend(ave.lf)

  deltafirst <- 0.00001
  delta <- Fdelta
  scales<-c(0, deltafirst, seq(0.01, Flimit, by=delta))
  yields<-numeric(length(scales))
  ssbs<-numeric(length(scales))
  for(i in 1:length(scales)){
    scale<-scales[i]
    F<-extend(ave.sl*scale)
    Z<-M+F
    for(a in 2:length(N)){
      N[a]<-N[a-1]*exp(-Z[a-1])  #there is no plus group in here and pe is ignored
    }
    C<-F/Z*(1-exp(-Z))*N*lf
    Y<-sum(C*cw)
    yields[i]<-Y  #ypr
    ssbs[i]<-sum(N*pm*sw)  #spr
  }

  fmaxidx<-which.max(yields)
  fmax<-scales[fmaxidx]

  deltaY<-diff(yields)
  f01idx<-which.min((deltaY/delta-0.1*deltaY[1]/deltafirst)^2)+1
  f01<-scales[f01idx]

  f30spridx<-which.min((ssbs-0.30*ssbs[1])^2)+1
  f30<-scales[f30spridx]
  f35spridx<-which.min((ssbs-0.35*ssbs[1])^2)+1
  f35<-scales[f35spridx]
  f40spridx<-which.min((ssbs-0.40*ssbs[1])^2)+1
  f40<-scales[f40spridx]

  rec <- rectable(fit)
  meanrec <- mean(rec[which(rownames(rec) %in% rec.years),1])
  f30ssb <- ssbs[which.min((scales-f30)^2)]*meanrec
  f35ssb <- ssbs[which.min((scales-f35)^2)]*meanrec
  f40ssb <- ssbs[which.min((scales-f40)^2)]*meanrec
  f30yield <- yields[which.min((scales-f30)^2)]*meanrec
  f35yield <- yields[which.min((scales-f35)^2)]*meanrec
  f40yield <- yields[which.min((scales-f40)^2)]*meanrec
  f30U <- f30yield/f30ssb
  f35U <- f35yield/f35ssb
  f40U <- f40yield/f40ssb

  fbarlab <- substitute(bar(F)[X - Y], list(X = fit$conf$minAge, Y = fit$conf$maxAge))
  ret<-list(fbar=scales, ssb=ssbs, yield=yields, fbarlab=fbarlab, f30=f30, f35=f35, f40=f40, f01=f01, fmax=fmax,f30Idx=f30spridx,
            f35Idx=f35spridx, f40Idx=f40spridx,f01Idx=f01idx, fmaxIdx=fmaxidx, f30ssb=f30ssb, f35ssb=f35ssb, f40ssb=f40ssb,
            f30yield=f30yield, f35yield=f35yield, f40yield=f40yield, f30U=f30U, f35U=f35U, f40U=f40U)
  class(ret)<-"ccamypr"
  return(ret)
}

##' @rdname ypr
##' @method ypr ccamset
##' @export
ypr.ccamset <- function(fit,what=names(ypr(fit[[1]])),...){
    ret <- tableset(fit, fun=ypr,...)
    if(!identical(what,names(ypr(fit[[1]])))){
        elem <- names(ypr(fit[[1]]))
        take <- which(elem %in% what)
        take <- as.vector(sapply(take, function(x) seq(from=x,by=length(elem),length.out = length(fit))))

        ret <- ret[take]
        names(ret) <- rep(elem,length(fit))[take]
    }
    return(ret)
}

