##' continuous ratio logit transformation of catch-at-age matrix
##' @param x catch-at-age vector, dataframe or matrix (rows=ages, columns=years).
##' @details (c) NOEL CADIGAN, adapted by EVB for use on vectors
##' @export

crlTransform = function(x){

    trans=function(x){
    na=length(x)
    ctot=sum(x)
    ctot.vector=rep(ctot,na)
    catnump=x
    ind= x==0
    catnump[ind]=0.5
    ctotp=sum(catnump)
    ctotp.vector=rep(ctotp,na)

    catnump[!ind] = catnump[!ind]*ctot.vector[!ind]/ctotp.vector[!ind]
    catch.prop = catnump/ctot.vector
    cum_catch.prop = cumsum(catch.prop)
    cond_catch.prop = catch.prop[1:(na-1)]
    for(i in 2:(na-1)){cond_catch.prop[i] = catch.prop[i]/(1 - cum_catch.prop[i-1])}
    crl = log(cond_catch.prop/(1-cond_catch.prop))
    return(crl)
    }

    if(is.vector(x)) {
        if(length(x)<2) stop('x should be a vector, matrix of data frame')
        crl=trans(x)
    }
    if(is.matrix(x)) crl=apply(x,2,trans)
    if(is.data.frame(x)) crl=apply(as.matrix(x),2,trans)

    return(crl)
}

##' continuous ratio logit transformed matrix of catches-at-age to matrix of proportions
##' @param x catch-at-age matrix (crl transformed, rows=ages, columns=years).
##' @details ...
##' @export

crlInverse = function(x){

  trans=function(x){
     na=length(x)
     mnew=exp(x)/(exp(x)+1)
     mprop=mnew
     mprop[2] = mprop[2]*(1-mprop[1])
     for(i in 3:na){mprop[i] = mnew[i]*(1 - sum(mprop[1:(i-1)]))}
     mcum=cumsum(mprop)
     mcum=c(mcum,1)
     mprop= c(mcum[1],diff(mcum))
     return(mprop)
  }

  if(is.vector(x)){
      if(length(x)<2) stop('x should be a vector, matrix of data frame')
      mprop=trans(x)
  }
  if(is.matrix(x)) mprop=apply(x,2,trans)
  if(is.data.frame(x)) mprop=apply(as.matrix(x),2,trans)

  return(mprop)
}

##' logit function
##' @param x a vector
##' @details ...
##' @export
logit = function (x)
{
    log(x/(1 - x))
}

##' inverse logit function
##' @param x a vector
##' @details ...
##' @export

invlogit = function (x)
{
    1/(1 + exp(-x))
}

##' Standardized proportions at age
##' @param x matrix to be standardized (years in columns, ages in rows)
##' @rdname spya
##' @details They show the age composition of the catch each year, and can be used to track cohorts trends. They show more clearly when catches are above or below average. Negative values are below average, positive values are above average.
##' @export
spya <- function(tdat){
    stdat <- apply(tdat,1,sum,na.rm=T)
    na <- length(stdat)
    ny <- ncol(tdat)
    mstdat <- matrix(stdat,nrow=na,ncol=ny,byrow=F)
    #mtdat <- matrix(tdat,nrow=na,ncol=ny,byrow=F)
    ptdat <- tdat/mstdat
    ptdat.dev <- ptdat - matrix(apply(ptdat,1,mean,na.rm=T),nrow=na,ncol=ny,byrow=T)
    ptdat.std <- sqrt(apply(ptdat.dev^2,1,mean,na.rm=T))
    ztdat <- ptdat.dev/matrix(ptdat.std,nrow=na,ncol=ny,byrow=F)
    return(ztdat)
}

##' Standardized proportions at year
##' @param x matrix to be standardized (years in columns, ages in rows)
##' @rdname spay
##' @details Standardized proportions at year can show cohort patterns more clearly.
##' @export
spay <- function(tdat){
    stdat <- apply(tdat,2,sum,na.rm=T)
    ny <- length(stdat)
    na <- nrow(tdat)
    mstdat <- matrix(stdat,nrow=na,ncol=ny,byrow=T)
    #mtdat <- matrix(tdat,nrow=na,ncol=ny,byrow=F)
    ptdat <- tdat/mstdat
    ptdat.dev <- ptdat - matrix(apply(ptdat,1,mean,na.rm=T),nrow=na,ncol=ny,byrow=F)
    ptdat.std <- sqrt(apply(ptdat.dev^2,1,mean,na.rm=T))
    ztdat <- ptdat.dev/matrix(ptdat.std,nrow=na,ncol=ny,byrow=F)
    return(ztdat)
}


