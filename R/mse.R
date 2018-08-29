##' fleet 3
##' @param x simulations of states
##' @param fit object from ccam.fit
##' @param f fleet (numeric)
##' @param deterministic logistic
##' @details generates observations of fleet 3 (keybiomasstreat 4 not yet implemented, 3 needs to be verified)
##' @return returns an array of dim (na-1) x 2 (observation columns) x number of simulations
##' @export
fleet3 <- function(x,fit,f,simpara,deterministic,ULpool=NULL, nm=nm, sw=sw, mo=mo, pm=pm, pf=pf, cw=cw){
    if(fit$conf$keyBiomassTreat[f]==0) obs <- ssb(x, nm, sw, mo, pm, pf)
    if(fit$conf$keyBiomassTreat[f] %in% c(1,3)) obs <- catch(x,nm, cw)
    if(fit$conf$keyBiomassTreat[f]==2) obs <- fsb(x, cw, nm)
    if(fit$conf$keyBiomassTreat[f]<3){
        q <- simpara[,which(colnames(simpara)=="logFpar"),drop=FALSE]
        q <- q[,fit$conf$keyLogFpar[f,1]+1]
        obs <- obs*exp(q)
    }

    obs <- cbind(obs,NA)
    if(fit$conf$obsLikelihoodFlag[f]=='CE'){
        if(deterministic){
            obs[,2] <- obs[,1]
        }else{
            obs[,1] <- mean(ULpool[,1],1)*obs[,1]
            obs[,2] <- mean(ULpool[,2],1)*obs[,1]
        }
    }
    obs[obs[,2]<3,2] <- 5 #at least 5 fish (so lower limit does not get to zero)
    obs <- log(obs)
    return(t(obs))
}

##' fleet 6
##' @param x simulations of states
##' @param nm simulation of natural mortality
##' @details generates observations of fleet 6 (deterministic)
##' @return returns an array of dim (na-1) x 2 (observation columns) x number of simulations
##' @export
fleet6 <- function(x,nm){
    caasim <- caa(x, nm)
    obs <- crlTransform(t(caasim))
    obs <- lapply(split(t(obs),1:ncol(obs)),cbind,NA) #results to list with second column for upper limit (NA)
    obs <- sapply(1:length(obs), function(y) obs[[y]], simplify = 'array') #list to array
    return(obs)
}

##' create new observations (logobs), only for one new year
##' @param fit returned from ccam.fit
##' @param sim simulation of states
##' @param ULpool matrix with ratios of predicted values vs the upper or lower bound
##' @param deterministic logical
##' @param nm vector of nm (natural mortality)
##' @param sw vector of sw (stock mean weight)
##' @param mo vector of mo (propotion mature)
##' @param pm vector of pm (proportion natural mortality)
##' @param pf vector of pf (proportion fishing mortality)
##' @param cw vector of cw (catch weight)
##' @details create new observations to apply model based management procedures
##' @return returns and array of size (n observations,2 (lower and upper column), number of simulations)
##' @export
newobs <- function(fit, sim, simpara, ULpool,deterministic, nm, sw=sw, mo=mo, pm=pm, pf=pf, cw=cw){
    aux <- unique(fit$data$aux[,2:3])
    fleetsim <- array(dim=c(nrow(aux),2,nrow(sim)))
    for(f in unique(aux[,1])){
        ft <- fit$data$fleetTypes[f]
        if(ft==3) {
            mysim <- fleet3(x=sim,fit,f,simpara,deterministic,ULpool, nm=nm, sw=sw, mo=mo, pm=pm, pf=pf, cw=cw)
            fleetsim[which(aux[,1]==f),,] <- mysim
        }
        if(ft==6) {
            mysim <- fleet6(x=sim,nm)
            fleetsim[which(aux[,1]==f),,] <- mysim
        }
    }
    return(fleetsim)
}

##' add one year to parameters
##' @param x fit returned from ccam.fit
##' @export
steppar <- function(x){
    x$logFy <- c(x$logFy,tail(x$logFy,1))
    x$logN <- cbind(x$logN,x$logN[,ncol(x$logN)])
    return(x)
}

##' add one year to configuration
##' @param x fit returned from ccam.fit
##' @export
stepconf <- function(x){
    x$keySel <- rbind(x$keySel,tail(x$keySel,1))
    return(x)
}

##' add one year to data
##' @param x fit returned from ccam.fit
##' @param ... next year data objects
##' @export
stepdat <- function(x,obs,aux,mo,sw,cw,nm,lf,dw,lw,pm,pf,en){
    newY <- tail(x$years,1)+1
    x$noYears <- x$noYears+1
    x$years <- c(x$years,newY)
    x$aux <- rbind(x$aux,aux)
    x$logobs <- rbind(x$logobs,obs)
    x$nobs <- nrow(x$logobs)
    newyear <- min(as.numeric(x$aux[,1])):max(as.numeric(x$aux[,1]))
    newfleet <- min(as.numeric(x$aux[,2])):max(as.numeric(x$aux[,2]))
    mmfun<-function(f,y, ff){idx <- which(x$aux[,1]==y & x$aux[,2]==f); ifelse(length(idx)==0, NA, ff(idx-1))}
    x$idx1 <- outer(newfleet, newyear, Vectorize(mmfun,c("f","y")), ff=min)
    x$idx2 <- outer(newfleet, newyear, Vectorize(mmfun,c("f","y")), ff=max)
    x$weight <- c(x$weight, tail(x$weight, nrow(aux)))
    x$propMat <- rbind(x$propMat,mo)
    rownames(x$propMat)[length(rownames(x$propMat))] <- newY
    x$stockMeanWeight <- rbind(x$stockMeanWeight ,sw)
    rownames(x$stockMeanWeight)[length(rownames(x$stockMeanWeight))] <- newY
    x$catchMeanWeight <- rbind(x$catchMeanWeight ,cw)
    rownames(x$catchMeanWeight )[length(rownames(x$catchMeanWeight ))] <- newY
    x$natMor <- rbind(x$natMor ,nm)
    rownames(x$natMor )[length(rownames(x$natMor ))] <- newY
    x$landFrac <- rbind(x$landFrac ,lf)
    rownames(x$landFrac )[length(rownames(x$landFrac))] <- newY
    x$disMeanWeight <- rbind(x$disMeanWeight ,dw)
    rownames(x$disMeanWeight )[length(rownames(x$disMeanWeight ))] <- newY
    x$landMeanWeight <- rbind(x$landMeanWeight ,lw)
    rownames(x$landMeanWeight )[length(rownames(x$landMeanWeight ))] <- newY
    x$propM <- rbind(x$propM ,pm)
    rownames(x$propM )[length(rownames(x$propM))] <- newY
    x$propF <- rbind(x$propF ,newY=pf)
    rownames(x$propF)[length(rownames(x$propF))] <- newY
    x$env <- rbind(x$env,en)
    rownames(x$env)[length(rownames(x$env))] <- newY
    return(x)
}
