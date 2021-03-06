##' Fit CCAM model
##' @param data data for the ccam model as returned from the setup.ccam.data function
##' @param conf model configuration which can be set up using the \code{\link{defcon}} function and then modified either directly in R or by saving it to a text file using the function \code{\link{saveConf}}, modifying the text file, and then reading the configuration from the textfile using the function \code{\link{loadConf}}. For more details about the configuration see details.
##' @param parameters initial values which can be set up using the \code{\link{defpar}} function and then modified.
##' @param newtonsteps optional extra true newton steps
##' @param rm.unidentified option to eliminate unidentified model parameters based on gradient in initial value (somewhat experimental)
##' @param run if FALSE return AD object without running the optimization
##' @param lower named list with lower bounds for optimization (only met before extra newton steps)
##' @param upper named list with upper bounds for optimization (only met before extra newton steps)
##' @param sim.condRE logical with default \code{TRUE}. Simulated observations will be conditional on estimated values of F and N, rather than also simulating F and N forward from their initial values.
##' @param paracheck logical. Check if initial parameters are consistent?
##' @param debug logical. print from cpp
##' @param ... extra arguments to MakeADFun
##' @return an object of class \code{ccam}
##' @details The model configuration object \code{conf} is a list of different objects defining different parts of the model. The different elements of the list are:
##' \describe{
##' \item{$minAge:}{A single integer defining the the lowest age class in the assessment.}
##' \item{$maxAge:}{A single integer defining the the highest age class in the assessment.}
##' \item{$maxAgePlusGroup:}{Is last age group considered a plus group (1 yes, or 0 no).}
##' \item{$keyLogFsta:}{A matrix of integers. The number of rows is equal to the number of fleets and the number of columns is equal to the number of age classes. The matrix describes the coupling of the fishing mortality states (normally only first row is used). '-1' is used for entries where no fishing mortality applies (e.g. age groups in survey fleets, or unobserved age groups). For the valid entries consecutive integers starting at zero must be used, because they are used as indices in the corresponding state vector. If the same number is used for two age classes, then the fishing mortality for those age classes are assumed equal (linked to the same state).}
##' \item{$corFlag:}{A single integer to specify the correlation structure of log-scale fishing mortality increments (0 independent, 1 compound symmetry, or 2 AR(1)).}
##' \item{$keyLogFpar:}{A matrix of integers. The number of rows is equal to the number of fleets and the number of columns is equal to the number of age classes. The matrix describes the coupling of survey catchability parameters (so only used for survey fleets). '-1' is used for entries where catchability should not be specified (e.g. fleet - age groups combinations where fishing mortality is specified above, or unobserved fleet - age group combinations).  For the valid entries consecutive integers starting at zero must be used, because they are used as indices in the corresponding parameter vector. If the same number is used for two age classes, then the catchability for those age classes are assumed equal (linked to the same parameter).}
##' \item{$keyQpow:}{A matrix of integers. The number of rows is equal to the number of fleets and the number of columns is equal to the number of age classes. The matrix describes the coupling of density dependent catchability power parameters. This can only be applied to fleets - age combinations where a catchability is defined. '-1' is used for entries where this cannot be applied (e.g. fleet - age groups combinations where fishing mortality is specified above, or unobserved fleet - age group combinations). '-1' is also used to specify that density dependent catchability power parameters is turned off (the most common setup). For entries where density dependent catchability power parameter is to be estimates entries consecutive integers starting at zero must be used. If the same number is used for two age classes, then the density dependent catchability power parameter for those age classes are assumed equal (linked to the same parameter).}
##' \item{$keyVarF:}{A matrix of integers. The number of rows is equal to the number of fleets and the number of columns is equal to the number of age classes. The matrix describes the coupling of variance parameters for the different states in the log-scale fishing mortality random walk process. '-1' should be used for entries where no fishing mortality state is defined in \code{keyLogFsta} above. For the valid entries consecutive integers starting at zero must be used, because they are used as indices in the corresponding parameter vector. If the same number is used for two age classes, then the catchability for those age classes are assumed equal (linked to the same parameter). ((a curiosity of this setup is that it is possible to set different variance parameter indices for F-states that are coupled in \code{keyLogFsta}. This is ignored and the index corresponding to the lowest F-state number is used)).}
##' \item{$keyVarLogN:}{A vector of integers. The length of the vector is equal to the number of age classes. The vector describes the coupling of variance parameters for the log(N)-process. Consecutive integers starting at zero must be used, because they are used as indices in the corresponding parameter vector. If the same number is used for two age classes, then the catchability for those age classes are assumed equal. A typical setup is to use a unique index for the first age group, because that corresponds to the variance in the (stock-)recruitment, which is often not similar to the variance in the survival process from year to year.}
##' \item{$keyVarObs:}{A matrix of integers. The number of rows is equal to the number of fleets and the number of columns is equal to the number of age classes. The matrix describes the coupling of observation variance parameters. '-1' should be used for entries where no observations are available. For the valid entries consecutive integers starting at zero must be used, because they are used as indices in the corresponding parameter vector. If the same number is used for two age classes, then the observation variance for those age classes are assumed equal (linked to the same parameter).}
##' \item{$obsCorStruct:}{A factor specifying the covariance structure used across ages for each fleet. The length of the factor is equal to the number of fleets. The possible options are: ("ID" independent, "AR" AR(1), or "US" for unstructured).}
##' \item{$keyCorObs:}{A matrix of integers. The number of rows is equal to the number of fleets and the number of columns is equal to the number of age classes _minus_ _one_. The matrix describes the coupling AR correlations between age classes, and hence is only meaningful for fleets where the "AR" observation correlation structure is chosen. '-1' should be used for entries where no observations are available. Notice that the matrix has one column less than the number of age classes, which is because the correlation between age classes is described. Consecutive integers starting at zero must be used. If the same number is used for a given fleet it means that a normal AR(1) structure is used. If different numbers are used for a fleet it means that the correlation parameter changes where the numbers differ. If the "AR" structure is specified above, then the corresponding row in this matrix must have valid non-negative entries.}
##' \item{$stockRecruitmentModelCode:}{A single integer to specify the stock recruitment connection to use (0 for plain random walk on log recruitment, 1 for Ricker, and 2 for Beverton-Holt).}
##' \item{$noScaledYears:}{A single integer specifying the number of years where catch scaling is to be estimated (most often 0, as this is a somewhat exotic option).}
##' \item{$keyScaledYears:}{A vector of the years where catch scaling is applied (length should match \code{noScaledYears}) (most often empty, as this is a somewhat exotic option).}
##' \item{$keyParScaledYA:}{A matrix of integers specifying the couplings of scale parameters (nrow = \code{noScaledYears}, ncols = no ages) (most often empty, as this is a somewhat exotic option).}
##' \item{$fbarRange:}{An integer vector of length 2 specifying lowest and highest age included in Fbar (average fishing mortality summary).}
##' \item{$keyBiomassTreat:}{A vector of integers with length equal to the number of fleets. '-1' should be used for entries where the corresponding fleet is not a mass index. A the corresponding fleet is a mass index, then three options are available (0 SSB index, 1 catch index, and 2 FSB index).}
##' \item{$obsLikelihoodFlag:}{A factor specifying the type of likelihood to use for each fleet. The length of the factor is equal to the number of fleets. The possible options are: ("LN" for log-normal and "ALN" Additive logistic normal).}
##' \item{$fixVarToWeight:}{A single integer. If weight attribute is supplied for observations this option defines how it is treated (0 as relative weight, 1 as a fixed variance = weight).}
##' }
##' @importFrom TMB MakeADFun sdreport
##' @importFrom stats nlminb optimHess
##' @importFrom utils relist packageDescription
##' @useDynLib CCAM
##' @export
##' @examples
##' data(canmackData)
##' data(canmackConf)
##' data(canmackParameters)
##' fit <- ccam.fit(canmackData, canmacknscodConf, canmackParameters)
ccam.fit <- function(data, conf, parameters, newtonsteps=3, rm.unidentified=FALSE, run=TRUE, lower=getLowerBounds(parameters), upper=getUpperBounds(parameters), sim.condRE=TRUE, paracheck=TRUE, debug=FALSE,...){
  definit <- defpar(data, conf)
  if(paracheck){
      if(!identical(parameters,relist(unlist(parameters), skeleton=definit))){
          warning("Initial values are not consistent, so running with default init values from defpar()")
          parameters<-definit
      }
  }
  censored <- 'CE' %in% conf$obsLikelihoodFlag
  if(censored){
          TrueObs <- data$logobs
          data$logobs[!is.na(data$logobs[,2]),1] <- log(apply(exp(data$logobs[!is.na(data$logobs[,2]),]),1,mean))
          conf$obsLikelihoodFlag[1] <- 'LN'
          if(!all(conf$keyVarObs[1,]<0)){
              conf$keyVarObs <- conf$keyVarObs-1
              conf$keyVarObs[conf$keyVarObs<0] <- -1
              parameters$logSdLogObs <- parameters$logSdLogObs[-(unique(conf$keyVarObs[1,])+1)]
          }
  }

  data <- clean.void.catches(data,conf)
  tmball <- c(data, conf, simFlag=as.numeric(sim.condRE))
  if(is.null(tmball$resFlag)){tmball$resFlag <- 0}
  tmball$debug <- as.numeric(debug)
  nmissing <- sum(is.na(data$logobs[,1]))
  parameters$missing <- numeric(nmissing)
  ran <- c("logN", "logFy", "missing")
  obj <- MakeADFun(tmball, parameters, random=ran, DLL="CCAM",...)
  if(rm.unidentified){
    skel <- parameters[!names(parameters)%in%ran]
    gr <- obj$gr()
    safemap <- relist(gr,skel)
    safemap <- lapply(safemap, function(x)factor(ifelse(abs(x)>1.0e-15,1:length(x),NA)))
    obj <- MakeADFun(tmball, parameters, random=ran, map=safemap, DLL="CCAM", ...)
  }

  lower2<-rep(-Inf,length(obj$par))
  upper2<-rep(Inf,length(obj$par))
  for(nn in names(lower)) lower2[names(obj$par)==nn]=lower[[nn]]
  for(nn in names(upper)) upper2[names(obj$par)==nn]=upper[[nn]]

  if(!run) return( list(sdrep=NA, pl=parameters, plsd=NA, data=data, conf=conf, opt=NA, obj=obj) )
  opt <- nlminb(obj$par, obj$fn,obj$gr ,control=list(trace=1, eval.max=2000, iter.max=1000),lower=lower2,upper=upper2)

  if(censored){
      parameters <- obj$env$parList(opt$par)
      data$logobs <- TrueObs
      conf$obsLikelihoodFlag[1]='CE'

      tmball <- c(data, conf, simFlag=as.numeric(sim.condRE))
      if(is.null(tmball$resFlag)){tmball$resFlag <- 0}
      tmball$debug <- as.numeric(debug)
      obj <- MakeADFun(tmball, parameters, random=ran, DLL="CCAM",...)

          if(rm.unidentified){
              skel <- parameters[!names(parameters)%in%ran]
              gr <- obj$gr()
              safemap <- relist(gr,skel)
              safemap <- lapply(safemap, function(x)factor(ifelse(abs(x)>1.0e-15,1:length(x),NA)))
              obj <- MakeADFun(tmball, parameters, random=ran, map=safemap, DLL="CCAM", ...)
          }
      opt <- nlminb(obj$par, obj$fn,obj$gr ,control=list(trace=1, eval.max=2000, iter.max=1000),lower=lower2,upper=upper2)
  }
  for(i in seq_len(newtonsteps)) { # Take a few extra newton steps
      g <- as.numeric( obj$gr(opt$par) )
      h <- optimHess(opt$par, obj$fn, obj$gr)
      opt$par <- opt$par - solve(h, g)
      opt$objective <- obj$fn(opt$par)
  }
  rep <- obj$report()
  sdrep <- sdreport(obj,opt$par)

  # Last two states
  parnames <- paste('par',unique(names(sdrep$par.fixed)),sep=".")
  idx <- which(names(sdrep$value) %in% c("lastLogN","lastLogF","logitSel",parnames))
  sdrep$estY <- sdrep$value[idx]
  sdrep$covY <- sdrep$cov[idx,idx]

  idx <- which(names(sdrep$value) %in% c("beforeLastLogN","beforeLastLogF","logitSel",parnames))
  sdrep$estYm1 <- sdrep$value[idx]
  sdrep$covYm1 <- sdrep$cov[idx,idx]

  pl <- as.list(sdrep,"Est")
  plsd <- as.list(sdrep,"Std")

  sdrep$cov<-NULL # save memory

  ret <- list(sdrep=sdrep, pl=pl, plsd=plsd, data=data, conf=conf, opt=opt, obj=obj, rep=rep, low=lower2, hig=upper2)
  attr(ret, "RemoteSha") <- substr(packageDescription("ccam")$RemoteSha, 1, 12)
  attr(ret, "Version") <- packageDescription("ccam")$Version
  class(ret)<-"ccam"
  return(ret)
}

##' Bounds
##' @param parameters initial values for the model in a format similar to what is returned from the defpar function
##' @return a named list
getLowerBounds<-function(parameters){
    list(sigmaObsParUS=rep(-10,length(parameters$sigmaObsParUS)))
}

##' Bounds
##' @param parameters initial values for the model in a format similar to what is returned from the defpar function
##' @return a named list
getUpperBounds<-function(parameters){
    list(sigmaObsParUS=rep(10,length(parameters$sigmaObsParUS)))
}

##' remove void catches
##' @param dat data for the ccam model as returned from the setup.ccam.data function
##' @param conf model configuration which can be set up using the \code{\link{defcon}} function and then modified
##' @return an updated dataset without the catches where F is fixed to zero
clean.void.catches<-function(dat, conf){
    caa=ifelse("Catch-at-age proportions" %in% attributes(dat)$fleetNames,3,1)
    rmidx <- ((dat$aux[,3]%in%(conf$minAge:conf$maxAge)[which(conf$keyLogFsta[1,]==(-1))])&dat$aux[,2]==caa)
    dat$aux <- dat$aux[!rmidx,]
    dat$logobs <- dat$logobs[!rmidx,]
    dat$weight <- dat$weight[!rmidx]
    dat$nobs<-sum(!rmidx)
    dat$minAgePerFleet<-as.integer(tapply(dat$aux[,"age"], INDEX=dat$aux[,"fleet"], FUN=min))
    dat$maxAgePerFleet<-as.integer(tapply(dat$aux[,"age"], INDEX=dat$aux[,"fleet"], FUN=max))
    newyear<-min(as.numeric(dat$aux[,"year"])):max(as.numeric(dat$aux[,"year"]))
    newfleet<-min(as.numeric(dat$aux[,"fleet"])):max(as.numeric(dat$aux[,"fleet"]))
    mmfun<-function(f,y, ff){idx<-which(dat$aux[,"year"]==y & dat$aux[,"fleet"]==f); ifelse(length(idx)==0, NA, ff(idx)-1)}
    dat$idx1<-outer(newfleet, newyear, Vectorize(mmfun,c("f","y")), ff=min)
    dat$idx2<-outer(newfleet, newyear, Vectorize(mmfun,c("f","y")), ff=max)
    dat
}


##' Jitter runs
##' @param data data to fit a model
##' @param conf model configuration
##' @param nojit a list of vectors. Each element in the list specifies a run where the fleets mentioned are omitted
##' @param par initial values to jitter around. The defaule ones are returned from the defpar function
##' @param sd the standard deviation used to jitter the initial values (most parameters are on a log scale, so similar to cv)
##' @param ncores the number of cores to attemp to use
##' @return A "ccamset" object, which is basically a list of ccam fits
##' @details ...
##' @importFrom parallel detectCores makeCluster clusterEvalQ parLapply stopCluster
##' @importFrom stats rnorm
##' @export
jit <- function(data, conf, par=defpar(data, conf), nojit=10,  sd=.25, ncores=detectCores(),parallell=TRUE,...){
  parv <- unlist(par)
  pars <- lapply(1:nojit, function(i)relist(parv+rnorm(length(parv),sd=sd), par))

  if(parallell){
      cl <- makeCluster(ncores) #set up nodes
      clusterEvalQ(cl, {library(CCAM)}) #load the package to each node
      clusterExport(cl, c('conf','data'), envir=environment()) #load the data to each node
      fits <- parLapply(cl, pars, function(p) ccam.fit(data, conf, p, silent = TRUE, paracheck=FALSE,...))
      stopCluster(cl) #shut it down
  }else{
      fits <- lapply(pars, function(p) ccam.fit(data, conf, p, silent = TRUE, paracheck=FALSE,...))
  }

  attr(fits,"fit") <- fit
  attr(fits,"jitflag") <- 1
  class(fits) <- c("ccamset")
  fits
}

##' Force fit
##' @param data data to fit a model
##' @param conf model configuration
##' @param parameters initial parameter values
##' @param nojit maximum number of tries with different initial values
##' @param silent logical
##' @param ... agruments to ccam.fit
##' @return ccam object
##' @details Try fitting the model with different initial parameters untill it fits (or untill it is clear it most likely never will)
##' @export
force.fit <- function(data,conf, parameters,nojit=10,silent=TRUE,...){
    # try to run the model
    fit <- tryCatch(ccam.fit(data, #if(u==1){x}else{attr(x,paste0('alternative',u-1))},
                             conf,
                             parameters,
                             silent=silent,
                             paracheck=FALSE,...),error=identity)
    # if it doesn't work, try with default parameters
    if(inherits(fit, "error")){  #if there is an error
        if(!identical(parameters,defpar(data,conf))){
            fit <- tryCatch(ccam.fit(data, #if(u==1){x}else{attr(x,paste0('alternative',u-1))},
                                     conf,
                                     defpar(data,conf),
                                     silent=silent,...),error=identity)
        }
        # and if it still doesn't work, try 10 times with different initial parameters
        if(inherits(fit, "error")){
            for(i in 1:nojit){
                fit <- tryCatch(jit(data,conf,parameters,nojit=1,parallell = FALSE,...)[[1]],error=identity)
                if(!inherits(fit, "error")){
                    if(fit$opt$convergence==0) break
                }else{
                    if(i==nojit){
                        fit <- paste0('CCAM model: error after about',nojit+1,' trials: ', fit)
                    }
                }
            }

        }
    }
    return(fit)
}

