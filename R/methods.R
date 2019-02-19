##' Plot ccam object
##' @method plot ccam
##' @param  x ...
##' @param  ... extra arguments (not possible to use add=TRUE --- please collect to a list of fits using e.g the c(...), and then plot that collected object)
##' @importFrom gridExtra gtable_rbind
##' @import ggplot2 grid
##' @details ...
##' @export
plot.ccam<-function(x, ...){
    gA <- ggplotGrob(ssbplot(x,...))
    gB <- ggplotGrob(fbarplot(x,...))
    gC <- ggplotGrob(recplot(x,...))
    grid.newpage()
    grid.draw(gtable_rbind(gA, gB, gC))
}

##' Plot ccam object
##' @method plot ccamset
##' @param  x ...
##' @param  ... extra arguments
##' @importFrom gridExtra grid.arrange arrangeGrob gtable_rbind
##' @import ggplot2 grid
##' @details ...
##' @export
plot.ccamset<-function(x, ...){
    extractLegend<-function(p){
        tmp <- ggplot_gtable(ggplot_build(p))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
    }

    mylegend<-extractLegend(ssbplot(x,...)+theme(legend.position='right'))
    lheight <- sum(mylegend$height)
    lwidth <- sum(mylegend$width)

    gA <- ggplotGrob(ssbplot(x,...)+ theme(legend.position="none"))
    gB <- ggplotGrob(fbarplot(x,...)+ theme(legend.position="none"))
    gC <- ggplotGrob(recplot(x,...)+ theme(legend.position="none"))

    grid.arrange(gtable_rbind(gA, gB, gC),
                 mylegend,
                 ncol = 2,
                 widths = unit.c(unit(1, "npc") - lwidth, lwidth))
}

##' Plot ccamforecast object
##' @method plot ccamforecast
##' @param  x ...
##' @param  ... extra arguments
##' @importFrom gridExtra gtable_rbind
##' @import ggplot2 grid
##' @details ...
##' @export
plot.ccamforecast<-function(x, ...){
    gA <- ggplotGrob(ssbplot(x,...))
    gB <- ggplotGrob(fbarplot(x,...))
    gC <- ggplotGrob(recplot(x,...))
    grid.newpage()
    grid.draw(gtable_rbind(gA, gB, gC))
}

##' Plot ccamforecast object
##' @method plot forecastset
##' @param  x ...
##' @param  ... extra arguments
##' @importFrom grid unit.pmax
##' @import ggplot2
##' @details ...
##' @export
plot.forecastset<-function(x, ...){
    extractLegend<-function(p){
        tmp <- ggplot_gtable(ggplot_build(p))
        leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
        legend <- tmp$grobs[[leg]]
        return(legend)
    }

    mylegend<-extractLegend(ssbplot(x,...)+theme(legend.position='right'))
    lheight <- sum(mylegend$height)
    lwidth <- sum(mylegend$width)

    gA <- ggplotGrob(ssbplot(x,...)+ theme(legend.position="none"))
    gB <- ggplotGrob(fbarplot(x,...)+ theme(legend.position="none"))
    gC <- ggplotGrob(recplot(x,...)+ theme(legend.position="none"))

    grid.arrange(gtable_rbind(gA, gB, gC),
                 mylegend,
                 ncol = 2,
                 widths = unit.c(unit(1, "npc") - lwidth, lwidth))
}

##' Collect ccam objects
##' @method c ccam
##' @param  ... ccam fits to be combined
##' @details ...
##' @export
c.ccam<-function(...){
  ret<-list(...)
  class(ret)<-"ccamset"
  ret
}
##' Collect ccamforecast objects
##' @method c ccamforecast
##' @param  ... ccam forecasts to be combined
##' @details ...
##' @export
c.ccamforecast<-function(...){
    ret<-list(...)
    class(ret)<-"forecastset"
    ret
}

##' Compute process residuals (single joint sample)
##' @param fit the fitted object as returned from the ccam.fit function
##' @param ... extra arguments (not currently used)
##' @return an object of class \code{ccamres}
##' @details ...
##' @import TMB
##' @export
procres <- function(fit, ...){
  pp<-fit$pl
  attr(pp,"what") <- NULL
  pp$missing <- NULL
  fit.co<-ccam.fit(fit$data, fit$conf, pp, run=FALSE)
  fit.co$obj$env$data$resFlag<-1
  fit.co$obj$retape()
  sdrep <- sdreport(fit.co$obj,fit$opt$par)
  ages <- as.integer(colnames(fit.co$data$natMor))
  iF<-fit.co$conf$keySel[1,]
  if (exists(".Random.seed")){
    oldseed <- get(".Random.seed", .GlobalEnv)
    oldRNGkind <- RNGkind()
  }
  set.seed(123456)
  idx <- which(names(sdrep$value)=="resN")
  resN <- rmvnorm(1,mu=sdrep$value[idx], Sigma=sdrep$cov[idx,idx])
  resN <- matrix(resN, nrow=nrow(fit.co$pl$logN))
  resN <- data.frame(year=fit.co$data$years[as.vector(col(resN))],
                     fleet=1,
                     age=ages[as.vector(row(resN))],
                     residual=as.vector(resN))
  ret <- resN
  attr(ret, "fleetNames") <- "Joint sample residuals log(N)"
  class(ret) <- "ccamres"
  if (exists("oldseed")){
    do.call("RNGkind",as.list(oldRNGkind))
    assign(".Random.seed", oldseed, .GlobalEnv)
  }
  return(ret)
}

##' Plot ccam residuals
##' @method plot ccamres
##' @param  x ...
##' @details ...
##' @import ggplot2
##' @importFrom plyr ddply
##' @importFrom gridExtra grid.arrange
##' @export
##' @examples
##' \dontrun{
##' data(canmackData)
##' data(canmackConf)
##' data(canmackParameters)
##' fit <- ccam.fit(canmackData, canmackConf, canmackParameters)
##' plot(residuals(fit))
##' }
plot.ccamres<-function(x,qq=TRUE,fleet=NULL){

  slope <- function(x){diff(quantile(x[!is.na(x)], c(0.25, 0.75)))/diff(qnorm(c(0.25, 0.75)))}
  intercept <- function(x){quantile(x[!is.na(x)], c(0.25, 0.75))[1L] - diff(quantile(x[!is.na(x)], c(0.25, 0.75)))/diff(qnorm(c(0.25, 0.75))) * qnorm(c(0.25, 0.75))[1L]}

  neg.age <- (x$age < -1.0e-6)
  x$age[neg.age] <- mean(x$age[!neg.age],na.rm=TRUE)

  df <- data.frame(do.call('cbind',x))
  df$by <- attr(x,"fleetNames")[x$fleet]
  df[df$by=='Total catch','age'] <- rep(c(0,1),each=length(df[df$by=='Total catch','age'])/2)
  if('observation' %in% colnames(df)) df <- df[!is.na(df$observation),]

  df$negpos <- ifelse(df$residual<0,'-','+')
  df <- ddply(df,c('by'),transform,slope=slope(residual),int=intercept(residual) )

  if(!is.null(fleet)){df <- df[df$fleet %in% fleet,]}

  p1 <- ggplot(df,aes(x=year,y=age,size=abs(residual),col=negpos))+geom_point(alpha=0.6)+
      scale_size(range = c(1, 9)) +ylab('Age') + xlab('Year')+
      scale_color_manual(values=c('darkred','darkgreen'))+
      facet_grid(by~.,scales = 'free_y')+
      guides(col=FALSE)+
      theme(legend.position='left',
            strip.background = element_blank(),
            strip.text.y = element_blank())

  if(isTRUE(qq)){
      p2 <- ggplot(df, aes(sample = residual)) + stat_qq() +
          facet_grid(by~.)+
          geom_abline(aes(slope = slope, intercept = int))
      grid.arrange(p1,p2,ncol=2,widths=c(4,1))
  }else{
      p1
  }
}

##' Print ccam object
##' @method print ccam
##' @param  x ...
##' @param  ... extra arguments
##' @details ...
##' @export
print.ccam<-function(x, ...){
  cat("CCAM model: log likelihood is", logLik.ccam(x,...),"Convergence", ifelse(0==x$opt$convergence, "OK\n", "failed\n"))
}

##' Print ccamres object
##' @method print ccamres
##' @param  x ...
##' @param  ... extra arguments
##' @details ...
##' @export
print.ccamres<-function(x, ...){
  class(x)<-NULL
  print(as.data.frame(x))
}

##' Log likelihood of ccam object
##' @method logLik ccam
##' @param  object ccam fitted object (result from ccam.fit)
##' @param  ... extra arguments
##' @details ...
##' @export
logLik.ccam<-function(object, ...){
  ret<- -object$opt$objective
  attr(ret,"df")<-length(object$opt$par)
  class(ret)<-"logLik"
  ret
}

##' Extract fixed coefficients of ccam object
##' @method coef ccam
##' @param  object ccam fitted object (result from ccam.fit)
##' @param  ... extra arguments
##' @details ...
##' @importFrom stats coef
##' @export
coef.ccam <- function(object, ...){
  ret <- object$sdrep$par.fixed
  attr(ret,"cov") <- object$sdrep$cov.fixed
  attr(ret,"sd") <- sqrt(diag(object$sdrep$cov.fixed))
  class(ret)<-"ccamcoef"
  ret
}

##' Print ccamcoef object
##' @method print ccamcoef
##' @param  x ...
##' @param  ... extra arguments
##' @details ...
##' @export
print.ccamcoef<-function(x, ...){
  y<-as.vector(x)
  names(y)<-names(x)
  print(y)
}


##' Extract number of observations from ccam object
##' @method nobs ccam
##' @param object ccam fitted object (result from ccam.fit)
##' @param ... extra arguments
##' @importFrom stats nobs
##' @details ...
##' @export
nobs.ccam<-function(object, ...){
  as.integer(object$data$nobs)
}

##' Extract residuals from ccam object
##' @method residuals ccam
##' @param object ccam fitted object (result from ccam.fit)
##' @param discrete logical if model contain discrete observations
##' @param ... extra arguments for TMB's oneStepPredict
##' @importFrom stats residuals
##' @importFrom TMB oneStepPredict
##' @details ...
##' @export
residuals.ccam<-function(object, discrete=FALSE, ...){
  cat("One-observation-ahead residuals. Total number of observations: ", nobs(object), "\n")
  res <- oneStepPredict(object$obj, observation.name="logobs", data.term.indicator="keep", discrete=discrete,...)
  cat("One-observation-ahead residuals. Done\n")
  ret <- cbind(object$data$aux, res)
  attr(ret,"fleetNames") <- attr(object$data, "fleetNames")
  class(ret)<-"ccamres"
  ret
}

##' Summary of ccam object
##' @method summary ccam
##' @param object ccam fitted object (result from ccam.fit)
##' @param ... extra arguments
##' @details ...
##' @export
summary.ccam<-function(object, ...){
  ret <- cbind(round(rectable(object)[,-4]), round(ssbtable(object)[,-4]), round(fbartable(object),3)[,-4])
  colnames(ret)[1] <- paste("R(age ", object$conf$minAge, ")", sep="")
  colnames(ret)[4] <- "SSB"
  colnames(ret)[7] <- paste("Fbar(",object$conf$fbarRange[1], "-", object$conf$fbarRange[2], ")", sep="")
  ret
}

##' Simulate from a ccam object
##' @method simulate ccam
##' @param object ccam fitted object (result from ccam.fit)
##' @param nsim number of response lists to simulate. Defaults to 1.
##' @param seed random number seed
##' @param full.data logical, should each inner list contain a full list of data. Defaults to TRUE
##' @param ... extra arguments
##' @importFrom stats simulate
##' @details ...
##' @return returns a list of lists. The outer list has length \code{nsim}. Each inner list contains simulated values of \code{logF}, \code{logN}, and \code{obs} with dimensions equal to those parameters.
##' @export
simulate.ccam<-function(object, nsim=1, seed=NULL, full.data=TRUE, ...){
  if(!is.null(seed)) set.seed(seed)
  est <- unlist(object$pl)
  if(full.data){
    ret <- replicate(nsim,
    	c(object$data[names(object$data)!="logobs"],#all the old data
    	object$obj$simulate(est)["logobs"])#simulated observations
    	, simplify=FALSE)
    ret<-lapply(ret, function(x){attr(x,"fleetNames") <- attr(object$data,"fleetNames");x})
  }else{
  	ret <- replicate(nsim, object$obj$simulate(est), simplify=FALSE)
  }
  ret
}


##' Plot ccam object
##' @method plot ccamypr
##' @param  x ...
##' @param  ... extra arguments
##' @importFrom gridExtra grid.arrange
##' @import ggplot2
##' @details ...
##' @export
plot.ccamypr<-function(x, ...){
    df <- data.frame(YPR=x$yield,Fbar=x$fbar, SPR=x$ssb)

    yend <- df[which(abs(df$Fbar-x$f01)==min(abs(df$Fbar-x$f01))),'YPR']
    yend.max <- df[which(abs(df$Fbar-x$fmax)==min(abs(df$Fbar-x$fmax))),'YPR']
    yend.spr30=df[which(abs(df$Fbar-x$f30)==min(abs(df$Fbar-x$f30))),'SPR']
    yend.spr40=df[which(abs(df$Fbar-x$f40)==min(abs(df$Fbar-x$f40))),'SPR']
    maxf <- max(df$Fbar)
    maxypr <- max(df$YPR)

    a <- ggplot(df,aes(x=Fbar,y=YPR))+geom_line(size=1.5)+xlab(x$fbarlab)+ylab('Yield per Recruit')+
        geom_segment(aes(x=x$f01 , y = 0, xend = x$f01, yend =yend), linetype = "dotted")+
        geom_text(aes(x=x$f01*1.1,y=yend*0.9,label=paste("F0.1 =",x$f01)),vjust=0,hjust=0)+
        scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0),limits=c(0,max(df$YPR)*1.1))+
        theme_bw()

    if(x$fmax>max(df$Fbar)){
        a <- a+geom_segment(aes(x=maxf*0.9,y=maxypr*0.9,xend=maxf,yend=maxypr*0.9),arrow = arrow(length = unit(0.2,"cm")))+
            geom_text(aes(x=maxf*0.8,y=maxypr*0.9),label='Fmax')
    }else{
        a <- a+geom_segment(aes(x=x$fmax, y = 0, xend = x$fmax, yend = yend.max), linetype = "dotted")+
            geom_text(aes(x=x$fmax*1.1,y=yend.max*0.9,label=paste("Fmax =",x$fmax)),vjust=0,hjust=0)
    }


    b <- ggplot(df,aes(x=Fbar,y=SPR),environment=environment())+geom_line(size=1.5)+xlab(x$fbarlab)+ylab('Spawners per Recruit')+
        geom_segment(aes(x=x$f30, y = 0, xend = x$f30, yend = yend.spr30), linetype = "dotted")+
        geom_segment(aes(x=x$f40, y = 0, xend = x$f40, yend = yend.spr40), linetype = "dotted")+
        geom_text(aes(x=x$f30*1.05,y=yend.spr30*1.05,label=paste("F30% =",x$f30)),hjust=0,vjust=0.2)+
        geom_text(aes(x=x$f40*1.05,y=yend.spr40*1.05,label=paste("F40% =",x$f40)),hjust=0,vjust=0.2)+
        scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+
        theme_bw()

    grid.arrange(a,b)
}

##' Print ccamypr object
##' @method print ccamypr
##' @param  x an object as returned from the ypr function
##' @param  ... extra arguments
##' @details ...
##' @export
print.ccamypr <- function(x, ...){
  idx <- c(x$fmaxIdx, x$f01Idx, x$f30Idx,x$f35Idx,x$f40Idx)
  ret <- cbind(x$fbar[idx],x$ssb[idx],x$yield[idx])
  rownames(ret) <- c("Fmax", "F01","F30", "F35","F40")
  colnames(ret) <- c("Fbar", "SSB", "Yield")
  print(ret)
}

##' Print ccamforecast object
##' @method print ccamforecast
##' @param  x an object as returned from the forecast function
##' @param  ... extra arguments
##' @details ...
##' @export
print.ccamforecast<-function(x, ...){
  print(attr(x,"tab"))
}

##' Print ccamset object
##' @method print ccamset
##' @param  x a list of ccam models
##' @param  ... extra arguments
##' @details ...
##' @importFrom stats logLik
##' @export
print.ccamset<-function(x,...){
  if("jitflag"%in%names(attributes(x))){
    fit<-attr(x,"fit")
    maxabsdiff <- apply(abs(do.call(cbind, lapply(x, function(f)unlist(f$pl)-unlist(fit$pl)))),1,max)
    maxlist <- relist(maxabsdiff, fit$pl)
    ret <- as.data.frame(unlist(lapply(maxlist,function(x)if(length(x)>0)max(x) else NULL)))
    fbar <- max(unlist(lapply(x, function(f)abs(fbartable(f)[,1]-fbartable(fit)[,1]))))
    ssb <- max(unlist(lapply(x, function(f)abs(ssbtable(f)[,1]-ssbtable(fit)[,1]))))
    rec <- max(unlist(lapply(x, function(f)abs(rectable(f)[,1]-rectable(fit)[,1]))))
    catch <- max(unlist(lapply(x, function(f)abs(catchtable(f)[,1]-catchtable(fit)[,1]))))
    logLik <- max(abs(unlist(lapply(x, logLik))-logLik(fit)))
    ret <- rbind(ret, ssb=ssb,  fbar=fbar, rec=rec, catch=catch, logLik=logLik)
    names(ret) <- "max(|delta|)"
    print(ret)
  }else{
    print.default(x,...)
  }
}


