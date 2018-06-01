##' Plot helper
##' @param x dataframe from table()
##' @param ylab y label
##' @param xlab x label
##' @details The basic plotting used by many of the plotting functions (e.g. ssbplot, fbarplot ...)
plotit <-function (x,ylab='Estimate',xlab='Year',ci=TRUE){
    UseMethod("plotit")
}

##' @rdname plotit
##' @method plotit dfccam
##' @export
plotit.dfccam <- function(x,ylab='Estimate',xlab='Year',ci=TRUE,col="black"){
 p <- ggplot(x,aes(x=year,y=Estimate))+geom_line(size=1.5,col=col)+
    scale_y_continuous(expand=c(0,0))+
    scale_x_continuous(expand=c(0,0))+
    theme_bw()+ylab(ylab)+xlab(xlab)
 if(ci){
     p <- p+geom_ribbon(aes(ymax=High, ymin=Low), alpha=0.2,fill=col)
 }
 return(p)
}

##' @rdname plotit
##' @method plotit dfccamset
##' @export
plotit.dfccamset <- function(x,ylab='Estimate',xlab='Year', ci=TRUE,col=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")){
    coll <- c(rep(col,length(unique(x$fit)) %/% length(col)),col[1:length(unique(x$fit)) %% length(col)])
    p <- ggplot(x,aes(x=year,y=Estimate,group=fit))+geom_line(aes(col=fit),size=1.5)+
        scale_y_continuous(expand=c(0,0))+
        scale_x_continuous(expand=c(0,0))+
        scale_fill_manual(values=coll)+
        scale_color_manual(values=coll)+
        labs(col='',fill='')+
        theme_bw()+ylab(ylab)+xlab(xlab)+
        theme(legend.position = c(0.8,0.8),
              legend.background = element_blank())
    if(ci){
        p <- p+geom_ribbon(aes(ymax=High, ymin=Low, fill=fit), alpha=0.2)
    }
    return(p)
}

##' @rdname plotit
##' @method plotit dfccamforecast
##' @export
plotit.dfccamforecast <- function(x,ylab='Estimate',xlab='Year',ci=TRUE,col='black'){
    p <- ggplot(x,aes(x=year,y=Estimate))+geom_line(size=1.5,col=col)+
        scale_x_continuous(expand=c(0,0))+
        theme_bw()+ylab(ylab)+xlab(xlab)+
        geom_errorbar(data=x[x$period=='Future',],aes(ymax=High, ymin=Low),col=col,alpha=0.2)+
        scale_y_continuous(expand=c(0,0),limits=c(0,max(df$High)*1.05))
    if(ci){
        p <- p+geom_ribbon(data=x[x$period=='Passed',],aes(ymax=High, ymin=Low), alpha=0.2,fill=col)
    }
    return(p)
}

##' @rdname plotit
##' @method plotit dfforecastset
##' @export
plotit.dfforecastset <- function(x,ylab='Estimate',xlab='Year', ci=TRUE, col=c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")){
    coll <- c(rep(col,length(unique(x$fit)) %/% length(col)),col[1:length(unique(x$fit)) %% length(col)])
    ggplot(x,aes(x=year,y=Estimate,group=fit))+geom_line(size=1.5,aes(col=fit))+
        geom_errorbar(data=x[x$period=='Future',],aes(ymax=High, ymin=Low,col=fit))+
        scale_y_continuous(expand=c(0,0.5),limits=c(0,max(x$High)*1.1))+
        scale_x_continuous(expand=c(0,0))+
        theme_bw()+ylab(ylab)+xlab(xlab)+
        scale_fill_manual(values=coll)+
        scale_color_manual(values=coll)+
        labs(col='',fill='')+
        theme(legend.position = c(0.8,0.8),
              legend.background = element_blank())
    if(ci){
        p <- p+geom_ribbon(data=x[x$period=='Passed',],aes(ymax=High, ymin=Low,fill=fit), alpha=0.2)
    }
    return(p)
}

##' CCAM SSB plot
##' @param x the object(s) returned from ccam.fit or forecast
##' @param ... extra arguments transferred to plotit
##' @details Plot of spawning stock biomass
##' @export
ssbplot<-function(x, ...){
    UseMethod("ssbplot")
}
##' @rdname ssbplot
##' @method ssbplot default
##' @export
ssbplot.default <- function(x,...){
    df <- ssbtable(x)
    plotit(df,ylab='SSB',...)
}

##' CCAM TSB plot
##' @param x the object(s) returned from ccam.fit
##' @param ... extra arguments transferred to plotit
##' @details Plot of total stock biomass
##' @export
tsbplot<-function(x, ...){
    UseMethod("tsbplot")
}
##' @rdname ssbplot
##' @method ssbplot default
##' @export
tsbplot.default <- function(x,...){
    df <- tsbtable(x)
    plotit(df,ylab='TSB',...)
}

##' CCAM exploitation rate plot
##' @param x the object(s) returned from ccam.fit
##' @param ... extra arguments transferred to plotit
##' @details Plot of total stock biomass
##' @export
expplot<-function(x, ...){
    UseMethod("expplot")
}
##' @rdname expplot
##' @method expplot default
##' @export
expplot.default <- function(x,...){
    df <- exptable(x)
    plotit(df,ylab='Exploitation rate',...)
}

##' CCAM recruitment plot
##' @param x the object(s) returned from ccam.fit or forecast
##' @param ... extra arguments transferred to plotit
##' @details Plot of recruitment
##' @export
recplot<-function(x, ...){
    UseMethod("recplot")
}
##' @rdname recplot
##' @method recplot default
##' @export
recplot.default <- function(x,...){
    df <- rectable(x)
    plotit(df,ylab='Recruitment',...)
}

##' CCAM fbar plot
##' @param x the object(s) returned from ccam.fit or forecast
##' @param ... extra arguments transferred to plotit
##' @details Plot of mean F
##' @export
fbarplot<-function(x, ...){
    UseMethod("fbarplot")
}
##' @rdname fbarplot
##' @method fbarplot default
##' @export
fbarplot.default <- function(x,...){
    df <- fbartable(x)
    plotit(df,ylab='Fbar',...)
}

##' CCAM selectivity plot
##' @param x the object(s) returned from ccam.fit or forecast
##' @param ... extra arguments transferred to plotit
##' @details Plot of selectivity
##' @export
selplot<-function(x, ...){
    UseMethod("selplot")
}
##' @rdname selplot
##' @method selplot default
##' @export
selplot.default <- function(x,...){
    df <- seltable(x)
    plotit(df,ylab='Selectivity',...)
}

##' CCAM catch plot
##' @param x the object(s) returned from ccam.fit or forecast
##' @param fleet add observations from this fleet
##' @param ... extra arguments transferred to plotit
##' @details Plot of spawning stock biomass
##' @export
catchplot<-function(x,fleet=NULL,...){
    UseMethod("catchplot")
}
##' @rdname catchplot
##' @method catchplot default
##' @export
catchplot.default <- function(x,fleet=NULL,...){
    df <- catchtable(x,fleet=fleet)
    a <- list(...)
    a$x <- df
    p <- do.call(plotit, a)
    if(!is.null(fleet)){
        if(class(x) %in% c('ccam','ccamforecast')) p <- p+geom_line(aes(y=aux1))+geom_line(aes(y=aux2))
        if(class(x) %in% c('ccamset','forecastset')) p <- p+geom_line(aes(y=aux1,col=fit))+geom_line(aes(y=aux2,col=fit))
    }
    p
}

##' CCAM model plot (ssb, F, rec)
##' @param x the object(s) returned from ccam.fit or forecast
##' @details Plot of spawning stock biomass
##' @export
modelplot<-function(x){
    UseMethod("modelplot")
}
##' @rdname modelplot
##' @method modelplot default
##' @importFrom gridExtra grid.arrange
##' @export
modelplot.default <- function(x){
    grid.arrange(ssbplot(x),
                 fbarplot(x),
                 recplot(x),
                 ncol=1)
}

# srplot
##' Plots the stock recruitment
##' @param fit the object returned from ccam.fit
##' @param textcol color of years on plot
##' @param linecol color of lines
##' @param curve add SR curve
##' @import ggplot2
##' @export
srplot<-function(fit,textcol="red",linecol='black',curve=FALSE){
    UseMethod("srplot")
}

##' @rdname srplot
##' @method srplot ccam
##' @export
srplot.ccam <- function(fit,...){
    fit <- c(fit)
    srplot(fit,...)
}


##' @rdname srplot
##' @method srplot ccamset
##' @export
srplot.ccamset <- function(fit,textcol="red",linecol='black',curve=FALSE){
    li <- lapply(1:length(fit),function(x){
        X <- summary(fit[[x]])
        n<-nrow(X)
        lag <- fit[[x]]$conf$minAge
        idxR <- (lag+1):n
        idxS <- 1:(n-lag)
        R<-X[idxR,1]
        S<-X[idxS,4]
        Rnam<<-colnames(X)[1]
        Snam<<-colnames(X)[4]
        y<-rownames(X)[idxR]
        Fnam <- names(fit)[x]
        if(is.null(Fnam)) Fnam <- x
        df<- data.frame(R=R,S=S, Rnam=Rnam,Snam=Snam,y=y,fit=Fnam)
        if(curve){
            a <- fit[[x]]$opt$par['rec_loga']
            b <- fit[[x]]$opt$par['rec_logb']
            e <- fit[[x]]$opt$par['rec_e']
            if(is.na(e)) e <- 0
            env <-fit[[x]]$data$env[idxS,1]
            df$SRenv <- exp(a+log(S)-log(1.0+exp(b)*S)+e*env)
        }
        return(df)
    })
    df <- do.call('rbind',li)

    p <- ggplot(df,aes(x=S,y=R))+geom_point()+
        geom_path()+
        ylab(Rnam)+xlab(Snam)+
        theme_bw()+
        geom_text(aes(label=y),col=textcol,hjust=0,vjust=0)

    if(curve){
        p <- p+geom_line(aes(y=SRenv),size=1.5)
    }
    if(length(fit)>1){
        p<- p+ facet_wrap(~fit)+
            theme(strip.background = element_blank(),
                  panel.border = element_rect(colour = "black"))
    }
    p
}

##' CCAM PE plot
##' @param x the object(s) returned from ccam.fit or forecast
##' @param heat logical. If true a heat plot is used, if false a bubble plot
##' @param ... parameters used in heat() or bubble()
##' @details Plot of proces error
##' @import ggplot2
##' @importFrom reshape2 melt
##' @export
peplot<-function(x,heat=TRUE,...){
    UseMethod("peplot")
}

##' @rdname peplot
##' @method peplot ccam
##' @export
peplot.ccam <- function(x,...){
    x <- c(x)
    peplot(x,...)
}

##' @rdname peplot
##' @method peplot ccamset
##' @export
peplot.ccamset <- function(x,heat=TRUE,...){
    pe <- lapply(x,function(y){pe <- t(y$rep$resN)
                               pe <- rbind(0,pe)
                               dimnames(pe) <- dimnames(y$data$natMor)
                               pe
    })
    if(heat){
     heat(pe,...)
    }else{
        bubble(pe,...)
    }
}

##' Plots residuals
##' @param fit the object returned from ccam.fit
##' @param fleets an integer vector of fleets to plot. Default is all of them
##' @param type Type of residual plot
##' @param std logical; standardize residuals?
##' @param out number of outliers to which to add identifier
##' @param ... extra arguments to plot
##' @details types are residuals ~ year (1); residuals ~ predicted (2); predicted ~ observed (3)
##' @export
resplot <- function(fit, log=TRUE, ...){
    UseMethod("resplot")
}
##' @rdname resplot
##' @method resplot ccam
##' @import ggplot2
##' @importFrom reshape2 dcast
##' @export
resplot.ccam <- function(fit, trans=function(x) x,fleets=unique(fit$data$aux[,"fleet"]), type=1, std=TRUE,out=0,...){
    idx<-fit$data$aux[,"fleet"]%in%fleets
    p <- trans(fit$obj$report(c(fit$sdrep$par.fixed,fit$sdrep$par.random))$predObs[idx])
    o <- trans(fit$data$logobs[idx,1])
    res <- p-o
    aa <- fit$data$aux[idx,"age"]
    neg.age <- (aa < -1.0e-6)
    aa[neg.age] <- NA
    Year <- fit$data$aux[idx,"year"]
    sds <- exp(fit$obj$par[which(names(fit$obj$par)=='logSdLogObs')])
    keysd <- fit$conf$keyVarObs[fleets,]
    keysd <- keysd[keysd>-1]
    sd <- sds[keysd+1]
    if(std) res <- res/sd[ifelse(is.na(aa),1,aa)]
    df <- data.frame(year=Year, p=p, o=o, res=res, age=aa)

    if(out!=0){
        res <- res[!is.na(res)]
        e <- sort(abs(res),decreasing = T)[as.numeric(1:out)]
        e <- which(res %in% c(e,-e))
        if(all(is.na(aa[e]))) aprint='' else aprint=paste0(".",aa[e])
        lab <- paste0(Year[e],aprint)
    }

    if(all(is.na(aa))){  #not age structured
        switch (type,
                'one' = {
                    p <- ggplot(df,aes(x=year,y=res))+geom_point()+theme_bw()+ylab('Residuals')+xlab('Year')+
                        geom_hline(yintercept = 0,lty='dashed')
                },
                'two' = {
                    p <- ggplot(df,aes(x=p,y=res))+geom_point()+theme_bw()+ylab('Residuals')+xlab('Predicted')+
                        geom_hline(yintercept = 0,lty='dashed')
                },
                'three' = {
                    p <- ggplot(df,aes(x=o,y=p))+geom_point()+theme_bw()+ylab('Predicted')+xlab('Observed')+
                        geom_abline(intercept=0, slope=1, lty='dashed')
                },
                'four' = {
                    p <- ggplot(df,aes(x=year))+geom_line(aes(y=p))+geom_point(aes(y=o))+
                    theme_bw()+
                    ylab('Value')+xlab('Year')
                },
                {
                    stop('type should be between 1 and 4 (without age level)')
                }
        )
    }else{
        switch (type,
                'one' = {
                    m <- as.matrix(dcast(data.frame(x=Year,y=aa,value=res),x~y)[,-1])
                    rownames(m) <- unique(Year)
                    p <- heat(m,posneg=TRUE)
                },
                'two' = {
                    p <- ggplot(df,aes(x=p,y=res,col=age))+geom_point()+theme_bw()+ylab('Residuals')+xlab('Predicted')+
                        geom_hline(yintercept = 0,lty='dashed')+
                        labs(col='Age')+scale_color_gradient(low='violetred',high='orange')
                },
                'three' = {
                    p <- ggplot(df,aes(x=p,y=o,col=age))+geom_point()+theme_bw()+ylab('Residuals')+xlab('Observed')+
                        geom_abline(intercept=0, slope=1, lty='dashed')+
                        labs(col='Age')+scale_color_gradient(low='violetred',high='orange')
                },
                'four' = {
                    p <- ggplot(df,aes(x=year))+geom_line(aes(y=p))+geom_point(aes(y=o))+
                        facet_wrap(~age,ncol=3)+
                        theme_bw()+
                        theme(panel.border = element_rect(colour = "black", fill=NA),
                              strip.background = element_blank())+
                        ylab('Value')+xlab('Year')
                },
                'five' = {
                    p <- ggplot(df,aes(x=year,y=res))+geom_point()+
                        facet_wrap(~age,ncol=3)+
                        geom_hline(yintercept = 0, lty='dashed')+
                        geom_smooth(method = "loess",col='black')+
                        theme_bw()+
                        theme(panel.border = element_rect(colour = "black", fill=NA),
                              strip.background = element_blank())+
                        ylab('Residuals')+xlab('Year')
                },
                {
                    stop('type should be between 1 and 5 (with age level)')
                }
        )
    }
    if(out!=0) p <- p + geom_text(data=df[e,],aes(label=lab),hjust=0,vjust=0)
    print(p)

}


##' CCAM parameter plot
##' @param fit the object returned from ccam.fit
##' @param cor.report.limit correlations with absolute value > this number is reported in the plot
##' @param col colors
##' @param ... extra arguments transferred to plot
##' @details Plot of all estimated model parameters (fixed effects). Shown with confidence interval.
##' @export
##' @importFrom stats cov2cor
parplot<-function(fit, cor.report.limit=0.95, col=NULL){
    UseMethod("parplot")
}

##' @rdname parplot
##' @method parplot default
##' @export
parplot.default <- function(fit, cor.report.limit=0.95, col=NULL){
    if(class(fit)=='ccam') fit=c(fit)
    if(is.null(col)){
        col <- if(length(fit)==1) 'black' else c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")
    }
    coll <- c(rep(col,length(fit) %/% length(col)),col[1:length(fit) %% length(col)])
    lab <- names(fit)
    if(is.null(lab)) lab <- 1:length(fit)
    param <- lapply(fit, coef)
    nam <- names(param[[1]])
    dup <- duplicated(nam)
    namadd <- rep(0,length(nam))
    for(i in 2:length(dup)){
        if(dup[i])namadd[i] <- namadd[i-1]+1
    }
    nam <- paste(nam, namadd, sep="_")

    mat <- lapply(1:length(param),function(i){
        m <- param[[i]]+t(c(-2,0,2)%o%attr(param[[i]],"sd"))
        rownames(m) <- 1:nrow(m)
        m <- data.frame(m)
        colnames(m) <- c('low','est','high')
        m$nam <- nam
        m$fit <- lab[i]
        m}
    )
    mat <- do.call('rbind',mat)

    corrs <- lapply(1:length(param), function(x) {corrs <- cov2cor(attr(param[[x]], "cov"))-diag(length(param[[x]]))
                                        rownames(corrs)<-nam
                                        colnames(corrs)<-nam
                                        co <- data.frame(melt(corrs))
                                        names(co) <- c('nam','label','value')
                                        co$hc <- ifelse(co$value>cor.report.limit,round(co$value*100),NA)
                                        co$lc <- ifelse(co$value< -cor.report.limit,round(co$value*100),NA)
                                        co <- co[!(is.na(co$hc) & is.na(co$lc)),]
                                        co$fit <- lab[x]
                                        co
    })
    print(do.call('rbind',corrs))

    p <- ggplot(mat,aes(x=nam,y=est,col=as.factor(fit),group=as.factor(fit)))+
        geom_point(position=position_dodge(width=0.5))+
        geom_errorbar(aes(ymin=low,ymax=high),position=position_dodge(width=0.5),width=0)+
        xlab('Estimate')+ylab('Parameter')+
        theme_bw()+labs(col='')+
        coord_flip()+
        scale_color_manual(values=coll)

    if(length(fit)==1) p <- p + theme(legend.position = 'none')
    p
}

##' Plot surveys
##' @param x surveys from read.ices
##' @rdname surveyplot
##' @import ggplot2
##' @importFrom reshape2 melt
##' @export
surveyplot <- function(x){
     ages <- do.call('c',lapply(x,colnames))

     x2 <- lapply(1:length(x),function(i){
        y <- x[[i]]
        m <- matrix(NA,nrow=nrow(y),ncol=length(ages),dimnames = list(rownames(y),ages))
        m[,colnames(y)] <- y
        m <- as.data.frame(m)
        m$year <- as.numeric(rownames(m))
        n <- names(x)[i]
        if(is.null(n)) n <- i
        m$survey <- n
        m
     })
     x2 <- do.call('rbind',x2)
     df <- melt(x2,id=c('year','survey'),variable.name = 'age')
     p <- ggplot(df,aes(x=year,y=value,col=age))+geom_line(size=1.5)+
        ylab('Value')+xlab('Year')+labs(col='Age')+
        theme_bw()+
        facet_wrap(~survey,scale='free')+
         theme(strip.background = element_blank(),
               panel.border = element_rect(colour = "black"))
     if(length(unique(df$survey)==1) & length(unique(df$age))){
         p <- p+ theme(legend.position = 'none')+
             scale_color_manual(values='black')
     }
     p
}


##' Plot implementation error
##' @param IEmeans list of mean values
##' @param IEsds list of sd values
##' @param col colors
##' @rdname IEplot
##' @import ggplot2
##' @export
IEplot <- function(IEmeans,IEsds, col=c("#000000", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#882255", "#AA4499")){
    myIE <-cbind(melt(do.call('rbind',IEmeans)),sd=melt(do.call('rbind',IEsds))[,3])
    p<- ggplot(myIE,aes(x=Var2-1,y=value))+
        geom_ribbon(aes(ymin=value-sd,ymax=value+sd,fill=Var1),alpha=0.2)+
        geom_line(aes(col=Var1),size=1.5)+
        geom_point(aes(col=Var1))+
        theme(legend.position="none")+
        ylab('Undeclared catch (t)')+xlab('Year')+labs(col='',fill='')+
        theme_bw()+
        scale_fill_manual(values=col[1:length(unique(myIE$Var1))])+
        scale_color_manual(values=col[1:length(unique(myIE$Var1))])
    return(p)
}


##' Bubble plot
##' @param x matrix
##' @param scale scaler for the bubble size (scale_range)
##' @param col color for the bubbles. if 2 colors are given, they indicate positive and negative values
##' @param alpha transparency of the bubbles
##' @param xlab x label
##' @param ylab y label
##' @rdname bubble
##' @import ggplot2
##' @export
bubble <- function(x,xlab='Year',ylab='Age',scale=15,col=c('black','darkgreen'),alpha=0.8){
    if(class(x)!='list') x=list(x)
    n <- names(x)
    if(is.null(n)) n <- 1:length(x)
    xx <- lapply(n,function(y){
        xx <- melt(x[[y]])
        xx$id <- y
        xx
    })
    xx <- do.call('rbind',xx)

    xx$posneg <- ifelse(xx$value<0,'#8B000099','#00640099')
    coll <- length(unique(xx$posneg))
    p <- ggplot(xx,aes(x=Var1,y=Var2,size=value,col=posneg))+geom_point(alpha=alpha)+
        scale_size(range = c(1,scale)) +
        theme_bw()+ylab(ylab)+xlab(xlab)+
        labs(size="")+
        scale_color_manual(values=col[1:coll])+
        guides(col=FALSE)

    if(length(unique(xx$id))>1){
        p <- p+facet_wrap(~id)+
            theme(strip.background = element_blank(),
                  panel.border = element_rect(colour = "black"))

    }
    p
}

##' Heat plot
##' @param x matrix
##' @param col colors that indicate the color scale limit
##' @param ncol number of colors to use in plot
##' @param leg logical (plot legend?)
##' @param leground number of digits to round legend
##' @rdname bubble
##' @importFrom reshape2 melt
##' @import ggplot2
##' @export
heat <- function(x,high=c('grey','darkgreen'), low='darkred',xlab='Year',ylab='Age',posneg=FALSE){
    if(class(x)!='list') x=list(x)
    n <- names(x)
    if(is.null(n)) n <- 1:length(x)
    xx <- lapply(n,function(y){
        xx <- melt(x[[y]])
        xx$id <- y
        xx
    })
    xx <- do.call('rbind',xx)


    xx$posneg <- ifelse(xx$value<0,'-','+')
    p <- ggplot(xx,aes(x=Var1,y=Var2,fill=value))+geom_tile()+
        theme_bw()+ylab(ylab)+xlab(xlab)+
        scale_x_continuous(expand=c(0,0))+
        scale_y_continuous(expand=c(0,0))+
        labs(fill='')+
        scale_fill_gradient(low=low,high=high)
    if(posneg){
        p <- p+geom_text(aes(label=posneg))
    }
    if(length(unique(xx$id))>1){
        p <- p+facet_wrap(~id)+
            theme(strip.background = element_blank(),
                  panel.border = element_rect(colour = "black"))

    }
    p
}

##' Plots fit to data
##' @param x the object(s) returned from ccam.forecast
##' @param what.y statistic returned from ccam.forecast
##' @param what.x statistic returned from ccam.forecast
##' @param ylab y label
##' @param xlab x label
##' @param rect color zone above this y value
##' @param ci logical (ci are quantile ranges)
##' @param vline plot vertical line(s)
##' @param hline plot horizontal line(s)
##' @param by IE, MP, OM
##' @param text logical. plot year labels
##' @param ... extra arguments to plot
##' @details Plot of probability of what over time
##' @importFrom graphics points
##' @export
foreplot <- function(x, what.y,what.x=NULL, ylab=what.y,xlab='Year',rect=NULL,ci=TRUE,vline=NULL,hline=NULL,by=NULL,text=FALSE){
    UseMethod("foreplot")
}
##' @rdname foreplot
##' @method foreplot ccamforecast
##' @export
foreplot.ccamforecast <- function(x,...){
    x=c(x)
    foreplot(x,...)
}

##' @rdname foreplot
##' @method foreplot forecastset
##' @export
foreplot.forecastset <- function(x, what.y,what.x=NULL, ylab=what.y,xlab='Year',rect=NULL,ci=TRUE,vline=NULL,hline=NULL,by=NULL, text=FALSE){

    byopt <- c('OM','IE','MP')
    if(!is.null(by) & !all(by %in% byopt)){stop('"by" may only include OM, IE or MP')}

    # create data frame
    df <- extract(x,what.y,add=TRUE)
    colnames(df)[1]='y'
    df[is.na(df$IE),'IE'] <- 'IE0'
    df[is.na(df$MP),'MP'] <- 'MP0'
    df[is.na(df$IE),'OM'] <- 'OM0'

    inter <- ifelse(ncol(df)==6,FALSE,TRUE)
    if(inter){
        colnames(df)[2:3]=c('ylow','yhigh')
    }

    if(!is.null(what.x)){
       df$x <- extract(x,what.x)[,1]
       xlab <- what.x
    }else{
       colnames(df)[which(colnames(df)=='year')]='x'
    }

    # aestethics
    if(length(by)==0){
        p <- ggplot(df,aes(x=x,y=y,col=as.factor(id)))
    }
    if(length(by)==1){
        int <- paste0("interaction(", paste0(byopt[byopt!=by], collapse =  ", "), ")")
        p <- ggplot(df,aes_string(x='x',y='y',col=int))+facet_wrap(reformulate(by) )
    }
    if(length(by)==2){
        p <- ggplot(df,aes_string(x='x',y='y',col=byopt[!byopt %in% by]))+facet_grid(reformulate(by[1],by[2]) )
    }
    if(length(by)==3){
        p <- ggplot(df,aes(x=x,y=yEstimate))+facet_grid(reformulate(by[1],by[2:3])  )
    }
    # rectangle
    if(!is.null(rect)) p <- p+ geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=rect, ymax=Inf),fill='lightgrey',alpha=0.3,col='white')

    # basis
    p <- p+ geom_point()+
            geom_path()+
            theme_bw()+ylab(ylab)+xlab(xlab)+
            labs(col='')+
            theme(strip.background = element_blank(),
                  panel.border = element_rect(colour = "black"))

    # vline
    if(!is.null(vline)){
        for(i in 1:length(vline)){
            p <- p+geom_vline(xintercept = vline[i],lty='dashed',col='darkgrey')
        }
    }

    # hline
    if(!is.null(hline)){
        for(i in 1:length(hline)){
            p <- p+geom_hline(yintercept = hline[i],lty='dashed',col='darkgrey')
        }
    }

    # text
    if(text){p <- p+geom_text(aes(label=year),hjust=-0.5,vjust=-0.5)}

    # quantiles
    if(inter & ci){
        colnames(df)[2:3]=c('low','high')
        p <- p+geom_errorbar(data = df,aes(ymin=low,ymax=high),size=1,width=0.5,alpha=0.2)+
            scale_y_continuous(expand=c(0,0),limits=c(0,max(df$high)*1.05))
    }else{
        p <- p+scale_x_continuous(expand=c(0,0))
    }

    # color
    if(length(unique(df$id))==1){
        p <- p+scale_color_manual(values ='black')+theme(legend.position = '')
    }

    return(p)
}

##' Plot MSE criteria
##' @param x the object(s) returned from forecast
##' @details types are residuals ~ year (1); residuals ~ predicted (2); predicted ~ observed (3)
##' @export
MSEplot <- function(x){
    UseMethod("MSEplot")
}

##' @rdname MSEplot
##' @method MSEplot ccamforecast
##' @export
MSEplot.ccamforecast <- function(x,...){
    x <- list(x)
    MSEplot(x,...)
}

##' @rdname MSEplot
##' @method MSEplot forecastset
##' @import ggplot2
##' @importFrom plyr ddply
##' @importFrom gtable gtable_filter gtable_matrix gtable_add_grob
##' @importFrom grid grid.newpage unit grid.draw
##' @details plots a table with the result for each objective (as TRUE/FALSE or a number) by OM, IE and MP. Objectives need to be adapted inside the plot code.
##' @export
MSEplot.forecastset <- function(x){

    meltfacet <- function(x,n1,n2){
        z <- ggplotGrob(x)
        locations <- grep("strip-t", z$layout$name)
        strip <- gtable_filter(z, "strip-t", trim = FALSE)
        ix <- 1:nrow(strip$layout)
        top <- strip$layout$t[1]
        l   <- strip$layout$l[seq(1,max(ix),by=n2) ]
        r   <- strip$layout$r[seq(n2,max(ix),by=n2) ]
        mat   <- matrix(vector("list", length = (n1*2-1)*2), nrow = 2)
        mat[] <- list(zeroGrob())
        res <- gtable_matrix("toprow", mat, unit(c(rep(c(1,0),n1-1),1), "null"), unit(c(1, 1), "null"))
        se <- seq(1,length(locations),by=n2)
        for(i in se){
            zz. <- gtable_add_grob(res,z$grobs[[locations[i]]]$grobs[[1]], 1, 1, 1, (n2)*2-1)
            z <- gtable_add_grob(z, zz., t = top,  l = l[which(se==i)],  b = top,  r = r[which(se==i)], name = c("add-strip"))
        }
        grid.newpage()
        print(grid.draw(z))
    }

    # 1) by next 5 year 75 % probCZ
    x=runlist
    df <- extract(x,'probCZ',add=TRUE)
    df <- df[df$year==(min(df$year)+5),]
    df$ylog <- ifelse(df$var>=0.75,TRUE,FALSE)
    df$ytext <- NA
    df$performance <- "in 5 year 75 % >CZ"
    out <- df
    # 2) by next 10 years 75% probHZ
    df <- extract(x,'probHZ',add=TRUE)
    df <- df[df$year==(min(df$year)+5),]
    df$ylog <- ifelse(df$var>=0.75,TRUE,FALSE)
    df$performance <- "in 10 year 75 % >HZ"
    df$ytext <- NA
    out <- rbind(out,df)
    # 3) never below 50% chance of 20% growth when in the critical zone (95%)
    df <- extract(x,'probgrowth20',add=TRUE)
    df2 <- extract(x,'probCZ',add=TRUE)
    ys <- ddply(df2, c('id'),summarise,year=year[var<0.95])  ##ddply!!!
    sel <- merge(df,ys)
    ylog <- ddply(sel,c('id'),summarise,ylog=any(var<0.50))
    sel <- ddply(sel,c('id'),summarise,year=max(year))
    sel <- merge(merge(ylog,sel),df)
    sel$performance <- "at least 20% growth (>50%) when in CZ (95%)"
    sel$ytext=NA
    out <- rbind(out,sel)
    # 4) show number of years TAC at or above 8000t
    df <- extract(x,'TAC',add=TRUE)
    colnames(df)[1] <- 'var'
    sel <- ddply(df[!df$year==min(df$year),],c('id','OM','MP','IE'),summarise,ytext=length(var[var>=8000]))
    sel$year <- NA
    sel$var <- NA
    sel$ylog <- NA
    sel$performance <- 'n years TAC >= 8000t'
    out <- rbind(out,sel)
    # 5) range of change
    #df <- extract(x,'TACrel',add=TRUE)
    #colnames(df)[1] <- 'var'
    #sel <- ddply(df[!df$year==min(df$year),],c('id','OM','MP','IE'),summarise,y=paste(range(var,collapse = '-')))
    #sel$year <- NA
    #sel$var <- NA
    #sel$performance <- 'range of change'
    #out <- rbind(out,sel)


    out$colo <- ifelse(is.logical(out$y),out$y,'black')

    p <- ggplot(out,aes(x=1,y=performance))+geom_point(aes(col=ylog),size=2)+
        geom_text(aes(label=ytext))+
        facet_grid(OM~MP+as.numeric(as.factor(IE)))+
        scale_x_continuous(expand=c(0,0))+
        theme_classic()+
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position="")+
        ylab('')+xlab('')+
        scale_color_manual(values=c('red','green'))

    meltfacet(p,length(unique(out$MP)),length(unique(out$IE)))

}

