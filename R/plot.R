##' Plot helper
##' @param x dataframe from table()
##' @param ylab y label
##' @param xlab x label
##' @param ci logical (plot confidence interval?)
##' @param years years to plot
##' @param linesize linesize
##' @param linetype linetype
##' @param scale values by which to devide the estimate
##' @import ggplot2 viridis
##' @details The basic plotting used by many of the plotting functions (e.g. ssbplot, fbarplot ...)
plotit <-function (x,ylab='Estimate',xlab='Year',ci=TRUE,years=unique(x$year),linetype=1,linesize=2,scale=1,...){
    UseMethod("plotit")
}

##' @rdname plotit
##' @export
plotit.dfccam <- function(x,ylab='Estimate',xlab='Year',ci=TRUE,years=unique(x$year),col="black",linetype=1,scale=1){
    x <- x[x$year %in% years,]
    x$Estimate <- x$Estimate/scale
    if(ci){x$High <- x$High/scale;x$Low <- x$Low/scale}
     p <- ggplot(x,aes(x=year,y=Estimate))+geom_line(col=col,linetype=linetype)+
    scale_y_continuous(expand=c(0,0))+
    scale_x_continuous(expand=c(0,0))+
    ylab(ylab)+xlab(xlab)
 if(ci){
     p <- p+geom_ribbon(aes(ymax=High, ymin=Low), alpha=0.2,fill=col)
 }
 return(p)
}

##' @rdname plotit
##' @method plotit dfccamset
##' @export
plotit.dfccamset <- function(x,ylab='Estimate',xlab='Year', ci=TRUE,years=unique(x$year),col=NULL,linetype=1,legendnames=NULL,scale=1){
    x <- x[x$year %in% years,]
    if(!is.null(legendnames)){
        x$fit <- as.factor(x$fit)
        levels(x$fit) <- legendnames
    }
    x$Estimate <- x$Estimate/scale
    if(ci){x$High <- x$High/scale;x$Low <- x$Low/scale}
    p <- ggplot(x,aes(x=year,y=Estimate,group=fit))+geom_line(aes(col=fit),linetype=linetype)+
        scale_y_continuous(expand=c(0,0))+
        scale_x_continuous(expand=c(0,0))+

        labs(col='',fill='')+
        ylab(ylab)+xlab(xlab)
    if(is.null(col)){
        p<- p+scale_color_viridis(discrete = TRUE)+
            scale_fill_viridis(discrete = TRUE)
    }else{
        coll <- c(rep(col,length(unique(x$fit)) %/% length(col)),col[1:length(unique(x$fit)) %% length(col)])
        p <- p+scale_fill_manual(values=coll)+
            scale_color_manual(values=coll)
    }
    if(ci){
        p <- p+geom_ribbon(aes(ymax=High, ymin=Low, fill=fit), alpha=0.2)
    }
    return(p)
}

##' @rdname plotit
##' @method plotit dfccamforecast
##' @export
plotit.dfccamforecast <- function(x,ylab='Estimate',xlab='Year',ci=TRUE,years=unique(x$year),col='black',linetype=1,final=TRUE,scale=1){
    x <- x[x$year %in% years,]
    x$Estimate <- x$Estimate/scale
    if(ci){x$High <- x$High/scale;x$Low <- x$Low/scale}
    if(!isTRUE(final))x[x$period=='Future','Estimate'][1]<-tail(x[x$period=='Passed','Estimate'],1)
    p <- ggplot(x,aes(x=year,y=Estimate))+geom_line(col=col,linetype=linetype)+
        scale_x_continuous(expand=c(0,0))+
        ylab(ylab)+xlab(xlab)+
        geom_errorbar(data=x[x$period=='Future',],aes(ymax=High, ymin=Low),col=col,alpha=0.2)+
        scale_y_continuous(expand=c(0,0),limits=c(0,max(x[x$year %in% years,]$High)*1.05))
    if(ci){
        p <- p+geom_ribbon(data=x[x$period=='Passed',],aes(ymax=High, ymin=Low), alpha=0.2,fill=col)
    }
    return(p)
}

##' @rdname plotit
##' @method plotit dfforecastset
##' @export
plotit.dfforecastset <- function(x,ylab='Estimate',xlab='Year', ci=TRUE, years=unique(x$year),col=NULL,linetype=1,legendnames=NULL,final=TRUE,scale=1){
    x <- x[x$year %in% years,]
    x$Estimate <- x$Estimate/scale
    if(ci){x$High <- x$High/scale;x$Low <- x$Low/scale}
    if(!isTRUE(final))x[x$period=='Future' & x$year==2018,'Estimate']<-x[x$period=='Passed' & x$year==2018,'Estimate']
    if(!is.null(legendnames)){
        x$fit <- as.factor(x$fit)
        levels(x$fit) <- legendnames
    }
    p <- ggplot(x[x$year %in% years,],aes(x=year,y=Estimate,group=fit))+geom_line(aes(col=fit),linetype=linetype)+
        ylab(ylab)+xlab(xlab)+
        labs(col='',fill='')
    if(is.null(col)){
        p<- p+scale_color_viridis(discrete = !is.numeric(x$fit))+
                scale_fill_viridis(discrete = !is.numeric(x$fit))
    }else{
        coll <- c(rep(col,length(unique(x$fit)) %/% length(col)),col[1:length(unique(x$fit)) %% length(col)])
        p <- p+scale_fill_manual(values=coll)+
            scale_color_manual(values=coll)
    }
    if(ci){
        p <- p+geom_ribbon(data=x[x$period=='Passed'& x$year %in% years,],aes(ymax=High, ymin=Low,fill=fit), alpha=0.2)+
            geom_errorbar(data=x[x$period=='Future',],aes(ymax=High, ymin=Low,col=fit))+
            scale_y_continuous(expand=c(0,0.5),limits=c(0,max(x$High)*1.1))+
            scale_x_continuous(expand=c(0,0))
    }else{
        ymax <- max(x$Estimate)
        p <- p+scale_y_continuous(expand=c(0,0.5),limits=c(0,ymax*1.1))+
            scale_x_continuous(expand=c(0,0))
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

##' CCAM SSB0 plot
##' @param x the object(s) returned from ccam.fit or forecast
##' @param ... extra arguments transferred to plotit
##' @details Plot of spawning stock biomass
##' @export
ssb0plot<-function(x, ...){
    UseMethod("ssb0plot")
}
##' @rdname ssb0plot
##' @method ssb0plot default
##' @export
ssb0plot.default <- function(x,...){
    df <- ssb0table(x)
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
    plotit(df,ylab='Exploitation rate',...)+scale_y_continuous(limits=c(0,1),expand=c(0,0))
}

##' CCAM recruitment plot
##' @param x the object(s) returned from ccam.fit or forecast
##' @param trans exp by default
##' @param ... extra arguments transferred to plotit
##' @details Plot of recruitment
##' @export
recplot<-function(x,trans=exp, ...){
    UseMethod("recplot")
}
##' @rdname recplot
##' @method recplot default
##' @export
recplot.default <- function(x,trans=exp,...){
    df <- rectable(x,trans=trans)
    plotit(df,ylab='Recruitment',...)
}

##' CCAM fbar plot
##' @param x the object(s) returned from ccam.fit or forecast
##' @param trans exp by default
##' @param ... extra arguments transferred to plotit
##' @details Plot of mean F
##' @export
fbarplot<-function(x, trans=exp,...){
    UseMethod("fbarplot")
}
##' @rdname fbarplot
##' @method fbarplot default
##' @export
fbarplot.default <- function(x,trans=exp,...){
    df <- fbartable(x,trans=trans)
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
    plotit(df,ylab='Selectivity',xlab='Age',...)+scale_y_continuous(limits=c(0,1.1),expand=c(0,0))
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
    p <- do.call(plotit, c(a, ylab='Catch'))
    if(!is.null(fleet)){
        if(class(x) %in% c('ccam','ccamforecast')) p <- p+geom_line(aes(y=aux1),col='darkgrey')+geom_line(aes(y=aux2),col='darkgrey')
        if(class(x) %in% c('ccamset','forecastset')) p <- p+geom_line(aes(y=aux1,col=fit))+geom_line(aes(y=aux2,col=fit))
    }
    p
}

##' CCAM ssb and catch
##' @param x the object(s) returned from ccam.fit or forecast
##' @param ci logical (confidence interval)
##' @param ... extra arguments transferred to plotit
##' @details Plot of spawning stock biomass
##' @export
scplot<-function(x,ci=TRUE,...){
    UseMethod("scplot")
}
##' @rdname scplot
##' @method scplot default
##' @export
scplot.default <- function(x,ci=TRUE,...){
    df1 <- catchtable(x)
    df1$var<-'catch'
    df2 <- ssbtable(x)
    df2$var<-'ssb'

    refs <- ypr(x)

    df<-rbind(df1,df2)
    p <- ggplot(df,aes(x=year,y=Estimate,group=var))+
        geom_hline(yintercept=refs$f40ssb,col='darkgrey',linetype='dashed')+
        geom_hline(yintercept=refs$f40ssb*0.8,col='darkgrey',linetype='dashed')+
        geom_hline(yintercept=refs$f40ssb*0.4,col='darkgrey',linetype='dashed')+
        geom_line(aes(col=var))+
        scale_color_manual(values=c('darkred','black'))+
        labs(col='')+
        geom_text(x=max(df$year),y=refs$f40ssb,label=expression('SSB'['F40%']),hjust=1,vjust=-0.1,col='darkgrey')+
        geom_text(x=max(df$year),y=refs$f40ssb*0.8,label=expression('SSB'['upp']),hjust=1,vjust=-0.1,col='darkgrey')+
        geom_text(x=max(df$year),y=refs$f40ssb*0.4,label=expression('SSB'['lim']),hjust=1,vjust=-0.1,col='darkgrey')+
        scale_x_continuous(expand=c(0,0))+
        xlab('Year')
    if(ci){
        p <- p+geom_ribbon(aes(ymax=High, ymin=Low,fill=var), alpha=0.2)+
            scale_fill_manual(values=c('darkred','black'))+labs(fill='')+
            scale_y_continuous(expand=c(0,0),limits=c(0,max(df$High)*1.1))
    }else{
        p <-p+scale_y_continuous(expand=c(0,0),limits=c(0,max(df$Estimate)*1.1))
    }
    return(p)
}

# srplot
##' Plots the stock recruitment
##' @param fit the object returned from ccam.fit
##' @param textcol color of years on plot
##' @param linecol color of lines
##' @param curve add SR curve
##' @import ggplot2 viridis
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
srplot.ccamset <- function(fit,text=TRUE,linecol='black',curve=FALSE){
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
            df$SR <- exp(a+log(S)-log(1.0+exp(b)*S)+e*env)
        }
        return(df)
    })
    df <- do.call('rbind',li)

    p <- ggplot(df,aes(x=S,y=R))+geom_point()+
        geom_path()+
        ylab(Rnam)+xlab(Snam)

    if(text){
        p <- p+geom_text(aes(label=y,col=y),hjust=0,vjust=0)+
            scale_color_viridis(direction=-1,discrete=TRUE)+
            theme(legend.position = 'none')
    }

    if(curve){
        p <- p+geom_line(aes(y=SR))
    }
    if(length(fit)>1){
        p<- p+ facet_wrap(~fit)
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
    p <- fit$rep$predObs[idx]
    o <- fit$data$logobs[idx,1]
    res <- o-p
    aa <- fit$data$aux[idx,"age"]

    neg.age <- (aa < -1.0e-6)
    aa[neg.age] <- NA
    Year <- fit$data$aux[idx,"year"]
    sds <- exp(fit$obj$par[which(names(fit$obj$par)=='logSdLogObs')])
    keysd <- fit$conf$keyVarObs[fleets,]
    keysd <- keysd[keysd>-1]
    sd <- sds[keysd+1]
    if(std) res <- res/sd[ifelse(is.na(aa),1,aa)]

    df <- data.frame(year=Year, p=p,o=o,res=res, age=aa)
    if(type %in% c(2,3,4) & !identical(trans,function(x) x)){
        if(type==2) stop('Residuals are on original scale, predicted values should also be')
        if(identical(trans,invcrl)){
            ix <- which(aa == fit$data$minAgePerFleet[fleets])
            ix <- rep(ix,each=length(unique(aa)))
            df <- data.frame(year=rep(unique(Year),each=length(unique(aa))+1),
                             p=unlist(by(p,ix,trans)),
                             o=unlist(by(o,ix,trans)),
                             age=rep(c(unique(aa),max(unique(aa))+1),length(unique(Year))))
        }else{
            df$p <- trans(p)
            df$o <- trans(o)
        }
    }

    if(out!=0){
        e <- sort(abs(res),decreasing = T)[as.numeric(1:out)]
        e <- which(res %in% c(e,-e))
        if(all(is.na(aa[e]))) aprint='' else aprint=paste0(".",aa[e])
        lab <- paste0(Year[e],aprint)
    }

    if(all(is.na(aa))){  #not age structured
        switch (type,
                'one' = {
                    p <- ggplot(df,aes(x=year,y=res))+geom_point()+ylab('Residuals')+xlab('Year')+
                        geom_hline(yintercept = 0,lty='dashed')
                },
                'two' = {
                    p <- ggplot(df,aes(x=p,y=res))+geom_point()+ylab('Residuals')+xlab('Predicted')+
                        geom_hline(yintercept = 0,lty='dashed')
                },
                'three' = {
                    p <- ggplot(df,aes(x=o,y=p))+geom_point()+ylab('Predicted')+xlab('Observed')+
                        geom_abline(intercept=0, slope=1, lty='dashed')
                },
                'four' = {
                    p <- ggplot(df,aes(x=year))+geom_line(aes(y=p))+geom_point(aes(y=o))+
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
                    p <- heat(m,posneg=TRUE,...)
                },
                'two' = {
                    p <- ggplot(df,aes(x=p,y=res,col=age))+geom_point()+ylab('Residuals')+xlab('Predicted')+
                        geom_hline(yintercept = 0,lty='dashed')+
                        labs(col='Age')+scale_color_viridis()
                },
                'three' = {
                    p <- ggplot(df,aes(x=o,y=p,col=age))+geom_point()+ylab('Predicted')+xlab('Observed')+
                        geom_abline(intercept=0, slope=1, lty='dashed')+
                        labs(col='Age')+scale_color_viridis()
                },
                'four' = {
                    p <- ggplot(df,aes(x=year))+geom_line(aes(y=p))+geom_point(aes(y=o),col='darkgrey')+
                        facet_wrap(~age,ncol=3)+
                        ylab('Value')+xlab('Year')
                },
                'five' = {
                    p <- ggplot(df,aes(x=year,y=res))+geom_point()+
                        facet_wrap(~age,ncol=3)+
                        geom_hline(yintercept = 0, lty='dashed')+
                        #geom_smooth(method = "loess",col='black')+
                        ylab('Residuals')+xlab('Year')
                },
                'six' = {
                    p <- ggplot(df,aes(x=year,y=res))+geom_text(aes(label=age))+
                        geom_hline(yintercept = 0, lty='dashed')+
                        ylab('Residuals')+xlab('Year')
                },
                'seven' = {
                    p <- ggplot(df,aes(x=age,y=res))+geom_point()+
                        geom_hline(yintercept = 0, lty='dashed')+
                        ylab('Residuals')+xlab('Age')
                },
                {
                    stop('type should be between 1 and 7 (with age level)')
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
##' @importFrom reshape2 melt
parplot<-function(fit, cor.report.limit=0.95, col=NULL){
    UseMethod("parplot")
}

##' @rdname parplot
##' @method parplot default
##' @export
parplot.default <- function(fit, cor.report.limit=0.95, col=NULL){
    if(class(fit)=='ccam') fit=c(fit)
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

    # corrs <- lapply(1:length(param), function(x) {corrs <- cov2cor(attr(param[[x]], "cov"))-diag(length(param[[x]]))
    #                                     rownames(corrs)<-nam
    #                                     colnames(corrs)<-nam
    #                                     co <- data.frame(melt(corrs))
    #                                     names(co) <- c('nam','label','value')
    #                                     co$hc <- ifelse(co$value>cor.report.limit,round(co$value*100),NA)
    #                                     co$lc <- ifelse(co$value< -cor.report.limit,round(co$value*100),NA)
    #                                     co <- co[!(is.na(co$hc) & is.na(co$lc)),]
    #                                     co$fit <- lab[x]
    #                                     co
    # })
    # print(do.call('rbind',corrs))

    p <- ggplot(mat,aes(x=nam,y=est,col=as.factor(fit),group=as.factor(fit)))+
        geom_point(position=position_dodge(width=0.5))+
        geom_errorbar(aes(ymin=low,ymax=high),position=position_dodge(width=0.5),width=0)+
        xlab('Estimate')+ylab('Parameter')+
        labs(col='')+
        coord_flip()

    if(is.null(col)){
        p<- p+scale_color_viridis(discrete = TRUE)
    }else{
        coll <- c(rep(col,length(fit) %/% length(col)),col[1:length(fit) %% length(col)])
        p <- p+ scale_color_manual(values=coll)
    }

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

     p <- ggplot(df,aes(x=year,y=value,col=as.factor(age)))+geom_line()+geom_point(col='grey30')+
        ylab('Value')+xlab('Year')+labs(col='Age')+
        facet_wrap(~survey,scale='free')
     if(length(unique(df$survey)==1) & length(unique(df$age))){
         p <- p+ theme(legend.position = 'none')+
             scale_color_manual(values='black')
     }else{
         p <- p+scale_color_viridis(discrete = TRUE)
     }
     p
}

##' Plot surplus production vs SSB
##' @param x fit from ccam
##' @import ggplot2 viridis
##' @export
prodplot<-function(x){
    UseMethod("prodplot")
}

##' @rdname prodplot
##' @method prodplot ccam
##' @export
prodplot.ccam <- function(x){
    B <- ssbtable(x)[,1]
    C <- catchtable(x)[,1]
    P <- diff(B) - C[1:(length(C)-1)]
    Endy <- x$data$noYears

    df <- data.frame(year=x$data$years[-Endy],B=B[-Endy],P=P)
    ggplot(df,aes(x=B,y=P))+
        geom_hline(yintercept = 0,linetype='dashed',col='darkgrey')+
        geom_point()+
        geom_text(aes(label=year,col=year),vjust=0,hjust=0)+
        geom_path()+
        geom_text(y=0,x=max(df$B),hjust=1,vjust=-1,label='Production > Catch',col='darkgrey')+
        geom_text(y=0,x=max(df$B),hjust=1,vjust=1.5,label='Production < Catch',col='darkgrey')+
        scale_color_viridis(direction = -1)+
        ylab('Surplus Production')+
        xlab('SSB')+
        theme(legend.position = 'none')
}

##' Kobe plot
##' @param x fit from ccam
##' @import ggplot2
##' @export
##' @details Instead of msy F40% is used as a proxy
kobeplot<-function(x){
    UseMethod("kobeplot")
}

##' @rdname kobeplot
##' @method kobeplot ccam
##' @export
kobeplot.ccam <- function(x,textsize=2,limit=NULL,legend=FALSE){
    refs <- ypr(x)
    Fref <- refs$f40
    SSBref <- refs$f40ssb

    status <- data.frame(x=ssbtable(x)[,1]/SSBref,y=fbartable(x)[,1]/Fref,lab=x$data$years)

    if(is.null(limit)) limit <- ceiling(max(status[,c(1,2)]))
    status[status$y>limit,"y"] <- limit
    status[status$x>limit,"x"] <- limit

    zones <- data.frame(zone=rep(c("Recovery","Overfishing","Fishery reduction","Lightly exploited"),each=4),x=c(0,1,1,0,0,1,1,0,1,limit,limit,1,1,limit,limit,1),y=c(0,0,1,1,1,1,limit,limit,1,1,limit,limit,0,0,1,1))

    p <- ggplot(status,aes(x=x,y=y),environment=environment())+
        geom_polygon(data=zones,aes(x=x,y=y,fill=zone))+
        geom_text(aes(label=lab),size=textsize)+
        geom_path(linetype='dotted')+
        geom_hline(yintercept=1)+geom_vline(xintercept = 1)+
        ylab('F/Fmsy')+xlab('SSB/SSBmsy')+
        scale_fill_manual(values=c("khaki3","olivedrab4","indianred3","khaki1"))+
        scale_y_continuous(breaks=1:limit,labels=c(1:(limit-1),paste0('>',limit)),expand=c(0,0))+
        scale_x_continuous(breaks=1:limit,labels=c(1:(limit-1),paste0('>',limit)),expand=c(0,0))

    if(!legend) p <- p+theme(legend.position = 'none')
}

##' Plot implementation error
##' @param IEmeans list of mean values
##' @param IEsds list of sd values
##' @param col colors
##' @param ylab y label
##' @param linesize linesize
##' @rdname IEplot
##' @import ggplot2
##' @importFrom reshape2 melt
##' @export
IEplot <- function(IEmeans,IEsds, col=NULL,ylab='Undeclared catch (t)',linesize=2){
    myIE <-cbind(melt(do.call('rbind',IEmeans)),sd=melt(do.call('rbind',IEsds))[,3])
    p<- ggplot(myIE,aes(x=Var2-1,y=value))+
        geom_ribbon(aes(ymin=value-sd,ymax=value+sd,fill=Var1),alpha=0.2)+
        geom_line(aes(col=Var1),size=linesize)+
        #geom_point(aes(col=Var1))+
        theme(legend.position="none")+
        ylab(ylab)+xlab('Year')+labs(col='',fill='')

    if(is.null(col)){
        p<- p+scale_color_viridis(discrete = TRUE)+
            scale_fill_viridis(discrete = TRUE)
    }else{
        p <- p+  scale_fill_manual(values=col)+
            scale_color_manual(values=col)
    }

    return(p)
}

##' Observation plot
##' @param x data object
##' @param trans logical. transform data?
##' @param years years to plot
##' @param fleets fleets to plot
##' @param type line, bar or area
##' @import ggplot2 viridis plyr
##' @details plot observations (logobs) by fleet
##' @rdname plotobs
##' @export
plotobs <- function(x,trans=TRUE,years=unique(x$aux[,1]),fleets=unique(x$aux[,2]),type='line',size=2){
    df <- as.data.frame(cbind(x$aux,x$logobs))
    df <- df[df$year %in% years,]
    df <- df[df$fleet %in% fleets,]
    #df[df$age<0,'age'] <- NA

    nam <- attr(x,'fleetNames')
    transcrl <- 'Catch-at-age proportions' %in% nam[fleets]

    if(trans){
        df[df$fleet != 2,c('aux1','aux2')] <- exp(df[df$fleet != 2,c('aux1','aux2')])
        if(transcrl){
            amin <- unname(x$minAgePerFleet[2])
            amax <- unname(x$maxAgePerFleet[2])
            new <- ddply(df[df$fleet == 2,],c('year','fleet','aux2'),summarise,aux1=invcrl(aux1))
            ny <- nrow(df[df$fleet == 2,])/length(amin:amax)
            new$age <- rep(amin:(amax+1),ny)
            new <- new[,colnames(df)]
            df[df$fleet == 2,] <- new[new$age %in% amin:amax,]
            df <- rbind(df,new[new$age==c(amax+1),])
        }
    }
    if(!all(is.na(df$age))) df$age <- factor(df$age, levels = rev(unique(df$age)))
    df$fleet <- factor(df$fleet)
    levels(df$fleet) <- nam[fleets]

    if(type=='line'){
        p <- ggplot(df,aes(x=year,y=aux1,col=age))+geom_line(size=size)+labs(col='')+scale_color_viridis(direction = -1,discrete = TRUE)
    }
    if(type=='area'){
        p <- ggplot(df,aes(x=year,y=aux1,fill=age))+geom_area()+labs(fill='')+scale_fill_viridis(direction = -1,discrete = TRUE)
    }
    if(type=='bar'){
        p <- ggplot(df,aes(x=year,y=aux1,fill=age))+geom_bar(stat='identity')+labs(fill='')+scale_fill_viridis(direction = -1,discrete = TRUE)
    }
    p <- p+facet_wrap(~fleet,scale='free_y',ncol=1)+ylab('Value')+xlab('Year')+ geom_line(data=df,aes(x=year,y=aux2),size=size)
    p
}

##' fitplot
##' @param x fit from ccam
##' @param type 'AIC' (default) or 'nll'
##' @param n logical. Add the number of parameters.
##' @import ggplot2
##' @export
##' @details plots the AIC or nll for every model with the number of parameters on top
fitplot<-function(x,type='AIC',n=TRUE){
    UseMethod("fitplot")
}

##' @rdname fitplot
##' @method fitplot default
##' @export
fitplot.default <- function(x,type=c('AIC','nll'),n=TRUE){
    tab <- as.data.frame(modeltable(x))
    tab$fit <- rownames(tab)
    names(tab)[2] <-'par'

    type <- match.arg(type)
    if(tolower(type)=='aic') names(tab)[3] <- 'y'
    if(tolower(type)=='nll') names(tab)[1] <- 'y'

    p <- ggplot(tab,aes(x=fit,y=y))+geom_point()+
        xlab('Model')+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))+
        ylab(type)
    if(isTRUE(n)){
        p <- p+geom_text(aes(label=par,x=fit,y=y),vjust=-0.3)
    }
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
##' @importFrom reshape2 melt
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

    xx$posneg <- ifelse(xx$value<=0,'#8B000099','#00640099')
    coll <- length(unique(xx$posneg))
    p <- ggplot(xx,aes(x=Var1,y=Var2,size=value,col=posneg))+geom_point(alpha=alpha)+
        scale_size(range = c(1,scale)) +
        ylab(ylab)+xlab(xlab)+
        labs(size="")+
        scale_color_manual(values=col[1:coll])+
        guides(col=FALSE)

    if(length(unique(xx$id))>1){
        p <- p+facet_wrap(~id)
    }
    p
}

##' Heat plot
##' @param x matrix
##' @param high colors that indicate the higher spectrum
##' @param low colors that indicate the lower spectrum
##' @param ncol number of colors to use in plot
##' @param leg logical (plot legend?)
##' @param leground number of digits to round legend
##' @rdname heat
##' @importFrom reshape2 melt
##' @import ggplot2
##' @import viridis
##' @export
heat <- function(x,high=NULL,low=NULL,xlab='Year',ylab='Age',posneg=FALSE){
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
        ylab(ylab)+xlab(xlab)+
        scale_x_continuous(expand=c(0,0))+
        scale_y_continuous(expand=c(0,0))+
        labs(fill='')

    if(is.null(high) & is.null(low)){
      p <- p+scale_fill_viridis()
    }else{
      p <- p+scale_fill_gradient(low=low,high=high)
    }

    if(posneg){
        p <- p+geom_text(aes(label=posneg))
    }
    if(length(unique(xx$id))>1){
        p <- p+facet_wrap(~id)
    }
    p
}

##' prettymatplot plot
##' @param x matrix
##' @param ylab y label
##' @param xlab x label
##' @param legend logical
##' @param col colors
##' @param lwd line with
##' @importFrom reshape2 melt
##' @import ggplot2
##' @import viridis
##' @export
prettymatplot <- function(x,ylab=ifelse(is.null(names(dimnames(x))[1]),'y',names(dimnames(x))[1]),
                          xlab=ifelse(is.null(names(dimnames(x))[1]),'x',names(dimnames(x))[1]),
                          legend=guide_legend(),
                          col=NULL,lwd=1){
    dimnames(x) <- lapply(1:2,function(k){
        nam <- dimnames(x)[[k]]
        if(is.null(nam)|all(duplicated(nam)[2:length(nam)])) {
            warning(paste('No dimnames for dim',k))
            nam <- 1:length(nam)
        }
        return(nam)
    })
    xm <- type.convert(melt(x))
    l <- list(x=xm[,1],y=xm[,3])
    if(length(unique(xm[,2]))>1) l <- c(l,list(col=xm[,2],group=xm[,2]))
    p <- ggplot()+
        geom_line(do.call(aes,l),lwd=lwd)+
        labs(col='',y=ylab,x=xlab)
    if(!is.null(col)) p <- p+scale_color_manual(values=col,guide = legend) else p <- p+scale_color_viridis(guide= legend)
    p
}

##' Plots forecasting results
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
##' @param IE numeric. which IE to use for naming if vector of IEs is used.
##' @param legend names vector of new legend names
##' @param data logical data returned?
##' @param ratio logical. devide by LRP?
##' @details Plot of probability of what over time
##' @importFrom graphics points
##' @import ggplot2 viridis
##' @export
foreplot <- function(x, what.y,what.x=NULL, ylab=what.y,xlab='Year',rect=NULL,ci=TRUE,vline=NULL,hline=NULL,by=NULL,text=FALSE,IE=NULL,year=NULL,legendnames=NULL,data=FALSE,ratio=FALSE){
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
foreplot.forecastset <- function(x, what.y,what.x=NULL, ylab=what.y,xlab='Year',rect=NULL,ci=TRUE,vline=NULL,hline=NULL,by=NULL,text=FALSE,IE=NULL,year=NULL,legendnames=NULL,data=FALSE,ratio=FALSE){

    byopt <- c('OM','IE','MP')
    if(!is.null(by) & !all(by %in% byopt)){stop('"by" may only include OM, IE or MP')}

    # create data frame
    df <- extract(x,what.y,add=TRUE)
    colnames(df)[1]='y'
    if(!is.null(IE)){
        df$IE <- unlist(lapply(strsplit(df$IE,'[.]'),'[[',IE))
    }else{
        byopt <- byopt[-2]
    }
    if(ratio){
        LRP=unlist(lapply(x,function(y) ypr(attr(y,'fit'))$f40ssb*0.4))
        df$y=df$y/LRP
    }
    if(!is.null(year)){
        df <- df[df$year %in% year,]
    }

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

    if(data) return(df)

    # aestethics
    if(length(by)==0){
        p <- ggplot(df,aes(x=x,y=y,col=as.factor(id)))
    }
    if(length(by)==1){
        int <- paste0("interaction(", paste0(byopt[byopt!=by], collapse =  ", "), ")")
        p <- ggplot(df,aes_string(x='x',y='y',col=int))+facet_wrap(reformulate(by) )
    }
    if(length(by)==2){
        int <- byopt[!byopt %in% by]
        p <- ggplot(df,aes_string(x='x',y='y',col=int))+facet_grid(reformulate(by[1],by[2]) )
    }
    if(length(by)==3){
        p <- ggplot(df,aes(x=x,y=y))+facet_grid(reformulate(by[1],by[2:3])  )
    }
    # rectangle
    if(!is.null(rect)) p <- p+ geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=rect, ymax=Inf),fill='lightgrey',alpha=0.3,col='white')

    # basis
    p <- p+ geom_point()+
            geom_path()+
            ylab(ylab)+xlab(xlab)+
            labs(col='')

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
    }else{
        if(length(by)<3) {
            if(!is.null(legendnames)) p <- p+scale_color_viridis(discrete = TRUE,labels=legendnames) else
             p <- p+scale_color_viridis(discrete = TRUE)
        }
    }
    return(p)
}

##' Diamond plots of MSE output
##' @param x the object(s) returned from ccam.forecast
##' @param what statistic returned from ccam.forecast
##' @param ylab y label
##' @param xlab x label
##' @param hline plot horizontal line(s)
##' @param area color plot underneath this yintercept
##' @param year year for which the statistics needs to be plotted
##' \itemize{
##'   \item{NULL}{ Takes median all years}
##'   \item{'numeric vector'}{ Takes median of specified years}
##'   \item{'threshold'}{ Takes year at which threshold is reached (see threshold param)}
##'   \item{'percyear'}{ Takes percentage of years where threshold is reached (see threshold param)}
##'   \item{'zone'}{ Makes a facet grid in which the statistic is grouped by zone (critical, healthy or cautious)}
##' }
##' @param IE Which IE to plot (if vector)? Defaults to collapsed IE. Can be numeric or string (grep).
##' @param IEnames character string containing alternative names for the IEs
##' @param legnames character string containing alternative names for the legend (OM types)
##' @param threshold if year=TRUE than a threshold can be specified and the year at which it is reached becomes the y value for 'what'
##' @param data logical. return data instead of plot.
##' @param ratio logical. devide by LRP?
##' @param OMtype logical. OMbase/core/whatever present?
##' @details performance (what) by management procedure, IE and OM type. OM type is defined by the OM name (numbers and 'OM' are removed).
##' @importFrom plyr ddply
##' @import ggplot2 viridis
##' @export
diamondplot <- function(x, what, ylab=what, xlab='Management Procedure',hline=NULL, area=NULL,year, IE=NULL,IEnames=NULL,legnames=NULL,threshold=NULL,data=FALSE,ratio=FALSE,OMtype=FALSE){
    UseMethod("diamondplot")
}
##' @rdname diamondplot
##' @method diamondplot ccamforecast
##' @export
diamondplot.ccamforecast <- function(x,...){
    x=c(x)
    diamondplot(x,...)
}

##' @rdname diamondplot
##' @method diamondplot forecastset
##' @export
diamondplot.forecastset <- function(x, what, ylab=what, xlab='Management Procedure',hline=NULL, area=NULL,year='all',IE=NULL,IEnames=NULL,legnames=NULL,threshold=NULL,data=FALSE,ratio=FALSE,OMtype=FALSE){

    # create data frame
    df <- extract(x,what,add=TRUE)
    colnames(df)[1]='y'
    if(!is.null(IE)){
        if(is.numeric(IE)) df$IE <- unlist(lapply(strsplit(df$IE,'[.]'),'[[',IE))
        if(is.character(IE)) df$IE <-  unlist(lapply(strsplit(df$IE,'[.]'),function(o) {
                                                    n <- o[grepl(IE,o)]
                                                     if(length(n)==0) n <- 'IE0'
                                                     return(n)
                                                     }))
    }
    if(ratio){
        LRP=unlist(lapply(x,function(y) ypr(attr(y,'fit'))$LRP))
        df$y=df$y/LRP
    }

    #clean up
    MPnum <- as.numeric(unlist(lapply(strsplit(df$MP, "\\D+"),'[[',2)))
    if(length(MPnum)>0) df$MP <- MPnum
    df$IE <- as.factor(df$IE)
    if(!is.null(IEnames)) {
         levels(df$IE) <- IEnames
        if(!is.null(attr(IEnames,'order'))) {
           df$IE <- factor(df$IE,levels(df$IE)[attr(IEnames,'order')])
        }
    }

    inter <- ifelse(ncol(df)==6,FALSE,TRUE)
    if(inter){
        colnames(df)[2:3]=c('ylow','yhigh')
    }

    df$type <- gsub("OM","",df$OM)
    if(OMtype) df$type <- gsub('[0-9]+', '', df$type)
    df$type <- as.factor(df$type)
    if(!is.null(legnames)) levels(df$type) <- legnames


    colnames(df)[which(colnames(df)=='year')]='x'
    df <- df[df$x!=min(df$x),] #remove the first year because not part of the future

    if(is.numeric(year)){ #take average specified years
        df <- df[df$x %in% year,]
        yearlab <- range(df$x)
        if(yearlab[1]==yearlab[2]){yearlab=as.character(yearlab[1])}else{yearlab <- paste(yearlab,collapse = '-')}
        df <- ddply(df,c('OM','MP','IE','id','type'),summarise,y=median(y))
    }else{
    if(year=='all'){ # take average all years
        yearlab <- paste(range(df$x),collapse = '-')
        df <- ddply(df,c('OM','MP','IE','id','type'),summarise,y=median(y))
        message('Taking the average value')
    }
    if(year=='zone'){
        cz <-  extract(x,'probCZ',add=TRUE)
        cz <- cz[cz$year!=min(cz$year),]
        hz <- extract(x,'probHZ',add=TRUE)
        hz <- hz[hz$yeqr!=min(hz$year),]
        df$zones <- ifelse(cz[,1]<0.75,'Critical Zone',ifelse(hz[,1]>0.75,'Healthy Zone','Cautious Zone'))
        perc <- table(df$zones)/sum(table(df$zones))
        df <- df[!df$zones %in% names(perc)[which(perc<0.1)],] #remove zones if it doesn't get there
    }
    if(year=='threshold'){ # if TRUE take year that meets threshold
        if(is.null(threshold)){stop('If year=threshold a threshold should be defined')}
        m <- min(as.numeric(df$x))-1
        df$x <- as.numeric(df$x)
        df <- ddply(df,c('OM','MP','IE','id','type'),summarise,y=ifelse(any(y>=threshold),min(x[which(y>=threshold)]),max(x)))
        df$y <- df$y-m
        yearlab <-''
    }
    if(year=='yearperc'){ # if TRUE take percentage of year that are above threshold
        if(is.null(threshold)){stop('If year=percyear a threshold should be defined')}
        df <- ddply(df,c('OM','MP','IE','id','type'),summarise,y=length(y[y>=threshold])/length(y)*100)
        yearlab <-''
    }
    if(year=='zone'){
        df <- ddply(df,c('OM','MP','IE','id','type','zones'),summarise,y=median(y))
        yearlab <-''
    }
    }

    if(data) return(df)

    p <- ggplot(df,aes(x=MP,y=y,col=type,shape=type))

    if(!is.null(hline)) p <- p+ geom_hline(yintercept=hline,linetype='dashed',col='darkgreen')
    if(!is.null(area)) p <- p+ geom_area(aes(y=area),fill='darkgreen',alpha=0.5,col='darkgreen')

    p <- p+geom_point(size=1.5)+labs(col='OM',shape='OM')+
        ylab(ylab)+xlab(xlab)+
        ggtitle(yearlab)+
        scale_x_continuous(labels = as.character(unique(df$MP)), breaks = unique(df$MP))

    if(!is.null(IE)){
        if(year=='zone'){ p <- p+facet_grid(zones~IE)}else{p <- p+facet_wrap(~IE)}
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
        strip <- gtable::gtable_filter(z, "strip-t", trim = FALSE)
        ix <- 1:nrow(strip$layout)
        top <- strip$layout$t[1]
        l   <- strip$layout$l[seq(1,max(ix),by=n2) ]
        r   <- strip$layout$r[seq(n2,max(ix),by=n2) ]
        mat   <- matrix(vector("list", length = (n1*2-1)*2), nrow = 2)
        mat[] <- list(zeroGrob())
        res <- gtable::gtable_matrix("toprow", mat, unit(c(rep(c(1,0),n1-1),1), "null"), unit(c(1, 1), "null"))
        se <- seq(1,length(locations),by=n2)
        for(i in se){
            zz. <- gtable::gtable_add_grob(res,z$grobs[[locations[i]]]$grobs[[1]], 1, 1, 1, (n1)*2-1)
            z <- gtable::gtable_add_grob(z, zz., t = top,  l = l[which(se==i)],  b = top,  r = r[which(se==i)], name = c("add-strip"))
        }
        grid::grid.newpage()
        print(grid::grid.draw(z))
    }
    # 1) by next 5 year 75 % probCZ
    df <- extract(x,'probCZ',add=TRUE)
    df <- df[df$year==(min(df$year)+5),]
    df$ylog <- ifelse(df$var>=0.75,TRUE,FALSE)
    df <- df[,!(colnames(df) %in% c('var','year'))]
    df$ytext <- NA
    df$performance <- "in 5 year 75 % >CZ"
    out <- df
    # 2) by next 10 years 75% probHZ
    df <- extract(x,'probHZ',add=TRUE)
    df <- df[df$year==(min(df$year)+10),]
    df$ylog <- ifelse(df$var>=0.75,TRUE,FALSE)
    df <- df[,!(colnames(df) %in% c('var','year'))]
    df$performance <- "in 10 year 75 % >HZ"
    df$ytext <- NA
    out <- rbind(out,df)
    # 3) never below 50% chance of 30% growth when in the critical zone (95%)
    df <- extract(x,'probgrowth20',add=TRUE)
    names(df)[1] <- 'probgrowth20'
    df$probgrowth <- extract(x,'probgrowth',add=TRUE)[,1]
    df$probCZ <- extract(x,'probCZ',add=TRUE)[,1]
    df$probHZ <- extract(x,'probHZ',add=TRUE)[,1]

    df <- df[df$year != min(df$year),]
    df$zone <- 'cautious'
    df[df$probCZ<0.25,'zone'] <- 'critical'
    df[df$probHZ>0.75,'zone'] <- 'healthy'

    df$ylog <- TRUE
    df[df$zone=='cautious' & df$probgrowth<0.50,'ylog'] <- FALSE
    df[df$zone=='critical' & df$probgrowth<0.50,'ylog'] <- FALSE

    df <- ddply(df,c('OM','MP','IE','id'),summarise,ylog=all(ylog==TRUE)) #all years need to adhere to the growth rule

    df$performance <- "ANY growth (>50%) when in CZ or MZ (75%)"
    df$ytext=NA
    out <- rbind(out,df)
    # 4) show number of years TAC at or above 8000t
    df <- extract(x,'TAC',add=TRUE)
    colnames(df)[1] <- 'var'
    df <- ddply(df[!df$year==min(df$year),],c('id','OM','MP','IE'),summarise,ytext=length(var[var>=8000]))
    df$ylog <- NA
    df$performance <- 'n years TAC >= 8000t'
    out <- rbind(out,df)
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


