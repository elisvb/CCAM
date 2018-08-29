###########################################################################################################
## Functions for different management procedures                                                         ##
##                                                                                                       ##
## Function names should start with 'MP'                                                                 ##
## They should be off class "MP" so they can be found by the forecast function                           ##
## There are different types here, developped for Atlantic mackerel:                                  ##
###########################################################################################################

##' MPeggsimple
##' @param data list with data objects
##' @param TAC.base TAC of previous year
##' @param i year in the future
##' @rdname MPeggsimple
##' @details simply changes TAC in function of change in I (last year/ geometric mean 3 years before). Change can not be bigger smaller than 0.5x or 2x.
##' @export
MPeggsimple<-function(data,TAC.base){
        minbound <- 0.75
        maxbound <- 1.25
        change <- unlist(lapply(data,function(x){ Ix <- tail(x$logobs[which(x$aux[,2]==3),1],4)
                                                change <- Ix[4]/gmean(Ix[1:3])
                                                change <- min(max(change,minbound),maxbound)
                                                return(change)}))
        TAC <- TAC.base*change
    return(TAC)
}
class(MPeggsimple) <- append(class(MPeggsimple),"MP")
attr(MPeggsimple,'model')=FALSE

##' MPeggcomplex
##' @param data list with data objects
##' @param TAC.base TAC of previous year
##' @param i year in the future
##' @rdname MPeggcomplex
##' @details Based on the Gerromont and Butterworth rule.
##' @export
MPeggcomplex<-function(data,TAC.base){
    Ihist <- cbind(I=exp(data[[1]]$logobs[,1]),data[[1]]$aux)
    Ihist <- Ihist[which(Ihist[,3]==3),c(1,2)]

    Ytarget <- c(1979,1983,1984,1987,1988,1989,1990,1992,1993,1994,1995)
    Y0 <- c(1996:2000,2004:2008)
    Ycurrent <- unname(tail(Ihist[,2],3))

    Itarget <- gmean(Ihist[which(Ihist[,2] %in% Ytarget),1])
    I0 <- gmean(Ihist[which(Ihist[,2] %in% Y0),1])
    TACtarget <- 25000
    TACbase <- 8000
    TACmin <- 8000
    TACmax <- 25000
    w = TACbase/TACtarget
    Irecent <- exp(unlist(lapply(data,function(x){gmean(tail(x$logobs[which(x$aux[,2]==3),1],3))})))

    ItoTac <- function(x){
        if(x<I0){
            TAC <- w*TACtarget*(x/I0)^2
        }else{
            TAC <- TACtarget*(w+(1-w)*((x-I0)/(Itarget-I0)))
        }
        TAC <- ifelse(TAC<TACmin,0,TAC)
        TAC <- ifelse(TAC>TACmax,TACmax,TAC)
        return(TAC)
    }

    TAC <- sapply(Irecent, ItoTac)
    if(length(TAC)==1) TAC <- rep(TAC,length(TAC.base))
    return(TAC)

    # MPdf=data.frame(Ihist)
    # MPdf[MPdf$year %in% Ytarget,'class'] <- 'Itarget'
    # MPdf[MPdf$year %in% Y0,'class'] <- 'I0'
    # MPdf[MPdf$year %in% Ycurrent,'class'] <- 'Iy'
    #
    # p1 <- ggplot(MPdf,aes(x=year,y=I))+geom_line()+geom_point(aes(col=class),size=2)+labs(col='',y=expression(SSB["egg"] (t)),x='Year')+
    #     theme(legend.position = c(0.8,0.8),
    #           legend.background = element_rect(fill=alpha('white', 0)))+
    #     geom_hline(yintercept = I0,col='red',linetype='dashed')+
    #     geom_hline(yintercept = Itarget,col='green',linetype='dashed')+
    #     geom_hline(yintercept = Irecent,col='blue',linetype='dotted')+
    #     scale_color_manual(values=c('red','green','blue'),breaks=c("I0", "Itarget",'Iy'))+
    #     scale_y_continuous(limits=c(0,max(MPdf$I)*1.05),expand=c(0,0))+
    #     geom_text(label='A)',y=Inf,x=-Inf,hjust=-0.3,vjust=1.3)
    #
    # MPdf2 <- data.frame(x=seq(0,Itarget*1.2,length.out = 500),y=sapply(seq(0,Itarget*1.2,length.out = 500), ItoTac))
    #
    # p2 <- ggplot(MPdf2,aes(x=x,y=y))+geom_line(size=2)+xlab(expression(SSB["egg"] (t)))+ylab('TAC')+
    #     geom_hline(yintercept = TACtarget)+geom_text(x=(Itarget-I0)/2,y=TACtarget,label='TAC target',vjust=-0.4,hjust=0.2)+
    #     geom_hline(yintercept = w*TACtarget)+geom_text(x=(Itarget-I0)/2,y=w*TACtarget,label='TAC base',vjust=-0.4,hjust=0.2)+
    #     geom_vline(xintercept=Itarget,col='green',linetype='dashed',size=1.5)+geom_text(x=Itarget,y=0,label='Itarget',vjust=-0.4,hjust=-0.2,col='green')+
    #     geom_vline(xintercept=I0,col='red',linetype='dashed',size=1.5)+  geom_text(x=I0,y=0,label='I0',vjust=-0.4,hjust=-0.4,col='red')+
    #     geom_vline(xintercept=Irecent,col='blue',linetype='dotted',size=1.5)+geom_text(x=Irecent,y=0,label='Iy',vjust=-0.4,hjust=-0.2,col='blue')+
    #     scale_y_continuous(expand=c(0,0),limits=c(0,30000))+scale_x_continuous(expand=c(0,0))+
    #     geom_text(label='B)',y=Inf,x=-Inf,hjust=-0.3,vjust=1.3)
    #
    # library(gridExtra)
    # savepng(grid.arrange(p1,p2,ncol=2),wdIMG,"/HCR/HCR_target",c(30,7))

}
class(MPeggcomplex) <- append(class(MPeggcomplex),"MP")
attr(MPeggcomplex,'model')=FALSE

##' MPeggcomplexinterim
##' @param data list with data objects
##' @param TAC.base TAC of previous year
##' @param i year in the future
##' @rdname MPeggcomplexinterim
##' @details Based on the Gerromont and Butterworth rule.
##' @export
MPeggcomplexinterim<-function(data,TAC.base){
    if(!i %% 2 == 1){return(TAC.base)}

    Ihist <- cbind(I=exp(data[[1]]$logobs[,1]),data[[1]]$aux)
    Ihist <- Ihist[which(Ihist[,3]==3),c(1,2)]

    Ytarget <- c(1979,1983,1984,1987,1988,1989,1990,1992,1993,1994,1995)
    Y0 <- c(1996:2000,2004:2008)
    Ycurrent <- unname(tail(Ihist[,2],3))

    Itarget <- gmean(Ihist[which(Ihist[,2] %in% Ytarget),1])
    I0 <- gmean(Ihist[which(Ihist[,2] %in% Y0),1])
    TACtarget <- 25000
    TACbase <- 8000
    TACmin <- 8000
    TACmax <- 25000
    w = TACbase/TACtarget
    Irecent <- exp(unlist(lapply(data,function(x){gmean(tail(x$logobs[which(x$aux[,2]==3),1],3))})))

    ItoTac <- function(x){
        if(x<I0){
            TAC <- w*TACtarget*(x/I0)^2
        }else{
            TAC <- TACtarget*(w+(1-w)*((x-I0)/(Itarget-I0)))
        }
        TAC <- ifelse(TAC<TACmin,0,TAC)
        TAC <- ifelse(TAC>TACmax,TACmax,TAC)
        return(TAC)
    }

    TAC <- sapply(Irecent, ItoTac)
    if(length(TAC)==1) TAC <- rep(TAC,length(TAC.base))
    return(TAC)
}
class(MPeggcomplexinterim) <- append(class(MPeggcomplexinterim),"MP")
attr(MPeggcomplexinterim,'model')=FALSE

##' MPeggcomplexramp
##' @param data list with data objects
##' @param TAC.base TAC of previous year
##' @param i year in the future
##' @rdname MPeggcomplexramp
##' @details Based on the Gerromont and Butterworth rule.
##' @export
MPeggcomplexramp<-function(data,TAC.base){
    Ihist <- cbind(I=exp(data[[1]]$logobs[,1]),data[[1]]$aux)
    Ihist <- Ihist[which(Ihist[,3]==3),c(1,2)]

    Ytarget <- c(1979,1983,1984,1987,1988,1989,1990,1992,1993,1994,1995)
    Y0 <- c(1996:2000,2004:2008)
    Ycurrent <- unname(tail(Ihist[,2],3))

    Itarget <- gmean(Ihist[which(Ihist[,2] %in% Ytarget),1])
    I0 <- gmean(Ihist[which(Ihist[,2] %in% Y0),1])
    TACtarget <- 25000
    TACbase <- 8000
    TACmin <- 0
    TACmax <- 25000
    w = TACbase/TACtarget
    Irecent <- exp(unlist(lapply(data,function(x){gmean(tail(x$logobs[which(x$aux[,2]==3),1],3))})))

    ItoTac <- function(x){
        if(x<I0){
            TAC <- w*TACtarget*(x/I0)^3
        }else{
            TAC <- TACtarget*(w+(1-w)*((x-I0)/(Itarget-I0)))
        }
        TAC <- ifelse(TAC<TACmin,0,TAC)
        TAC <- ifelse(TAC>TACmax,TACmax,TAC)
        return(TAC)
    }

    TAC <- sapply(Irecent, ItoTac)
    if(length(TAC)==1) TAC <- rep(TAC,length(TAC.base))
    return(TAC)

    # MPdf=data.frame(Ihist)
    # MPdf[MPdf$year %in% Ytarget,'class'] <- 'Itarget'
    # MPdf[MPdf$year %in% Y0,'class'] <- 'I0'
    # MPdf[MPdf$year %in% Ycurrent,'class'] <- 'Iy'
    #
    # p1 <- ggplot(MPdf,aes(x=year,y=I))+geom_line()+geom_point(aes(col=class),size=2)+labs(col='',y=expression(SSB["egg"] (t)),x='Year')+
    #     theme(legend.position = c(0.8,0.8),
    #           legend.background = element_rect(fill=alpha('white', 0)))+
    #     geom_hline(yintercept = I0,col='red',linetype='dashed')+
    #     geom_hline(yintercept = Itarget,col='green',linetype='dashed')+
    #     geom_hline(yintercept = Irecent,col='blue',linetype='dotted')+
    #     scale_color_manual(values=c('red','green','blue'),breaks=c("I0", "Itarget",'Iy'))+
    #     scale_y_continuous(limits=c(0,max(MPdf$I)*1.05),expand=c(0,0))+
    #     geom_text(label='A)',y=Inf,x=-Inf,hjust=-0.3,vjust=1.3)
    #
    # MPdf2 <- data.frame(x=seq(0,Itarget*1.2,length.out = 500),y=sapply(seq(0,Itarget*1.2,length.out = 500), ItoTac))
    #
    # p2 <- ggplot(MPdf2,aes(x=x,y=y))+geom_line(size=2)+xlab(expression(SSB["egg"] (t)))+ylab('TAC')+
    #     geom_hline(yintercept = TACtarget)+geom_text(x=(Itarget-I0)/2,y=TACtarget,label='TAC target',vjust=-0.4,hjust=0.2)+
    #     geom_hline(yintercept = w*TACtarget)+geom_text(x=(Itarget-I0)/2,y=w*TACtarget,label='TAC base',vjust=-0.4,hjust=0.2)+
    #     geom_vline(xintercept=Itarget,col='green',linetype='dashed',size=1.5)+geom_text(x=Itarget,y=0,label='Itarget',vjust=-0.4,hjust=-0.2,col='green')+
    #     geom_vline(xintercept=I0,col='red',linetype='dashed',size=1.5)+  geom_text(x=I0,y=0,label='I0',vjust=-0.4,hjust=-0.4,col='red')+
    #     geom_vline(xintercept=Irecent,col='blue',linetype='dotted',size=1.5)+geom_text(x=Irecent,y=0,label='Iy',vjust=-0.4,hjust=-0.2,col='blue')+
    #     scale_y_continuous(expand=c(0,0),limits=c(0,30000))+scale_x_continuous(expand=c(0,0))+
    #     geom_text(label='B)',y=Inf,x=-Inf,hjust=-0.3,vjust=1.3)
    #
    # library(gridExtra)
    # savepng(grid.arrange(p1,p2,ncol=2),wdIMG,"/HCR/HCR_target_ramp",c(30,7))

}
class(MPeggcomplexramp) <- append(class(MPeggcomplexramp),"MP")
attr(MPeggcomplexramp,'model')=FALSE

##' MPeggcomplex0
##' @param data list with data objects
##' @param TAC.base TAC of previous year
##' @param i year in the future
##' @rdname MPeggcomplex0
##' @details Based on the Gerromont and Butterworth rule.
##' @export
MPeggcomplex0<-function(data,TAC.base){
    Ihist <- cbind(I=exp(data[[1]]$logobs[,1]),data[[1]]$aux)
    Ihist <- Ihist[which(Ihist[,3]==3),c(1,2)]

    Ytarget <- c(1979,1983,1984,1987,1988,1989,1990,1992,1993,1994,1995)
    Y0 <- c(1996:2000,2004:2008)
    Ycurrent <- unname(tail(Ihist[,2],3))

    Itarget <- gmean(Ihist[which(Ihist[,2] %in% Ytarget),1])
    I0 <- gmean(Ihist[which(Ihist[,2] %in% Y0),1])
    TACtarget <- 25000
    TACbase <- 0
    TACmin <- 0
    TACmax <- 25000
    w = TACbase/TACtarget
    Irecent <- exp(unlist(lapply(data,function(x){gmean(tail(x$logobs[which(x$aux[,2]==3),1],3))})))

    ItoTac <- function(x){
        if(x<I0){
            TAC <- w*TACtarget*(x/I0)^2
        }else{
            TAC <- TACtarget*(w+(1-w)*((x-I0)/(Itarget-I0)))
        }
        TAC <- ifelse(TAC<TACmin,0,TAC)
        TAC <- ifelse(TAC>TACmax,TACmax,TAC)
        return(TAC)
    }

    TAC <- sapply(Irecent, ItoTac)
    if(length(TAC)==1) TAC <- rep(TAC,length(TAC.base))
    return(TAC)

    # MPdf=data.frame(Ihist)
    # MPdf[MPdf$year %in% Ytarget,'class'] <- 'Itarget'
    # MPdf[MPdf$year %in% Y0,'class'] <- 'I0'
    # MPdf[MPdf$year %in% Ycurrent,'class'] <- 'Iy'
    #
    # p1 <- ggplot(MPdf,aes(x=year,y=I))+geom_line()+geom_point(aes(col=class),size=2)+labs(col='',y=expression(SSB["egg"] (t)),x='Year')+
    #     theme(legend.position = c(0.8,0.8),
    #           legend.background = element_rect(fill=alpha('white', 0)))+
    #     geom_hline(yintercept = I0,col='red',linetype='dashed')+
    #     geom_hline(yintercept = Itarget,col='green',linetype='dashed')+
    #     geom_hline(yintercept = Irecent,col='blue',linetype='dotted')+
    #     scale_color_manual(values=c('red','green','blue'),breaks=c("I0", "Itarget",'Iy'))+
    #     scale_y_continuous(limits=c(0,max(MPdf$I)*1.05),expand=c(0,0))+
    #     geom_text(label='A)',y=Inf,x=-Inf,hjust=-0.3,vjust=1.3)
    #
    # MPdf2 <- data.frame(x=seq(0,Itarget*1.2,length.out = 500),y=sapply(seq(0,Itarget*1.2,length.out = 500), ItoTac))
    #
    # p2 <- ggplot(MPdf2,aes(x=x,y=y))+geom_line(size=2)+xlab(expression(SSB["egg"] (t)))+ylab('TAC')+
    #     geom_hline(yintercept = TACtarget)+geom_text(x=(Itarget-I0)/2,y=TACtarget,label='TAC target',vjust=-0.4,hjust=0.2)+
    #     geom_hline(yintercept = w*TACtarget)+geom_text(x=(Itarget-I0)/2,y=w*TACtarget,label='TAC base',vjust=-0.4,hjust=0.2)+
    #     geom_vline(xintercept=Itarget,col='green',linetype='dashed',size=1.5)+geom_text(x=Itarget,y=0,label='Itarget',vjust=-0.4,hjust=-0.2,col='green')+
    #     geom_vline(xintercept=I0,col='red',linetype='dashed',size=1.5)+  geom_text(x=I0,y=0,label='I0',vjust=-0.4,hjust=-0.4,col='red')+
    #     geom_vline(xintercept=Irecent,col='blue',linetype='dotted',size=1.5)+geom_text(x=Irecent,y=0,label='Iy',vjust=-0.4,hjust=-0.2,col='blue')+
    #     scale_y_continuous(expand=c(0,0),limits=c(0,30000))+scale_x_continuous(expand=c(0,0))+
    #     geom_text(label='B)',y=Inf,x=-Inf,hjust=-0.3,vjust=1.3)
    #
    # library(gridExtra)
    # savepng(grid.arrange(p1,p2,ncol=2),wdIMG,"/HCR/HCR_target_0",c(30,7))

}
class(MPeggcomplex0) <- append(class(MPeggcomplex0),"MP")
attr(MPeggcomplex0,'model')=FALSE

##' MPeggcomplex8000
##' @param data list with data objects
##' @param TAC.base TAC of previous year
##' @param i year in the future
##' @rdname MPeggcomplex8000
##' @details Based on the Gerromont and Butterworth rule.
##' @export
MPeggcomplex8000<-function(data,TAC.base){
    Ihist <- cbind(I=exp(data[[1]]$logobs[,1]),data[[1]]$aux)
    Ihist <- Ihist[which(Ihist[,3]==3),c(1,2)]

    Ytarget <- c(1979,1983,1984,1987,1988,1989,1990,1992,1993,1994,1995)
    Y0 <- c(1996:2000,2004:2008)
    Ycurrent <- unname(tail(Ihist[,2],3))

    Itarget <- gmean(Ihist[which(Ihist[,2] %in% Ytarget),1])
    I0 <- gmean(Ihist[which(Ihist[,2] %in% Y0),1])
    TACtarget <- 25000
    TACbase <- 8000
    TACmin <- 8000
    TACfloor <- 8000
    TACmax <- 25000
    w = TACbase/TACtarget
    Irecent <- exp(unlist(lapply(data,function(x){gmean(tail(x$logobs[which(x$aux[,2]==3),1],3))})))

    ItoTac <- function(x){
        if(x<I0){
            TAC <- w*TACtarget*(x/I0)^2
        }else{
            TAC <- TACtarget*(w+(1-w)*((x-I0)/(Itarget-I0)))
        }
        TAC <- ifelse(TAC<TACmin,TACfloor,TAC)
        TAC <- ifelse(TAC>TACmax,TACmax,TAC)
        return(TAC)
    }

    TAC <- sapply(Irecent, ItoTac)
    if(length(TAC)==1) TAC <- rep(TAC,length(TAC.base))
    return(TAC)

    # MPdf=data.frame(Ihist)
    # MPdf[MPdf$year %in% Ytarget,'class'] <- 'Itarget'
    # MPdf[MPdf$year %in% Y0,'class'] <- 'I0'
    # MPdf[MPdf$year %in% Ycurrent,'class'] <- 'Iy'
    #
    # p1 <- ggplot(MPdf,aes(x=year,y=I))+geom_line()+geom_point(aes(col=class),size=2)+labs(col='',y=expression(SSB["egg"] (t)),x='Year')+
    #     theme(legend.position = c(0.8,0.8),
    #           legend.background = element_rect(fill=alpha('white', 0)))+
    #     geom_hline(yintercept = I0,col='red',linetype='dashed')+
    #     geom_hline(yintercept = Itarget,col='green',linetype='dashed')+
    #     geom_hline(yintercept = Irecent,col='blue',linetype='dotted')+
    #     scale_color_manual(values=c('red','green','blue'),breaks=c("I0", "Itarget",'Iy'))+
    #     scale_y_continuous(limits=c(0,max(MPdf$I)*1.05),expand=c(0,0))+
    #     geom_text(label='A)',y=Inf,x=-Inf,hjust=-0.3,vjust=1.3)
    #
    # MPdf2 <- data.frame(x=seq(0,Itarget*1.2,length.out = 500),y=sapply(seq(0,Itarget*1.2,length.out = 500), ItoTac))
    #
    # p2 <- ggplot(MPdf2,aes(x=x,y=y))+geom_line(size=2)+xlab(expression(SSB["egg"] (t)))+ylab('TAC')+
    #     geom_hline(yintercept = TACtarget)+geom_text(x=(Itarget-I0)/2,y=TACtarget,label='TAC target',vjust=-0.4,hjust=0.2)+
    #     geom_hline(yintercept = w*TACtarget)+geom_text(x=(Itarget-I0)/2,y=w*TACtarget,label='TAC base',vjust=-0.4,hjust=0.2)+
    #     geom_vline(xintercept=Itarget,col='green',linetype='dashed',size=1.5)+geom_text(x=Itarget,y=0,label='Itarget',vjust=-0.4,hjust=-0.2,col='green')+
    #     geom_vline(xintercept=I0,col='red',linetype='dashed',size=1.5)+  geom_text(x=I0,y=0,label='I0',vjust=-0.4,hjust=-0.4,col='red')+
    #     geom_vline(xintercept=Irecent,col='blue',linetype='dotted',size=1.5)+geom_text(x=Irecent,y=0,label='Iy',vjust=-0.4,hjust=-0.2,col='blue')+
    #     scale_y_continuous(expand=c(0,0),limits=c(0,30000))+scale_x_continuous(expand=c(0,0))+
    #     geom_text(label='B)',y=Inf,x=-Inf,hjust=-0.3,vjust=1.3)
    #
    # library(gridExtra)
    # savepng(grid.arrange(p1,p2,ncol=2),wdIMG,"/HCR/HCR_target_8000",c(30,7))

}
class(MPeggcomplex8000) <- append(class(MPeggcomplex8000),"MP")
attr(MPeggcomplex8000,'model')=FALSE

##' MPeggcomplex10000
##' @param data list with data objects
##' @param TAC.base TAC of previous year
##' @param i year in the future
##' @rdname MPeggcomplex10000
##' @details Based on the Gerromont and Butterworth rule.
##' @export
MPeggcomplex10000<-function(data,TAC.base){
    Ihist <- cbind(I=exp(data[[1]]$logobs[,1]),data[[1]]$aux)
    Ihist <- Ihist[which(Ihist[,3]==3),c(1,2)]

    Ytarget <- c(1979,1983,1984,1987,1988,1989,1990,1992,1993,1994,1995)
    Y0 <- c(1996:2000,2004:2008)
    Ycurrent <- unname(tail(Ihist[,2],3))

    Itarget <- gmean(Ihist[which(Ihist[,2] %in% Ytarget),1])
    I0 <- gmean(Ihist[which(Ihist[,2] %in% Y0),1])
    TACtarget <- 25000
    TACbase <- 10000
    TACmin <- 10000
    TACfloor <- 10000
    TACmax <- 25000
    w = TACbase/TACtarget
    Irecent <- exp(unlist(lapply(data,function(x){gmean(tail(x$logobs[which(x$aux[,2]==3),1],3))})))

    ItoTac <- function(x){
        if(x<I0){
            TAC <- w*TACtarget*(x/I0)^2
        }else{
            TAC <- TACtarget*(w+(1-w)*((x-I0)/(Itarget-I0)))
        }
        TAC <- ifelse(TAC<TACmin,TACfloor,TAC)
        TAC <- ifelse(TAC>TACmax,TACmax,TAC)
        return(TAC)
    }

    TAC <- sapply(Irecent, ItoTac)
    if(length(TAC)==1) TAC <- rep(TAC,length(TAC.base))
    return(TAC)

    # MPdf=data.frame(Ihist)
    # MPdf[MPdf$year %in% Ytarget,'class'] <- 'Itarget'
    # MPdf[MPdf$year %in% Y0,'class'] <- 'I0'
    # MPdf[MPdf$year %in% Ycurrent,'class'] <- 'Iy'
    #
    # p1 <- ggplot(MPdf,aes(x=year,y=I))+geom_line()+geom_point(aes(col=class),size=2)+labs(col='',y=expression(SSB["egg"] (t)),x='Year')+
    #     theme(legend.position = c(0.8,0.8),
    #           legend.background = element_rect(fill=alpha('white', 0)))+
    #     geom_hline(yintercept = I0,col='red',linetype='dashed')+
    #     geom_hline(yintercept = Itarget,col='green',linetype='dashed')+
    #     geom_hline(yintercept = Irecent,col='blue',linetype='dotted')+
    #     scale_color_manual(values=c('red','green','blue'),breaks=c("I0", "Itarget",'Iy'))+
    #     scale_y_continuous(limits=c(0,max(MPdf$I)*1.05),expand=c(0,0))+
    #     geom_text(label='A)',y=Inf,x=-Inf,hjust=-0.3,vjust=1.3)
    #
    # MPdf2 <- data.frame(x=seq(0,Itarget*1.2,length.out = 500),y=sapply(seq(0,Itarget*1.2,length.out = 500), ItoTac))
    #
    # p2 <- ggplot(MPdf2,aes(x=x,y=y))+geom_line(size=2)+xlab(expression(SSB["egg"] (t)))+ylab('TAC')+
    #     geom_hline(yintercept = TACtarget)+geom_text(x=(Itarget-I0)/2,y=TACtarget,label='TAC target',vjust=-0.4,hjust=0.2)+
    #     geom_hline(yintercept = w*TACtarget)+geom_text(x=(Itarget-I0)/2,y=w*TACtarget,label='TAC base',vjust=-0.4,hjust=0.2)+
    #     geom_vline(xintercept=Itarget,col='green',linetype='dashed',size=1.5)+geom_text(x=Itarget,y=0,label='Itarget',vjust=-0.4,hjust=-0.2,col='green')+
    #     geom_vline(xintercept=I0,col='red',linetype='dashed',size=1.5)+  geom_text(x=I0,y=0,label='I0',vjust=-0.4,hjust=-0.4,col='red')+
    #     geom_vline(xintercept=Irecent,col='blue',linetype='dotted',size=1.5)+geom_text(x=Irecent,y=0,label='Iy',vjust=-0.4,hjust=-0.2,col='blue')+
    #     scale_y_continuous(expand=c(0,0),limits=c(0,30000))+scale_x_continuous(expand=c(0,0))+
    #     geom_text(label='B)',y=Inf,x=-Inf,hjust=-0.3,vjust=1.3)
    #
    # library(gridExtra)
    # savepng(grid.arrange(p1,p2,ncol=2),wdIMG,"/HCR/HCR_target_10000",c(30,7))

}
class(MPeggcomplex10000) <- append(class(MPeggcomplex10000),"MP")
attr(MPeggcomplex10000,'model')=FALSE

##' MPeggcomplex2000
##' @param data list with data objects
##' @param TAC.base TAC of previous year
##' @param i year in the future
##' @rdname MPeggcomplex2000
##' @details Based on the Gerromont and Butterworth rule.
##' @export
MPeggcomplex2000<-function(data,TAC.base){
    Ihist <- cbind(I=exp(data[[1]]$logobs[,1]),data[[1]]$aux)
    Ihist <- Ihist[which(Ihist[,3]==3),c(1,2)]

    Ytarget <- c(1979,1983,1984,1987,1988,1989,1990,1992,1993,1994,1995)
    Y0 <- c(1996:2000,2004:2008)
    Ycurrent <- unname(tail(Ihist[,2],3))

    Itarget <- gmean(Ihist[which(Ihist[,2] %in% Ytarget),1])
    I0 <- gmean(Ihist[which(Ihist[,2] %in% Y0),1])
    TACtarget <- 25000
    TACbase <- 2000
    TACmin <- 2000
    TACfloor <- 2000
    TACmax <- 25000
    w = TACbase/TACtarget
    Irecent <- exp(unlist(lapply(data,function(x){gmean(tail(x$logobs[which(x$aux[,2]==3),1],3))})))

    ItoTac <- function(x){
        if(x<I0){
            TAC <- w*TACtarget*(x/I0)^2
        }else{
            TAC <- TACtarget*(w+(1-w)*((x-I0)/(Itarget-I0)))
        }
        TAC <- ifelse(TAC<TACmin,TACfloor,TAC)
        TAC <- ifelse(TAC>TACmax,TACmax,TAC)
        return(TAC)
    }

    TAC <- sapply(Irecent, ItoTac)
    if(length(TAC)==1) TAC <- rep(TAC,length(TAC.base))
    return(TAC)

    # MPdf=data.frame(Ihist)
    # MPdf[MPdf$year %in% Ytarget,'class'] <- 'Itarget'
    # MPdf[MPdf$year %in% Y0,'class'] <- 'I0'
    # MPdf[MPdf$year %in% Ycurrent,'class'] <- 'Iy'
    #
    # p1 <- ggplot(MPdf,aes(x=year,y=I))+geom_line()+geom_point(aes(col=class),size=2)+labs(col='',y=expression(SSB["egg"] (t)),x='Year')+
    #     theme(legend.position = c(0.8,0.8),
    #           legend.background = element_rect(fill=alpha('white', 0)))+
    #     geom_hline(yintercept = I0,col='red',linetype='dashed')+
    #     geom_hline(yintercept = Itarget,col='green',linetype='dashed')+
    #     geom_hline(yintercept = Irecent,col='blue',linetype='dotted')+
    #     scale_color_manual(values=c('red','green','blue'),breaks=c("I0", "Itarget",'Iy'))+
    #     scale_y_continuous(limits=c(0,max(MPdf$I)*1.05),expand=c(0,0))+
    #     geom_text(label='A)',y=Inf,x=-Inf,hjust=-0.3,vjust=1.3)
    #
    # MPdf2 <- data.frame(x=seq(0,Itarget*1.2,length.out = 500),y=sapply(seq(0,Itarget*1.2,length.out = 500), ItoTac))
    #
    # p2 <- ggplot(MPdf2,aes(x=x,y=y))+geom_line(size=2)+xlab(expression(SSB["egg"] (t)))+ylab('TAC')+
    #     geom_hline(yintercept = TACtarget)+geom_text(x=(Itarget-I0)/2,y=TACtarget,label='TAC target',vjust=-0.4,hjust=0.2)+
    #     geom_hline(yintercept = w*TACtarget)+geom_text(x=(Itarget-I0)/2,y=w*TACtarget,label='TAC base',vjust=-0.4,hjust=0.2)+
    #     geom_vline(xintercept=Itarget,col='green',linetype='dashed',size=1.5)+geom_text(x=Itarget,y=0,label='Itarget',vjust=-0.4,hjust=-0.2,col='green')+
    #     geom_vline(xintercept=I0,col='red',linetype='dashed',size=1.5)+  geom_text(x=I0,y=0,label='I0',vjust=-0.4,hjust=-0.4,col='red')+
    #     geom_vline(xintercept=Irecent,col='blue',linetype='dotted',size=1.5)+geom_text(x=Irecent,y=0,label='Iy',vjust=-0.4,hjust=-0.2,col='blue')+
    #     scale_y_continuous(expand=c(0,0),limits=c(0,30000))+scale_x_continuous(expand=c(0,0))+
    #     geom_text(label='B)',y=Inf,x=-Inf,hjust=-0.3,vjust=1.3)
    #
    # library(gridExtra)
    # savepng(grid.arrange(p1,p2,ncol=2),wdIMG,"/HCR/HCR_target_2000",c(30,7))

}
class(MPeggcomplex2000) <- append(class(MPeggcomplex2000),"MP")
attr(MPeggcomplex2000,'model')=FALSE

##' MPeggcomplex4000
##' @param data list with data objects
##' @param TAC.base TAC of previous year
##' @param i year in the future
##' @rdname MPeggcomplex4000
##' @details Based on the Gerromont and Butterworth rule.
##' @export
MPeggcomplex4000<-function(data,TAC.base){
    Ihist <- cbind(I=exp(data[[1]]$logobs[,1]),data[[1]]$aux)
    Ihist <- Ihist[which(Ihist[,3]==3),c(1,2)]

    Ytarget <- c(1979,1983,1984,1987,1988,1989,1990,1992,1993,1994,1995)
    Y0 <- c(1996:2000,2004:2008)
    Ycurrent <- unname(tail(Ihist[,2],3))

    Itarget <- gmean(Ihist[which(Ihist[,2] %in% Ytarget),1])
    I0 <- gmean(Ihist[which(Ihist[,2] %in% Y0),1])
    TACtarget <- 25000
    TACbase <- 4000
    TACmin <- 4000
    TACfloor <- 4000
    TACmax <- 25000
    w = TACbase/TACtarget
    Irecent <- exp(unlist(lapply(data,function(x){gmean(tail(x$logobs[which(x$aux[,2]==3),1],3))})))

    ItoTac <- function(x){
        if(x<I0){
            TAC <- w*TACtarget*(x/I0)^2
        }else{
            TAC <- TACtarget*(w+(1-w)*((x-I0)/(Itarget-I0)))
        }
        TAC <- ifelse(TAC<TACmin,TACfloor,TAC)
        TAC <- ifelse(TAC>TACmax,TACmax,TAC)
        return(TAC)
    }

    TAC <- sapply(Irecent, ItoTac)
    if(length(TAC)==1) TAC <- rep(TAC,length(TAC.base))
    return(TAC)

    # MPdf=data.frame(Ihist)
    # MPdf[MPdf$year %in% Ytarget,'class'] <- 'Itarget'
    # MPdf[MPdf$year %in% Y0,'class'] <- 'I0'
    # MPdf[MPdf$year %in% Ycurrent,'class'] <- 'Iy'
    #
    # p1 <- ggplot(MPdf,aes(x=year,y=I))+geom_line()+geom_point(aes(col=class),size=2)+labs(col='',y=expression(SSB["egg"] (t)),x='Year')+
    #     theme(legend.position = c(0.8,0.8),
    #           legend.background = element_rect(fill=alpha('white', 0)))+
    #     geom_hline(yintercept = I0,col='red',linetype='dashed')+
    #     geom_hline(yintercept = Itarget,col='green',linetype='dashed')+
    #     geom_hline(yintercept = Irecent,col='blue',linetype='dotted')+
    #     scale_color_manual(values=c('red','green','blue'),breaks=c("I0", "Itarget",'Iy'))+
    #     scale_y_continuous(limits=c(0,max(MPdf$I)*1.05),expand=c(0,0))+
    #     geom_text(label='A)',y=Inf,x=-Inf,hjust=-0.3,vjust=1.3)
    #
    # MPdf2 <- data.frame(x=seq(0,Itarget*1.2,length.out = 500),y=sapply(seq(0,Itarget*1.2,length.out = 500), ItoTac))
    #
    # p2 <- ggplot(MPdf2,aes(x=x,y=y))+geom_line(size=2)+xlab(expression(SSB["egg"] (t)))+ylab('TAC')+
    #     geom_hline(yintercept = TACtarget)+geom_text(x=(Itarget-I0)/2,y=TACtarget,label='TAC target',vjust=-0.4,hjust=0.2)+
    #     geom_hline(yintercept = w*TACtarget)+geom_text(x=(Itarget-I0)/2,y=w*TACtarget,label='TAC base',vjust=-0.4,hjust=0.2)+
    #     geom_vline(xintercept=Itarget,col='green',linetype='dashed',size=1.5)+geom_text(x=Itarget,y=0,label='Itarget',vjust=-0.4,hjust=-0.2,col='green')+
    #     geom_vline(xintercept=I0,col='red',linetype='dashed',size=1.5)+  geom_text(x=I0,y=0,label='I0',vjust=-0.4,hjust=-0.4,col='red')+
    #     geom_vline(xintercept=Irecent,col='blue',linetype='dotted',size=1.5)+geom_text(x=Irecent,y=0,label='Iy',vjust=-0.4,hjust=-0.2,col='blue')+
    #     scale_y_continuous(expand=c(0,0),limits=c(0,30000))+scale_x_continuous(expand=c(0,0))+
    #     geom_text(label='B)',y=Inf,x=-Inf,hjust=-0.3,vjust=1.3)
    #
    # library(gridExtra)
    # savepng(grid.arrange(p1,p2,ncol=2),wdIMG,"/HCR/HCR_target_4000",c(30,7))

}
class(MPeggcomplex4000) <- append(class(MPeggcomplex4000),"MP")
attr(MPeggcomplex4000,'model')=FALSE

##' MPeggcomplex6000
##' @param data list with data objects
##' @param TAC.base TAC of previous year
##' @param i year in the future
##' @rdname MPeggcomplex6000
##' @details Based on the Gerromont and Butterworth rule.
##' @export
MPeggcomplex6000<-function(data,TAC.base){
    Ihist <- cbind(I=exp(data[[1]]$logobs[,1]),data[[1]]$aux)
    Ihist <- Ihist[which(Ihist[,3]==3),c(1,2)]

    Ytarget <- c(1979,1983,1984,1987,1988,1989,1990,1992,1993,1994,1995)
    Y0 <- c(1996:2000,2004:2008)
    Ycurrent <- unname(tail(Ihist[,2],3))

    Itarget <- gmean(Ihist[which(Ihist[,2] %in% Ytarget),1])
    I0 <- gmean(Ihist[which(Ihist[,2] %in% Y0),1])
    TACtarget <- 25000
    TACbase <- 6000
    TACmin <- 6000
    TACfloor <- 6000
    TACmax <- 25000
    w = TACbase/TACtarget
    Irecent <- exp(unlist(lapply(data,function(x){gmean(tail(x$logobs[which(x$aux[,2]==3),1],3))})))

    ItoTac <- function(x){
        if(x<I0){
            TAC <- w*TACtarget*(x/I0)^2
        }else{
            TAC <- TACtarget*(w+(1-w)*((x-I0)/(Itarget-I0)))
        }
        TAC <- ifelse(TAC<TACmin,TACfloor,TAC)
        TAC <- ifelse(TAC>TACmax,TACmax,TAC)
        return(TAC)
    }

    TAC <- sapply(Irecent, ItoTac)
    if(length(TAC)==1) TAC <- rep(TAC,length(TAC.base))
    return(TAC)

    # MPdf=data.frame(Ihist)
    # MPdf[MPdf$year %in% Ytarget,'class'] <- 'Itarget'
    # MPdf[MPdf$year %in% Y0,'class'] <- 'I0'
    # MPdf[MPdf$year %in% Ycurrent,'class'] <- 'Iy'
    #
    # p1 <- ggplot(MPdf,aes(x=year,y=I))+geom_line()+geom_point(aes(col=class),size=2)+labs(col='',y=expression(SSB["egg"] (t)),x='Year')+
    #     theme(legend.position = c(0.8,0.8),
    #           legend.background = element_rect(fill=alpha('white', 0)))+
    #     geom_hline(yintercept = I0,col='red',linetype='dashed')+
    #     geom_hline(yintercept = Itarget,col='green',linetype='dashed')+
    #     geom_hline(yintercept = Irecent,col='blue',linetype='dotted')+
    #     scale_color_manual(values=c('red','green','blue'),breaks=c("I0", "Itarget",'Iy'))+
    #     scale_y_continuous(limits=c(0,max(MPdf$I)*1.05),expand=c(0,0))+
    #     geom_text(label='A)',y=Inf,x=-Inf,hjust=-0.3,vjust=1.3)
    #
    # MPdf2 <- data.frame(x=seq(0,Itarget*1.2,length.out = 500),y=sapply(seq(0,Itarget*1.2,length.out = 500), ItoTac))
    #
    # p2 <- ggplot(MPdf2,aes(x=x,y=y))+geom_line(size=2)+xlab(expression(SSB["egg"] (t)))+ylab('TAC')+
    #     geom_hline(yintercept = TACtarget)+geom_text(x=(Itarget-I0)/2,y=TACtarget,label='TAC target',vjust=-0.4,hjust=0.2)+
    #     geom_hline(yintercept = w*TACtarget)+geom_text(x=(Itarget-I0)/2,y=w*TACtarget,label='TAC base',vjust=-0.4,hjust=0.2)+
    #     geom_vline(xintercept=Itarget,col='green',linetype='dashed',size=1.5)+geom_text(x=Itarget,y=0,label='Itarget',vjust=-0.4,hjust=-0.2,col='green')+
    #     geom_vline(xintercept=I0,col='red',linetype='dashed',size=1.5)+  geom_text(x=I0,y=0,label='I0',vjust=-0.4,hjust=-0.4,col='red')+
    #     geom_vline(xintercept=Irecent,col='blue',linetype='dotted',size=1.5)+geom_text(x=Irecent,y=0,label='Iy',vjust=-0.4,hjust=-0.2,col='blue')+
    #     scale_y_continuous(expand=c(0,0),limits=c(0,30000))+scale_x_continuous(expand=c(0,0))+
    #     geom_text(label='B)',y=Inf,x=-Inf,hjust=-0.3,vjust=1.3)
    #
    # library(gridExtra)
    # savepng(grid.arrange(p1,p2,ncol=2),wdIMG,"/HCR/HCR_target_6000",c(30,7))

}
class(MPeggcomplex6000) <- append(class(MPeggcomplex6000),"MP")
attr(MPeggcomplex6000,'model')=FALSE


##' MPeggcomplex15000
##' @param data list with data objects
##' @param TAC.base TAC of previous year
##' @param i year in the future
##' @rdname MPeggcomplex15000
##' @details Based on the Gerromont and Butterworth rule.
##' @export
MPeggcomplex15000<-function(data,TAC.base){
    Ihist <- cbind(I=exp(data[[1]]$logobs[,1]),data[[1]]$aux)
    Ihist <- Ihist[which(Ihist[,3]==3),c(1,2)]

    Ytarget <- c(1979,1983,1984,1987,1988,1989,1990,1992,1993,1994,1995)
    Y0 <- c(1996:2000,2004:2008)
    Ycurrent <- unname(tail(Ihist[,2],3))

    Itarget <- gmean(Ihist[which(Ihist[,2] %in% Ytarget),1])
    I0 <- gmean(Ihist[which(Ihist[,2] %in% Y0),1])
    TACtarget <- 25000
    TACbase <- 15000
    TACmin <- 15000
    TACfloor <- 15000
    TACmax <- 25000
    w = TACbase/TACtarget
    Irecent <- exp(unlist(lapply(data,function(x){gmean(tail(x$logobs[which(x$aux[,2]==3),1],3))})))

    ItoTac <- function(x){
        if(x<I0){
            TAC <- w*TACtarget*(x/I0)^2
        }else{
            TAC <- TACtarget*(w+(1-w)*((x-I0)/(Itarget-I0)))
        }
        TAC <- ifelse(TAC<TACmin,TACfloor,TAC)
        TAC <- ifelse(TAC>TACmax,TACmax,TAC)
        return(TAC)
    }

    TAC <- sapply(Irecent, ItoTac)
    if(length(TAC)==1) TAC <- rep(TAC,length(TAC.base))
    return(TAC)

    # MPdf=data.frame(Ihist)
    # MPdf[MPdf$year %in% Ytarget,'class'] <- 'Itarget'
    # MPdf[MPdf$year %in% Y0,'class'] <- 'I0'
    # MPdf[MPdf$year %in% Ycurrent,'class'] <- 'Iy'
    #
    # p1 <- ggplot(MPdf,aes(x=year,y=I))+geom_line()+geom_point(aes(col=class),size=2)+labs(col='',y=expression(SSB["egg"] (t)),x='Year')+
    #     theme(legend.position = c(0.8,0.8),
    #           legend.background = element_rect(fill=alpha('white', 0)))+
    #     geom_hline(yintercept = I0,col='red',linetype='dashed')+
    #     geom_hline(yintercept = Itarget,col='green',linetype='dashed')+
    #     geom_hline(yintercept = Irecent,col='blue',linetype='dotted')+
    #     scale_color_manual(values=c('red','green','blue'),breaks=c("I0", "Itarget",'Iy'))+
    #     scale_y_continuous(limits=c(0,max(MPdf$I)*1.05),expand=c(0,0))+
    #     geom_text(label='A)',y=Inf,x=-Inf,hjust=-0.3,vjust=1.3)
    #
    # MPdf2 <- data.frame(x=seq(0,Itarget*1.2,length.out = 500),y=sapply(seq(0,Itarget*1.2,length.out = 500), ItoTac))
    #
    # p2 <- ggplot(MPdf2,aes(x=x,y=y))+geom_line(size=2)+xlab(expression(SSB["egg"] (t)))+ylab('TAC')+
    #     geom_hline(yintercept = TACtarget)+geom_text(x=(Itarget-I0)/2,y=TACtarget,label='TAC target',vjust=-0.4,hjust=0.2)+
    #     geom_hline(yintercept = w*TACtarget)+geom_text(x=(Itarget-I0)/2,y=w*TACtarget,label='TAC base',vjust=-0.4,hjust=0.2)+
    #     geom_vline(xintercept=Itarget,col='green',linetype='dashed',size=1.5)+geom_text(x=Itarget,y=0,label='Itarget',vjust=-0.4,hjust=-0.2,col='green')+
    #     geom_vline(xintercept=I0,col='red',linetype='dashed',size=1.5)+  geom_text(x=I0,y=0,label='I0',vjust=-0.4,hjust=-0.4,col='red')+
    #     geom_vline(xintercept=Irecent,col='blue',linetype='dotted',size=1.5)+geom_text(x=Irecent,y=0,label='Iy',vjust=-0.4,hjust=-0.2,col='blue')+
    #     scale_y_continuous(expand=c(0,0),limits=c(0,30000))+scale_x_continuous(expand=c(0,0))+
    #     geom_text(label='B)',y=Inf,x=-Inf,hjust=-0.3,vjust=1.3)
    #
    # library(gridExtra)
    # savepng(grid.arrange(p1,p2,ncol=2),wdIMG,"/HCR/HCR_target_15000",c(30,7))
    #
}
class(MPeggcomplex15000) <- append(class(MPeggcomplex15000),"MP")
attr(MPeggcomplex15000,'model')=FALSE

##' MPccam
##' @param data lists with data lists
##' @param parameters one list of parameters
##' @param conf one list with configurations
##' @param TAC.base TAC of previous year
##' @param i year in the future
##' @rdname MPccam
##' @details Determine TAC by fiting censored model, determining F40% and health zones and applying F based on the latter. Attributed alternative data  can be used in case simulated data does not lead to convergence.
##' @export
##' @import parallel
MPccam <- function(data, parameters, conf, TAC.base, i){
    if(!i %% 2 == 1){return(TAC.base)} # first future year there is an evaluation, than there's not
    #ncores <- detectCores()
    #cl <- makeCluster(ncores-1) #set up nodes, outfile="debug_MPccam.txt"
    #clusterExport(cl, varlist=c("conf","parameters"), envir=environment())
    #clusterCall(cl, function() library(CCAM))
    #TAC <- parLapply(cl,data,function(x){
    counter <-0 # both for debugging and later on if superassigment

    penv <- dynGet('myenv',inherits = TRUE)
    nfit <- dynGet('nfit',inherits = TRUE)
    nfail <- dynGet('nfail',inherits = TRUE)

    #myds <<- data #just for debugging
    TAC <- lapply(data,function(x){
        counter <<- counter+1
        nfit <<- nfit+1

        #make sure fleet 2 (caa) goes not extremely low because of observation error (mackerel specific...)
        for(i in x$noYears:ncol(x$idx1)){
            n <- x$logobs[(x$idx1[2,i]+1):(x$idx2[2,i]+1),1]
            n <- crlInverse(n)
            n[n<0.0001] <- 0.0001
            x$logobs[(x$idx1[2,i]+1):(x$idx2[2,i]+1),1] <- crlTransform(n)
        }

        fit <- force.fit(x,conf, parameters,20,silent = TRUE)
        if(class(fit)!='ccam'){  #if there is an error or no model convergence
                nfail <<- nfail+1
                return(TAC.base[counter])
        }

        ssbp <- tail(ssbtable(fit),1)[,1]
        ref <- ypr(fit,what=c('f40ssb','f40'))
        ssbref <- ref$f40ssb
        fref <- ref$f40
        Bupper <- 0.8*ssbref
        Blim <- 0.4*ssbref

            if(ssbp<Blim){
                TAC <- 1
            }else{
                if(ssbp<Bupper){
                    fref <- fref*((ssbp-Blim)/(Bupper-Blim)) #scale if in cautious zone
                }
                fut <- tryCatch(forecast(fit,fval=fref,nosim = 300,
                                         ave.years = max(fit$data$years)+(-9:0),
                                         rec.years = 1969:max(fit$data$years),
                                         rec.meth = 2,
                                         verbose=FALSE,
                                         Flim=2.5,
                                         deadzone=1000,
                                         determini=TRUE), error = identity)
                if (inherits(fut, "error")) {nfail <<- nfail+1;return(TAC.base[counter])}
                Clim <- tail(attr(fut,"tab")[,'catch:median'],1)
                declared <- exp(fit$data$logobs[which(!is.na(fit$data$logobs[,2])),1])
                predicted <- exp(fit$rep$predObs[which(!is.na(fit$data$logobs[,2]))])
                U <- predicted-declared
                U <- max(c(0,mean(tail(U,3)))) # undeclared catch 3 last years estimated by model, never lower than 0
                TAC <-  max(Clim-U,0) #still gives catches that include US catch...
                TAC <- TAC*0.8 # presuming 20% is for the US #plot(ct[,1]/(0.25*ctUSA[-c(1:8),1]+ct[,1]))
            }
        return(TAC)
    })
    TAC <- unlist(TAC)
    if(length(TAC)==1) TAC <- rep(TAC, length(TAC.base))
    #stopCluster(cl)
    assign('nfit',nfit,envir = penv)
    assign('nfail',nfail,envir = penv)
    return(TAC)
}
class(MPccam) <- append(class(MPccam),"MP")
attr(MPccam,'model')=TRUE

##' MPspm
##' @param data list with data objects
##' @param TAC.base TAC of previous state
##' @param i year
##' @rdname MPspm
##' @details determine TAC by fiting surplus production model from SPiCT package.
##' @export
##' @importFrom  spict fit.spict get.par manage prop.F check.inp make.datin
MPspm <- function(data,TAC.base,i){
    if(!i %% 2 == 1){return(TAC.base)}
    #ncores <- detectCores()
    #cl <- makeCluster(ncores-1) #set up nodes
    #clusterCall(cl, function() library(spict))
    #obsfit <- parLapply(cl,data,function(x){
    fit.spict.rob <- function (inp, datin, pl, dbg = 0, sd=0.25){
        rep <- NULL
        for (i in 1:inp$nphases) {
            if (inp$nphases > 1) cat(paste("Estimating - phase", i, "\n"))
            for(b in 1:60){
                mysd <- ifelse(b==1,0,sd)
                parv <- unlist(pl)
                pl2 <- relist(parv+rnorm(length(parv),sd=mysd), pl)
                obj <- make.obj(datin, pl2, inp, phase = i)
                opt <- try(nlminb(obj$par, obj$fn, obj$gr, control = inp$optimiser.control))

                if (class(opt) != "try-error") {
                    pl <- obj$env$parList(opt$par)
                }
                if(opt$convergence==0){
                    msy <- try(obj$report()$Bmsy)
                    if(msy>0) break
                }
              }
        }
        if (dbg < 1) {
            optfailflag <- class(opt) == "try-error"
            sdfailflag <- FALSE
            if (optfailflag) {
                cat("obj$par:\n")
                print(obj$par)
                cat("obj$fn:\n")
                print(obj$fn())
                cat("obj$gr:\n")
                print(obj$gr())
                stop("Could not fit model. Error msg:", opt)
            }
            else {
                if (inp$do.sd.report) {
                    verflag <- as.numeric(gsub("[.]", "", as.character(packageVersion("TMB")))) >=
                        171
                    if (verflag) {
                        rep <- try(TMB::sdreport(obj, getJointPrecision = inp$getJointPrecision,
                                                 bias.correct = inp$bias.correct, bias.correct.control = inp$bias.correct.control,
                                                 getReportCovariance = inp$getReportCovariance))
                    }
                    else {
                        rep <- try(TMB::sdreport(obj, getJointPrecision = inp$getJointPrecision,
                                                 bias.correct = inp$bias.correct, bias.correct.control = inp$bias.correct.control))
                    }
                    sdfailflag <- class(rep) == "try-error"
                    if (sdfailflag) {
                        warning("Could not calculate sdreport.\n")
                        rep <- NULL
                    }
                }
                if (is.null(rep)) {
                    rep <- list()
                    if (sdfailflag) {
                        rep <- list()
                        rep$sderr <- 1
                        rep$par.fixed <- opt$par
                        rep$cov.fixed <- matrix(NA, length(opt$par),
                                                length(opt$par))
                    }
                }
                rep$inp <- inp
                rep$obj <- obj
                rep$opt <- opt
                rep$opt$gr <- rep$obj$gr(rep$opt$par)
                rep$pl <- obj$env$parList(opt$par)
                obj$fn()
                rep$Cp <- obj$report()$Cp
                rep$report <- obj$report()
                if (!sdfailflag & inp$reportall) {
                    if (!inp$osar.method == "none") {
                        reposar <- try(calc.osa.resid(rep))
                        if (class(reposar) != "try-error") {
                            rep <- reposar
                        }
                    }
                }
            }
        }
        if (!is.null(rep)) {
            class(rep) <- "spictcls"
        }
        return(rep)
    }

    myds <<- data
    counter <- 0

    penv <- dynGet('myenv',inherits = TRUE)
    nfit <- dynGet('nfit',inherits = TRUE)
    nfail <- dynGet('nfail',inherits = TRUE)

    obsfit <- lapply(data,function(x){
        counter <<- counter+1
        nfit <<- nfit+1
            o <- x$logobs
            aux <- x$aux
            d <- list(obsC = unname(exp(o[which(aux[,2]==1),1])), #simulated catch is simulated Can+ 0.25% US catch
                      timeC = unname(aux[which(aux[,2]==1),1]),
                      obsI = unname(exp(o[which(aux[,2]==3),1])),
                      timeI = unname(aux[which(aux[,2]==3),1]))
            d$obsC[d$obsC<1] <- 1
            d$obsI[d$obsI<1] <- 1

            inp <- check.inp(d)
            datin <- make.datin(inp, 0)
            pl <- inp$parlist
            if(i!=1) pl[c(1:19,26:27)] <- mypl[c(1:19,26:27)] #improve initial values
            res <- tryCatch(fit.spict.rob(inp,datin,pl),error=identity)
            myres <<- res
            if('error' %in% class(res)){
                nfail <<- nfail+1
               return(TAC.base[counter])
            }
            msy <- get.par('Bmsy', res)[2]
            if(is.na(msy))  msy <- -1
            if(0 !=res$opt$convergence | msy<0){
                nfail <<- nfail+1
              return(TAC.base[counter])
            }
        Bmsy <- get.par('Bmsy', res)[2]
        Blast <- get.par('logBl', res, exp=TRUE)[2]
        Blim <- Bmsy*0.4
        Bupper <- Bmsy*0.8
        if(Blast<Blim){TAC=1}else{
            # maximal fishing
            if(Blast>Bupper){
                res <- tryCatch(manage(res, scenarios = 3),error=function(e)print(e)) #fish at msy
                if('error' %in% class(res)){
                    nfail <<- nfail+1
                    return(TAC.base[counter])
                }
            }else{
                # in between fishing
                Fmsy <- get.par("logFmsyd", res, exp = TRUE)[2]
                Fref <- (Fmsy*(Blast-Blim))/(Bupper-Blim) # catious zone
                Flast <- get.par("logF", res, exp = TRUE)[res$inp$indpred[1], 2]
                Fchange <- Fref/Flast
                res <- tryCatch(prop.F(Fchange, res$inp, res, which(res$inp$time >= res$inp$manstart)),error=function(e)print(e))
                if('error' %in% class(res)){
                    nfail <<- nfail+1
                    return(TAC.base[counter])
                }
            }
            TAC <- get.par("Cp", res)[2]
            TAC <- TAC*0.8 #Because the observed catches are Canadian declared catches + 0.25 US landings
        }
        if(i==1) mypl <<- res$pl
        return(TAC)
       })
    #stopCluster(cl)
    TAC <- do.call('c',obsfit)
    if(length(TAC)==1) TAC <- rep(TAC, length(TAC.base))
    assign('nfit',nfit,envir = penv)
    assign('nfail',nfail,envir = penv)
    return(TAC)
}
class(MPspm) <- append(class(MPspm),"MP")
attr(MPspm,'model')=TRUE
