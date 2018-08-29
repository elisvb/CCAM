library(CCAM)

ct <- read.ices("ct.dat")
cn<-read.ices("cn.dat")
cw<-read.ices("cw.dat")
dw<-read.ices("dw.dat")
lf<-read.ices("lf.dat")
lw<-read.ices("lw.dat")
mo<-read.ices("mo.dat")
nm<-read.ices("nm.dat")
pf<-read.ices("pf.dat")
pm<-read.ices("pm.dat")
sw<-read.ices("sw.dat")
surveys<-read.ices("survey.dat")

dat<-setup.ccam.data(surveys=surveys,
                    residual.fleet=cn,
                    #total.catch = ct,
                    prop.mature=mo,
                    stock.mean.weight=sw,
                    catch.mean.weight=cw,
                    dis.mean.weight=dw,
                    land.mean.weight=lw,
                    prop.f=pf,
                    prop.m=pm,
                    natural.mortality=nm,
                    land.frac=lf)

conf<-defcon(dat)
conf$fbarRange<-c(2,6)
conf$keySel<- matrix(c(0,1,2,3,4,4,4,4,4),ncol=ncol(conf$keySel),nrow=nrow(conf$keySel),byrow = TRUE)

par<-defpar(dat,conf)
fit<-ccam.fit(dat,conf,par,silent=TRUE)
fit

