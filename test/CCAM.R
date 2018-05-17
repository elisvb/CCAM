setwd("C:/Users/VANBE/Desktop/post-doc/DATA/CCAM/")
setwd("C:/Users/VANBE/Desktop/post-doc/DATA/CCAM/tests/canmack")

library(CCAM)

cn <- read.ices("cn.dat")
ct <- read.ices("ct.dat")
cw <- read.ices("cw.dat")
dw <- read.ices("dw.dat")
lf <- read.ices("lf.dat")
lw <- read.ices("lw.dat")
mo <- read.ices("mo.dat")
nm <- read.ices("nm.dat")
pf <- read.ices("pf.dat")
pm <- read.ices("pm.dat")
sw <- read.ices("sw.dat")
surveys <- read.ices("survey.dat")
env <- read.ices("env.dat")
env[,1]=0

dat <- setup.ccam.data(surveys=surveys,
                      residual.fleet=cn, # add argument called split.catch and see that residual.fleet does nothing if null
                      total.catch=ct,
                      prop.mature=mo,
                      stock.mean.weight=sw,
                      catch.mean.weight=cw,
                      dis.mean.weight=dw,
                      land.mean.weight=lw,
                      prop.f=pf,
                      prop.m=pm,
                      natural.mortality=nm,
                      land.frac=lf,
                      env=env)
dat$fleetTypes

rm=is.na(dat$logobs[,1])
dat$logobs=dat$logobs[!rm,]
dat$aux=dat$aux[!rm,]
#head(cbind(dat$aux,dat$logobs),20)
conf <- defcon(dat) # most keys should be the same for 1 and 2
conf$keyLogFsta[1,]=c(0,1,2,3,3,3,3,3,3,3)
conf$keyVarObs[1,]=-1                       # sd around total catch always estimated as nearly zero, and problem calculating sd
conf$keyVarObs[2,1:9]=c(0,1,2,2,2,2,2,2,1)
conf$keyVarObs[3,1]=3
conf$stockRecruitmentModelCode=0

#conf$stockRecruitmentModelCode=2

par <- defpar(dat,conf)


data=dat;parameters=par;sim.condRE=TRUE
definit <- defpar(data, conf)
if(!identical(parameters,relist(unlist(parameters), skeleton=definit))){
    warning("Initial values are not consistent, so running with default init values from defpar()")
    parameters<-definit
}

conf$obsLikelihoodFlag[1]='CE'
data$logobs[!is.na(data$logobs[,2]),1] =log(apply(exp(data$logobs[!is.na(data$logobs[,2]),]),1,mean))
rem=unique(conf$keyVarObs[1,])+1
if(rem!=0) parameters$logSdLogObs=parameters$logSdLogObs[-rem]
conf$obsLikelihoodFlag[1]='LN'
if(!all(conf$keyVarObs[1,]<0)) conf$keyVarObs = conf$keyVarObs-1
conf$keyVarObs[conf$keyVarObs<0]=-1

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
data<-clean.void.catches(data,conf)
tmball <- c(data, conf, simFlag=as.numeric(sim.condRE))
if(is.null(tmball$resFlag)){tmball$resFlag <- 0}
nmissing <- sum(is.na(data$logobs[,1]))
parameters$missing <- numeric(nmissing)
ran <- c("logN", "logFy","missing") ## change in function

TMB::compile(paste0("C:/Users/VANBE/Desktop/post-doc/DATA/CCAM/test/CCAM.cpp"),"-O1 -g",DLLFLAGS ="")
dyn.load(TMB::dynlib("C:/Users/VANBE/Desktop/post-doc/DATA/CCAM/test/CCAM"))   #dyn.unload(dynlib(file))
#dyn.load(TMB::dynlib("C:/Users/VANBE/Desktop/post-doc/DATA/CCAM/src-x64/CCAM"))
obj <- TMB::MakeADFun(tmball,parameters,random=ran,DLL="CCAM",silent=FALSE)

lower2<-rep(-Inf,length(obj$par))
upper2<-rep(Inf,length(obj$par))

opt <- nlminb(obj$par, obj$fn,obj$gr ,control=list(trace=1, eval.max=2000, iter.max=1000),lower=lower2,upper=upper2)




rep = obj$rep()
opt$par
for(i in seq_len(3)) { # Take a few extra newton steps
    g <- as.numeric( obj$gr(opt$par) )
    h <- optimHess(opt$par, obj$fn, obj$gr)
    opt$par <- opt$par - solve(h, g)
    opt$objective <- obj$fn(opt$par)
}
rep <- obj$report()
sdrep <- sdreport(obj,opt$par)
sdrep

idx <- c(which(names(sdrep$value)=="lastLogN"),which(names(sdrep$value)=="lastLogF"))
sdrep$estY <- sdrep$value[idx]
sdrep$covY <- sdrep$cov[idx,idx]


### test package
 fit1 <- ccam.fit(dat,conf,par,silent=TRUE)              # run uncensored with lognormal error dist Ctot (sd = parameter)

 conf$obsLikelihoodFlag[1]='CE'

 fit2 <- ccam.fit(dat,conf,par,phase=1)     # run uncensored with lognormal error dist Ctot (sd = 0.01) and Clower becomes Cmean
 fit3 <- ccam.fit(dat,conf,par,phase=2)     # run phase 1 + censored

 conf$stockRecruitmentModelCode=2
 fit4 <- ccam.fit(dat,conf,par,phase=2)

 srplot(fit3)
 ypr(fit3)
 catchplot(fit3,ci=FALSE)


 fits=c(fit1,fit2,fit3)
 names(fits)=c('normal','uncens','cens')

 modeltable(fits)
 partable(fits)
 ssbplot(fits)
 ssbplot(fits,ci=FALSE,addCI=FALSE)
 fbarplot(fits,ci=FALSE)
 recplot(fits,ci=FALSE)
 catchplot(fits,ylim=c(0,70000))

 catchtable(fit1)

 # plot still to clean up
 fitplot(fit1,fleets = 2)
 fitplot(fit1,fleets = 2,log=F) # this is crl transformed and not in code. Anyway, ugly plots
 res <- residuals(fit1)
 plot(res)
 resp <- procres(fit)
 plot(resp)
 sims=simulate(fit)


 ### MSE

##different OMS
 fit=fit1; fscale=NULL; catchval=NULL; fval=NULL; nosim=2; year.base=max(fit$data$years);
 ave.years=max(fit$data$years)+(-9:0); rec.years=max(fit$data$years)+(-39:0); label=NULL;
 overwriteSelYears=NULL; deterministic=FALSE;threshold=103000;IE=FALSE;IEpar1=0;IEpar2=0;
 capType='abs';capLower=0;capUpper=NULL;MP = rep('MPf40base',3);UL.years=max(fit$data$years)+(-4:0)

 #uniroot error because eggsurvey estimates always go trough the roof, leading to insane catches
test1a=forecast(fit=fit1,MP = rep('MPeggsurvey',10),nosim = 200,ave.years = max(fit1$data$years)+(-9:0),rec.years = max(fit1$data$years)+(-39:0),
               threshold=103000,capType='abs',capLower = 0,deterministic = FALSE)
test2a=forecast(fit3,MP = rep('MPeggsurvey',10),nosim = 3,ave.years = max(fit1$data$years)+(-9:0),rec.years = max(fit1$data$years)+(-39:0),
               threshold=103000,capType='abs',capLower = 0)
### different MPs
test1b=forecast(fit=fit1,MP = rep('MPf40base',3),nosim = 2,ave.years = max(fit1$data$years)+(-9:0),rec.years = max(fit1$data$years)+(-39:0),
               threshold=103000,capType='abs',capLower = 0)

test1c=forecast(fit=fit1,MP = c(14000,16000,rep(1,10)),nosim = 200,ave.years = max(fit1$data$years)+(-9:0),rec.years = max(fit1$data$years)+(-39:0),
                threshold=103000,capType='abs',capLower = 0)
# why the hell are catches so high

### different capstuff
ssbplot(test1a)
ssbplot(c(test1a,test1c),col=c('black','red'))

catchplot(test1a)
catchplot(c(test1a,test1c),col=c('black','red'))

recplot(test1a)
recplot(c(test1a,test1c),col=c('black','red'))

fbarplot(test1c)
fbarplot(c(test1a,test1c),col=c('black','red'))

tacplot(test1c)
tacplot(c(test1a,test1c),col=c('black','red'))

probplot(test1a,what='ULR')  #always says no applicable method
probplot(test1a,what='Growth',add=T,col='red')
probplot(test1a,what='Double',add=T,col='green')
abline(h=75,lty=2)

probplot(c(test1a,test1c),what='ULR',leg=c('egg','moratorium'),zone = 75,main='ULR')
probplot(c(test1a,test1c),what='Growth',leg=c('egg','moratorium'),zone = 75,main='Growth')

tradeplot(test1a,what.x='ULR', what.y='catch',zone=c(75,8000),ci=TRUE,lab=TRUE)
tradeplot(test1a,what.x='ULR', what.y='catchcumul',zone=c(75),ci=TRUE)

tradeplot(c(test1a,test1c),what.x='ULR', what.y='ssb',zone=c(75),ci=FALSE,leg=c(1,2))
tradeplot(c(test1a,test1c),what.x='ULR', what.y='catchcumul',zone=c(75),ci=FALSE,leg=c(1,2))





