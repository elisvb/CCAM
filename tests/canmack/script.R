setwd("C:/Users/VANBE/Desktop/post-doc/DATA/CCAM/")
setwd("C:/Users/VANBE/Desktop/post-doc/DATA/CCAM/tests/canmack")
library(CCAM)
library(plyr)

#################################################################################################################
########### READ IN DATA ########################################################################################
#################################################################################################################

cn <- read.ices("cn.dat")
ct <- read.ices("ct.dat")
ctUSA <- read.ices("ctUSA.dat")
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

rm=is.na(dat$logobs[,1])
dat$logobs=dat$logobs[!rm,]
dat$aux=dat$aux[!rm,]
dat$weight=dat$weight[!rm]

conf <- defcon(dat)
conf$keyLogFsta[1,]=c(0,1,2,3,3,3,3,3,3,3)
conf$keyVarObs[1,]=-1                       # sd around total catch always estimated as nearly zero, and problem calculating sd
conf$keyVarObs[2,1:9]=c(0,1,2,2,2,2,2,2,1)
conf$keyVarObs[3,1]=3
conf$stockRecruitmentModelCode=2 #0: RW, 1: ricker, 2: BH
conf$fbarRange=c(4,10)

par <- defpar(dat,conf)

            # check input data
            # bubble(t(cn),scale=0.03) # raw data
            # bubble(spay(t(cn)))      # to check cohorts
            # bubble(spya(t(cn)))      # to compare between years
            # heat(mo,ncol=100)
            # heat(cw,col=c('black','darkgreen'))
            # matplot(mo,type='l',ylab='Proportion mature',xlab='Year')
            # matplot(cw,type='l',ylab='Catch weight',xlab='Year')
            # matplot(ct,type='l',ylab='Catch (t)', xlab='Year')
            # surveyplot(surveys)

#################################################################################################################
########### fit model ###########################################################################################
#################################################################################################################

fit1 <- ccam.fit(dat,conf,par,silent=TRUE)              # run uncensored with lognormal error dist Ctot (sd = parameter)

conf$obsLikelihoodFlag[1]='CE'

fit2 <- ccam.fit(dat,conf,par,phase=1)     # run uncensored with lognormal error dist Ctot (sd = 0.01) and Clower becomes Cmean
fit3 <- ccam.fit(dat,conf,par,phase=2)     # run phase 1 + censored

fits=c(fit1,fit2,fit3)
names(fits)=c('normal','uncens','cens')


#################################################################################################################
########### PLOTS & TABLES ######################################################################################
#################################################################################################################

srplot(fit3,curve=TRUE)
ypr(fit3,rec.years=1980:2016)
catchplot(fit3,ci=FALSE)

ssbplot(fits)
ssbplot(fits,ci=FALSE,addCI=FALSE)
fbarplot(fits,ci=FALSE)
recplot(fits,ci=FALSE)
catchplot(fits,ylim=c(0,70000))
plot(fits)

modeltable(fits)
partable(fits)
catchtable(fits)
ssbtable(fits)

# Residuals
par(mfrow=c(3,2))
resplot(fit1,fleets = 3,type=1)
resplot(fit1,fleets = 3,type=2,out=1)
resplot(fit1,fleets = 3,type=3)

resplot(fit1,fleets = 2,type=1)
resplot(fit1,fleets = 2,type=2,out=3)
resplot(fit1,fleets = 2,type=3,out=3)
resplot(fit1,fleets = 2,type=4)

fitplot(fit1,fleets = 3,log = FALSE,pch=16)
fitplot(fit1,fleets = 2,log = FALSE,pch=16)
fitplot(fit1,fleets = 1,log = FALSE,pch=16)

#unfinished
# res <- residuals(fit1)  #problems because logobs has become a matrix (instead of vector)!!
# plot(res)
# resp <- procres(fit) #needs to be adapted because of separable F
# plot(resp)
# sims=simulate(fit)  #simulate data (not yet ok for lower and upper catch limits)

#retrospective analysis
retro <- retro(fit1,year=7)  #maybe make plot with relative change
plot(retro)
mohn(retro)

#################################################################################################################
########### MSE #################################################################################################
#################################################################################################################

##### debugging
fit=fit1; fscale=NULL; catchval=NULL; fval=NULL; nosim=2; year.base=max(fit$data$years);
ave.years=max(fit$data$years)+(-9:0); rec.years=max(fit$data$years)+(-39:0); MPlabel=NULL;OMlabel=NULL;
overwriteSelYears=NULL; deterministic=FALSE;CZ=103000;HZ=257500;IE=NULL;rec.meth=3
capLower=0;capUpper=NULL;MP = rep('MPspm',3);UL.years=max(fit$data$years)+(-4:0);TAC.base=8000

ny=3
nosim=5

#***************************************************************************
#************* define Operating Models *************************************
#***************************************************************************

# base model
OMbase <- list(fit=fit1,
               nosim=nosim,
               OMlabel='OMbase',
               year.base=2016,
               ave.years=tail(fit1$data$years,10),
               rec.years=1977:2016,
               rec.meth=2,
               UL.years=tail(fit1$data$years,10))

copy(x=OMbase,n=c(5,3),name=c('OMcore','OMstress'))

# recruitment
OMcore1$rec.meth=6  # trailing sampling recruitment
OMcore2$rec.meth=1  # BH AC
OMstress1$rec.scale=0.7 # mean with AC and reduced by 0.75

# different M
newdat1 <- dat
newdat1$natMor[,] <- 0.15

fitMunc <- ccam.fit(newdat1,conf,par,phase=1)     # run uncensored with lognormal error dist Ctot (sd = 0.01) and Clower becomes Cmean
fitM <- ccam.fit(newdat1,conf,par,phase=2)     # run phase 1 + censored

OMcore3$fit=fitM

OMstress2$bio.scale=list('nm'=0.8)
OMstress3$bio.scale=list('nm'=1.2)

# different Upper limit
newdat2 <- dat
oldUpper <- newdat2$logobs[which(!is.na(newdat2$logobs[,2])),1]
newUpper <- log(exp(oldUpper) + ctUSA[-c(1:8),1])
newdat2$logobs[which(!is.na(newdat2$logobs[,2])),2] <- newUpper

fitCunc <- ccam.fit(newdat2,conf,par,phase=1)     # run uncensored with lognormal error dist Ctot (sd = 0.01) and Clower becomes Cmean
fitC <- ccam.fit(newdat2,conf,par,phase=2)     # run phase 1 + censored

OMcore4$fit=fitC

#*****************************************************************************
#************* define Harvest Control Rules **********************************
#*****************************************************************************

nMP=5

MP1 <- list(MP=rep(0,ny),
            MPlabel='MP1',
            CZ=103000,  #this should have uncertainty as well...
            HZ=257500,
            IE=NULL,
            capLower=0,
            TAC.base=8000)

copy(x=MP1,n=c(nMP-1),name=c('MP'))

avail('MP')

MP2$MP <- rep('MPeggsurvey',ny)
MP3$MP <- rep('MPeggsurveytrail3interim',ny)
MP4$MP <- rep('MPeggsurveytrail3',ny)
MP5$MP <- rep('MPeggsurveytargetinterim',ny)
MP6$MP <- rep('MPeggsurveytarget',ny)
MP7$MP <- rep('MPf40base',ny)
MP8$MP <- rep('MPspm',ny)


# Plot different IEs
IEmeans=list(IE1 = rep(6000,100),
             IE2 = rep(7200,100),
             IE3 = c(Reduce(function(v, x) .8*v , x=numeric(3),  init=6000, accumulate=TRUE)[-1],rep(3000,97)),
             IE4 = c(Reduce(function(v, x) .75*v , x=numeric(6),  init=6000, accumulate=TRUE)[-1],rep(1000,94)),
             IE5 = rep(6000*0.8,100)
)
IEsds=lapply(IEmeans,'/',3)

z <- array( c( do.call('rbind',IEmeans) , do.call('rbind',IEsds)  ) , dim = c( 5 , 100 , 2 ) )

png(file=paste0("./IMG/Ubars.png"),units="cm",res=300,width=18,height=13)
ggplot(U,aes(x=time,y=value))+
    geom_ribbon(aes(ymin=value-sd,ymax=value+sd,fill=variable),alpha=0.3)+
    geom_line(aes(col=variable))+
    geom_point(aes(col=variable))+
    theme(legend.position="none")+
    ylab('Undeclared catch (t)')+xlab('Year')
dev.off()

#******************************************************************************
#************* forecast for each combination **********************************
#******************************************************************************
myOMs <- 'OMbase'
myMPs <- 'MP1'

myOMs <- grep('OM',ls(),value = TRUE)
myMPs <- paste0('MP',1:nMP)

scenmat <- expand.grid(OM=myOMs, MP=myMPs)
scennames <- apply(scenmat,1,paste,collapse = ".")

# create a list with all scenarios to test (combos MP/OM)
scenlist <- lapply(split(scenmat,1:nrow(scenmat)),function(x){
    c(get(as.character(x[1,1])),get(as.character(x[1,2])))
})
names(scenlist) <- scennames

# forecast each scenario (combos MP/OM)
wdRdata <- "C:/Users/vanbe/Desktop/Post-Doc/DATA/MSE/Rdata/"

runlist <- lapply(scenlist,function(x){
    RUN <- do.call(forecast, x)
    save(RUN,file=paste0(wdRdata,names(x),'.Rdata'))
})

files <- split(paste0(wdRdata,scennames ,'.Rdata'),1:nrow(scenmat))
runlist <- lapply(files, function(x) get(load(x)))

names(runlist)=scennames
class(runlist)='forecastset'

#******************************************************************************
#************* plot each forecast and compare**********************************
#******************************************************************************

ssbplot(RUNbase.1)
ssbplot(RUNbase.list) #no applicable method for 'plotit' applied to an object of class "list"

catchplot(RUNbase.1)
catchplot(RUNbase.list)

recplot(RUNbase.1)
recplot(RUNbase.list)

fbarplot(RUNbase.1)
fbarplot(RUNbase.list)

foreplot(RUNbase.1,what='probCZ',zone=0.75)
foreplot(RUNbase.1,what='probHZ',zone=0.75)
foreplot(RUNbase.1,what='probGrowth')
foreplot(RUNbase.1,what='TAC')
foreplot(RUNbase.1,what='ssbmsyratio',ylab='ssb/ssbmsy')
foreplot(RUNbase.1,what='fmsyratio',ylab='F/Fmsy')
foreplot(RUNbase.1,what='Umsyratio',ylab='U/Umsy')


### Objective 1: rebuild out of critical zone and into healthy zone with 75% prob
foreplot(RUNbase.list,what='probCZ',zone=0.75,ylab='Probability out of the CZ')
foreplot(RUNbase.list,what='probHZ',zone=0.75,ylab='Probability into the HZ')

## Objective 2: maintain a positive growth trajectory
foreplot(RUNbase.list,what='probGrowth')


# stuff should be added here (Numberof years prob of growth below x%)


## Objective 3: Maximize annual catches
foreplot(RUNbase.list,what='catch')

## Obective 4: maximise fishery stability



## tradeplots
tradeplot(RUNbase.1,what.x='ULR', what.y='catch',zone=c(75,8000),ci=TRUE,lab=TRUE)
tradeplot(RUNbase.1,what.x='ULR', what.y='catchcumul',zone=c(75),ci=TRUE)

tradeplot(RUNbase.list,what.x='ULR', what.y='ssb',zone=c(75),ci=FALSE,leg=c(1,2))
tradeplot(RUNbase.list,what.x='ULR', what.y='catchcumul',zone=c(75),ci=FALSE,leg=c(1,2))

#-----------------------------------------------------------------------------------------
#-- compare recruitment methods ----------------------------------------------------------
#-----------------------------------------------------------------------------------------
wdIMG <- "C:/Users/vanbe/Desktop/Post-Doc/DATA/MSE/IMG/"

OMrec1 <- list(fit=fit3,
               nosim=1000,
               OMlabel='OMrec1',
               year.base=2016,
               ave.years=tail(fit1$data$years,10),
               rec.years=1977:2016,
               rec.meth=1,
               UL.years=tail(fit1$data$years,10))

MP1 <- list(MP=c(rep(0,40)),
            MPlabel='MP1',
            CZ=103000,  #this should have uncertainty as well...
            HZ=257500,
            IE=NULL,
            capLower=0,
            TAC.base=8000)
copy(x=OMrec1,n=c(9),name=c('OMrec'))
OMrec2$rec.meth=2
OMrec3$rec.meth=3
OMrec4$rec.meth=4
OMrec5$rec.meth=5
OMrec6$rec.meth=6
OMrec7$rec.meth=7
OMrec8$rec.meth=8

MSErecs <- lapply(as.list(1:8),function(x) c(get(paste0("OMrec",x)),MP1))
RUNrec.list <- lapply(MSErecs,function(x) do.call(forecast, x))
class(RUNrec.list)<-"forecastset"
names(RUNrec.list) <- c('BH.AC','mean.AC','BH.noAC','mean.noAC','samp','sampBack','sampRPS','sampRPSBack')

savepng(ssbplot(RUNrec.list,ci=FALSE),wdIMG,"recruitment.methods.ssb.all",c(22,14))
savepng(recplot(RUNrec.list,ci=FALSE),wdIMG,"recruitment.methods.rec.all",c(22,14))

RUNrec.list.reduced=RUNrec.list
RUNrec.list.reduced=RUNrec.list.reduced[c(1:6)]
class(RUNrec.list.reduced)<-"forecastset"

savepng(ssbplot(RUNrec.list.reduced,ci=FALSE),wdIMG,"recruitment.methods.ssb.6",c(22,14))
savepng(recplot(RUNrec.list.reduced,ci=FALSE,ylim=c(0,400000)),wdIMG,"recruitment.methods.rec.6",c(22,14))

savepng(foreplot(RUNrec.list.reduced,what='probCZ',zone=0.75)+scale_x_continuous(limits=c(2017,2030)),wdIMG,"recruitment.methods.CZ.6",c(22,14))
CZ=extract(RUNrec.list.reduced,'probCZ')
names(CZ)[1]='prob'
reasonable=ddply(CZ,'id',summarise,length(prob[prob<0.75]))
reasonable$id=names(RUNrec.list.reduced)
reasonable
