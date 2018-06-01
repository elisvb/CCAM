setwd("C:/Users/VANBE/Desktop/post-doc/DATA/CCAM/")
setwd("C:/Users/VANBE/Desktop/post-doc/DATA/CCAM/tests/canmack")

wdIMG <- "C:/Users/vanbe/Desktop/Post-Doc/DATA/MSE/IMG/"
wdRdata <- "C:/Users/vanbe/Desktop/Post-Doc/DATA/MSE/Rdata/"

library(CCAM)
library(plyr)

#################################################################################################################
########### READ IN DATA ########################################################################################
#################################################################################################################

cn <- read.ices("cn.dat")
ct <- read.ices("ct.dat")
ctUSA <- read.ices("ctUSA.dat")
ctForeign <- read.ices("ctForeign.dat")
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
            # bubble(cn,col='darkblue') # raw data
            # bubble(t(spay(t(cn))))   # to check cohorts
            # bubble(t(spya(t(cn))))   # to compare between years
            # heat(mo)
            # heat(cw)
            # matplot(mo,type='l',ylab='Proportion mature',xlab='Year')
            # matplot(cw,type='l',ylab='Catch weight',xlab='Year')
            # matplot(ct,type='l',ylab='Catch (t)', xlab='Year')
            # surveyplot(surveys)
            # allC <- cbind(ct[,1],ctUSA[-c(1:8),1],ctForeign[-c(1:8),1])
            # allC <- cbind(allC, rowSums(allC))
            # matplot(rownames(allC),allC,type='l',lty = c(1,1,1,2),ylab='Catch',xlab='Year',bty="l")
            # legend('topright', c('CAN','USA','Foreign','Total'),fill=seq_len(ncol(allC)),box.col = 'white',cex=1.2)

#################################################################################################################
########### fit model ###########################################################################################
#################################################################################################################

conf$obsLikelihoodFlag[1]='LN'
fit1 <- ccam.fit(dat,conf,par,silent=TRUE)              # lognormal distribution (sd 0.000001)

conf$obsLikelihoodFlag[1]='CE'
fit2 <- ccam.fit(dat,conf,par,silent=TRUE)              # censored

fitBase=fit2
    #save(fitBase, file=paste0(wdRdata,'fitBase.Rdata'))
    #load(file=paste0(wdRdata,'fitBase.Rdata'))

fits=c(fit1,fit2)
names(fits)=c('uncens','cens')


#################################################################################################################
########### PLOTS & TABLES ######################################################################################
#################################################################################################################

srplot(fit2,curve=TRUE)

catchplot(fit2,fleet = 1)
catchplot(fits)

ssbplot(fit2)
ssbplot(fits,ci=TRUE)

fbarplot(fit2)
fbarplot(fits)
expplot(fit2)
expplot(fits)

recplot(fit2)
recplot(fits)

peplot(fit2)
peplot(fits)

y <- ypr(fit2,rec.years=1980:2016)
y
plot(y)

plot(fit2)

modeltable(fits)
partable(fits)

# Residuals (not one step ahead!!!)
resplot(fit2,fleets = 3,type=1)
resplot(fit2,fleets = 3,type=2,out=1)
resplot(fit2,fleets = 3,type=3)
resplot(fit2,fleets = 3,type=4)

resplot(fit2,fleets = 2,type=1)
resplot(fit2,fleets = 2,type=2,out=3)
resplot(fit2,fleets = 2,type=3,out=3)
resplot(fit2,fleets = 2,type=4)
resplot(fit2,fleets = 2,type=5)

#unfinished
# res <- residuals(fit1)  #problems because logobs has become a matrix (instead of vector)!!
# plot(res)
# resp <- procres(fit) #needs to be adapted because of separable F
# plot(resp)
# sims=simulate(fit)  #simulate data (not yet ok for lower and upper catch limits)

#retrospective analysis
retro <- retro(fit2,year=7)  #maybe make plot with relative change
plot(retro,ci=FALSE)
mohn(retro)

#################################################################################################################
########### MSE #################################################################################################
#################################################################################################################

##### debugging
fit=fitBase; fscale=NULL; catchval=NULL; fval=NULL; nosim=2; year.base=max(fit$data$years);
ave.years=max(fit$data$years)+(-9:0); rec.years=max(fit$data$years)+(-39:0); MPlabel=NULL;OMlabel=NULL;
overwriteSelYears=NULL; deterministic=FALSE;CZ=103000;HZ=257500;IE=NULL;rec.meth=3;bio.scale=NULL;rec.scale=1
capLower=0;capUpper=NULL;MP = rep('MPspm',5);UL.years=max(fit$data$years)+(-4:0);TAC.base=8000

ny=15
nosim=100

#***************************************************************************
#************* define Operating Models *************************************
#***************************************************************************

#--------------------- base model ------------------------------------------
OMbase <- list(fit=fitBase,
               nosim=nosim,
               OMlabel='OMbase',
               year.base=2016,
               ave.years=tail(fit1$data$years,10),
               rec.years=1977:2016,
               rec.meth=2,
               UL.years=tail(fit1$data$years,10))

copy(x=OMbase,n=c(5,3),name=c('OMcore','OMstress'))

# --------------------- recruitment ------------------------------------------
OMcore1$rec.meth=6  # trailing sampling recruitment
OMcore2$rec.meth=1  # BH AC
OMstress1$rec.scale=0.7 # mean with AC and reduced by 0.75

# ---------------------  M ----------------------------------------------------
newdat1 <- dat
newdat1$natMor[,] <- 0.15

fitM <- ccam.fit(newdat1,conf,par,phase=2)     # run phase 1 + censored

    #save(fitM, file=paste0(wdRdata,'fitM.Rdata'))
    #load(file=paste0(wdRdata,'fitM.Rdata'))

OMcore3$fit=fitM

OMstress2$bio.scale=list('nm'=0.8)
OMstress3$bio.scale=list('nm'=1.2)

# --------------------- Upper limit ------------------------------------------
newdat2 <- dat
oldUpper <- newdat2$logobs[which(!is.na(newdat2$logobs[,2])),1]
newUpper <- log(exp(oldUpper) + ctUSA[-c(1:8),1] + ctForeign[-c(1:8),1])
newdat2$logobs[which(!is.na(newdat2$logobs[,2])),2] <- newUpper

fitC <- ccam.fit(newdat2,conf,par,phase=2)     # run phase 1 + censored

    #save(fitC, file=paste0(wdRdata,'fitC.Rdata'))
    #load(file=paste0(wdRdata,'fitC.Rdata'))

OMcore4$fit=fitC

# --------------------- OM list ------------------------------------------

OM.list=list(OMbase=OMbase,
             OMcore1=OMcore1,
             OMcore2=OMcore2,
             OMcore3=OMcore3,
             OMcore4=OMcore4,
             OMstress1=OMstress1,
             OMstress2=OMstress2,
             OMstress3=OMstress3)


OMfits=c(fitBase=fitBase,FitM=fitM,fitC=fitC)
ssbplot(OMfits,ci=FALSE)
catchplot(OMfits)

#*****************************************************************************
#************* define Harvest Control Rules **********************************
#*****************************************************************************

# --------------------- base MPs ----------------------------------------------

nMP=8

MP1 <- list(MP=rep(0,ny),
            MPlabel='MP1',
            CZ=103000,  #this should have uncertainty as well...
            HZ=257500,
            IE=NULL,
            capLower=0,
            TAC.base=10000)

copy(x=MP1,n=c(nMP),name=c('MP'))

avail('MP')

MP2$MP <- rep('MPeggsurvey',ny)
MP3$MP <- rep('MPeggsurveytrail3interim',ny)
MP4$MP <- rep('MPeggsurveytrail3',ny)
MP5$MP <- rep('MPeggsurveytargetinterim',ny)
MP6$MP <- rep('MPeggsurveytarget',ny)
MP7$MP <- rep('MPf40base',ny)
MP8$MP <- rep('MPspm',ny)

# --------------------- MP list (include different IEs) ------------------------------------------

nIE=5

MPmat=expand.grid(MP=paste0('MP',1:nMP), IE=c(paste0('IEnormgamma',1:nIE),'IEnothing'))
MP.list <- lapply(split(MPmat,1:nrow(MPmat)),function(x){
    MPx <- get(as.character(x[1,1]))
    if(!is.na(x[1,2])) MPx$IE <- rep(as.character(x[1,2]),ny)
    return(MPx)
})
names(MP.list) <- paste(MPmat[,1],MPmat[,2],sep=".")

# --------------------- plot IE  ------------------------------------------

IEmeans=list(IE1 = rep(6000,100),
             IE2 = rep(7200,100),
             IE3 = c(Reduce(function(v, x) .8*v , x=numeric(3),  init=6000, accumulate=TRUE)[-1],rep(3000,97)),
             IE4 = c(Reduce(function(v, x) .75*v , x=numeric(6),  init=6000, accumulate=TRUE)[-1],rep(1000,94)),
             IE5 = rep(6000*0.8,100),
             IE6 = rep(0,100)
)
IEmeans=lapply(IEmeans,function(x) c(6000,x))
IEsds=lapply(IEmeans,'/',4)
savepng(IEplot(IEmeans,IEsds)+ scale_x_continuous(limits = c(0,ny+1),expand = c(0,0)),wdIMG,"/HCR/IE",c(9,5))

#******************************************************************************
#************* forecast for each combination **********************************
#******************************************************************************
#OM.list <- list(OMbase=OMbase)
#MP.list <- list(MP1=MP1)

scenmat <- expand.grid(OM=names(OM.list), MP=names(MP.list))
scennames <- apply(scenmat,1,paste,collapse = ".")

# create a list with all scenarios to test (combos MP/OM)
scen.list <- lapply(split(scenmat,1:nrow(scenmat)),function(x){
    c(OM.list[[as.character(x[1,1])]],MP.list[[as.character(x[1,2])]])
})
names(scen.list) <- scennames

length(scen.list)

# forecast each scenario (combos MP/OM)
Date = "2018-05-18"
Date = Sys.Date()
DateDir = paste0(wdRdata,Date,"/")
dir.create(DateDir)

    #save(scen.list, file=paste0(wdRdata,'scen.list.Rdata'))
    #load(file=paste0(wdRdata,'scen.list.Rdata'))

x=names(scen.list)[48]


scen.list=scen.list[-grep('MP7',names(scen.list))]

# 48 is missing
runlist <- lapply(names(scen.list)[-grep('MP8',names(scen.list))],function(x){
    y <- scen.list[[x]]
    RUN <- do.call(forecast, y)
    save(RUN,file=paste0(DateDir,x,'.Rdata'))
    return(RUN)
})

runlist <- lapply(names(scen.list)[grep('MP8',names(scen.list))],function(x){
    y <- scen.list[[x]]
    RUN <- do.call(forecast, y)
    save(RUN,file=paste0(DateDir,x,'.Rdata'))
    return(RUN)
})

filenames <- dir(DateDir, pattern = ".Rdata")
files <- paste0(DateDir,filenames)
  #files <- split(paste0(DateDir,names(scen.list) ,'.Rdata'),1:length(names(scen.list)))
runlist <- lapply(files, function(x) {print(x);get(load(x))})
names(runlist)=gsub(pattern = ".Rdata",replacement = "",x = filenames)
class(runlist)='forecastset'

#******************************************************************************
#************* plot each forecast and compare**********************************
#******************************************************************************

# examples of plots for 1 forecast
foreplot(runlist[[1]],what.y='probCZ',rect=0.75)
foreplot(runlist[[1]],what.y='probHZ',rect=0.75)
foreplot(runlist[[1]],what.y='probGrowth')
foreplot(runlist[[1]],what.y='TAC')
foreplot(runlist[[1]],what.y='ssbmsyratio',ylab='ssb/ssbmsy')
foreplot(runlist[[1]],what.y='fmsyratio',ylab='F/Fmsy')
foreplot(runlist[[1]],what.y='Umsyratio',ylab='U/Umsy')

# trade off plots
foreplot(runlist,what.y='catchcumul',what.x='probCZ',by='OM',ci=FALSE)
foreplot(runlist,what.y='Umsyratio',what.x='ssbmsyratio',by='OM',ci=FALSE,hline=1,vline=1)

### Objective 1: rebuild out of critical zone and into healthy zone with 75% prob
foreplot(runlist,what.y='probCZ',ylab='Probability out of the CZ',by=c('OM','MP'),vline=c(5,10)+2016,rect=0.75)

## Objective 2: maintain a positive growth trajectory
foreplot(runlist,what.y='probgrowth20',by=c('OM','MP'))
foreplot(runlist,what.y='probgrowth20',by='OM')

# stuff should be added here (Numberof years prob of growth below x%)


## Objective 3: Maximize annual catches
foreplot(runlist,what='catch',by=c('OM','MP'))

## Obective 4: maximise fishery stability

## all objectives
MSEplot(runlist)


#-----------------------------------------------------------------------------------------
#-- compare recruitment methods ----------------------------------------------------------
#-----------------------------------------------------------------------------------------

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

### test with just one scenario


ny=5

OMbase <- list(fit=fitBase,
               nosim=5,
               OMlabel='OMbase',
               year.base=2016,
               ave.years=tail(fit1$data$years,10),
               rec.years=1977:2016,
               rec.meth=3,
               UL.years=tail(fit1$data$years,10))

MP1 <- list(fval=rep(0.01,ny),
            MPlabel='MP1',
            CZ=103000,  #this should have uncertainty as well...
            HZ=257500,
            IE=NULL,
            capLower=0,
            TAC.base=8000)


MSEbase.1 <- c(OMbase,MP1)

RUNbase.1 <- do.call(forecast, MSEbase.1)

