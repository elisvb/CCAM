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
ct[,1]=ct[,1]*1.10
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
conf$obsLikelihoodFlag[1]='CE'

par <- defpar(dat,conf)

            # check input data
            savepng(bubble(cn,col='darkblue',scale = 10),wdIMG,"/input/cn.raw",c(15,8)) # raw data
            savepng(bubble(t(spay(t(cn))),scale = 8),wdIMG,"/input/cn.spay.compCohorts",c(15,8))   # to check cohorts
            savepng(bubble(t(spya(t(cn))),scale = 8),wdIMG,"/input/cn.spya.comyYears",c(15,8))   # to compare between years
            savepng(heat(mo),wdIMG,"/input/pm.heat",c(15,8))
            savepng(heat(cw),wdIMG,"/input/cw.heat",c(15,8))
            savepng(matplot2(mo,ylab='Proportion mature',xlab='Year'),wdIMG,"/input/pm.line",c(15,8))
            savepng(matplot2(cw,ylab='Catch weight',xlab='Year'),wdIMG,"/input/cw.line",c(15,8))
            savepng(matplot2(ct,ylab='Catch (t)', xlab='Year',col=c('black','darkgrey')),wdIMG,"/input/ct",c(15,8))
            savepng(surveyplot(surveys),wdIMG,"/input/survey",c(15,8))
            allC <- cbind(Canada=ct[,1],USA=ctUSA[-c(1:8),1],Foreign=ctForeign[-c(1:8),1])
            allC <- cbind(allC, Total=rowSums(allC))
            savepng(matplot2(allC,ylab='Catch (t)', xlab='Year',col=c('orange','yellowgreen','mediumorchid','black')),wdIMG,"/input/ct.all",c(15,8))
            plotobs(dat,type='line')

#################################################################################################################
########### fit model ###########################################################################################
#################################################################################################################

fitBase <- ccam.fit(dat,conf,par,silent=TRUE)            # censored

    #save(fitBase, file=paste0(wdRdata,'fitBase.Rdata'))
    #load(file=paste0(wdRdata,'fitBase.Rdata'))

fits=c(fitBase,fitBase)
names(fits)=c('cens','cens')

#################################################################################################################
########### PLOTS & TABLES ######################################################################################
#################################################################################################################

saveallplots(fitBase,wd=wdIMG,name='test')

srplot(fitBase,curve=TRUE)
catchplot(fits)
ssbplot(fits,ci=TRUE)
expplot(fitBase)
recplot(fitBase,years=1969:2016)
recplot(fits)
peplot(fits)
y <- ypr(fits,rec.years=1980:2016)
plot(refBase)
plot(fitBase)

modeltable(fits)
partable(fits)

# Residuals (not one step ahead!!!)
resplot(fitBase,fleets = 3,type=1)
resplot(fitBase,fleets = 3,type=2,out=1)
resplot(fitBase,fleets = 3,type=3)
resplot(fitBase,fleets = 3,type=4)
resplot(fitBase,fleets = 3,type=4,trans = exp)

resplot(fitBase,fleets = 2,type=1,low=c('red','orange'),high=c('grey','green','darkgreen'))
resplot(fitBase,fleets = 2,type=2,out=3)
resplot(fitBase,fleets = 2,type=3)
resplot(fitBase,fleets = 2,type=4)
resplot(fite,fleets = 2,type=5)
resplot(fitBase,fleets = 2,type=6)
resplot(fitBase,fleets = 2,type=7)

#unfinished
# res <- residuals(fit1)  #problems because logobs has become a matrix (instead of vector)!!
# plot(res)
# resp <- procres(fit) #needs to be adapted because of separable F
# plot(resp)
# sims=simulate(fit)  #simulate data (not yet ok for lower and upper catch limits)

#retrospective analysis
retro <- retro(fitBase,year=7)  #maybe make plot with relative change
savepng(plot(retro,ci=FALSE),wdIMG,"/fitBase/retro",c(25,20))
m <- round(mohn(retro),2)
write.table(m,paste0(wdIMG,"/fitBase/mohn.txt"))

#################################################################################################################
########### MSE #################################################################################################
#################################################################################################################

##### debugging
fit=fitBase; fscale=NULL; catchval=NULL; fval=NULL; nosim=2; year.base=max(fit$data$years);
ave.years=max(fit$data$years)+(-9:0); rec.years=max(fit$data$years)+(-39:0); MPlabel=NULL;OMlabel=NULL;
overwriteSelYears=NULL; deterministic=FALSE;IE=NULL;rec.meth=3;bio.scale=NULL;rec.scale=1
capLower=0;capUpper=NULL;MP = rep('MPf40base',5);UL.years=max(fit$data$years)+(-4:0);TAC.base=8000

ny=15
nosim=50

#***************************************************************************
#************* define Operating Models *************************************
#***************************************************************************

#--------------------- base model ------------------------------------------
OMbase <- list(fit=fitBase,
               nosim=nosim,
               OMlabel='OMbase',
               year.base=2016,
               ave.years=tail(fitBase$data$years,10),
               rec.years=1977:2016,
               rec.meth=2,
               UL.years=tail(fitBase$data$years,10))

copy(x=OMbase,n=c(5,3),name=c('OMcore','OMstress'))

# --------------------- recruitment ------------------------------------------
OMcore1$rec.meth=6  # trailing sampling recruitment
OMcore2$rec.meth=1  # BH AC
OMstress1$rec.scale=0.7 # mean with AC and reduced by 0.75

# ---------------------  M ----------------------------------------------------
newdat1 <- dat
newdat1$natMor[,] <- 0.15

fitM <- ccam.fit(newdat1,conf,par)     # run phase 1 + censored

    #save(fitM, file=paste0(wdRdata,'fitM.Rdata'))
    #load(file=paste0(wdRdata,'fitM.Rdata'))

OMcore3$fit=fitM

OMstress2$bio.scale=list('nm'=0.8)
OMstress3$bio.scale=list('nm'=1.2)

# --------------------- Upper limit ------------------------------------------
newdat2 <- dat
oldUpper <- newdat2$logobs[which(!is.na(newdat2$logobs[,2])),1]
newUpper <- log(exp(oldUpper) + 0.75*(ctUSA[-c(1:8),1] + ctForeign[-c(1:8),1]))
newdat2$logobs[which(!is.na(newdat2$logobs[,2])),2] <- newUpper

fitC <- ccam.fit(newdat2,conf,par)     # run phase 1 + censored

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
catchplot(OMfits,ci=FALSE)

savepng(ssbplot(OMfits,ci=FALSE)+scale_y_continuous(limits=c(0,5e5),expand=c(0,0)),wdIMG,"/OMs/ssb",c(17,10))
savepng(catchplot(OMfits,ci=FALSE)+scale_y_continuous(limits=c(0,1e5),expand=c(0,0))+ylab('Catch'),wdIMG,"/OMs/catch",c(17,10))

#*****************************************************************************
#************* define Harvest Control Rules **********************************
#*****************************************************************************

# --------------------- base MPs ----------------------------------------------

nMP=11

MP1 <- list(MPlabel='MP1',
            IE=NULL,
            capLower=0,
            TAC.base=10000)

copy(x=MP1,n=nMP,name=c('MP'))

avail('MP')
MP1$catchval <- rep(0,ny)
MP2$catchval <- rep(8000,ny)
MP3$catchval <- rep(10000,ny)
MP4$catchval <- rep(30000,ny)
MP5$MP <- rep('MPeggsurvey',ny)
MP6$MP <- rep('MPeggsurveytrail3interim',ny)
MP7$MP <- rep('MPeggsurveytrail3',ny)
MP8$MP <- rep('MPeggsurveytargetinterim',ny)
MP9$MP <- rep('MPeggsurveytarget',ny)
MP10$MP <- rep('MPccam',ny)
MP11$MP <- rep('MPspm',ny)


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
Date = "2018-06-12"
Date = Sys.Date()
DateDir = paste0(wdRdata,Date,"/")
dir.create(DateDir)

    #save(scen.list, file=paste0(wdRdata,'scen.list.2018.06.12.Rdata'))
    #load(file=paste0(wdRdata,'scen.list.Rdata'))

x=names(scen.list)[49] #ccam
x=names(scen.list)[1]
scen.list=scen.list[-grep('MP7',names(scen.list))]

# 48 is missing
x=names(scen.list)[grep('MP11',names(scen.list))][1]

sublist <- names(scen.list)[-grep(c('MP11|MP10'),names(scen.list))]
runlist <- lapply(sublist,function(x){ #ccam
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

### with no fishing at all
ssbplot(runlist[[6]])
foreplot(runlist[[6]],what.y='ssb')
foreplot(runlist[[6]],what.y='probCZ',rect=0.75)


# examples of plots for 1 forecast
foreplot(runlist[[6]],what.y='probCZ',rect=0.75)
foreplot(runlist[[6]],what.y='probHZ',rect=0.75)
foreplot(runlist[[6]],what.y='probGrowth')
foreplot(runlist[[6]],what.y='TAC')
foreplot(runlist[[6]],what.y='ssbmsyratio',ylab='ssb/ssbmsy')
foreplot(runlist[[6]],what.y='fmsyratio',ylab='F/Fmsy')
foreplot(runlist[[6]],what.y='Umsyratio',ylab='U/Umsy')

# trade off plots
foreplot(runlist,what.y='catchcumul',what.x='probCZ',by='OM',ci=FALSE)
foreplot(runlist,what.y='Umsyratio',what.x='ssbmsyratio',by='OM',ci=FALSE,hline=1,vline=1)

### Objective 1: rebuild out of critical zone and into healthy zone with 75% prob
foreplot(runlist,what.y='probCZ',ylab='Probability out of the CZ',by=c('OM','MP'),vline=c(5,10)+2016,rect=0.75)

## Objective 2: maintain a positive growth trajectory
foreplot(runlist,what.y='probgrowth30',by=c('OM','MP'))
foreplot(runlist,what.y='probgrowth30',by='OM')

# stuff should be added here (Numberof years prob of growth below x%)


## Objective 3: Maximize annual catches
foreplot(runlist,what='catch',by=c('OM','MP'))

## Obective 4: maximise fishery stability
foreplot(runlist,what.y='TACrel',by=c('OM','MP'))

## all objectives
MSEplot(runlist)


