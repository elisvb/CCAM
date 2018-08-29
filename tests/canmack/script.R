#setwd("C:/Users/VANBE/Desktop/post-doc/DATA/CCAM/") #for development
setwd("C:/Users/VANBE/Desktop/post-doc/DATA/CCAM/tests/canmack")

wdIMG <- "C:/Users/vanbe/Desktop/Post-Doc/DATA/MSE/IMG/"
wdRdata <- "C:/Users/vanbe/Desktop/Post-Doc/DATA/MSE/Rdata/"

library(CCAM)

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
nm[,] <-0.15
pf <- read.ices("pf.dat")
pm <- read.ices("pm.dat")
sw <- read.ices("sw.dat")
surveys <- read.ices("survey.dat")
surveys[[1]] <- surveys[[1]][!is.na(surveys[[1]]),1,drop=FALSE]
attr(surveys[[1]],'time') <- c(0.5)
env <- read.ices("env.dat")
env[,1]=0

# define catch limits
ct[,1] <- ct[,1]*1.10 + ctUSA[-c(1:8),1]*0.25
ct[,2] <- ct[,2] + ctUSA[-c(1:8),1]*0.25

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

conf <- defcon(dat)
conf$keySel <- matrix(c(0,1,2,3,4,4,4,4,4,4), nrow=49, ncol=10,byrow = T)

#conf$keySel <- rbind(matrix(c(0,1,2,3,4,4,4,4,4,4), nrow=28, ncol=10,byrow = T),
#                     matrix(c(5,6,7,8,9,9,9,9,9,9), nrow=49-28, ncol=10,byrow = T))

conf$keyVarObs[1,]=-1                       # sd around total catch always estimated as nearly zero, and problem calculating sd
conf$keyVarObs[2,1:9]=c(0,1,2,2,2,2,2,2,3)
conf$keyVarObs[3,1]=4
conf$stockRecruitmentModelCode=2 #0: RW, 1: ricker, 2: BH, 3:mean
conf$fbarRange=c(5,10)
conf$obsLikelihoodFlag[1]='CE'

par <- defpar(dat,conf)

# plotInput script

save(dat,'')

#################################################################################################################
########### fit model ###########################################################################################
#################################################################################################################

fitBase2 <- ccam.fit(dat,conf,par,silent=FALSE,paracheck = FALSE,debug=TRUE)            # censored
fitBase
ssbplot(fitBase)
    #save(fitBase, file=paste0(wdRdata,'fitBase.Rdata'))
    #load(file=paste0(wdRdata,'fitBase.Rdata'))

check=fitBase
fits=c(fitBase,fitBase2)
names(fits)=c('check','cens')

fits=c(fitBase31,fitBase32,fitbase33,fitBase34,fitBase)


#################################################################################################################
########### PLOTS & TABLES ######################################################################################
#################################################################################################################

saveallplots(fitBase1969,wd=wdIMG,name='fitBase')

srplot(fitBase,curve=TRUE)
catchplot(fits)
ssbplot(fits,ci=TRUE)
expplot(fitBase)
recplot(fitBase,years=1969:2016)
recplot(fits)
peplot(fits)
refBase <- ypr(fitBase,rec.years=1980:2016)
plot(refBase)
plot(fitBase)
plotprod(fitBase)
selplot(fitBase)

modeltable(fits)
partable(fits)

# Residuals (not one step ahead!!!)
resplot(fitBase,fleets = 3,type=1)
resplot(fitBase,fleets = 3,type=2,out=1)
resplot(fitBase,fleets = 3,type=3)
resplot(fitBase,fleets = 3,type=4)
resplot(fitBase,fleets = 3,type=4,trans = exp)

resplot(fitBase,fleets = 2,type=1)
resplot(fitBase,fleets = 2,type=2,out=3)
resplot(fitBase,fleets = 2,type=3)
resplot(fitBase,fleets = 2,type=4,trans=crlInverse)
resplot(fitBase,fleets = 2,type=4)
resplot(fitBase,fleets = 2,type=5,std=TRUE)
resplot(fitBase,fleets = 2,type=6)
resplot(fitBase,fleets = 2,type=7)

#run mackerel in SAM and run in CCAM with cn
#unfinished
 res <- residuals(fitBase)  #problems because logobs has become a matrix (instead of vector)!!
 plot(res)
# resp <- procres(fitBase) #needs to be adapted because of separable F
# plot(resp)
# sims=simulate(fit)  #simulate data (not yet ok for lower and upper catch limits)

#retrospective analysis
retro <- retro(fitBase,year=7)  #maybe make plot with relative change
savepng(plot(retro,ci=FALSE),wdIMG,"/fitBase/retro",c(25,20))
m <- round(mohn(retro),2)
write.table(m,paste0(wdIMG,"/fitBase/mohn.txt"))

plot(unlist(lapply(lapply(retro,ypr),'[','f40ssb')),type='l',ylab='SSBF40%',xlab='peel')
plot(unlist(lapply(lapply(retro,ypr),'[','f40')),type='l',ylab='F40%',xlab='peel')

#################################################################################################################
########### MSE #################################################################################################
#################################################################################################################

##### debugging
fit=fitBase; fscale=NULL; catchval=NULL; fval=NULL; nosim=2; year.base=max(fit$data$years);
ave.years=max(fit$data$years)+(-9:0); rec.years=max(fit$data$years)+(-39:0); MPlabel=NULL;OMlabel=NULL;
overwriteSelYears=NULL; deterministic=FALSE;IE=NULL;rec.meth=3;bio.scale=NULL;rec.scale=1
capLower=0;capUpper=NULL;MP = rep('MPccam',5);UL.years=max(fit$data$years)+(-4:0);TAC.base=8000;deadzone=1000

ny=5
nosim=2

test=forecast(fit=fitBase, MP=rep('MPspm',7), nosim=200,
                     ave.years=max(fitBase$data$years)+(-9:0), rec.years=max(fitBase$data$years)+(-39:0),rec.meth=1,
                     UL.years=max(fitBase$data$years)+(-4:0),
                     deadzone=1000,verbose = TRUE)

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

MPmat=expand.grid(MP=paste0('MP',1:nMP), IE=c(paste0('IEnorm',1:nIE),'IEnothing'))
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
    #load(file=paste0(wdRdata,'scen.list.2018.06.12.Rdata'))

x=names(scen.list)[49] #ccam
x=names(scen.list)[1]
scen.list=scen.list[-grep('MP7',names(scen.list))]

# 48 is missing
x=names(scen.list)[grep('MP11',names(scen.list))][1]

sublist <- names(scen.list)[grep(c('MP10'),names(scen.list))]
x <- sublist[1]
runlist <- lapply(sublist,function(x){ #ccam
    y <- scen.list[[x]]
    #y$nosim <- 10
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

runlist=runlist[-grep(c('MP10|MP11'),names(runlist))]
class(runlist)='forecastset'


runlist=runlist[-grep(c('MP10|MP11'),names(runlist))]
class(runlist)='forecastset'



#******************************************************************************
#************* plot each forecast and compare**********************************
#******************************************************************************

### with no fishing at all
savepng(ssbplot(runlist[[6]]),wdIMG,'/MSE/SSB_base_0fish',c(17,10))
savepng(fbarplot(runlist[[6]]),wdIMG,'/MSE/Fbar_base_0fish',c(17,10))
savepng(foreplot(runlist[[6]],what.y='probCZ',rect=0.75),wdIMG,'/MSE/CZ_base_0fish',c(17,10))
savepng(recplot(runlist[[6]]),wdIMG,'/MSE/rec_base_0fish',c(17,10))
savepng(foreplot(runlist[[6]],what.y='probHZ',rect=0.75),wdIMG,'/MSE/HZ_base_0fish',c(17,10))
savepng(foreplot(runlist[[6]],what.y='ssbmsyratio',ylab='ssb/ssbmsy'),wdIMG,'/MSE/SSBmsy_base_0fish',c(17,10))
savepng(foreplot(runlist[[6]],what.y='fmsyratio',ylab='F/Fmsy'),wdIMG,'/MSE/Fmsy_base_0fish',c(17,10))
savepng(foreplot(runlist[[6]],what.y='Umsyratio',ylab='U/Umsy'),wdIMG,'/MSE/Umsy_base_0fish',c(17,10))


### with 8000 t and 6000 undeclared
usual=runlist[grep(c('.MP2.IEnormgamma1'),names(runlist))]
class(usual)='forecastset'

savepng(ssbplot(usual,ci=FALSE)+scale_y_continuous(limits=c(0,800000)),wdIMG,'/MSE/SSB_TAC8000',c(25,15))
savepng(fbarplot(usual,ci=FALSE),wdIMG,'/MSE/Fbar_TAC8000',c(25,15))
savepng(foreplot(usual,what.y='probCZ',rect=0.75),wdIMG,'/MSE/CZ_TAC8000',c(17,10))
savepng(recplot(usual,ci=FALSE)+scale_y_continuous(limits=c(0,1000000)),wdIMG,'/MSE/rec_TAC8000',c(25,15))
savepng(foreplot(usual,what.y='probHZ',rect=0.75),wdIMG,'/MSE/HZ_TAC8000',c(17,10))
savepng(foreplot(usual,what.y='ssbmsyratio',ylab='ssb/ssbmsy',hline=1),wdIMG,'/MSE/SSBmsy_TAC8000',c(17,10))
savepng(foreplot(usual,what.y='fmsyratio',ylab='F/Fmsy',hline=1),wdIMG,'/MSE/Fmsy_TAC8000',c(17,10))
savepng(foreplot(usual,what.y='Umsyratio',ylab='U/Umsy',hline=1),wdIMG,'/MSE/Umsy_TAC8000',c(17,10))


# trade off plots
savepng(foreplot(runlist,what.y='ssb',by=c('OM','MP'),ci=FALSE),wdIMG,'/MSE/SSB_ALL',c(40,25))

foreplot(runlist,what.y='catchcumul',what.x='probCZ',by='OM',ci=FALSE)
foreplot(runlist,what.y='Umsyratio',what.x='ssbmsyratio',by='OM',ci=FALSE,hline=1,vline=1)

### Objective 1: rebuild out of critical zone and into healthy zone with 75% prob
savepng(foreplot(runlist,what.y='probCZ',ylab='Probability out of the CZ',by=c('OM','MP'),vline=c(5,10)+2016,rect=0.75),wdIMG,'/MSE/CZ_ALL',c(40,25))

savepng(diamondplot(runlist,what='probCZ',ylab='Probability out of the CZ',year=2016+5,hline=0.75)+ scale_color_manual(values=c('orange','grey25','grey60')),wdIMG,'/MSE/CZ_diamond',c(25,17))
savepng(diamondplot(runlist,what='probHZ',ylab='Probability into the HZ',year=2016+10,hline=0.75)+ scale_color_manual(values=c('orange','grey25','grey60')),wdIMG,'/MSE/HZ_diamond',c(25,17))
savepng(diamondplot(runlist,what='catch',ylab='Average catch')+ scale_color_manual(values=c('orange','grey25','grey60')),wdIMG,'/MSE/Catchmean_diamond',c(25,17))

## Objective 2: maintain a positive growth trajectory
foreplot(runlist,what.y='probgrowth30',by=c('OM','MP'))
foreplot(runlist,what.y='probgrowth30',by='OM')

# stuff should be added here (Numberof years prob of growth below x%)


## Objective 3: Maximize annual catches
foreplot(runlist,what='catch',by=c('OM','MP'))

## Obective 4: maximise fishery stability
foreplot(runlist,what.y='TACrel',by=c('OM','MP'))

## all objectives
savepng(MSEplot(runlist),wdIMG,'/MSE/objectives_ALL',c(40,25))




