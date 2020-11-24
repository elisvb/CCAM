library(CCAM)

#################################################################################################################
########### READ IN DATA ########################################################################################
#################################################################################################################

cn <- read.ices('cn.dat')
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
sw0 <- read.ices("sw0.dat")
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
                      stock.start.weight=sw0,
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
conf$keyVarObs[2,1:9]=c(0,1,2,2,2,2,2,2,1)
conf$keyVarObs[3,1]=3
conf$stockRecruitmentModelCode=2 #0: RW, 1: ricker, 2: BH, 3:mean
conf$fbarRange=c(5,10)
conf$obsLikelihoodFlag[1]='CE'

par <- defpar(dat,conf)

#################################################################################################################
########### fit model ###########################################################################################
#################################################################################################################

fitBase <- ccam.fit(dat,conf,par,paracheck = FALSE)            # censored

plot(fitBase)

