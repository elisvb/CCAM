#library(testthat)
#devtools::install_github("iml-mackerel/CCAM")
#usethis::use_testthat(3)
#usethis::use_test()
library(TMB)
compile('src/CCAM.cpp');dyn.load(dynlib("src/CCAM")) # check c++ code pour grammar error


#library(gert)
#git_branch_list()
#dans le terminal
#  git diff origin/master #pour comparer avec la version d'orignien tous les commits depuis
#git diff origin/master..lmPeriod > diff_complet.txt four faire un fichier text mais sans culeur c,est difficile ligne précédente est mieux

library(devtools)
load_all()
#library(CCAM)

test_2024_fit <-  function(){
dir <- "example/2024/"

##########FULL READ########

cn <- read.ices(paste0(dir,'cn.dat'))
ct <- read.ices(paste0(dir,'ct.dat'))
ctUSA <- read.ices(paste0(dir,'ctUSA.dat'))

ctForeign <- read.ices(paste0(dir,'ctForeign.dat'))
ctunac <- read.ices(paste0(dir,'ctUnaccounted.dat'))
cw <- read.ices(paste0(dir,'cw.dat'))
mo <- read.ices(paste0(dir,'mo.dat'))
nm <- read.ices(paste0(dir,'nm.dat'))

pf <- read.ices(paste0(dir,'pf.dat'))
pm <- read.ices(paste0(dir,'pm.dat'))
sw <- read.ices(paste0(dir,'sw.dat'))
sw0 <- read.ices(paste0(dir,'sw0.dat'))
survey <- read.ices(paste0(dir,'tep.dat'))
survey[[1]] <- survey[[1]][!is.na(survey[[1]]),1,drop=FALSE]
survey[[1]][,1] <- survey[[1]][,1]*10^9  #to verify this in reality should be 10^12, but in the model conversion from kg to t not done, so corrected here instead.
attr(survey[[1]],'time') <- c(0.47)
pfem <- read.ices(paste0(dir,'propFemale.dat'))
fec <- read.ices(paste0(dir,'fec.dat'))

# redefine catch limits (add 20-80% US catch)
ctwusa <- ct                                   # reported landings
ctwusa[,2] <- ctwusa[,1] + ctunac[,1]          # add max missing in Canada
ctwusa[,1] <- ctwusa[,1]*1.10 + ctUSA[,1]*0.2     # lower bound: increase and add US
ctwusa[,2] <- ctwusa[,2] + ctUSA[,1]*0.8          # add us upper


dat <- setup.ccam.data(surveys=survey,
                       residual.fleet=cn[-1,],
                       total.catch=ctwusa[-1,],
                       prop.mature=mo[-1,],
                       stock.mean.weight=sw[-1,],
                       stock.start.weight=sw0[-1,],
                       catch.mean.weight=cw[-1,],
                       prop.f=pf[-1,],
                       prop.m=pm[-1,],
                       natural.mortality=nm[-1,],
                       prop.fem=pfem[-1,],
                       fec=fec[-1,])

conf <- defcon(dat)
conf$keySel <- matrix(c(0,1,2,3,4,4,4,4,4,4), nrow=nrow(conf$keySel), ncol=ncol(conf$keySel),byrow = T)
conf$keyVarObs[1,]=-1
conf$keyVarObs[2,1:9]=c(0,1,2,2,2,2,1,1,1) # 2024
#conf$keyVarObs[2,1:9]=c(0,1,2,2,2,2,2,2,1) # 2022


conf$keyVarObs[3,1]=3
conf$stockRecruitmentModelCode=2 #0: RW, 1: ricker, 2: BH, 3:mean
conf$obsLikelihoodFlag[1]='CE'
conf$keyBiomassTreat[3]=5
conf$fbarRange=c(5,10) #fully recruited fish
conf$do_lmPeriod= NULL
conf$lmPeriod= NULL
#conf$iyr_lmPeriod = (length(dat$years) - (conf$lmPeriod)) : (length(dat$years)-1)
conf$iyr_lmPeriod = NULL
#length(conf$iyr_lmPeriod)

par <- defpar(dat,conf)
par$alpha_lmPeriod=NULL
par$beta_lmPeriod=NULL
par$log_sig_lmPeriod=NULL
fit <- ccam.fit(dat,conf,par,debug=T)

##########END OF FULL READ########
load(paste0(dir, "fit.Rdata"))
load(paste0(dir, "conf.Rdata"))
load(paste0(dir, "dat.Rdata"))
load(paste0(dir, "par.Rdata"))


fit <- ccam.fit(data=dat,conf,par,debug=T)
fit

return(fit)
}

#load_all()
#test_that("get_fit ", {
#    expect_snapshot_output(cat(get_fit()))
#    })
