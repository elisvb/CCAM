#library(testthat)
#devtools::install_github("iml-mackerel/CCAM")
#usethis::use_testthat(3)
#usethis::use_test()
library(TMB)
#compile('src/CCAM.cpp');dyn.load(dynlib("src/CCAM")) # check c++ code pour grammar error

#library(gert)
#git_branch_list()
#dans le terminal
#  git diff origin/master #pour comparer avec la version d'orignien tous les commits depuis
#git diff origin/master..lmPeriod > diff_complet.txt four faire un fichier text mais sans culeur c,est difficile ligne précédente est mieux

library(devtools)
load_all() # works better of using install in the build tab
#library(CCAM)

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


datnew <- setup.ccam.data(surveys=survey,
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

confnew <- defcon(datnew)
confnew$keySel <- matrix(c(0,1,2,3,4,4,4,4,4,4), nrow=nrow(confnew$keySel), ncol=ncol(confnew$keySel),byrow = T)
confnew$keyVarObs[1,]=-1
confnew$keyVarObs[2,1:9]=c(0,1,2,2,2,2,1,1,1) # 2024
#confnew$keyVarObs[2,1:9]=c(0,1,2,2,2,2,2,2,1) # 2022


confnew$keyVarObs[3,1]=3
confnew$stockRecruitmentModelCode=2 #0: RW, 1: ricker, 2: BH, 3:mean
confnew$obsLikelihoodFlag[1]='CE'
confnew$keyBiomassTreat[3]=5
confnew$fbarRange=c(5,10) #fully recruited fish
confnew$do_lmPeriod= 1
confnew$lmPeriod= 15
confnew$iyr_lmPeriod = (length(datnew$years) - (confnew$lmPeriod)) : (length(datnew$years)-1)
length(confnew$iyr_lmPeriod)

parnew <- defpar(datnew,confnew)

newfit2 <- ccam.fit(datnew,confnew,parnew,debug=T, rm.unidentified = F, run=T)


out<- newfit2$sdrep$value
ssn_trend<- data.frame(year=datnew$years[confnew$iyr_lmPeriod+1],
pred_ssn=exp(out[grep(names(out), pattern="log_ssn_pred")]),
ssn=exp(out[grep(names(out), pattern="logssn")])[confnew$iyr_lmPeriod+1])


 ggplot()+geom_point(data=ssn_trend,aes(x=year, y=ssn)) +
     geom_line(data=ssn_trend,aes(x=year, y=pred_ssn), col="red") +
    ylab("Spawning stock abundance") +
    geom_text(aes(x=2020, y=max(ssn_trend$ssn),
                  label=paste("% change =",as.numeric(round(out[grep(names(out), pattern="pchange")])))))
ggsave("example/2024/DeclineSSN2010_2024.png", width=5, height=3, dpi=300, units="in")

##########END OF FULL READ########

load(paste0(dir, "dat.Rdata"))
load(paste0(dir, "conf.Rdata"))
conf$do_lmPeriod= 0
conf$lmPeriod= 15
conf$iyr_lmPeriod = (length(dat$years) - (conf$lmPeriod)) : (length(dat$years)-1)
length(confnew$iyr_lmPeriod)
par <- defpar(dat,conf)

oldfit <- ccam.fit(data=dat,conf,par,debug=T, rm.unidentified = T) # use rm.unidentified to remove unused parameters from the linear regression.
oldfit

########return results
library(tidyverse)

newfit2$sdrep$sd[grep(x=names(newfit2$sdrep$value), pattern="pchange")]
newfit2$sdrep$sd[grep(x=names(newfit2$sdrep$value), pattern="pchange")]

checkfit<- full_join(ssbtable(oldfit) |>  as_tibble() |>  dplyr::mutate(type="oldfit"),
ssbtable(newfit2) |>  as_tibble() |>  mutate(type="newfit"))
ggplot(checkfit, aes(x=year, y=Estimate, col=type)) +geom_line() + ylab("SSB")
ggsave("example/2024/SSBchangewithlm2010_2024.png", width=5, height=3, dpi=300, units="in")

pt<- partable(c(oldfit,newfit2)) |> as_tibble(.name_repair = "unique") |>  dplyr::select(par...1, `par...2`, fit) |>
    pivot_wider(names_from=fit, values_from=par...2) |>
    mutate(newfit=round(`1`,3), oldfit=round(`2`,3)) |>  dplyr::select(-`1`, -`2`)
write.table(pt, "example/2024/changeinparameterswithlm.txt", sep="\t", dec=".", row.names=F)



