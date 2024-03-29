##' Function to supress incomplete final line warning
##' @param ... arguments
##' @importFrom utils read.table
##' @details ...
read.table.nowarn<-function(...){
  tryCatch.W.E <- function(expr)
  {
    W <- NULL
    w.handler <- function(w){ # warning handler
      if(!grepl('incomplete final line',w))W<<-w
      invokeRestart("muffleWarning")
    }
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),warning = w.handler),warning = W)
  }
  lis<-tryCatch.W.E(read.table(...))
  if(!is.null(lis$warning))warning(lis$warning)
  lis$value
}

##' Function to test if x is ...
##' @param x number
##' @param tol precision
##' @details ...
is.whole.positive.number <- function(x, tol = .Machine$double.eps^0.5){
  (abs(x - round(x)) < tol)&(x>=0)
}

##' Function to read ices survey format
##' @param filen the file
##' @details ...
read.surveys<-function(filen){
  # Function to read ices survey file
  lin<-readLines(filen,warn=FALSE)[-c(1:2)]
  empty<-which(lapply(lapply(strsplit(lin, split='[[:space:]]+'),
               paste, collapse=''), nchar)==0)
  if(length(empty)>0){
    lin<-lin[-empty]
  }
  lin<-sub("^\\s+", "",lin)
  idx1<-grep('^[A-Z#]', lin, ignore.case=TRUE)
  idx2<-c(idx1[-1]-1,length(lin))
  names<-lin[idx1]
  years<-matrix(as.numeric(unlist(strsplit(lin[idx1+1], '[[:space:]]+'))), ncol=2, byrow=TRUE)
  times<-matrix(as.numeric(unlist(strsplit(lin[idx1+2], '[[:space:]]+'))), ncol=4, byrow=TRUE)[,3:4,drop=FALSE]
  ages<-matrix(as.numeric(unlist(lapply(strsplit(lin[idx1+3], '[[:space:]]+'), function(x)x[1:2]))), ncol=2, byrow=TRUE)
  for(i in 1:length(names)){
    # Check years
    if(!is.whole.positive.number(years[i,1])){
      stop(paste("In file",filen, ": Minimum year is expected to be a positive integer number for fleet number",i))
    }
    if(!is.whole.positive.number(years[i,2])){
      stop(paste("In file",filen, ": Maximum year is expected to be a positive integer number for fleet number",i))
    }
    if(years[i,1]>years[i,2]){
      stop(paste("In file",filen, ": Maximum year is expected to be greater than minimum year for fleet number",i))
    }
    # Check ages
    ##if(!is.whole.positive.number(ages[i,1])){
    ##  stop(paste("In file",filen, ": Minimum age is expected to be a positive integer number for fleet number",i))
    ##}
    ##if(!is.whole.positive.number(ages[i,2])){
    ##  stop(paste("In file",filen, ": Maximum age is expected to be a positive integer number for fleet number",i))
    ##}
    if(ages[i,1]>ages[i,2]){
      stop(paste("In file",filen, ": Maximum age is expected to be greater than minimum age for fleet number",i))
    }
    # Check times
    if((times[i,1]<0)|(times[i,1]>1)){
      stop(paste("In file",filen, ": Minimum survey time is expected to be within [0,1] for fleet number",i))
    }
    if((times[i,2]<0)|(times[i,2]>1)){
      stop(paste("In file",filen, ": Maximum survey time is expected to be within [0,1] for fleet number",i))
    }
    if(times[i,2]<times[i,1]){
      stop(paste("In file",filen, ": Maximum survey time is expected greater than minimum survey time for fleet number",i))
    }
  }

  as.num <- function(x, na.strings = "NA") {
    stopifnot(is.character(x))
    na = x %in% na.strings
    x[na] = 0
    x = as.numeric(x)
    x[na] = NA_real_
    x
  }

  onemat<-function(i){
    lin.local<-gsub('^[[:blank:]]*','',lin[(idx1[i]+4):idx2[i]])
    nr<-idx2[i]-idx1[i]-3
    ret<-matrix(as.num(unlist((strsplit(lin.local,'[[:space:]]+')))),nrow=nr, byrow=TRUE)[,,drop=FALSE]   #[,1:(2+ages[i,2]-ages[i,1]),drop=FALSE]
    if(nrow(ret)!=(years[i,2]-years[i,1]+1)){
      stop(paste("In file",filen, ": Year range specified does not match number of rows for survey fleet number",i))
    }
    if((ncol(ret)-1)<(ages[i,2]-ages[i,1]+1)){
      stop(paste("In file",filen, ": Fewer columns than indicated by age range for survey fleet number",i))
    }
    if(!is.numeric(ret)){
      stop(paste("In file",filen, ": Non numeric data values detected for survey fleet number",i))
    }
    ret<-as.matrix(ret[,-1]/ret[,1])
    rownames(ret)<-years[i,1]:years[i,2]
    ret<-ret[,1:length(ages[i,1]:ages[i,2]),drop=FALSE]
    colnames(ret)<-ages[i,1]:ages[i,2]
    attr(ret,'time')<-times[i,]
    ret[ret<0]<-NA
    ret
  }
  obs<-lapply(1:length(names),onemat)
  names(obs)<-names
  obs
}

##' Function to read ICES/CEFAS data files and validate if input makes sense
##' @param filen The filename
##' @importFrom TMB MakeADFun sdreport
##' @details
##' First two lines are ignored and can be used for comments.
##' Can read formats 1 full, 2 row, 3 scalar, and 5 column
##'
##' Tests:
##' Formatcode is valid, years and ages are pos. integers
##' minimum <= maximum for years and ages
##' number of rows and coulmns match year and age ranges
##' data contains only numbers.
##'
##' Returns: A validated data matrix.
##' @export
read.ices<-function(filen){
  if(grepl("^[0-9]", scan(filen, skip=2, n=1, quiet=TRUE, what=""))){ # is not a survey file

    head<-scan(filen, skip=2, n=5, quiet=TRUE)
    minY<-head[1]
    maxY<-head[2]
    minA<-head[3]
    maxA<-head[4]
    datatype<-head[5]

    if(!is.whole.positive.number(minY)){
      stop(paste("In file",filen, ": Minimum year is expected to be a positive integer number"))
    }
    if(!is.whole.positive.number(maxY)){
      stop(paste("In file",filen, ": Maximum year is expected to be a positive integer number"))
    }
    if(!is.whole.positive.number(minA)){
      stop(paste("In file",filen, ": Minimum age is expected to be a positive integer number"))
    }
    if(!is.whole.positive.number(maxA)){
      stop(paste("In file",filen, ": Maximum age is expected to be a positive integer number"))
    }

    if(!(datatype%in%c(1,2,3,5,6,7))){
      stop(paste("In file",filen, ": Datatype code is expected to be one of the numbers 1, 2, 3, 5, 6 or 7"))
    }

    if(minY>maxY){
      stop(paste("In file",filen, ": Minimum year is expected to be less than maximum year"))
    }
    if(minA>maxA){
      stop(paste("In file",filen, ": Minimum age is expected to be less than maximum age"))
    }

    C<-as.matrix(read.table.nowarn(filen, skip=5, header=FALSE))

    if(datatype==1){
      if((maxY-minY+1)!=nrow(C)){
        stop(paste("In file",filen, ": Number of rows does not match the year range given"))
      }
      if((maxA-minA+1)>ncol(C)){
        stop(paste("In file",filen, ": Fewer columns than the age range given"))
      }
    }

    if(datatype==2){
      if(1!=nrow(C)){
        stop(paste("In file",filen, ": For datatype 2 only one row of data is expected"))
      }
      if((maxA-minA+1)>ncol(C)){
        stop(paste("In file",filen, ": Fewer columns than the age range given"))
      }
      C<-C[rep(1,maxY-minY+1),]
    }

    if(datatype==3){
      if(1!=nrow(C)){
        stop(paste("In file",filen, ": For datatype 3 only one row of data is expected"))
      }
      if(1!=ncol(C)){
        stop(paste("In file",filen, ": For datatype 3 only one column of data is expected"))
      }
      C<-C[rep(1,maxY-minY+1),rep(1,maxA-minA+1)]
    }

    if(datatype==5){
      if((maxY-minY+1)!=nrow(C)){
        stop(paste("In file",filen, ": Number of rows does not match the year range given"))
      }
      if(1!=ncol(C)){
        stop(paste("In file",filen, ": For datatype 5 only one column of data is expected"))
      }
      C<-C[,rep(1,maxA-minA+1)]
    }

    if(datatype==6){
        if((maxY-minY+1)!=nrow(C)){
            stop(paste("In file",filen, ": Number of rows does not match the year range given"))
        }
        if(2!=ncol(C)){
            stop(paste("In file",filen, ": For datatype 6 two columns of data are expected"))
        }
    }

    if(datatype==7){
        if(is.vector(C)){
            C=cbind(C,0)
        }
        if((maxY-minY+1)!=nrow(C)){
            stop(paste("In file",filen, ": Number of rows does not match the year range given"))
        }

    }


    rownames(C)<-minY:maxY
    if(datatype %in% c(1,2,3,5)){
    C<-C[,1:length(minA:maxA)]
    colnames(C)<-minA:maxA
    }else{
        if(datatype==6){
            colnames(C)<-c('min','max')
        }else{
            colnames(C)<-1:ncol(C)
        }
    }

    if(!is.numeric(C)){
        stop(paste("In file",filen, ": Non numeric data values detected (could for instance be comma used as decimal operator)"))
    }
    return(C)
  }else{
    return(read.surveys(filen))
  }
}

##' Combine the data sources to CCAM readable object
##' @param fleets comm fleets vith effort (currently unimplemented)
##' @param surveys surveys
##' @param residual.fleet total catch minus commercial
##' @param prop.mature pm
##' @param stock.mean.weight sw
##' @param stock.start.weight sw0
##' @param catch.mean.weight cw
##' @param dis.mean.weight dw
##' @param land.mean.weight lw
##' @param natural.mortality nm
##' @param prop.f ...
##' @param prop.m ...
##' @param land.frac ...
##' @param recapture ...
##' @param prop.fem ....
##' @param fec ...
##' @param env ...
##' @importFrom stats complete.cases
##' @details ...
##' @export
setup.ccam.data <- function(fleets=NULL, surveys=NULL, residual.fleet=NULL,
                            total.catch=NULL,prop.mature=NULL, stock.mean.weight=NULL, stock.start.weight=NULL,catch.mean.weight=NULL,
                           dis.mean.weight=NULL, land.mean.weight=NULL,
                           natural.mortality=NULL, prop.f=NULL, prop.m=NULL, land.frac=NULL, recapture=NULL,
                           prop.fem=NULL,fec=NULL,env=NULL
                           ){
  # Function to write records in state-space assessment format and create
  # collected data object for future use
  fleet.idx<-0
  type<-NULL
  time<-NULL
  name<-NULL
  dat<-data.frame(year=NA,fleet=NA,age=NA,aux1=NA,aux2=NA)
  weight<-NULL
  doone<-function(m){
    year<-rownames(m)[row(m)]
    fleet.idx<<-fleet.idx+1
    fleet<-rep(fleet.idx,length(year))
    age<-as.integer(colnames(m)[col(m)])
    aux1<-(as.vector(m))
    aux2<-NA
    dat<<-rbind(dat,data.frame(year,fleet,age,aux1,aux2))
    if("weight"%in%names(attributes(m))){
      weight<<-c(weight,as.vector(attr(m,"weight")))
    }else{
      weight<<-c(weight,rep(NA,length(year)))
    }
  }
  if(!is.null(total.catch)){
      year<-rownames(total.catch)
      fleet.idx<-fleet.idx+1
      fleet<-rep(fleet.idx,length(year))
      age <- -1
      aux1<-log(total.catch[,1])
      aux2<-log(total.catch[,2])
      dat<-rbind(dat,data.frame(year,fleet,age,aux1,aux2))
      dat<-dat[!is.infinite(dat$aux1),]
      if("weight"%in%names(attributes(total.catch))){
          weight<-c(weight,as.vector(attr(total.catch,"weight")))
      }else{
          weight<-c(weight,rep(NA,length(year)))
      }
      type<-c(type,3)
      time<-c(time,0)
      name<-c(name,"Total catch")
  }
  if(!is.null(residual.fleet)){
    if(!is.null(total.catch)){
       rf=t(crl(t(residual.fleet)))
       type<-c(type,6)
       name<-c(name,"Catch-at-age proportions")
    }else{
       residual.fleet[residual.fleet==0]=1
       residual.fleet[residual.fleet<0]=NA
       rf=log(residual.fleet)
       type<-c(type,0)
       name<-c(name,"Residual catch")
    }
    doone(rf)
    time<-c(time,0)
  }
  if(!is.null(fleets)){
    if(is.data.frame(fleets)|is.matrix(fleets)){
      doone(log(fleets))
      type<-c(type,1)
      time<-c(time,0)
      name<-c(name,"Comm fleet")
    }else{
      dummy<-lapply(log(fleets),doone)
      type<-c(type,rep(1,length(fleets)))
      time<-c(time,rep(0,length(fleets)))
      name<-c(name,strtrim(gsub("\\s", "", names(dummy)), 50))
    }
  }
  if(!is.null(surveys)){
    if(is.data.frame(surveys)|is.matrix(surveys)){
      doone(log(surveys))
      thistype<-ifelse(min(as.integer(colnames(surveys)))<(-.5),3,2)
      type<-c(type,thistype)
      time<-c(time,mean(attr(surveys,'time')))
      name<-c(name,"Survey fleet")
    }else{
      dummy<-lapply(lapply(surveys,log),doone)
      type<-c(type,unlist(lapply(surveys, function(x)ifelse(min(as.integer(colnames(x)))<(-.5), 3, 2))))
      time<-c(time,unlist(lapply(surveys, function(x)mean(attr(x,'time')))))
      name<-c(name,strtrim(gsub("\\s", "", names(dummy)), 50))
    }
  }
  ny <- length(min(as.numeric(dat$year),na.rm=T):max(as.numeric(dat$year),na.rm=T))
  if(is.null(land.frac)){
    land.frac<-matrix(as.integer(1),nrow=ny, ncol=ncol(residual.fleet),dimnames=dimnames(catch.mean.weight)) # should be pure 1
  }
  if(is.null(stock.start.weight)){
    stock.start.weight<-stock.mean.weight
  }
  if(is.null(dis.mean.weight)){
    dis.mean.weight<-catch.mean.weight
  }
  if(is.null(land.mean.weight)){
    land.mean.weight<-catch.mean.weight
  }
  if(is.null(prop.f)){
    prop.f<-matrix(0,nrow=ny, ncol=ncol(residual.fleet))
  }
  if(is.null(prop.m)){
    prop.m<-matrix(0,nrow=ny, ncol=ncol(residual.fleet))
  }
  if(is.null(prop.fem)){
    prop.fem<-matrix(0.5,nrow=ny, ncol=ncol(residual.fleet))
  }
  if(is.null(fec)){
    fec<-matrix(1,nrow=ny, ncol=ncol(residual.fleet))
  }
  if(is.null(env)){
    env<-matrix(0,nrow=ny, ncol=1)
    rownames(env)=min(as.numeric(dat$year),na.rm=T):max(as.numeric(dat$year),na.rm=T)
  }
  dat<-dat[!is.na(dat$year),]

  if(!is.null(recapture)){
    tag<-data.frame(year=recapture$ReleaseY)
    fleet.idx <- fleet.idx+1
    tag$fleet <- fleet.idx
    tag$age <- recapture$ReleaseY-recapture$Yearclass
    tag$aux1 <- recapture$r
    tag <- cbind(tag, recapture[,c("RecaptureY", "Yearclass", "Nscan", "R", "Type")])
    dat[names(tag)[!names(tag)%in%names(dat)]]<-NA
    dat<-rbind(dat, tag)
    weight<-c(weight,rep(NA,nrow(tag)))
    type<-c(type,5)
    time<-c(time,0)
    name<-c(name,"Recaptures")
  }
  dat<-dat[complete.cases(dat[,1:3]),]
  dat<-type.convert(dat,as.is=TRUE)

  o<-order(dat$year,dat$fleet,dat$age)
  attr(dat,'type')<-type
  names(time)<-NULL
  attr(dat,'time')<-time
  names(name)<-NULL
  attr(dat,'name')<-name
  dat<-dat[o,]
  weight<-weight[o]
  newyear<-min(dat$year):max(dat$year)
  newfleet<-min(dat$fleet):max(dat$fleet)
  mmfun<-function(f,y, ff){idx<-which(dat$year==y & dat$fleet==f); ifelse(length(idx)==0, NA, ff(idx)-1)}
  idx1<-outer(newfleet, newyear, Vectorize(mmfun,c("f","y")), ff=min)
  idx2<-outer(newfleet, newyear, Vectorize(mmfun,c("f","y")), ff=max)
  attr(dat,'idx1')<-idx1
  attr(dat,'idx2')<-idx2
  attr(dat,"minAgePerFleet")<-as.vector(tapply(dat[,"age"], INDEX=dat[,"fleet"], FUN=min))
  attr(dat,"maxAgePerFleet")<-as.vector(tapply(dat[,"age"], INDEX=dat[,"fleet"], FUN=max))
  attr(dat,'year')<-newyear
  attr(dat,'nyear')<-max(as.numeric(dat$year))-min(as.numeric(dat$year))+1 ##length(unique(dat$year))
  cutY<-function(x)x[rownames(x)%in%newyear,,drop=FALSE]
  attr(dat,'prop.mature')<-cutY(prop.mature)
  attr(dat,'stock.mean.weight')<-cutY(stock.mean.weight)
  attr(dat,'stock.start.weight')<-cutY(stock.start.weight)
  attr(dat,'catch.mean.weight')<-cutY(catch.mean.weight)
  attr(dat,'dis.mean.weight')<-cutY(dis.mean.weight)
  attr(dat,'land.mean.weight')<-cutY(land.mean.weight)
  attr(dat,'natural.mortality')<-cutY(natural.mortality)
  attr(dat,'prop.f')<-cutY(prop.f)
  attr(dat,'prop.m')<-cutY(prop.m)
  attr(dat,'prop.fem')<-cutY(prop.fem)
  attr(dat,'fec')<-cutY(fec)
  attr(dat,'env')<-cutY(env)
  attr(dat,'land.frac')<-cutY(land.frac)
  ret<-list(
    noFleets=length(attr(dat,'type')),
    fleetTypes=as.integer(attr(dat,'type')),
    sampleTimes=attr(dat,'time'),
    noYears=attr(dat,'nyear'),
    years=attr(dat,'year'),
    minAgePerFleet=attr(dat,"minAgePerFleet"),
    maxAgePerFleet=attr(dat,"maxAgePerFleet"),
    nobs=nrow(dat),
    idx1=attr(dat,'idx1'),
    idx2=attr(dat,'idx2'),
    aux=data.matrix(dat[,-c(4:5)]),
    logobs=data.matrix(dat[,4:5]),
    weight=as.numeric(weight),
    propMat=attr(dat,'prop.mature'),
    stockMeanWeight=attr(dat,'stock.mean.weight'),
    stockStartWeight=attr(dat,'stock.start.weight'),
    catchMeanWeight=attr(dat,'catch.mean.weight'),
    natMor=attr(dat,'natural.mortality'),
    landFrac=attr(dat,'land.frac'),
    disMeanWeight=attr(dat,'dis.mean.weight'),
    landMeanWeight=attr(dat,'land.mean.weight'),
    propF=attr(dat,'prop.f'),
    propM=attr(dat,'prop.m'),
    propFemale=attr(dat,'prop.fem'),
    fec=attr(dat,'fec'),
    env=attr(dat,'env')
  )
  attr(ret,"fleetNames")<-attr(dat,"name")
  return(ret)
}


##' Read a fitted model from stockassessment.org
##' @param stockname The short-form name of a stock on stockassessment.org. This will (currently?) not work for stocks defined via the AD Model builder version of CCAM.
##' @param character.only a logical indicating whether 'stockname' can be assumed to be a character string
##' @details ...
##' @export
fitfromweb <- function(stockname, character.only=FALSE){
  if (!character.only) stockname <- as.character(substitute(stockname))
  fit<-NA
  load(url(sub("SN",stockname,"https://stockassessment.org/datadisk/stockassessment/userdirs/user3/SN/run/model.RData")))
  fit
}
