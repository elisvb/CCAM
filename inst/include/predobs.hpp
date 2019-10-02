template <class Type>
array<Type> crltransform(array<Type> x){   // ADDED BY EVB (continuation ratio logit transformation) ages in cols, years in rows
  matrix<Type> xt=x.matrix().transpose();
  int na=xt.rows();
  int ny=xt.cols();
  matrix<Type> xprop(na,ny);
  for(int j = 0 ; j < ny; j++){
    xprop.col(j) = xt.col(j)/xt.col(j).sum(); 
  }

  Type total;
  matrix<Type> xprop_cond(na-1,ny);  
  for(int j = 0 ;j <ny; j++){
   xprop_cond(0,j) = xprop(0,j); 
   for(int i = 1 ;i < na-1; i++){
     total=0.0;
     for(int k = 0 ;k < i; k++){total += xprop(k,j);}  
     xprop_cond(i,j) = xprop(i,j)/(Type(1) - total);
   }
 }

 array<Type> xcrl(ny,na-1);  
 for(int j = 0 ;j <ny; j++){
  for(int i = 0 ; i < na-1; i++){
    xcrl(j,i) = log(xprop_cond(i,j)/(Type(1) - xprop_cond(i,j)));
  }
 }  
 return xcrl;
}


template <class Type> // ADDED BY EVB (total egg production)
Type TEP(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, int f, int y){   
  int stateDimN=logN.dim[0];
  Type tep=0;
  Type logtep;
  Type ts = dat.sampleTimes(f-1);
  Type Zprop;
   for(int a=0; a<stateDimN; ++a){
    Zprop=dat.natMor(y,a);
    if(conf.keySel(y,a)<0) Zprop+=exp(logF(a,y));
    tep += exp(logN(a,y))*dat.propMat(y,a)*exp(-Zprop*ts)*dat.propFemale(y,a)*dat.fec(y,a);
   }
  logtep = log(tep);
 return logtep;
}


template <class Type>
vector<Type> predObsFun(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logN, array<Type> &logF, vector<Type> &logssb, vector<Type> &logfsb, vector<Type> &logCatch, array<Type> &catNr,vector<Type> &logLand){
  vector<Type> pred(dat.nobs);
  pred.setZero();


  // for tagging
  vector<Type> releaseSurvival(par.logitReleaseSurvival.size());
  vector<Type> releaseSurvivalVec(dat.nobs);
  if(par.logitReleaseSurvival.size()>0){
    releaseSurvival=invlogit(par.logitReleaseSurvival);
    for(int j=0; j<dat.nobs; ++j){
      if(!isNAINT(dat.aux(j,7))){
        releaseSurvivalVec(j)=releaseSurvival(dat.aux(j,7)-1);
      }
    }
  }
  
  // for crl transform
  array<Type> cprop = crltransform(catNr);

  // predict observations
  int f, ft, a, y, yy, scaleIdx;  // a is no longer just ages, but an attribute (e.g. age or length) 
  int minYear=dat.aux(0,0);
  Type zz=Type(0);
  Type logtep;
  for(int i=0;i<dat.nobs;i++){
    y=dat.aux(i,0)-minYear;
    f=dat.aux(i,1);
    ft=dat.fleetTypes(f-1);
    a=dat.aux(i,2)-conf.minAge;
    if(ft==3){a=0;}
    if(ft<3){ 
      zz = dat.natMor(y,a);
      if(conf.keySel(y,a)>(-1)){
        zz+=exp(logF(a,y));
      }
    }    
    switch(ft){
      case 0: // Catch-at-age in numbers
        pred(i)=logN(a,y)-log(zz)+log(1-exp(-zz));
        if(conf.keySel(y,a)>(-1)){
          pred(i)+=logF(a,y);
        }
        scaleIdx=-1;
        yy=dat.aux(i,0);
        for(int j=0; j<conf.noScaledYears; ++j){
          if(yy==conf.keyScaledYears(j)){
            scaleIdx=conf.keyParScaledYA(j,a);
            if(scaleIdx>=0){
              pred(i)-=par.logScale(scaleIdx);
            }
            break;
          }
        }
      break;
  
      case 1:
  	error("Unknown fleet code");
        return(0);
      break;
  
      case 2: // Survey-at-age in numbers
        pred(i)=logN(a,y)-zz*dat.sampleTimes(f-1);
        if(conf.keyQpow(f-1,a)>(-1)){
          pred(i)*=exp(par.logQpow(conf.keyQpow(f-1,a))); 
        }
        if(conf.keyLogFpar(f-1,a)>(-1)){
          pred(i)+=par.logFpar(conf.keyLogFpar(f-1,a));
        }
        
      break;
  
      case 3:// annual data (whether it be I, C, land, TEP, etc.)
        if(conf.keyBiomassTreat(f-1)==0){
          pred(i) = logssb(y)+par.logFpar(conf.keyLogFpar(f-1,a)); 
        }
        if(conf.keyBiomassTreat(f-1)==1){
          pred(i) = logCatch(y)+par.logFpar(conf.keyLogFpar(f-1,a));
        }
        if(conf.keyBiomassTreat(f-1)==2){
          pred(i) = logfsb(y)+par.logFpar(conf.keyLogFpar(f-1,a));
        }
        if(conf.keyBiomassTreat(f-1)==3){ 
          pred(i) = logCatch(y);
        }
        if(conf.keyBiomassTreat(f-1)==4){
          pred(i) = logLand(y);
        }
        if(conf.keyBiomassTreat(f-1)==5){
          logtep = TEP(dat, conf, logN, logF, f, y);
          pred(i) = logtep+par.logFpar(conf.keyLogFpar(f-1,a));
        }
	break;
  
      case 4: 
          error("Unknown fleet code");
        return 0;
      break;
  
      case 5:// tag data  
        if((a+conf.minAge)>conf.maxAge){a=conf.maxAge-conf.minAge;} 
	      pred(i)=exp(log(dat.aux(i,6))+log(dat.aux(i,5))-logN(a,y)-log(1000))*releaseSurvivalVec(i);
      break;
  
      case 6: // Catch-at-age in proportions (crl transformed)
        pred(i) = cprop(y,a);
      break;
  
      case 7: 
          error("Unknown fleet code");
        return 0;
      break;
  
      default:
  	error("Unknown fleet code");
        return 0 ;
      break;
    }    
  }
  return pred;
}
