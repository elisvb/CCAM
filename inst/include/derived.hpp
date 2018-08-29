template <class Type> // ADDED BY EVB for debugging
Type print(array<Type> m){
  int rows=m.dim[0];
  int cols=m.dim[1];

    for(int x=0;x<rows;x++)
    {
        for(int y=0;y<cols;y++) 
        {
            std::cout<<m(x,y) << " ";
        }
    std::cout<<std::endl; 
    }

  return 0;
}

template <class Type> // ADDED BY EVB for debugging
Type print(vector<Type> v){
  int len=v.size();

    for(int x=0;x<len;x++)
    {
      std::cout<< v(x) << " ";
    }
  std::cout<<std::endl; 
  return 0;
}


template <class Type> // ADDED BY EVB
array<Type> SelFun(confSet &conf, vector<Type> &logitSel){
  int timeSteps=conf.keySel.dim[0];
  int stateDimN=conf.keySel.dim[1];
  int selblocks=1;
  vector<int> ix(timeSteps-1);
  ix(0)=0;

   for(int j=1;j<timeSteps;++j){
    if(conf.keySel.matrix().row(j)!=conf.keySel.matrix().row(j-1)){
      selblocks+=1;
      ix(j)=ix(j-1)+1;
    }else{
      ix(j)=ix(j-1);
    } 
   }

  array<Type> Sel(stateDimN,selblocks);  

  for(int i=0;i<stateDimN;++i){ 
   for(int j=0;j<timeSteps;++j){
    if(conf.keySel(j,i)<0){
      Sel(i,ix(j))=Type(0.0000001);
    }else{
      Sel(i,ix(j))=invlogit(logitSel(conf.keySel(j,i)));
    }
   }
  } 
  return Sel;
}


template <class Type> // Changed BY EVB
array<Type> FFun(confSet &conf, vector<Type> &logFy, vector<Type> &logitSel,int tran){
  if(conf.debug==1) std::cout << "-- f check"  << std::endl;
  int timeSteps=conf.keySel.dim[0];
  int stateDimN=conf.keySel.dim[1];
  array<Type> F(stateDimN,timeSteps);  

  for(int i=0;i<stateDimN;++i){ 
   for(int j=0;j<timeSteps;++j){
    if(conf.keySel(j,i)<0){
      F(i,j)=Type(0.0000001);
    }else{
      F(i,j)=exp(logFy(j));
      if(conf.keySel(j,i)!=conf.keySel.matrix().maxCoeff())  F(i,j)*=invlogit(logitSel(conf.keySel(j,i)));
    }
   }
  } 
  if(conf.debug==1) std::cout << "-- f check2"  << std::endl;
  if(tran==1) F=log(F);
  return F;
}


template <class Type> // changed by EVB
Type ssbi(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, int i, int s){
  int stateDimN=logN.dim[0];
  Type ssb=0;
  Type Zprop;
  for(int j=0; j<stateDimN; ++j){
    if(s==1){
      Zprop=0; // if s(tart) = 0 ssb is calcuted at another time of year then the beginning
    }else{
      Zprop=dat.natMor(i,j)*dat.propM(i,j);
        if(conf.keySel(i,j)<0) Zprop+=exp(logF(j,i))*dat.propF(i,j);
    }
    Zprop = exp(-Zprop);
    ssb += exp(logN(j,i))*Zprop*dat.propMat(i,j)*dat.stockMeanWeight(i,j);
  }
  return ssb;
}

template <class Type>
vector<Type> ssbFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, int s){
  int timeSteps=logF.dim[1];
  vector<Type> ssb(timeSteps);
  ssb.setZero();
  for(int i=0;i<timeSteps;i++){
    ssb(i)=ssbi(dat,conf,logN,logF,i,s);
  }
  return ssb;
}

template <class Type> // ADDED BY EVB
array<Type> catchNrFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF){
  int timeSteps=dat.catchMeanWeight.dim(0);
  int stateDimN=dat.catchMeanWeight.dim(1); 
  array<Type> cat(timeSteps,stateDimN);
  cat.setZero();
  for(int y=0;y<timeSteps;y++){
    for(int a=0;a<=(conf.maxAge-conf.minAge);a++){  
      Type z=dat.natMor(y,a);
      if(conf.keySel(y,a)>(-1)){
        z+=exp(logF(a,y));
        cat(y,a)=exp(logF(a,y))/z*exp(logN(a,y))*(Type(1.0)-exp(-z));
      }
    }
  }
  return cat;
}

template <class Type> // changed by EVB
vector<Type> catchFun(dataSet<Type> &dat, confSet &conf, array<Type> &catNr){
  int timeSteps=dat.catchMeanWeight.dim(0);
  vector<Type> cat(timeSteps);
  cat.setZero();
  for(int y=0;y<timeSteps;y++){
    for(int a=0;a<=(conf.maxAge-conf.minAge);a++){  
        cat(y)+=catNr(y,a)*dat.catchMeanWeight(y,a);
    }
  }
  return cat;
}

template <class Type> // changed by EVB
vector<Type> varLogCatchFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, paraSet<Type> &par){
  int len=dat.catchMeanWeight.dim(0);
  vector<Type> cat(len);
  cat.setZero();
  vector<Type> varLogCat(len);
  if(conf.keyVarObs.matrix().row(0).maxCoeff()<0) return varLogCat.setConstant(0.0001);
  varLogCat.setZero();
  Type CW=0;
  Type Ca=0;
  for(int y=0;y<len;y++){
    for(int a=0;a<=(conf.maxAge-conf.minAge);a++){
      Type z=dat.natMor(y,a);
      if(conf.keySel(y,a)>(-1)){
        z+=exp(logF(a,y));
        CW=dat.catchMeanWeight(y,a);
        Ca=exp(logF(a,y))/z*exp(logN(a,y))*(Type(1.0)-exp(-z));
        cat(y)+=Ca*CW;
        varLogCat(y)+=exp(2.0*par.logSdLogObs(conf.keyVarObs(0,a)))*CW*CW*Ca*Ca; 
      }
    }
    varLogCat(y)/=cat(y)*cat(y);
  }
  return varLogCat;
}

template <class Type>
vector<Type> landFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF){
  int len=dat.landMeanWeight.dim(0);
  vector<Type> land(len);
  land.setZero();
  Type LW=0;
  Type LF=0;
  Type LWLF=0;
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){
      Type z=dat.natMor(y,a-conf.minAge);
      if(conf.keySel(y,a-conf.minAge)>(-1)){
        z+=exp(logF(a-conf.minAge,y));
        LW=dat.landMeanWeight(y,a-conf.minAge);
        if(LW<0){LW=0;}
	      LF=dat.landFrac(y,a-conf.minAge);
        if(LF<0){LF=0;}
        LWLF=LW*LF;
        land(y)+=exp(logF(a-conf.minAge,y))/z*exp(logN(a-conf.minAge,y))*(Type(1.0)-exp(-z))*LWLF;
      }
    }
  }
  return land;
}

template <class Type>
vector<Type> varLogLandFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, paraSet<Type> &par){
  int len=dat.landMeanWeight.dim(0);
  vector<Type> land(len);
  land.setZero();
  vector<Type> varLogLand(len);
  varLogLand.setZero();
  Type LW=0;
  Type LF=0;
  Type LWLF=0;
  Type Ca=0;
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){
       if(conf.keyVarObs(0,a-conf.minAge)<0) continue; 
      Type z=dat.natMor(y,a-conf.minAge);
      if(conf.keySel(y,a-conf.minAge)>(-1)){
        z+=exp(logF(a-conf.minAge,y));
        LW=dat.landMeanWeight(y,a-conf.minAge);
        if(LW<0){LW=0;}
        LF=dat.landFrac(y,a-conf.minAge);
        if(LF<0){LF=0;}
        LWLF=LW*LF;
        Ca=exp(logF(a-conf.minAge,y))/z*exp(logN(a-conf.minAge,y))*(Type(1.0)-exp(-z));
        land(y)+=Ca*LWLF;
        varLogLand(y)+=exp(2.0*par.logSdLogObs(conf.keyVarObs(0,a-conf.minAge)))*LWLF*LWLF*Ca*Ca;
      }
    }

    if(land(y)==0) continue; 
    varLogLand(y)/=land(y)*land(y);
  }
  return varLogLand;
}

template <class Type>
vector<Type> fsbFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF){
  int len=dat.catchMeanWeight.dim(0);
  vector<Type> fsb(len);
  fsb.setZero();
  Type sumF;
  for(int y=0;y<len;y++){  // calc logfsb
    sumF=Type(0);
    for(int a=conf.minAge;a<=conf.maxAge;a++){  
      if(conf.keySel(y,a-conf.minAge)>(-1)){
        sumF+=exp(logF(a-conf.minAge,y));
      }
    }
    for(int a=conf.minAge;a<=conf.maxAge;a++){  
      Type z=dat.natMor(y,a-conf.minAge);
      if(conf.keySel(y,a-conf.minAge)>(-1)){
        z+=exp(logF(a-conf.minAge,y));
        fsb(y)+=(exp(logF(a-conf.minAge,y))/sumF)*exp(logN(a-conf.minAge,y))*exp(-Type(0.5)*z)*dat.catchMeanWeight(y,a-conf.minAge);
      }
    }
  }
  return fsb;
}

template <class Type>
vector<Type> tsbFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN){
  int timeSteps=logN.dim[1];
  vector<Type> tsb(timeSteps);
  tsb.setZero();  
  for(int y=0;y<timeSteps;y++){  
    for(int a=conf.minAge;a<=conf.maxAge;a++){  
      tsb(y)+=exp(logN(a-conf.minAge,y))*dat.stockMeanWeight(y,a-conf.minAge);
    }
  }
  return tsb;
}

template <class Type>
vector<Type> rFun(array<Type> &logN){
  int timeSteps=logN.dim[1];
  vector<Type> R(timeSteps);
  R.setZero();
  for(int y=0;y<timeSteps;y++){  
    R(y)=exp(logN(0,y));
  }
  return R;
}

template <class Type>
vector<Type> fbarFun(confSet &conf, array<Type> &logF){
  int timeSteps=logF.dim[1];
  vector<Type> fbar(timeSteps);
  fbar.setZero();
  for(int y=0;y<timeSteps;y++){  
    for(int a=conf.fbarRange(0);a<=conf.fbarRange(1);a++){  
      fbar(y)+=exp(logF(a-conf.minAge,y));
    }
    fbar(y)/=Type(conf.fbarRange(1)-conf.fbarRange(0)+1);
  }
  return fbar;
}

template <class Type>  // ADDED BY EVB
vector<Type> exploitFun(vector<Type> &cat, vector<Type> &ssb){
  int timeSteps=cat.size();
  vector<Type> exploit(timeSteps);
  for(int y=0;y<timeSteps;y++){
   exploit(y)=cat(y)/ssb(y);
  }
  return exploit;
}
