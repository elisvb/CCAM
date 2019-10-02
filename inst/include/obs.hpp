template <class Type>
matrix<Type> setupVarCovMatrix(int minAge, int maxAge, int minAgeFleet, int maxAgeFleet, vector<int> rhoMap, vector<Type> rhoVec, vector<int> sdMap, vector<Type> sdVec){

  using CppAD::abs;
  int dim = maxAgeFleet-minAgeFleet+1;
  int offset = minAgeFleet-minAge;
  matrix<Type> ret(dim,dim);
  ret.setZero();

  Type rho0 = Type(0.5);
  vector<Type> xvec(dim);
  xvec(0)=Type(0);
  int maxrm=-1;
  if(rhoVec.size()>0){
    for(int i=1; i<xvec.size(); i++) { 
      if(rhoMap(i-1+offset)>=0)
	xvec(i) = xvec(i-1)+rhoVec(rhoMap(i-1+offset)); 
      if(rhoMap(i-1+offset)>maxrm) maxrm=rhoMap(i-1+offset);
    } 
  }
   
  for(int i=0; i<dim; i++)
    for(int j=0; j<dim; j++){
      if(i!=j && maxrm>=0){	
	Type dist = abs(xvec(i)-xvec(j));
     	ret(i,j)=pow( rho0,dist)*sdVec( sdMap(i+offset) )*sdVec( sdMap(j+offset));
      } else if(i==j) ret(i,j) = sdVec( sdMap(i+offset) )*sdVec( sdMap(j+offset));
    }
  return ret;
}

template <class Type> 
density::UNSTRUCTURED_CORR_t<Type> getCorrObj(vector<Type> params){
  density::UNSTRUCTURED_CORR_t<Type> ret(params);
  return ret;
}

template <class Type>
vector<Type> addLogratio(vector<Type> logx){
  int n = logx.size();
  vector<Type> res(n-1);
  for(int i = 0; i < res.size(); ++i)
    res(i) = logx(i) - logx(n-1);
  return res;//log(x.head(x.size()-1)/x.tail(1));
}

template<class Type>
vector<Type> multLogratio(vector<Type> logx){
  vector<Type> res(logx.size()-1);
  for(int i = 0; i < res.size(); ++i)
    res(i) = logx(i)-log(Type(1.0)-exp(logExpSum(logx.head(i+1))));
  return res;
}

template <class Type>
Type log2expsum(vector<Type> x){
  return exp(x).sum();
}

template<class Type>
Type logExpSum(vector<Type> x){
  Type m = max(x);
  return m + log(exp(x-m).sum());
}

template <class Type>
vector<Type> log2proportion(vector<Type> x){
  return exp(x) / log2expsum(x);
}


template<class Type>
matrix<Type> buildJac(vector<Type> x, vector<Type> w){
  matrix<Type> res(x.size(),x.size()); 
  Type xs = x.sum();
  Type xsp = pow(xs,Type(2));
  for(int i = 0; i < res.rows(); ++i){
    for(int j = 0; j < res.cols(); ++j){
      if(i == j){
	res(i,j) = Type(1.0)/xs-x(i)/xsp;
      }else{
	res(i,j) = -x(i)/xsp;
      }
    }
  }
  for(int j = 0; j < res.cols(); ++j){
    res(res.rows()-1,j) = w(j);
  }
  return res;
}


template <class Type>
Type jacobianDet(vector<Type> x){
  vector<Type> w(x.size());
  w.fill(Type(1.0));
  return buildJac(x,w).determinant();
}

template <class Type>
Type jacobianDet(vector<Type> x,vector<Type> w){
  return buildJac(x,w).determinant();
}

template <class Type>
Type nllObs(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logN, vector<Type> &logFy,
     data_indicator<vector<Type>,Type> &keep, objective_function<Type> *of){

  if(conf.debug==1) std::cout << "-- start"  << std::endl;
  using CppAD::abs;
  Type nll=0;

  if(conf.debug==1) std::cout << "-- selectivity"  << std::endl;
  array<Type> Sel = SelFun(conf,par.logitSel); 

  if(conf.debug==1)   std::cout << "-- F"  << std::endl;
  array<Type> logF = FFun(conf, par, logFy, 1);
  vector<Type> fbar = fbarFun(conf, logF);
  vector<Type> logfbar = log(fbar);
  if(conf.debug==1) print(fbar);

  if(conf.debug==1)   std::cout << "-- ssb stuff"  << std::endl;
  vector<Type> ssb = ssbFun(dat, conf, logN, logF, 0);
  vector<Type> logssb = log(ssb);
  vector<Type> ssb0 = ssbFun(dat, conf, logN, logF, 1);
  vector<Type> logssb0 = log(ssb0);
  vector<Type> tsb = tsbFun(dat, conf, logN);
  vector<Type> logtsb = log(tsb);
  vector<Type> fsb = fsbFun(dat, conf, logN, logF);
  vector<Type> logfsb = log(fsb);
  if(conf.debug==1) print(ssb);
  
  if(conf.debug==1)   std::cout << "-- catch/land"  << std::endl;
  array<Type> catNr = catchNrFun(dat, conf, logN, logF);
  vector<Type> cat = catchFun(dat, conf, catNr);
  vector<Type> logCatch = log(cat);
  if(conf.debug==1) print(cat);

  if(conf.debug==1)   std::cout << "-- land"  << std::endl;
  vector<Type> land = landFun(dat, conf, logN, logF);
  vector<Type> logLand = log(land);

  if(conf.debug==1)   std::cout << "-- var catch/land"  << std::endl;
  vector<Type> varLogCatch = varLogCatchFun(dat, conf, logN, logF, par);
  vector<Type> varLogLand = varLogLandFun(dat, conf, logN, logF, par);

  if(conf.debug==1)   std::cout << "-- rec"  << std::endl;
  vector<Type> R = rFun(logN);
  vector<Type> logR = log(R);  

  if(conf.debug==1)   std::cout << "-- exploit"  << std::endl;
  vector<Type> exploit = exploitFun(cat, ssb0);

  if(conf.debug==1)   std::cout << "-- predicted observations"  << std::endl;
  vector<Type> predObs = predObsFun(dat, conf, par, logN, logF, logssb, logfsb, logCatch, catNr, logLand);

  // setup obs likelihoods
   if(conf.debug==1)   std::cout << "-- obs likelihood stuff"  << std::endl;
  vector< density::MVNORM_t<Type> >  nllVec(dat.noFleets);
  vector< density::UNSTRUCTURED_CORR_t<Type> > neg_log_densityObsUnstruc(dat.noFleets);
  vector< vector<Type> > obsCovScaleVec(dat.noFleets);
  vector<Type> varLogObs=exp(par.logSdLogObs*Type(2.0));
  vector<Type> IRARdist(par.transfIRARdist.size()); //[ d_1, d_2, ...,d_N-1 ]
  if(par.transfIRARdist.size()>0) IRARdist=exp(par.transfIRARdist);
  vector< vector<Type> > sigmaObsParVec(dat.noFleets);
  int aidx;
  vector< matrix<Type> > obsCov(dat.noFleets); // for reporting
   if(conf.debug==1)   std::cout << ". "  << std::endl;
  vector<Type> recapturePhi(par.logitRecapturePhi.size());
  vector<Type> recapturePhiVec(dat.nobs);
  if(par.logitRecapturePhi.size()>0){
    recapturePhi=invlogit(par.logitRecapturePhi);
    for(int j=0; j<dat.nobs; ++j){
      if(!isNAINT(dat.aux(j,7))){
        recapturePhiVec(j)=recapturePhi(dat.aux(j,7)-1);
      }
    }
  }
   if(conf.debug==1)   std::cout << ". "  << std::endl;
  int nfleet = conf.maxAge-conf.minAge+1;
  int dn=nfleet*(nfleet-1)/2;
  int from=-dn, to=-1; 
  for(int f=0; f<dat.noFleets; f++){
    if(conf.obsCorStruct(f)!=2) continue; // skip if not US 
    nfleet = dat.maxAgePerFleet(f)-dat.minAgePerFleet(f)+1;
    if(conf.obsLikelihoodFlag(f) == 1) nfleet-=1; // ALN has dim-1
    dn = nfleet*(nfleet-1)/2;
    from=to+1;
    to=to+dn;
    vector<Type> tmp(dn);
    for(int i=from; i<=to; i++) tmp(i-from) = par.sigmaObsParUS(i);
    sigmaObsParVec(f) = tmp; 
  }
   if(conf.debug==1)   std::cout << ". "  << std::endl;
  for(int f=0; f<dat.noFleets; ++f){
  	if(conf.debug==1) std::cout << "- fleet: "  << f <<std::endl;
    if(!((dat.fleetTypes(f)==5)||(dat.fleetTypes(f)==3)||(dat.fleetTypes(f)==6))){  //catch- & survey-at-age in numbers
      int thisdim=dat.maxAgePerFleet(f)-dat.minAgePerFleet(f)+1;
      if(conf.obsLikelihoodFlag(f) == 1) thisdim-=1; // ALN has dim-1
      matrix<Type> cov(thisdim,thisdim);
      cov.setZero();
      if(conf.obsCorStruct(f)==0){//ID (independent)  
        for(int i=0; i<thisdim; ++i){
          aidx = i+dat.minAgePerFleet(f)-conf.minAge;
  	      cov(i,i)=varLogObs(conf.keyVarObs(f,aidx));
        }
      } else if(conf.obsCorStruct(f)==1){//(AR) irregular lattice AR
        cov = setupVarCovMatrix(conf.minAge, conf.maxAge, dat.minAgePerFleet(f), dat.maxAgePerFleet(f), conf.keyCorObs.transpose().col(f), IRARdist, conf.keyVarObs.transpose().col(f) , exp(par.logSdLogObs) );
	if(conf.obsLikelihoodFlag(f) == 1){ // ALN has dim-1
	  cov.conservativeResize(thisdim,thisdim); // resize, keep contents but drop last row/col
	}

      } else if(conf.obsCorStruct(f)==2){//(US) unstructured
        neg_log_densityObsUnstruc(f) = getCorrObj(sigmaObsParVec(f));  
        matrix<Type> tmp = neg_log_densityObsUnstruc(f).cov();
  
        tmp.setZero();
        int offset = dat.minAgePerFleet(f)-conf.minAge;
        obsCovScaleVec(f).resize(tmp.rows());
        for(int i=0; i<tmp.rows(); i++) {
	  tmp(i,i) = sqrt(varLogObs(conf.keyVarObs(f,i+offset)));
	  obsCovScaleVec(f)(i) = tmp(i,i);
        }
        cov  = tmp*matrix<Type>(neg_log_densityObsUnstruc(f).cov()*tmp);
      } else { error("Unkown obsCorStruct code"); }
        nllVec(f).setSigma(cov);
        obsCov(f) = cov;
    }else{
      matrix<Type> dummy(1,1);
      dummy(0,0) = R_NaReal;
      obsCov(f) = dummy;
    }
  }

  //eval likelihood 
  if(conf.debug==1)   std::cout << "!! eval likelihood"  << std::endl;
  Type sd,ZU,ZL,res;
  for(int y=0;y<dat.noYears;y++){
    int totalParKey = 0;
    for(int f=0;f<dat.noFleets;f++){
      if(isNAINT(dat.idx1(f,y))) continue;
      if(conf.debug==1) std::cout << "- year: "  << y << " / fleet:"<< f <<std::endl;
      int idxfrom=dat.idx1(f,y);
      int idxlength=dat.idx2(f,y)-dat.idx1(f,y)+1;

      switch(dat.fleetTypes(f)){
        case 0: case 2: 
        {
        vector<Type> currentVar=nllVec(f).cov().diagonal();
        vector<Type> sqrtW(currentVar.size());
          switch(conf.obsLikelihoodFlag(f)){
            case 0: // (LN) log-Normal distribution
                    
                    for(int idxV=0; idxV<currentVar.size(); ++idxV){
                      if(isNA(dat.weight(idxfrom+idxV))){
                        sqrtW(idxV)=Type(1.0);
                      }else{
                        if(conf.fixVarToWeight==1){ 
                          sqrtW(idxV)=sqrt(dat.weight(idxfrom+idxV)/currentVar(idxV));
                        }else{
                          sqrtW(idxV)=sqrt(Type(1)/dat.weight(idxfrom+idxV));
                        }
                      }
                    }

              nll += nllVec(f)((dat.logobs.block(idxfrom,0,idxlength,1)-predObs.segment(idxfrom,idxlength))/sqrtW,keep.segment(idxfrom,idxlength));
                    nll += (log(sqrtW)*keep.segment(idxfrom,idxlength)).sum();
              SIMULATE_F(of){
                dat.logobs.block(idxfrom,0,idxlength,1) = predObs.segment(idxfrom,idxlength) + (nllVec(f).simulate()*sqrtW);
              }
              break;
            case 1: // (ALN) Additive logistic-normal proportions + log-normal total numbers
              nll +=  nllVec(f)(addLogratio((vector<Type>)dat.logobs.block(idxfrom,0,idxlength,1))-addLogratio((vector<Type>)predObs.segment(idxfrom,idxlength)));
              nll += log(log2proportion((vector<Type>)dat.logobs.block(idxfrom,0,idxlength,1))).sum();
              nll -= dnorm(log(log2expsum((vector<Type>)dat.logobs.block(idxfrom,0,idxlength,1))),
                         log(log2expsum((vector<Type>)predObs.segment(idxfrom,idxlength))),
                       exp(par.logSdLogTotalObs(totalParKey++)),true);
              nll += log(log2expsum((vector<Type>)dat.logobs.block(idxfrom,0,idxlength,1)));
              nll -= log(abs(jacobianDet((vector<Type>)dat.logobs.block(idxfrom,0,idxlength,1).exp())));
                    nll -= dat.logobs.block(idxfrom,0,idxlength,1).sum();
              SIMULATE_F(of){
                vector<Type> logProb(idxlength);
                logProb.setZero();
                logProb.segment(0,idxlength-1) = addLogratio(((vector<Type>)predObs.segment(idxfrom,idxlength))) + nllVec(f).simulate();
                Type logDenom = logExpSum(logProb);
                logProb -= logDenom;
                Type logTotal = rnorm(log(log2expsum((vector<Type>)predObs.segment(idxfrom,idxlength))),
                    exp(par.logSdLogTotalObs(totalParKey++)));
                dat.logobs.block(idxfrom,0,idxlength,1) = logProb + logTotal; 
              }
              break;
            default:
              error("obsLikelihoodFlag for fleetTypes 1 and 2 should be either LN or ALN");
              break;
            }
        break;
        }
        case 3:
        {
          switch(conf.obsLikelihoodFlag(f)){
            case 0:      // lognormal distributions
              if(conf.keyBiomassTreat(f)==3){ 
                  sd = sqrt(varLogCatch(y));   // if catch not split, based on logSdLogObs first row
              }else{
                if(conf.keyBiomassTreat(f)==4){
                  sd = sqrt(varLogLand(y));
                }else{
                  sd = exp(par.logSdLogObs(conf.keyVarObs(f,0)));
                }
              }
 			if(conf.debug==1) std::cout << "pred "  <<  predObs(idxfrom) <<  " /obs  "  <<  dat.logobs(idxfrom,0) << std::endl;
              nll += -keep(idxfrom)*dnorm(dat.logobs(idxfrom,0),predObs(idxfrom),sd,true);
              SIMULATE_F(of){
               dat.logobs(idxfrom,0) = rnorm(predObs(idxfrom),sd);
              }
             break;
            case 2:  // censored   
              if(conf.debug==1) std::cout << "pred "  <<  predObs(idxfrom) <<  " /obs1 "  <<  dat.logobs(idxfrom,0)  << " /obs2 "  <<  dat.logobs(idxfrom,1) << std::endl;
                ZU = (dat.logobs(idxfrom,1) - predObs(idxfrom))/0.01;
                ZL = (dat.logobs(idxfrom,0) - predObs(idxfrom))/0.01;
                nll -= keep(idxfrom)*log(pnorm(ZU) - pnorm(ZL));
				  //SIMULATE_F(of){
				   // Not figured out yet
				   //dat.logobs(idxfrom,0)=
				   //dat.logobs(idxfrom,1)=
				  //}
             break;
            case 3:      // robuust distributions (mixture lognormal and t(df1))
              if(conf.keyBiomassTreat(f)==3){ 
                  sd = sqrt(varLogCatch(y));   // if catch not split, based on logSdLogObs first row
              }else{
                if(conf.keyBiomassTreat(f)==4){
                  sd = sqrt(varLogLand(y));
                }else{
                  sd = exp(par.logSdLogObs(conf.keyVarObs(f,0)));
                }
              }
              res =(dat.logobs(idxfrom,0)-predObs(idxfrom))/sd; 
              nll += -keep(idxfrom)*(Type(0.95)*dnorm(dat.logobs(idxfrom,0),predObs(idxfrom),sd,true))+(Type(0.05)*dt(res,Type(1),true));
              SIMULATE_F(of){
               dat.logobs(idxfrom,0) = rnorm(predObs(idxfrom),sd); // too change
              }
             break;
            default:
              error("obsLikelihoodFlag for fleetType 3 should be either LN, CE or RO");
             break;
          }   
        break;
        }
        case 5:
        {
        for(int i=idxfrom; i<(idxfrom+idxlength); ++i){
          nll += -keep(i)*dnbinom(dat.logobs(i,0),predObs(i)*recapturePhiVec(i)/(Type(1.0)-recapturePhiVec(i)),recapturePhiVec(i),true);
            SIMULATE_F(of){
             dat.logobs(i,0) = rnbinom(predObs(i)*recapturePhiVec(i)/(Type(1.0)-recapturePhiVec(i)),recapturePhiVec(i));
            }
        }
        break;
        }
        case 6:
        {
           for(int i=idxfrom; i<(idxfrom+idxlength); ++i){
             sd = exp(par.logSdLogObs(conf.keyVarObs(f,dat.aux(i,2)-conf.minAge)));
             nll += -keep(i)*dnorm(dat.logobs(i,0),predObs(i),sd,true);
              SIMULATE_F(of){
                dat.logobs(i,0) = rnorm(predObs(i),sd);
              }
           }
        break;
        }
        default:
        {
          error("FleetType should be 0,2,3,5 or 6");        
        break;
        }
      }
      if(conf.debug==1) std::cout << "nll: "<<  nll << std::endl;
    }    
  }

  if(conf.debug==1) std::cout << "The end"  << std::endl;

  SIMULATE_F(of) {
    REPORT_F(logFy,of);
    REPORT_F(logN,of);
    array<Type> logobs=dat.logobs; 
    REPORT_F(logobs,of);
  }

  REPORT_F(obsCov,of);
  REPORT_F(predObs,of);
  REPORT_F(logN,of);
  REPORT_F(logFy,of);
  ADREPORT_F(logssb,of);
  ADREPORT_F(logssb0,of);
  ADREPORT_F(Sel,of);
  ADREPORT_F(logfbar,of);
  ADREPORT_F(logCatch,of);
  ADREPORT_F(logLand,of);
  ADREPORT_F(logtsb,of);
  ADREPORT_F(logR,of);
  ADREPORT_F(exploit,of);

  int timeSteps=logF.dim[1];
  
  vector<Type> lastLogN = logN.col(timeSteps-1);
  ADREPORT_F(lastLogN,of);
  vector<Type> lastLogF = logF.col(timeSteps-1);
  ADREPORT_F(lastLogF,of);  

  vector<Type> beforeLastLogN = logN.col(timeSteps-2);
  ADREPORT_F(beforeLastLogN,of);
  vector<Type> beforeLastLogF = logF.col(timeSteps-2);
  ADREPORT_F(beforeLastLogF,of);  

  ADREPORT_F(par.logFpar,of); 
  ADREPORT_F(par.logQpow,of);  
  ADREPORT_F(par.logSdLogFsta,of);
  ADREPORT_F(par.logSdLogN,of); 
  ADREPORT_F(par.logSdLogObs,of);
  ADREPORT_F(par.logSdLogTotalObs,of);
  ADREPORT_F(par.transfIRARdist,of);
  ADREPORT_F(par.sigmaObsParUS,of);
  ADREPORT_F(par.rec_loga,of);
  ADREPORT_F(par.rec_logb,of); 
  ADREPORT_F(par.rec_e,of);
  ADREPORT_F(par.logScale,of);
  ADREPORT_F(par.logitReleaseSurvival,of);
  ADREPORT_F(par.logitRecapturePhi,of);
  ADREPORT_F(par.logitSel,of); 

  return nll;
}
