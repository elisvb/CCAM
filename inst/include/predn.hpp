
template <class Type>
vector<Type> predNFun(dataSet<Type> &dat, confSet &conf, paraSet<Type> &par, array<Type> &logN, array<Type> &logF, int i){
  int stateDimN=logN.dim[0];
  vector<Type> predN(stateDimN); 
  Type thisSSB=Type(0); 
    if(conf.stockRecruitmentModelCode==0){ // straight RW 
      predN(0)=logN(0,i-1);
    }else{
      if((i-conf.minAge)>=0){
        thisSSB=ssbi(dat,conf,logN,logF,i-conf.minAge,0);
      }else{
        thisSSB=ssbi(dat,conf,logN,logF,0,0); // use first in beginning       
      } 
      if(conf.stockRecruitmentModelCode==1){//ricker
        predN(0)=par.rec_loga(0)+log(thisSSB)-exp(par.rec_logb(0))*thisSSB;
      }else{
        if(conf.stockRecruitmentModelCode==2){//BH
          predN(0)=par.rec_loga(0)+log(thisSSB)-log(1.0+exp(par.rec_logb(0))*thisSSB); 
        }else{
          if(conf.stockRecruitmentModelCode==3){ // around mean value
           predN(0)=par.rec_loga(0);
          }else{
          error("SR model code not recognized");
         }
        }
      }
    }

    for(int j=0; j<dat.env.dim(1); ++j){
      if(dat.env.matrix().col(j).sum()!=0){
        predN(0) += par.rec_e(j)*dat.env(i-1,j);
      }
    } 
    //std::cout << "* older ages"  << std::endl;  
    for(int j=1; j<stateDimN; ++j){
      if(conf.keyLogFsta(0,j-1)>(-1)){
        predN(j)=logN(j-1,i-1)-exp(logF(j-1,i-1))-dat.natMor(i-1,j-1); 
      }else{
        predN(j)=logN(j-1,i-1)-dat.natMor(i-1,j-1); 
      }
    }  
    //std::cout << "* plus"  << std::endl;  
    if(conf.maxAgePlusGroup==1){
      predN(stateDimN-1)=log(exp(logN(stateDimN-2,i-1)-exp(logF(stateDimN-2,i-1))-dat.natMor(i-1,stateDimN-2))+
                             exp(logN(stateDimN-1,i-1)-exp(logF(stateDimN-1,i-1))-dat.natMor(i-1,stateDimN-1))); 
    }
  return predN;  
}
