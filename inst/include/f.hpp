template <class Type>
Type nllF(confSet &conf, paraSet<Type> &par, vector<Type> &logFy, vector<Type> &logitSel,data_indicator<vector<Type>,Type> &keep, objective_function<Type> *of){
  using CppAD::abs;

  Type nll=0; 
  Type sd= exp(par.logSdLogFsta);
 
  for(int j = 1 ; j < logFy.size(); ++j){
   nll -= dnorm(logFy(j), logFy(j-1), sd, true);  // random walk for Fy
       SIMULATE_F(of){
        if(conf.simFlag==0){
          logFy(j)=rnorm(logFy(j-1),sd);
        }
      }
  }
 
  if(CppAD::Variable(keep.sum())){ // add wide prior for first state, but _only_ when computing ooa residuals
    Type huge = 10;
    nll -= dnorm(logFy(0), Type(0), huge, true);  
  } 

  return nll;
}
