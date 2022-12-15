#include "../interface/MCSFitFunctions.h"

Double_t MCSFitFunction(Double_t *x,Double_t *par){
  Double_t arg;
  
  if(x[0] < par[3]) {
    arg = par[1]*par[1] + par[2]*par[2]*(x[0]-par[3])*(x[0]-par[3]);
  }
  else {
    arg = par[1]*par[1] + par[4]*par[4]*(x[0]-par[3])*(x[0]-par[3]);
  }
  Double_t fitval = par[0]+sqrt(arg);
  
  // arg = par[1]*par[1]+par[2]*par[2]*(x[0]-par[3])*(x[0]-par[3]);
  //double fitval = par[0]+sqrt(arg);
  return fitval;
}

