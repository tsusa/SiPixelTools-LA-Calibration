#include "../interface/FitFPixMuH.h"

double FitFPixMuH::calcChi2(double par_0){

  double tanlorpertesla = par_0;
  //double thick = 290.;
  double tlpt2 = tanlorpertesla*tanlorpertesla;
  double v[3], tanx, tany, xshift, yshift;
  int n = cmvar.size();
  //   printf("npar = %d, f = %lf, par[0] = %lf, flag = %d, n = %d \n", npar, f, par[0], flag, n);
  double chi2 = 0.0;
  for (int i=0; i<n; i++){
    v[0] = -(tanlorpertesla*cmvar[i].bfield[1]+tlpt2*cmvar[i].bfield[0]*cmvar[i].bfield[2]);
    v[1] = tanlorpertesla*cmvar[i].bfield[0]-tlpt2*cmvar[i].bfield[1]*cmvar[i].bfield[2];
    v[2] = -(1.+tlpt2*cmvar[i].bfield[2]*cmvar[i].bfield[2]);
    
    tanx = v[0]/v[2];
    tany = v[1]/v[2];
    //     xshift = -tanx*thick/2.;
    // yshift = -tany*thick/2.;
    
    xshift = tanx;
    yshift = tany;
    
    chi2 +=
      (xshift - cmvar[i].shiftx)*(xshift - cmvar[i].shiftx)/cmvar[i].shiftx_err/cmvar[i].shiftx_err 
      + (yshift - cmvar[i].shifty)*(yshift - cmvar[i].shifty)/cmvar[i].shifty_err/cmvar[i].shifty_err;
    //
  }
  return chi2;  
}

void fcn_func(Int_t& npar, Double_t* deriv, Double_t& f, Double_t *par, Int_t flag){ 
  

  //f = (dynamic_cast<FitFPixMuH*>((TVirtualFitter::GetFitter())->GetObjectFit()))->calcChi2(par[0]);
  f=((dynamic_cast<FitFPixMuH*>((TVirtualFitter::GetFitter())->GetObjectFit()))->calcChi2(par[0]));
  //   printf("xshift/yshift = %lf/%lf \n", xshift,     
}    

FitFPixMuH::FitFPixMuH():
  muH(0),
  muHErr(0)
{}
FitFPixMuH::~FitFPixMuH(){}
//-----------------------------------------------------------
void FitFPixMuH::Add(const FPixMuH &fmh){        
  cmvar.push_back(fmh);
}
//-----------------------------------------------------------
int FitFPixMuH::Fit(){
  
  std::cout << "size: "  << cmvar.size() << std::endl; 
    
  //TVirtualFitter *minuit = TVirtualFitter::Fitter(0,1);    
  minuit = TVirtualFitter::Fitter(this,1);    
    
  minuit->SetFCN(fcn_func);
  // starting values 
  double tanlorpertesla = 0.08;
  //double startA = tanlorpertesla;
  
  // if not limited (vhigh <= vlow) 
  //minuit->SetParameter(0,"muH",startA,0.1,0.,0.);
  minuit->SetParameter(0, "muH", tanlorpertesla, 0.1, 0., 0.);
  // if not limited (vhigh <= vlow)
  //    minuit->SetParameter(1,"rho",startA,0.1,0.,0.);
  // create Minimizer (default is Migrad)
  
  Double_t arglist[100];
  arglist[0] = 3.;
  minuit->ExecuteCommand("SET PRINT", arglist, 1);
  Double_t up = 1.0;
  minuit->SetErrorDef(up);
  arglist[0] = 100000.;
  minuit->ExecuteCommand("MIGRAD", arglist, 0);
  
  std::cout << "Fitted parameter: " <<  minuit->GetParameter(0) << " "
	    << minuit->GetParError(0) << std::endl;   
  muH = minuit->GetParameter(0);
  muHErr = minuit->GetParError(0);
  
  return 0;
}
//-----------------------------------------------------------
