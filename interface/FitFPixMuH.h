#ifndef FitFPixMuH_h
#define FitFPixMuH_h

#include <string>
#include <vector>
#include <TObject.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TVirtualFitter.h>
#include <cstdint>
//#include "TFitter.h"
//#include "MCSFitFunctions.h"

struct FPixMuH{
  double bfield[3];
  double shiftx;        // Shift in x direction
  double shifty;        // Shift in y direction
  double shiftx_err;        // Shift in x direction
  double shifty_err;        // Shift in y direction  
};


class FitFPixMuH : public TObject {
  
 public:  
  FitFPixMuH();
  //-----------------------------------------------------------
  //void Add(const FPixMuH &fmh){        
  //  cmvar.push_back(fmh);
  // }
  ~FitFPixMuH() override;
  
  //~FitFPixMuH();
 
  void Add(const FPixMuH &fmh);
  int Fit();
  double GetMuH(){return muH;}
  double GetMuHErr(){return muHErr;}
  /*
  void Initialize(){
    cmvar.clear();
    muH = 0 ;
    muHErr = 0; 
  }
  */
  
 private:
  friend void fcn_func(Int_t& npar, Double_t* deriv, Double_t& f, Double_t *par, Int_t flag); 
  
  TVirtualFitter* minuit;
  std::vector<FPixMuH> cmvar;
  double muH;
  double muHErr; 
  double calcChi2(double par);
  
};

#endif
