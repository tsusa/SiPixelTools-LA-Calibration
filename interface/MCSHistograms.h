#ifndef MCSHistograms_h
#define MCSHistograms_h

#include <string>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "MCSFitFunctions.h"

class MCSHistograms {
  
 public:  
 MCSHistograms():
  isBooked(0),
    isMeanCalculated(0),
    h_size_cotanAngle(0),
    h_mean_cotanAngle(0),
    h_slice(0),
    ff(0),
    minCotAngle(-999),
    minCotAngleError(-999)    
      {}
  
  MCSHistograms(edm::Service<TFileService> fs, 
		std::string name, std::string angle, std::string title,
		int n_x, int n_y){
    Book(fs, name, angle, title, n_x, n_y);
  }
  
  void Book(edm::Service<TFileService> fs, std::string name, std::string angle, std::string title,
	    int n_x, int n_y){    
    isBooked = 1;
    
    std::string axis = angle == "alpha" ? "x": "y";    
    std::string histo_name = "h_size" + axis + "_cotan" + angle 
      + "_" + name;    
    std::string axisTitles = ";cotan(#" + angle +") ; cluster size " 
      + axis + " [pixel]"; 
    
    h_size_cotanAngle = fs->make<TH2F>(histo_name.c_str(), 
				       axisTitles.c_str(),
				       n_x, -3., 3., 
				       n_y, 0.5, n_y + 0.5);
    
    histo_name = "h_mean_cotan" + angle + "_" + name;
    axisTitles = ";cotan(#" + angle +") ; Average cluster size " + 
      axis + " [pixel]";     
    
    h_mean_cotanAngle = fs->make<TH1F>(histo_name.c_str(), 
				       axisTitles.c_str(), n_x, -3., 3.);
    
    histo_name = "h_slice_" + name;
    h_slice = fs->make<TH1F>(histo_name.c_str(), 
			     ";cluster size; Entries", 
			     n_y, 0.5, n_y + 0.5);   
  }  
  
  void Fill(float cotAngle, float size){    
    h_size_cotanAngle->Fill(cotAngle, size);
  }  
  
  void CalcMean(){
    
    isMeanCalculated = 1;    
    int n_x = h_size_cotanAngle->GetXaxis()->GetNbins();
    int n_y = h_size_cotanAngle->GetYaxis()->GetNbins();
    
    for(int i = 1; i <= n_x; i++){
      h_slice->Reset("ICE");
      
      //loop over bins in size
      
      for( int j = 1; j<= n_y; j++){
        h_slice-> SetBinContent(j, h_size_cotanAngle->GetBinContent(i,j));
      }
      double mean = h_slice->GetMean(1);
      double error = h_slice->GetMeanError(1);
      
      h_mean_cotanAngle->SetBinContent(i, mean);
      h_mean_cotanAngle->SetBinError(i, error);
    } // end loop over bins in depth
  }
  
  int checkIfBinsAreFilled(double fitMin, double fitMax){
    int binIdMin = h_mean_cotanAngle->GetBin(fitMin);
    int binIdMax = h_mean_cotanAngle->GetBin(fitMax);
    if  (binIdMin <= 0)
      binIdMin = 1;
    if  (binIdMin > h_mean_cotanAngle->GetNbinsX())
      binIdMax = h_mean_cotanAngle->GetNbinsX();
    
    for (int i = binIdMin; i<= binIdMax; ++i){
      if (h_mean_cotanAngle->GetBinContent(i) == 0)
	return 1;
    }
    return 0;
  }
  
  int Fit(double fitRegion){    
    if (h_size_cotanAngle->GetEntries() == 0)
      return 1;
    
    ff = new TF1("ff", MCSFitFunction, -3, 3, 5);
     
    ff->SetParNames("Offset","RMS Constant","SlopeL",
		    "cot(alpha)_min","SlopeR");
    ff->SetParameters(1., 0.1, -1.6, 0, 1.6);             
    if (!isMeanCalculated){
      CalcMean();
    }
    
    int nFits = 0;
    minCotAngle = -999;
    minCotAngleError = -999;   
    while (nFits < 5){
      nFits++;
      double fitMin = ff->GetParameter(3) - fitRegion;
      double fitMax = ff->GetParameter(3) + fitRegion;
      if (checkIfBinsAreFilled(fitMin, fitMax))
	return 1;      
      
      TFitResultPtr r = h_mean_cotanAngle->Fit(ff,"ERSM", "", fitMin, fitMax);
      //std::cout << " " << ff->GetParameter(3) << " +- " 
      //		<< ff->GetParError(3) << std::endl;
      //std::cout << "Fitting done " << h_mean_cotanAngle->GetName() << std::endl;
      if (r->IsValid()){
	minCotAngle = ff->GetParameter(3);
	minCotAngleError = ff->GetParError(3);   
      }
    }
    return 0;  
  }
  
  void Write(TFile *f){
    h_size_cotanAngle->Write();
    if (h_mean_cotanAngle != 0)
      h_mean_cotanAngle->Write();
  }
  
  TH2F *GetSizeCotanAngleHisto(){return h_size_cotanAngle;}
  TH1F* &GetMeCotanAngleHisto(){return h_mean_cotanAngle;}

  double GetMinCotAngle(){return minCotAngle;}
  double GetMinCotAngleError(){return minCotAngleError;}
  
  
 private:
  int isBooked;
  int isMeanCalculated;
  TH2F *h_size_cotanAngle;
  TH1F *h_mean_cotanAngle;
  TH1F *h_slice;  
  TF1  *ff;
  double minCotAngle;
  double minCotAngleError;    
};

#endif
