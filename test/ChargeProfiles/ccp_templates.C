#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "TF1.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TCanvas.h"

#include <sstream>
#include <iostream>
#include <fstream>

#define layersAndDisks 8 //irene
#define Layers 4 //irene
#define maxpix 500 //irene
#define modulesPerLayer 8

#define hist_drift_  400 //int
#define hist_depth_ 100//int
#define min_drift_ -1000 //double
#define max_drift_ 1000//double
#define min_depth_ -100//double
#define max_depth_ 400//double



struct pixinfo{

  Int_t npix; 
  Float_t row[maxpix];
  Float_t col[maxpix];
  Float_t adc[maxpix];
  Float_t x[maxpix];  
  Float_t y[maxpix];  

} pixinfo_;

struct trackhit{

  Float_t x;
  Float_t y;
  Double_t alpha;
  Double_t beta;
  Double_t gamma;

} simhit_, trackhit_;   //irene: track hit is the measured hit

typedef std::map<std::pair<TString, TString>, TH2D *> MapTH2;
TH2D * h_drift_depth;

void MyFill(MapTH2 * map, TString histtype, TString layer, TString bias, pixinfo * pixinfo_, double cotbeta,float ylim1,float ylim2, float xlim1, trackhit * trackhit_,double drcor);

void ccp_templates(TString NAME = "templates_reweight_102X_bugfix", TString filename = "reweightPaul_bugfix/reweight", int filenumber = 999, bool data = false){
  
  //********************  The results depend upon temperature and Hall Factor ****************************
  bool print = true;
  
  //**********************************  Define Run number and range of subfiles ******************************************

  //  int runnumber = 305406;  
  int runnumber;  
  TString sRunnumber;
  
  runnumber = 314649;
  //  runnumber[1] = 314650;
  //  runnumber[2] = 314663;
  sRunnumber.Form("%d",runnumber);
 
  double minfit = 5.; //needed
  double maxfit = 280.; //needed

  const double thick_ = 0.0285; //########### I suppose this is the sensor thickness in cm
  double ypitch_ = 0.0150; // ############## pixel long pitch (in the y direction) in cm 
  double halfypitch_ = ypitch_/2.;
  double cotbeta_min = 4.*ypitch_/thick_;  //########### see below, beta_min is considering 1/4 of the thickness


  //  typedef std::map<std::pair<TString, TString>, TH2D *> MapTH2;
  std::map<std::pair<TString, TString>, TH2D *>::iterator it2;
  typedef std::map<std::pair<TString, TString>, TH1D *> MapTH1;
  std::map<std::pair<TString, TString>, TH1D *>::iterator it1;
  std::map<std::pair<TString, TString>, TH1D *>::iterator it1c;

  MapTH2 m_drift_depth_adc;
  //  TH2D *h_drift_depth_adc2[BiasVoltages][Layers];
  MapTH2 m_drift_depth_noadc;
  MapTH1 m_Layers_depth_adc_ccp;
  MapTH1 m_Layers_depth_noadc_ccp;
  MapTH1 m_pq_vs_depth; //charge profile hist

  TH2D * h_drift_depth_adc;
  //  TH2D * h_drift_depth_adc2;
  TH2D * h_drift_depth_noadc;
  TH1D * hLayers_depth_adc_ccp;
  TH1D * hLayers_depth_noadc_ccp;
  TH1D * pq_vs_depth; //charge profile hist


  TString sLayer[Layers];
  sLayer[0] = "1";
  sLayer[1] = "2";
  sLayer[2] = "3";
  sLayer[3] = "4";

  Int_t run_;
  Int_t lumiblock_;
  ULong64_t event_;
  Int_t module_;
  Int_t ladder_;
  Int_t layer_;
  Int_t isflipped_;
  Float_t pt_;
  Float_t eta_;
  Float_t phi_;
  Double_t chi2_;
  Double_t ndof_;

   

  struct {

    Int_t ncol;
    Int_t dcol[maxpix];
    Float_t adc[maxpix];
    Float_t depth[maxpix];

  } colinfo_;
      

  struct {

    Float_t x;
    Float_t y;
    Float_t charge;
    Int_t size_x;
    Int_t size_y;
    Int_t maxPixelCol;
    Int_t maxPixelRow;
    Int_t minPixelCol;
    Int_t minPixelRow;

  } clust_;

   

  struct {

    Float_t x;
    Float_t y;

  } rechit_;   //irene: hit reconstucted

  //cuts
  float pt_cut = 1.5;
  float clusterSizeY_cut = 4; 
  float residual_cut = 0.005; //not used
  float normChi2_cut = 4.0; //not used
  float clusterCharge_cut = 30.;  //charge in ke per unit thickness
  float trackQuality_cut = 2.0; //not used
  int highPurity_cut = 1; //not used

  //cuts to compare with https://twiki.cern.ch/twiki/bin/view/CMSPublic/PixelOfflinePlotsAugust2017#Pixel_Charge_Profiles
  //FIXME exclude edge pixels 
  pt_cut = 3; //irene, to compare with https://twiki.cern.ch/twiki/bin/view/CMSPublic/PixelOfflinePlotsAugust2017#Pixel_Charge_Profiles (before 1.5)
  residual_cut = 0.01;//Hit residuals < 100 um //now sync
  float clusterSizeYlow_cut = 4; // 4 <= y < 8 
  float clusterSizeYhigh_cut = 8; //4 <= y < 8 
  clusterCharge_cut = 1000.;  //charge in ke < 600 ke- 
  //FIX ME histos : inner ladder,outer ladder, b , f , b inner, b outer, f inner, f outer, all



  int sfile = 0;
  //  int efile = 2301;  
  TString infile;
  TString srun;
  TString snfile;
  TFile *fout = new TFile(NAME+"/ccp_summary_nobetasel_"+NAME+".root","RECREATE");
  

  for(int nfile = sfile; nfile <= filenumber; ++nfile)
    {
      snfile=Form("%d",nfile);
      //	  cout << " Run : " << sRunnumber[k] << endl;
      //	  infile = "/eos/cms/store/group/dpg_tracker_pixel/comm_pixel/LATrees/2018/collisions/"+sRunnumber[k]+"_express/LATree_Collisions18_"+sRunnumber[k]+"_"+snfile+"_c.root";
      //infile = filename+snfile+".root";
      infile = filename+"_"+snfile+".root";
      if(data)
	infile = filename+"_"+snfile+"_c.root";
      
      TFile * f = new TFile(infile);
      if ((!f) || f->IsZombie()) 
	{
	  cout << "Error opening file "<< infile  << endl;
	  delete f; continue;
	}
      else
	{
	  cout << " opening file "<< infile  << endl;    
	  f->cd();
	  // fill the histrograms with the ntpl
	  TTree * LATree = (TTree*)f->Get("SiPixelLorentzAngleTree_"); //irene: it is in file
	  int nentries = LATree->GetEntries();
	  LATree->SetBranchAddress("run", &run_); //irene ok
	  LATree->SetBranchAddress("lumiblock", &lumiblock_); //irene ok
	  LATree->SetBranchAddress("event", &event_); //irene ok
	  LATree->SetBranchAddress("module", &module_); //irene ok
	  LATree->SetBranchAddress("ladder", &ladder_); //irene ok
	  LATree->SetBranchAddress("layer", &layer_); //irene ok
	  LATree->SetBranchAddress("isflipped", &isflipped_); //irene ok
	  LATree->SetBranchAddress("pt", &pt_); //irene ok
	  LATree->SetBranchAddress("eta", &eta_); //irene ok
	  LATree->SetBranchAddress("phi", &phi_);//irene ok
	  LATree->SetBranchAddress("chi2", &chi2_);//irene ok
	  LATree->SetBranchAddress("ndof", &ndof_);//irene ok
	  LATree->SetBranchAddress("trackhit", &trackhit_); //irene ok, with leaves
	  LATree->SetBranchAddress("simhit", &simhit_); //irene ok, with leaves
	  LATree->SetBranchAddress("npix", &pixinfo_.npix); //irene ok
	  LATree->SetBranchAddress("rowpix", pixinfo_.row); //irene ok
	  LATree->SetBranchAddress("colpix", pixinfo_.col); //irene ok
	  LATree->SetBranchAddress("adc", pixinfo_.adc); //irene ok
	  LATree->SetBranchAddress("xpix", pixinfo_.x); //irene ok
	  LATree->SetBranchAddress("ypix", pixinfo_.y); //irene ok
	  LATree->SetBranchAddress("clust", &clust_); // charge is given in 1000 e //irene ok, with leaves
	  LATree->SetBranchAddress("rechit", &rechit_); //irene ok, with leaves
	
	  cout << "Running over " << nentries << " hits" << endl;
	  
	  for(int ientrie = 0 ; ientrie < nentries; ientrie++)
	    {
	      if(print && ientrie%1000==0) cout << "inside for cycle " << ientrie << endl; //irene debug
	      LATree->GetEntry(ientrie);
	      double chi2_ndof = chi2_/ndof_;
	      if(pt_< pt_cut)continue;
	      //	    if(chi2_ndof > 2.)continue;
	      if(pixinfo_.npix > 25) continue;
	      bool large_pix = false;
	      int hindex = 0;
	      double cotbeta = 1./TMath::Tan(trackhit_.beta);
	      //	      if(fabs(cotbeta) <= cotbeta_min) continue; //irene commented to check L4
	      double cotalpha = 1./TMath::Tan(trackhit_.alpha);
	      double drdz = sqrt(1.+cotalpha*cotalpha+cotbeta*cotbeta);
	      float ccc = clusterCharge_cut*drdz;
	      //	    cout << " ccc = clusterCharge_cut*drdz "<<  ccc << endl;
	      //	    cout<< " clust_.charge " << clust_.charge << endl;
	      double qclus = clust_.charge/drdz;
	      double drcor = drdz/fabs(cotbeta);
	      
	      // Find the track ends for later use
	      float ylim1 = trackhit_.y - thick_*cotbeta/2.;
	      float ylim2 = trackhit_.y + thick_*cotbeta/2.;
	      float xlim1 = trackhit_.x - thick_*cotalpha/2.;  
	      colinfo_.ncol = 0;
	    
	      TString module;
		  
	      //rejectEvent_BpixCollision
	      for (int j = 0; j <  pixinfo_.npix; j++)
		{
		  int colpos = static_cast<int>(pixinfo_.col[j]);
		  if (      pixinfo_.row[j] == 0 || pixinfo_.row[j] == 79 || pixinfo_.row[j] == 80 || pixinfo_.row[j] == 159 || colpos % 52 == 0 || colpos % 52 == 51 ){ large_pix = true;} //irene: the event is going over several chips
		}//for on pixinfo_.npix
	  	  
	      if(large_pix) continue;//large_pix sync


	      hindex = 2*(layer_-1);
	      if(module_ > 4) { hindex +=1;}
	      
	      double residual = TMath::Sqrt( (trackhit_.x - rechit_.x) * (trackhit_.x - rechit_.x) + (trackhit_.y - rechit_.y) * (trackhit_.y - rechit_.y) ); //residual sync
	      //	    cout << "residual " << residual << endl;
	      //	    if( (clust_.size_y >= clusterSizeY_cut) && (clust_.charge < ccc) )
	      if( 
		 (clust_.size_y >= clusterSizeYlow_cut) &&
		 //	       && (clust_.size_y<clusterSizeYhigh_cut)  		 && 
		 (clust_.charge < clusterCharge_cut) 
		 && (residual < residual_cut)
		  )//cuts to compare
		{
		  if(layer_==1)
		    {
		      //cout << " layer " << layer_ << " module " << module_ << endl;
		      MyFill(&m_drift_depth_adc, "adc","L1","all", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
		      MyFill(&m_drift_depth_noadc, "noadc","L1","all", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
		      //cout << "my fill adc done" <<endl;
		      
		      // if(module_ ==1
		      // 	 ||module_==2
		      // 	 ||module_==7
		      // 	 ||module_==8
		      // 	 ) 
		      // 	{
		      // 	  //cout << " layer " << layer_ << " module " << module_ << endl;
		      // 	  MyFill(&m_drift_depth_adc, "adc","L1","3001", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
		      // 	  MyFill(&m_drift_depth_noadc, "noadc","L1","3001", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
		      // 	  //cout << "my fill adc done" <<endl;
		      // 	}
		      // if(module_ ==4
		      // 	 ||module_==5
		      // 	 ||module_==6
		      // 	 ||module_==3
		      // 	 ) 
		      // 	{
		      // 	  //cout << " layer " << layer_ << " module " << module_ << endl;
		      // 	  MyFill(&m_drift_depth_adc, "adc","L1","3002", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
		      // 	  MyFill(&m_drift_depth_noadc, "noadc","L1","3002", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
		      // 	  //cout << "my fill adc done" <<endl;
		      // 	}

		      
		      module.Form("%i",1);
		      if(module_ == 1)
			{
			  //cout << " layer " << layer_ << " module " << module_ << endl;
			  MyFill(&m_drift_depth_adc, "adc","L1",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  MyFill(&m_drift_depth_noadc, "noadc","L1",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  //cout << "my fill adc done" <<endl;
			}

		      module.Form("%i",2);
		      if(module_ == 2)
			{
			  //cout << " layer " << layer_ << " module " << module_ << endl;
			  MyFill(&m_drift_depth_adc, "adc","L1",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  MyFill(&m_drift_depth_noadc, "noadc","L1",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  //cout << "my fill adc done" <<endl;
			}
		      module.Form("%i",3);
		      if(module_ == 3)
			{
			  //cout << " layer " << layer_ << " module " << module_ << endl;
			  MyFill(&m_drift_depth_adc, "adc","L1",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  MyFill(&m_drift_depth_noadc, "noadc","L1",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  //cout << "my fill adc done" <<endl;
			}
		      module.Form("%i",4);
		      if(module_ == 4)
			{
			  //cout << " layer " << layer_ << " module " << module_ << endl;
			  MyFill(&m_drift_depth_adc, "adc","L1",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  MyFill(&m_drift_depth_noadc, "noadc","L1",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  //cout << "my fill adc done" <<endl;
			}
		      module.Form("%i",5);
		      if(module_ == 5)
			{
			  //cout << " layer " << layer_ << " module " << module_ << endl;
			  MyFill(&m_drift_depth_adc, "adc","L1",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  MyFill(&m_drift_depth_noadc, "noadc","L1",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  //cout << "my fill adc done" <<endl;
			}
		      module.Form("%i",6);
		      if(module_ == 6)
			{
			  //cout << " layer " << layer_ << " module " << module_ << endl;
			  MyFill(&m_drift_depth_adc, "adc","L1",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  MyFill(&m_drift_depth_noadc, "noadc","L1",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  //cout << "my fill adc done" <<endl;
			}
		      module.Form("%i",7);
		      if(module_ == 7)
			{
			  //cout << " layer " << layer_ << " module " << module_ << endl;
			  MyFill(&m_drift_depth_adc, "adc","L1",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  MyFill(&m_drift_depth_noadc, "noadc","L1",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  //cout << "my fill adc done" <<endl;
			}
		      module.Form("%i",8);
		      if(module_ == 8)
			{
			  //cout << " layer " << layer_ << " module " << module_ << endl;
			  MyFill(&m_drift_depth_adc, "adc","L1",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  MyFill(&m_drift_depth_noadc, "noadc","L1",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  //cout << "my fill adc done" <<endl;
			}
		      module.Form("%i",0);
		      if(module_ == 0)
			{
			  //cout << " layer " << layer_ << " module " << module_ << endl;
			  MyFill(&m_drift_depth_adc, "adc","L1",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  MyFill(&m_drift_depth_noadc, "noadc","L1",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  //cout << "my fill adc done" <<endl;
			}
	

		      if(clust_.size_y >= 5 )
			{
			  MyFill(&m_drift_depth_adc, "adc","L1","size5m", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  MyFill(&m_drift_depth_noadc, "noadc","L1","size5m", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			}
		      if(clust_.size_y >= 6 )
			{
			  MyFill(&m_drift_depth_adc, "adc","L1","size6m", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  MyFill(&m_drift_depth_noadc, "noadc","L1","size6m", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			}
		      if(clust_.size_y >= 7 )
			{
			  MyFill(&m_drift_depth_adc, "adc","L1","size7m", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  MyFill(&m_drift_depth_noadc, "noadc","L1","size7m", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			}

		      if(clust_.size_y >= 8 )
			{
			  MyFill(&m_drift_depth_adc, "adc","L1","size8m", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  MyFill(&m_drift_depth_noadc, "noadc","L1","size8m", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			}

		      if(clust_.size_y >= 9 )
			{
			  MyFill(&m_drift_depth_adc, "adc","L1","size9m", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  MyFill(&m_drift_depth_noadc, "noadc","L1","size9m", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			}

		      if(clust_.size_y == 10 )
			{
			  MyFill(&m_drift_depth_adc, "adc","L1","size10", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  MyFill(&m_drift_depth_noadc, "noadc","L1","size10", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			}

			


		    }
		  if(layer_ == 2)
		    {
		      
		      //cout << " layer " << layer_ << " module " << module_ << endl;
		      MyFill(&m_drift_depth_adc, "adc","L2","all", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
		      MyFill(&m_drift_depth_noadc, "noadc","L2","all", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
		      //cout << "my fill adc done" <<endl;
		      
		      module.Form("%i",1);
		      if(module_ == 1)
			{
			  //cout << " layer " << layer_ << " module " << module_ << endl;
			  MyFill(&m_drift_depth_adc, "adc","L2",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  MyFill(&m_drift_depth_noadc, "noadc","L2",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  //cout << "my fill adc done" <<endl;
			}

		      module.Form("%i",2);
		      if(module_ == 2)
			{
			  //cout << " layer " << layer_ << " module " << module_ << endl;
			  MyFill(&m_drift_depth_adc, "adc","L2",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  MyFill(&m_drift_depth_noadc, "noadc","L2",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  //cout << "my fill adc done" <<endl;
			}
		      module.Form("%i",3);
		      if(module_ == 3)
			{
			  //cout << " layer " << layer_ << " module " << module_ << endl;
			  MyFill(&m_drift_depth_adc, "adc","L2",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  MyFill(&m_drift_depth_noadc, "noadc","L2",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  //cout << "my fill adc done" <<endl;
			}
		      module.Form("%i",4);
		      if(module_ == 4)
			{
			  //cout << " layer " << layer_ << " module " << module_ << endl;
			  MyFill(&m_drift_depth_adc, "adc","L2",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  MyFill(&m_drift_depth_noadc, "noadc","L2",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  //cout << "my fill adc done" <<endl;
			}
		      module.Form("%i",5);
		      if(module_ == 5)
			{
			  //cout << " layer " << layer_ << " module " << module_ << endl;
			  MyFill(&m_drift_depth_adc, "adc","L2",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  MyFill(&m_drift_depth_noadc, "noadc","L2",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  //cout << "my fill adc done" <<endl;
			}
		      module.Form("%i",6);
		      if(module_ == 6)
			{
			  //cout << " layer " << layer_ << " module " << module_ << endl;
			  MyFill(&m_drift_depth_adc, "adc","L2",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  MyFill(&m_drift_depth_noadc, "noadc","L2",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  //cout << "my fill adc done" <<endl;
			}
		      module.Form("%i",7);
		      if(module_ == 7)
			{
			  //cout << " layer " << layer_ << " module " << module_ << endl;
			  MyFill(&m_drift_depth_adc, "adc","L2",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  MyFill(&m_drift_depth_noadc, "noadc","L2",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  //cout << "my fill adc done" <<endl;
			}
		      module.Form("%i",8);
		      if(module_ == 8)
			{
			  //cout << " layer " << layer_ << " module " << module_ << endl;
			  MyFill(&m_drift_depth_adc, "adc","L2",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  MyFill(&m_drift_depth_noadc, "noadc","L2",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  //cout << "my fill adc done" <<endl;
			}
		      module.Form("%i",0);
		      if(module_ == 0)
			{
			  //cout << " layer " << layer_ << " module " << module_ << endl;
			  MyFill(&m_drift_depth_adc, "adc","L2",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  MyFill(&m_drift_depth_noadc, "noadc","L2",module, &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
			  //cout << "my fill adc done" <<endl;
			}



		    }
		  if(layer_==3)
		    {		
		      //cout << " layer " << layer_ << " module " << module_ << endl;
		      MyFill(&m_drift_depth_adc, "adc","L3","all", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
		      MyFill(&m_drift_depth_noadc, "noadc","L3","all", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
		      //cout << "my fill adc done" <<endl;
		      // if(module_ ==1
		      // 	 ||module_==2
		      // 	 ||module_==3
		      // 	 ||module_==4
		      // 	 ) 
		      // 	{
		      // 	  //cout << " layer " << layer_ << " module " << module_ << endl;
		      // 	  MyFill(&m_drift_depth_adc, "adc","L3","3322", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
		      // 	  MyFill(&m_drift_depth_noadc, "noadc","L3","3322", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
		      // 	  //cout << "my fill adc done" <<endl;
		      // 	}
		      // if(module_ ==7
		      // 	 ||module_==5
		      // 	 ||module_==6
		      // 	 ||module_==8
		      // 	 ) 
		      // 	{
		      // 	  //cout << " layer " << layer_ << " module " << module_ << endl;
		      // 	  MyFill(&m_drift_depth_adc, "adc","L3","3422", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
		      // 	  MyFill(&m_drift_depth_noadc, "noadc","L3","3422", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
		      // 	  //cout << "my fill adc done" <<endl;
		      // 	}
		    }
		  if(layer_==4)
		    {
		      //cout << " layer " << layer_ << " module " << module_ << endl;
		      MyFill(&m_drift_depth_adc, "adc","L4","all", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
		      MyFill(&m_drift_depth_noadc, "noadc","L4","all", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
		      //cout << "my fill adc done" <<endl;
		      // if(module_ ==1
		      // 	 ||module_==2
		      // 	 ||module_==3
		      // 	 ||module_==4
		      // 	 ) 
		      // 	{
		      // 	  //cout << " layer " << layer_ << " module " << module_ << endl;
		      // 	  MyFill(&m_drift_depth_adc, "adc","L4","3732", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
		      // 	  MyFill(&m_drift_depth_noadc, "noadc","L4","3732", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
		      // 	  //cout << "my fill adc done" <<endl;
		      // 	}
		      // if(module_ ==7
		      // 	 ||module_==5
		      // 	 ||module_==6
		      // 	 ||module_==8
		      // 	 ) 
		      // 	{
		      // 	  //cout << " layer " << layer_ << " module " << module_ << endl;
		      // 	  MyFill(&m_drift_depth_adc, "adc","L4","3832", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
		      // 	  MyFill(&m_drift_depth_noadc, "noadc","L4","3832", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
		      // 	  //cout << "my fill adc done" <<endl;
		      // 	}
		    }
		}//if( (clust_.size_y >= clusterSizeY_cut) && (clust_.charge <ccc) )

	      
	      if((clust_.size_y == 2)  && (residual < residual_cut) && (clust_.charge < clusterCharge_cut) )
                {
                  if(layer_==3)
                    {
                      MyFill(&m_drift_depth_adc, "adc","L3","all_clsize2", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                      MyFill(&m_drift_depth_noadc, "noadc","L3","all_clsize2", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                    }
                  if(layer_==4)
                    {
                      MyFill(&m_drift_depth_adc, "adc","L4","all_clsize2", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                      MyFill(&m_drift_depth_noadc, "noadc","L4","all_clsize2", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                    }
                }
              if((clust_.size_y == 3)  && (residual < residual_cut) && (clust_.charge < clusterCharge_cut)  )
                {
                  if(layer_==3)
                    {
                      MyFill(&m_drift_depth_adc, "adc","L3","all_clsize3", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                      MyFill(&m_drift_depth_noadc, "noadc","L3","all_clsize3", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                    }
                  if(layer_==4)
                    {
                      MyFill(&m_drift_depth_adc, "adc","L4","all_clsize3", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                      MyFill(&m_drift_depth_noadc, "noadc","L4","all_clsize3", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                    }
                }
              if((clust_.size_y == 4)  && (residual < residual_cut)&& (clust_.charge < clusterCharge_cut)  )
                {
                  if(layer_==3)
                    {
                      MyFill(&m_drift_depth_adc, "adc","L3","all_clsize4", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                      MyFill(&m_drift_depth_noadc, "noadc","L3","all_clsize4", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                    }
                  if(layer_==4)
                    {
                      MyFill(&m_drift_depth_adc, "adc","L4","all_clsize4", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                      MyFill(&m_drift_depth_noadc, "noadc","L4","all_clsize4", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                    }
		}
	      if((clust_.size_y < 4)  && (residual < residual_cut)&& (clust_.charge < clusterCharge_cut)  )
		{
                  if(layer_==3)
                    {
                      MyFill(&m_drift_depth_adc, "adc","L3","all_clsizeless4", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                      MyFill(&m_drift_depth_noadc, "noadc","L3","all_clsizeless4", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                    }
                  if(layer_==4)
                    {
                      MyFill(&m_drift_depth_adc, "adc","L4","all_clsizeless4", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                      MyFill(&m_drift_depth_noadc, "noadc","L4","all_clsizeless4", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                    }
                }
	      if((clust_.size_y == 5)  && (residual < residual_cut)&& (clust_.charge < clusterCharge_cut)  )
		{
                  if(layer_==3)
                    {
                      MyFill(&m_drift_depth_adc, "adc","L3","all_clsize5", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                      MyFill(&m_drift_depth_noadc, "noadc","L3","all_clsize5", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                    }
                  if(layer_==4)
                    {
                      MyFill(&m_drift_depth_adc, "adc","L4","all_clsize5", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                      MyFill(&m_drift_depth_noadc, "noadc","L4","all_clsize5", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                    }
                }
	      if((clust_.size_y == 6)  && (residual < residual_cut)&& (clust_.charge < clusterCharge_cut)  )
		{
                  if(layer_==3)
                    {
                      MyFill(&m_drift_depth_adc, "adc","L3","all_clsize6", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                      MyFill(&m_drift_depth_noadc, "noadc","L3","all_clsize6", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                    }
                  if(layer_==4)
                    {
                      MyFill(&m_drift_depth_adc, "adc","L4","all_clsize6", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                      MyFill(&m_drift_depth_noadc, "noadc","L4","all_clsize6", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                    }
                }
	      if((clust_.size_y == 7)  && (residual < residual_cut)&& (clust_.charge < clusterCharge_cut)  )
		{
                  if(layer_==3)
                    {
                      MyFill(&m_drift_depth_adc, "adc","L3","all_clsize7", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                      MyFill(&m_drift_depth_noadc, "noadc","L3","all_clsize7", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                    }
                  if(layer_==4)
                    {
                      MyFill(&m_drift_depth_adc, "adc","L4","all_clsize7", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                      MyFill(&m_drift_depth_noadc, "noadc","L4","all_clsize7", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                    }
                }
	      if((clust_.size_y == 8)  && (residual < residual_cut)&& (clust_.charge < clusterCharge_cut)  )
		{
                  if(layer_==3)
                    {
                      MyFill(&m_drift_depth_adc, "adc","L3","all_clsize8", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                      MyFill(&m_drift_depth_noadc, "noadc","L3","all_clsize8", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                    }
                  if(layer_==4)
                    {
                      MyFill(&m_drift_depth_adc, "adc","L4","all_clsize8", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                      MyFill(&m_drift_depth_noadc, "noadc","L4","all_clsize8", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                    }
                }
	      if((clust_.size_y > 8 )  && (residual < residual_cut)&& (clust_.charge < clusterCharge_cut)  )
		{
                  if(layer_==3)
                    {
                      MyFill(&m_drift_depth_adc, "adc","L3","all_clsizemore8", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                      MyFill(&m_drift_depth_noadc, "noadc","L3","all_clsizemore8", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                    }
                  if(layer_==4)
                    {
                      MyFill(&m_drift_depth_adc, "adc","L4","all_clsizemore8", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                      MyFill(&m_drift_depth_noadc, "noadc","L4","all_clsizemore8", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);
                    }
                }
	      
	    }//entries for
	}//else: file is good!
    }//for cicle on files
  
  cout << "preparing root file" << endl;   
  fout->cd();
  
  for ( it2 = m_drift_depth_adc.begin(); it2 != m_drift_depth_adc.end(); it2++ )
    {
      it2->second->Write();
      hLayers_depth_adc_ccp  = it2->second->ProjectionY();
      hLayers_depth_adc_ccp->Write();
      m_Layers_depth_adc_ccp.insert(std::make_pair(it2->first,hLayers_depth_adc_ccp ));
    }
  for ( it2 = m_drift_depth_noadc.begin(); it2 != m_drift_depth_noadc.end(); it2++ )
    {
      (it2->second)->Write();
      hLayers_depth_noadc_ccp  = (it2->second)->ProjectionY();
      hLayers_depth_noadc_ccp->Write();
      m_Layers_depth_noadc_ccp.insert(std::make_pair(it2->first,hLayers_depth_noadc_ccp ));
    }     
  //  it1c= m_Layers_depth_noadc_ccp.begin();
  for ( it1 = m_Layers_depth_adc_ccp.begin(); it1 != m_Layers_depth_adc_ccp.end(); it1++ )
    {
      cout << "looping on adc projection map" << endl;         
      pq_vs_depth = (TH1D*) it1->second->Clone();
      //      cout << "cloned ccp histo" << endl;         
      pq_vs_depth->SetName("h_ccp_layer_"+(it1->first).first+"_"+(it1->first).second);
      //      cout << "set name" << endl;         
      pq_vs_depth->SetTitle("h_ccp_layer_"+(it1->first).first+"_"+(it1->first).second);
      //      cout << "set title" << endl;         
         
      it1c = m_Layers_depth_noadc_ccp.find(std::make_pair((it1->first).first,(it1->first).second));
      cout << "no adc iterator " << (it1->first).first << " " << (it1->first).second << endl;
      if(it1c != m_Layers_depth_noadc_ccp.end())
	{
	  pq_vs_depth->Divide(it1c->second);     
	  cout << "divided" << endl;         
	}
      else
	pq_vs_depth->SetName("h_ccp_layer_NOTDIVIDED_"+(it1->first).first+"_"+(it1->first).second);
      pq_vs_depth->Write();
      cout << "writed" << endl;         
      cout << it1->first.first << " " << it1->first.second << endl;
      cout << pq_vs_depth->GetEntries() << endl;
      m_pq_vs_depth.insert(std::make_pair(it1->first,pq_vs_depth));	  
      //      cout << "inserted ccp histo" << endl;         
      //      it1c++;
      //      cout << "incremented noadc projection map iterator" << endl;         
    }     
      
  
  fout->Close();
  
      
}
				

void MyFill(MapTH2 * map, TString histtype, TString layer, TString bias,pixinfo * pixinfo_, double cotbeta,float ylim1,float ylim2, float xlim1, trackhit * trackhit_,double drcor)
{
  double ypitch_ = 0.0150; // ############## pixel long pitch (in the y direction) in cm 
  double halfypitch_ = ypitch_/2.;
  auto it = map->find(std::make_pair(layer,bias));
 
  if(it  != map->end())
    {
      //cout << " found map with " << layer << " and " << bias  << endl;
      //      //cout << "infact key " << it->first.first << " " << it->first.second << endl;
    }
  else
    {
      //      //cout << " not found map with " << layer << " and " << bias << endl;
      h_drift_depth = new TH2D("h_drift_depth_"+histtype+"_"+layer+"_"+bias,"h_drift_depth_"+histtype+"_"+layer+"_"+bias,
                               hist_drift_ , min_drift_, max_drift_,
                               hist_depth_, min_depth_, max_depth_);
      //      //cout << " initialized hist" << endl;
      map->insert(std::make_pair(std::make_pair(layer,bias),h_drift_depth));
      //      //cout << "map key " << it->first.first << " " << it->first.second << endl;

      //                                      MyFill(h_drift_depth_adc, "adc","L1","REFstart-350V", &pixinfo_,cotbeta,ylim1,ylim2,xlim1, &trackhit_,drcor);                                                  
      //      //cout << " inserted " << endl;
      //     h_drift_depth->Reset();
    }

  for (int j = 0; j <  pixinfo_->npix; j++)
    {    
      // hist = new TH2D("h_drift_depth_"+histtype+"_"+layer+"_"+bias,"h_drift_depth_"+histtype+"_"+layer+"_"+bias,
      // 		      hist_drift_ , min_drift_, max_drift_,
      // 		      hist_depth_, min_depth_, max_depth_);
      // Handle the end pixels using the track segment in those pixels
      float ypixlow = pixinfo_->y[j]-halfypitch_;
      float ypixhigh = pixinfo_->y[j]+halfypitch_;
      if(cotbeta > 0.) 
	{
	  if(ylim1 > ypixlow) ypixlow = ylim1;
	  if(ylim2 < ypixhigh) ypixhigh = ylim2;
	}
      else 
	{
	  if(ylim2 > ypixlow) ypixlow = ylim2;
	  if(ylim1 < ypixhigh) ypixhigh = ylim1;
	}
			  
      float ypixavg = 0.5*(ypixlow+ypixhigh);		   
      float dypix = fabs(ypixhigh-ypixlow);
      float dx = (pixinfo_->x[j] - xlim1) * 10000.;   //irene:micron?      
      float dy = (pixinfo_->y[j] - ylim1) * 10000.;         //Tanja    
      //		    float dy = (ypixavg - ylim1) * 10000.;         //Morris    
      float depth = dy * tan(trackhit_->beta);
      float drift = dx - dy * tan(trackhit_->gamma);
      double drpix = 10000. * dypix*drcor;	  
	
      it = map->find(std::make_pair(layer,bias));	  
      if(it != map->end())
	{
	  //	  //cout << "map is there" << endl;
	  if(histtype == "adc")
	    (it->second)->Fill(drift, depth, pixinfo_->adc[j]);
	  //        map[std::make_pair(layer,bias)]->Fill(drift, depth, pixinfo_->adc[j]);
	  //	hist->Fill(drift, depth, pixinfo_->adc[j]);
	  if(histtype == "noadc")
	    (it->second)->Fill(drift, depth);        
	  //      map[std::make_pair(layer,bias)]->Fill(drift, depth);
	  //	hist->Fill(drift, depth);
	}
      else cout << "no map, why?" << endl;
    }
}
