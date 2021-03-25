#include "TH1.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <math.h>
#include "TMath.h"
#include "TFile.h"
#include "TFile.h"
#include "TString.h"

#include "tdrstyle.C"
#include "CMS_lumi.C"


using namespace std;
bool print=true;

#define comparisons 3
#define layersAndDisks 4


TString output_folder = "/afs/cern.ch/work/i/izoi/private/EPR/CMSSW_9_4_1/src/comparisonPlots/";


void DrawNorm55(TH1F **h, bool drawing[comparisons], TString name, TString layer, TString labels[comparisons], TString templates, TString additional[comparisons],float xmin = 0., float xmax = 285.,  float ymin = 0., float ymax = 20000., bool isNorm = false, float rescale = 1.);
void DrawNorm1(TH1F **h, bool drawing[comparisons], TString name, TString layer, TString labels[comparisons], TString templates, TString additional[comparisons],float xmin = 0., float xmax = 285.,  float ymin = 0., float ymax = 20000., bool isNorm = false);
void DrawNormMC(TH1F **h, bool drawing[comparisons], TString name, TString layer, TString labels[comparisons], TString templates, TString additional[comparisons],float xmin = 0., float xmax = 285.,  float ymin = 0., float ymax = 20000., bool isNorm = false);
void DrawNorm55ratio(TH1F **h, bool drawing[comparisons], TString name, TString layer, TString labels[comparisons], TString templates, TString additional[comparisons],float xmin = 0., float xmax = 285.,  float ymin = 0., float ymax = 20000., bool isNorm = false);

void TDR(TString lumi = "");
void TDR2(TCanvas * c_all);
void norm55(TH1F *h);
void norm1(TH1F *h);

void compare_2018ULshaded_scaled(TString nomignolo = "UL_StartEnd_shaded_rescale")
{

  TString sLayer[layersAndDisks];
  sLayer[0] = "L1";
  sLayer[1] = "L2";
  sLayer[2] = "L3";
  sLayer[3] = "L4";

  TString labels[comparisons];
  //  labels[1] = "data2017_run304292_rereco_forUL";
  //labels[0] = "MC2017_forULprep";//dataBegin2018_run315784";
  //  labels[0] = "data2016_run281707_forUL";
  labels[0] = "data2018Start_ULbow";
  labels[1] = "data2018End_ULbow";
  labels[2] = "MC2018rew_ULbow"; //reweight2018MCUL";

  TString shortlabels[layersAndDisks];
  shortlabels[0] = "Layer 1"; // (Phase-0)";//"data2018_run317435";
  shortlabels[1] = "Layer 2"; // (Phase-0)";
  shortlabels[2] = "Layer 3"; // (Phase-0)";
  shortlabels[3] = "Layer 4";

  float rescale[layersAndDisks] = {1.,1./0.8,1./0.84,1./0.88}; 
  TString mc[comparisons];

  mc[0] = "Data beginning 2018, Run 315257 (Era A)"; 
  mc[1] = "Data end 2018, Run 325175 (Era D)";
  mc[2] = "Run-2 Legacy simulation 2018";
  //mc[2] = "Realistic simulation 2016";

  TString base_dir = "/afs/cern.ch/work/i/izoi/private/EPR/CMSSW_9_4_1/src/";

  TString sName[comparisons];
  for(int i=0; i<comparisons;i++)
    {
	sName[i]=labels[i]+"/ccp_summary_nobetasel_"+labels[i]+".root";
	if(print) cout << sName[i] << endl;
    }


  TString data_input[comparisons];

  TFile * f_input[comparisons];
  TH1F * h_chargeProfile[layersAndDisks][comparisons];
  // TProfile *pq_vs_depth[layersAndDisks];
  
  bool drawing[comparisons];

  for(int i=0; i<comparisons; i++)
    {
           drawing[i]=true;
      
      for(int j = 0; j < layersAndDisks; j++)
	{

	  data_input[i] = base_dir+sName[i];
	  if(print) std::cout << "data_input  " << data_input[i]  <<std::endl;
	  
	  f_input[i] = new TFile(data_input[i]);
	      
	  h_chargeProfile[j][i] = (TH1F*)f_input[i]->Get("h_ccp_layer_"+sLayer[j]+"_all");
	  if(print) std::cout << "histos found  "  <<std::endl;
	  
	}
    }

      for(int j = 0; j<layersAndDisks; j++)
	{
	  DrawNorm55(h_chargeProfile[j], drawing, nomignolo+"data2018ULbow_StartEnd", shortlabels[j], mc,sLayer[j], mc,5., 285.,0.,20000, false,rescale[j]);

	  //	  DrawNorm55(h_chargeProfile[j], drawing, nomignolo+"norm55_MC_data2017run304292rereco", shortlabels[j], mc,sLayer[j], mc,0., 285.,0.,1.2, true);
	  //DrawNorm55(h_chargeProfile[j], drawing, nomignolo+"norm55_data2016UL_StartEnd", shortlabels[j], mc,sLayer[j], mc,0., 285.,0.,1.2, true);
	  //	  DrawNormMC(h_chargeProfile[j], drawing, nomignolo+"normMC2_data2018ULbow_StartEnd", shortlabels[j], mc,sLayer[j], mc,0., 285.,0.,1.2, true);
	}



}//compare



                                                                                                                                                                                 
void DrawNormMC(TH1F **h,bool drawing[comparisons], TString name, TString layer, TString labels[comparisons], TString templates, TString additional[comparisons],float xmin, float xmax,  float ymin, float ymax, bool isNorm)
{
  if(print) std::cout << " drawing " << std::endl;
   TDR();

  TCanvas *c;//= new TCanvas("c3","gaussian hist",10,10,700,900);                                                                                                                                           
  //  c = new TCanvas("c","c",50,50,600,600);
   c = new TCanvas("c","c",700,700);
  //  c->SetLeftMargin(0.2);                                                                                                                                                                                

  c->SetTopMargin(-0.3);//-0.15                                                                                                                                                                              
  c->SetBottomMargin(0.15);
  // c->SetRightMargin(-0.05);
  // c->SetLeftMargin(-0.2);
  c->SetRightMargin(-0.01);//tdr -0.01
  c->SetLeftMargin(-1.);

  
  if(isNorm){
    /*    
	  float scale = h[2]->GetEntries();  
	  h[2]->Scale(1./scale);
	  for(int i =0; i < comparisons-1; i++) {
	  h[i]->Scale(scale/h[i]->GetEntries());
	  }
    */
    int binMC = h[2]->FindFixBin(55);
    float scale = h[2]->GetBinContent(binMC);
    h[2]->Scale(1/scale);
    if(print) cout << scale << endl;
    for(int i =0; i < comparisons-1; i++) {                                                                                                                                                                                             
      h[i]->Scale(scale/h[i]->GetEntries());                                                                                                                                                                                              
    }                                                                                                                                                                                                                                   
  }

  int Nbins=h[0]->GetNbinsX();
  float depth[Nbins];
  float deptherror[Nbins];
  float mc[Nbins];

  float start[Nbins];
  float end[Nbins];

  for(int i = 0; i < Nbins; i++){
    depth[i] = h[0]->GetBinCenter(i);
    deptherror[i] = 0;
    mc[i]=h[2]->GetBinContent(i+1);
    start[i] = h[0]->GetBinContent(i+1);
    end[i] = h[1]->GetBinContent(i+1); 
 }



  TGraph *  grgreen = new TGraph(2*Nbins);
  for(int i = 0; i < Nbins; i++){
    grgreen->SetPoint(i,depth[i],start[i]);
    grgreen->SetPoint(Nbins+i,depth[Nbins-i-1],end[Nbins-i-1]);
  }




  gStyle->SetOptStat(0);
  // gPad->SetTickx();
  // gPad->SetTicky();
  if(print) std::cout << " canvas " << std::endl;
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetOptStat(0);
  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.9);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.2);
  gStyle->SetOptTitle(0);

  //  cout << Xaxis << hRec->GetTitle() << endl;                                                                                                                                                          
  TPad *grid = new TPad("grid","",0,0,1,1); 
  grid->Draw();
  grid->cd();
  grid->SetGrid();
  grid->SetFillStyle(4000); 

  TGraph * gr = new TGraph(Nbins, depth, mc);//, errx, rms);                                                                                                                                                                       
 
  //  h[0]->Draw();
  gPad->Update();
  TString YaxisLabel = "Average pixel charge [a.u.]";
  if(isNorm) YaxisLabel = "Normalized average pixel charge [a.u.]";
  TGaxis::SetMaxDigits(3);  

  h[2]->GetYaxis()->SetTitle(YaxisLabel);
  h[2]->GetXaxis()->SetTitle("Depth [#mum]");
  h[2]->GetXaxis()->SetLabelFont(42);
  h[2]->GetXaxis()->SetLabelSize(0.03);
  h[2]->GetXaxis()->SetTitleSize(0.05);
  h[2]->GetXaxis()->SetTitleOffset(0.8);
  h[2]->GetXaxis()->SetTitleFont(42);
  h[2]->GetYaxis()->SetLabelFont(42);
  h[2]->GetYaxis()->SetRangeUser(ymin,ymax);
  h[2]->GetXaxis()->SetRangeUser(xmin,xmax);
  h[2]->GetYaxis()->SetTitleSize(0.05);
  h[2]->GetYaxis()->SetTitleOffset(0.9);//0.8
  h[2]->GetYaxis()->SetTitleFont(42);
  h[2]->GetYaxis()->SetLabelSize(0.03);
  h[2]->SetLineWidth(2);
  h[2]->SetMarkerStyle(21);
  h[2]->SetMarkerColor(kRed);
  h[2]->SetLineColor(kRed);



  grgreen->GetYaxis()->SetTitle(YaxisLabel);
  grgreen->GetXaxis()->SetTitle("Depth [#mum]");
  grgreen->GetXaxis()->SetLabelFont(42);
  grgreen->GetXaxis()->SetLabelSize(0.03);
  grgreen->GetXaxis()->SetTitleSize(0.05);
  grgreen->GetXaxis()->SetTitleOffset(0.8);
  grgreen->GetXaxis()->SetTitleFont(42);
  grgreen->GetYaxis()->SetLabelFont(42);
  grgreen->GetYaxis()->SetRangeUser(ymin,ymax);
  grgreen->GetXaxis()->SetRangeUser(xmin,xmax);
  grgreen->GetYaxis()->SetTitleSize(0.05);
  grgreen->GetYaxis()->SetTitleOffset(0.9);//0.8
  grgreen->GetYaxis()->SetTitleFont(42);
  grgreen->GetYaxis()->SetLabelSize(0.03);

  grgreen->SetFillStyle(1001);
  grgreen->SetFillColor(kGray);

  grgreen->Draw("AF");
  h[2]->GetXaxis()->SetRangeUser(xmin,xmax);
  h[2]->SetMinimum(ymin);
  h[2]->SetMaximum(ymax);

  h[2]->Draw("EPsame");



  TLegend* leg4 = new TLegend(0.12,0.12,0.89,0.2);
  leg4->SetLineColor(0);
  leg4->SetTextSize(0.034);


  leg4->AddEntry(gr,labels[2],"ep"); // +", "+additional[i], "ep");

  leg4->Draw();



  TPaveText pt(0.6,0.82,0.85,0.89,"NDC");
  pt.SetTextSize(0.06);
  pt.AddText(layer);
  pt.SetBorderSize(0);
  pt.SetFillColor(kWhite);
  pt.Draw();  

   TDR2(c);

  TString plotname = "chargeprofiles_"+templates+"_"+name;  
  c->SaveAs(output_folder+plotname+".eps");
  c->SaveAs(output_folder+plotname+".pdf");
  c->SaveAs(output_folder+plotname+".png");
  c->SaveAs(output_folder+plotname+".root");
  if(print) std::cout << " saved " <<  std::endl;
  

}//drawnormMC                                                                                                                                                                     
void DrawNorm1(TH1F **h,bool drawing[comparisons], TString name, TString layer, TString labels[comparisons], TString templates, TString additional[comparisons],float xmin, float xmax,  float ymin, float ymax, bool isNorm)
{
  if(print) std::cout << " drawing " << std::endl;
   TDR();

  TCanvas *c;//= new TCanvas("c3","gaussian hist",10,10,700,900);                                                                                                                                           
  //  c = new TCanvas("c","c",50,50,600,600);
   c = new TCanvas("c","c",700,700);
  //  c->SetLeftMargin(0.2);                                                                                                                                                                                

  c->SetTopMargin(-0.3);//-0.15                                                                                                                                                                              
  c->SetBottomMargin(0.15);
  // c->SetRightMargin(-0.05);
  // c->SetLeftMargin(-0.2);
  c->SetRightMargin(-0.01);//tdr -0.01
  c->SetLeftMargin(-1.);


  if(isNorm){
    for(int i =0; i < comparisons; i++) {
      norm1(h[i]);
    }
  }

  int Nbins=h[0]->GetNbinsX();
  float depth[Nbins];
  float deptherror[Nbins];
  float mc[Nbins];

  float start[Nbins];
  float end[Nbins];

  for(int i = 0; i < Nbins; i++){
    depth[i] = h[0]->GetBinCenter(i);
    deptherror[i] = 0;
    mc[i]=h[2]->GetBinContent(i+1);
    start[i] = h[0]->GetBinContent(i+1);
    end[i] = h[1]->GetBinContent(i+1); 
 }



  TGraph *  grgreen = new TGraph(2*Nbins);
  for(int i = 0; i < Nbins; i++){
    grgreen->SetPoint(i,depth[i],start[i]);
    grgreen->SetPoint(Nbins+i,depth[Nbins-i-1],end[Nbins-i-1]);
  }




  gStyle->SetOptStat(0);
  // gPad->SetTickx();
  // gPad->SetTicky();
  if(print) std::cout << " canvas " << std::endl;
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetOptStat(0);
  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.9);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.2);
  gStyle->SetOptTitle(0);

  //  cout << Xaxis << hRec->GetTitle() << endl;                                                                                                                                                          
  TPad *grid = new TPad("grid","",0,0,1,1); 
  grid->Draw();
  grid->cd();
  grid->SetGrid();
  grid->SetFillStyle(4000); 

  TGraph * gr = new TGraph(Nbins, depth, mc);//, errx, rms);                                                                                                                                                                       
 
  //  h[0]->Draw();
  gPad->Update();
  TString YaxisLabel = "Average pixel charge [a.u.]";
  if(isNorm) YaxisLabel = "Normalized average pixel charge [a.u.]";
  TGaxis::SetMaxDigits(3);  

  h[2]->GetYaxis()->SetTitle(YaxisLabel);
  h[2]->GetXaxis()->SetTitle("Depth [#mum]");
  h[2]->GetXaxis()->SetLabelFont(42);
  h[2]->GetXaxis()->SetLabelSize(0.03);
  h[2]->GetXaxis()->SetTitleSize(0.05);
  h[2]->GetXaxis()->SetTitleOffset(0.8);
  h[2]->GetXaxis()->SetTitleFont(42);
  h[2]->GetYaxis()->SetLabelFont(42);
  h[2]->GetYaxis()->SetRangeUser(ymin,ymax);
  h[2]->GetXaxis()->SetRangeUser(xmin,xmax);
  h[2]->GetYaxis()->SetTitleSize(0.05);
  h[2]->GetYaxis()->SetTitleOffset(0.9);//0.8
  h[2]->GetYaxis()->SetTitleFont(42);
  h[2]->GetYaxis()->SetLabelSize(0.03);
  h[2]->SetLineWidth(2);
  h[2]->SetMarkerStyle(21);
  h[2]->SetMarkerColor(kRed);
  h[2]->SetLineColor(kRed);



  grgreen->GetYaxis()->SetTitle(YaxisLabel);
  grgreen->GetXaxis()->SetTitle("Depth [#mum]");
  grgreen->GetXaxis()->SetLabelFont(42);
  grgreen->GetXaxis()->SetLabelSize(0.03);
  grgreen->GetXaxis()->SetTitleSize(0.05);
  grgreen->GetXaxis()->SetTitleOffset(0.8);
  grgreen->GetXaxis()->SetTitleFont(42);
  grgreen->GetYaxis()->SetLabelFont(42);
  grgreen->GetYaxis()->SetRangeUser(ymin,ymax);
  grgreen->GetXaxis()->SetRangeUser(xmin,xmax);
  grgreen->GetYaxis()->SetTitleSize(0.05);
  grgreen->GetYaxis()->SetTitleOffset(0.9);//0.8
  grgreen->GetYaxis()->SetTitleFont(42);
  grgreen->GetYaxis()->SetLabelSize(0.03);

  grgreen->SetFillStyle(1001);
  grgreen->SetFillColor(kGray);

  grgreen->Draw("AF");
  h[2]->GetXaxis()->SetRangeUser(xmin,xmax);
  h[2]->SetMinimum(ymin);
  h[2]->SetMaximum(ymax);

  h[2]->Draw("EPsame");



  TLegend* leg4 = new TLegend(0.12,0.12,0.89,0.2);
  leg4->SetLineColor(0);
  leg4->SetTextSize(0.034);


  leg4->AddEntry(gr,labels[2],"ep"); // +", "+additional[i], "ep");

  leg4->Draw();



  TPaveText pt(0.6,0.82,0.85,0.89,"NDC");
  pt.SetTextSize(0.06);
  pt.AddText(layer);
  pt.SetBorderSize(0);
  pt.SetFillColor(kWhite);
  pt.Draw();  

   TDR2(c);

  TString plotname = "chargeprofiles_"+templates+"_"+name;  
  c->SaveAs(output_folder+plotname+".eps");
  c->SaveAs(output_folder+plotname+".pdf");
  c->SaveAs(output_folder+plotname+".png");
  if(print) std::cout << " saved " <<  std::endl;
  

}//drawnorm1

void DrawNorm55(TH1F **h,bool drawing[comparisons], TString name, TString layer, TString labels[comparisons], TString templates, TString additional[comparisons],float xmin, float xmax,  float ymin, float ymax, bool isNorm, float rescale =1.)
{
  if(print) std::cout << " drawing " << std::endl;
   TDR();

  TCanvas *c;//= new TCanvas("c3","gaussian hist",10,10,700,900);                                                                                                                                           
  //  c = new TCanvas("c","c",50,50,600,600);
   c = new TCanvas("c","c",700,700);
  //  c->SetLeftMargin(0.2);                                                                                                                                                                                

  c->SetTopMargin(-0.3);//-0.15                                                                                                                                                                              
  c->SetBottomMargin(0.15);
  // c->SetRightMargin(-0.05);
  // c->SetLeftMargin(-0.2);
  c->SetRightMargin(-0.01);//tdr -0.01
  c->SetLeftMargin(-1.);


  if(isNorm){
    for(int i =0; i < comparisons; i++) {
      norm55(h[i]);
    }
  }

  int Nbins=h[0]->GetNbinsX();
  float depth[Nbins];
  float deptherror[Nbins];
  float mc[Nbins];

  float start[Nbins];
  float end[Nbins];

  for(int i = 0; i < Nbins; i++){
    depth[i] = h[0]->GetBinCenter(i);
    deptherror[i] = 0;
    mc[i]=h[2]->GetBinContent(i+1);
    start[i] = h[0]->GetBinContent(i+1)*rescale;
    end[i] = h[1]->GetBinContent(i+1)*rescale; 
 }



  TGraph *  grgreen = new TGraph(2*Nbins);
  for(int i = 0; i < Nbins; i++){
    grgreen->SetPoint(i,depth[i],start[i]);
    grgreen->SetPoint(Nbins+i,depth[Nbins-i-1],end[Nbins-i-1]);
  }




  gStyle->SetOptStat(0);
  // gPad->SetTickx();
  // gPad->SetTicky();
  if(print) std::cout << " canvas " << std::endl;
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetOptStat(0);
  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.9);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.2);
  gStyle->SetOptTitle(0);

  //  cout << Xaxis << hRec->GetTitle() << endl;                                                                                                                                                          
  TPad *grid = new TPad("grid","",0,0,1,1); 
  grid->Draw();
  grid->cd();
  grid->SetGrid();
  grid->SetFillStyle(4000); 

  TGraph * gr = new TGraph(Nbins, depth, mc);//, errx, rms);                                                                                                                                                                       
 
  //  h[0]->Draw();
  gPad->Update();
  TString YaxisLabel = "Average pixel charge [a.u.]";
  if(isNorm) YaxisLabel = "Normalized average pixel charge [a.u.]";
  TGaxis::SetMaxDigits(3);  

  h[2]->GetYaxis()->SetTitle(YaxisLabel);
  h[2]->GetXaxis()->SetTitle("Depth [#mum]");
  h[2]->GetXaxis()->SetLabelFont(42);
  h[2]->GetXaxis()->SetLabelSize(0.03);
  h[2]->GetXaxis()->SetTitleSize(0.05);
  h[2]->GetXaxis()->SetTitleOffset(0.8);
  h[2]->GetXaxis()->SetTitleFont(42);
  h[2]->GetYaxis()->SetLabelFont(42);
  h[2]->GetYaxis()->SetRangeUser(ymin,ymax);
  h[2]->GetXaxis()->SetRangeUser(xmin,xmax);
  h[2]->GetYaxis()->SetTitleSize(0.05);
  h[2]->GetYaxis()->SetTitleOffset(0.9);//0.8
  h[2]->GetYaxis()->SetTitleFont(42);
  h[2]->GetYaxis()->SetLabelSize(0.03);
  h[2]->SetLineWidth(2);
  h[2]->SetMarkerStyle(21);
  h[2]->SetMarkerColor(kRed);
  h[2]->SetLineColor(kRed);



  grgreen->GetYaxis()->SetTitle(YaxisLabel);
  grgreen->GetXaxis()->SetTitle("Depth [#mum]");
  grgreen->GetXaxis()->SetLabelFont(42);
  grgreen->GetXaxis()->SetLabelSize(0.03);
  grgreen->GetXaxis()->SetTitleSize(0.05);
  grgreen->GetXaxis()->SetTitleOffset(0.8);
  grgreen->GetXaxis()->SetTitleFont(42);
  grgreen->GetYaxis()->SetLabelFont(42);
  grgreen->GetYaxis()->SetRangeUser(ymin,ymax);
  grgreen->GetXaxis()->SetRangeUser(xmin,xmax);
  grgreen->GetYaxis()->SetTitleSize(0.05);
  grgreen->GetYaxis()->SetTitleOffset(0.9);//0.8
  grgreen->GetYaxis()->SetTitleFont(42);
  grgreen->GetYaxis()->SetLabelSize(0.03);

  grgreen->SetFillStyle(1001);
  grgreen->SetFillColor(kGray);

  grgreen->Draw("AF");
  h[2]->GetXaxis()->SetRangeUser(xmin,xmax);
  h[2]->SetMinimum(ymin);
  h[2]->SetMaximum(ymax);

  h[2]->Draw("EPsame");



  TLegend* leg4 = new TLegend(0.12,0.12,0.89,0.25);
  leg4->SetLineColor(0);
  leg4->SetTextSize(0.034);


  leg4->AddEntry(h[2],labels[2],"ep"); // +", "+additional[i], "ep");
  leg4->AddEntry(grgreen,"Data 2018","f"); // +", "+additional[i], "ep");
  leg4->AddEntry(grgreen,"run 315257 (Era A) - run 325175 (Era D)",""); // +", "+additional[i], "ep");

  leg4->Draw();



  TPaveText pt(0.6,0.82,0.85,0.89,"NDC");
  pt.SetTextSize(0.06);
  pt.AddText(layer);
  pt.SetBorderSize(0);
  pt.SetFillColor(kWhite);
  pt.Draw();  

   TDR2(c);

  TString plotname = "chargeprofiles_"+templates+"_"+name;  
  c->SaveAs(output_folder+plotname+".eps");
  c->SaveAs(output_folder+plotname+".pdf");
  c->SaveAs(output_folder+plotname+".png");
  if(print) std::cout << " saved " <<  std::endl;
  

}//drawnorm55
                                                                                                                                                                                            
void DrawNorm55ratio(TH1F **h,bool drawing[comparisons], TString name, TString layer, TString labels[comparisons], TString templates, TString additional[comparisons],float xmin, float xmax,  float ymin, float ymax, bool isNorm)
{
  if(print) std::cout << " drawing " << std::endl;
   TDR();

  TCanvas *c;
  c = new TCanvas("c","c",700,700);
  c->SetTopMargin(-0.3);//-0.15                                                                                                                                                                              
  c->SetBottomMargin(0.15);
  c->SetRightMargin(-0.01);//tdr -0.01
  c->SetLeftMargin(-1.);

  gStyle->SetOptStat(0);

  // Upper plot will be in pad1                                                                                                                                                                                                                                                
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.23, 1, 1.0);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined                                                                                                                                                                                                               
  pad1->Draw();             // Draw the upper pad: pad1                                                                                                                                                                                                                       
  pad1->cd();               // pad1 becomes the current pad                                                                                                                                                                                                                    
  pad1->SetTickx();
  pad1->SetTicky();

  if(print) std::cout << " canvas " << std::endl;
  gROOT->SetStyle("Plain");
  gStyle->SetPadGridX(0);
  gStyle->SetPadGridY(0);
  gStyle->SetOptStat(0);
  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.9);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.2);
  gStyle->SetOptTitle(0);

  //  cout << Xaxis << hRec->GetTitle() << endl;                                                                                                                                                          
  TPad *grid = new TPad("grid","",0,0,1,1); 
  grid->Draw();
  grid->cd();
  grid->SetGrid();
  grid->SetFillStyle(4000); 
   
  //  h[0]->Draw();
  gPad->Update();
  TString YaxisLabel = "Average pixel charge [a.u.]";
  if(isNorm) YaxisLabel = "Normalized average pixel charge [a.u.]";
  TGaxis::SetMaxDigits(3);  

  for(int i =0; i < comparisons; i++)
    {
      if(isNorm)	norm55(h[i]);
      h[i]->GetYaxis()->SetTitle(YaxisLabel);
      //      h[i]->GetXaxis()->SetTitle("Depth [#mum]");
      //      h[i]->GetXaxis()->SetLabelFont(42);
      //h[i]->GetXaxis()->SetLabelSize(0.03);
      //h[i]->GetXaxis()->SetTitleSize(0.05);
      //h[i]->GetXaxis()->SetTitleOffset(0.8);
      //h[i]->GetXaxis()->SetTitleFont(42);
      h[i]->GetYaxis()->SetLabelFont(42);
      h[i]->GetYaxis()->SetRangeUser(ymin,ymax);
      h[i]->GetXaxis()->SetRangeUser(xmin,xmax);
      h[i]->GetYaxis()->SetTitleSize(0.05);
      h[i]->GetYaxis()->SetTitleOffset(0.9);//0.8
      h[i]->GetYaxis()->SetTitleFont(42);
      h[i]->GetYaxis()->SetLabelSize(0.04);
      h[i]->SetLineWidth(2);
      h[i]->SetMarkerStyle(20+i);
      if(i==1)      h[i]->SetMarkerStyle(25);
      h[i]->SetMarkerColor(kRed-1+comparisons-i);
      h[i]->SetLineColor(kRed-1+comparisons-i);
      if(i==0)
        {
          h[i]->SetMarkerColor(kBlack);
          h[i]->SetLineColor(kBlack);
        }
  
    if(drawing[i])      h[i]->Draw("epsame");

    }

  TLegend* leg5 = new TLegend(0.2,0.41,0.75,0.47);
  leg5->SetLineColor(0);
  leg5->SetTextSize(0.034);
  leg5->SetMargin(0.);
  leg5->AddEntry(h[0],"#int L dt = 29.4 (Run-1) + 4.2 (2015) + 20.1 (2016) fb^{-1}",""); 
  //leg5->AddEntry(h[0],"#int L dt = 0.05 (2010) + 6.10 (2011) + 23.30 (2012) + 4.21 (2015) + 20.12 (2016) fb^{-1}",""); 
  //leg5->AddEntry(h[0],"#int L dt = 0.05 (2010) + 6.10 (2011) + 23.30 (2012) + 4.21 (2015) + 20.12 (2016) fb^{-1}",""); 
  leg5->Draw();
  TLegend* leg4 = new TLegend(0.2,0.2,0.75,0.4);
  leg4->SetLineColor(0);
  leg4->SetTextSize(0.034);
  leg4->SetMargin(0.15);
  //leg4->AddEntry(h[0],"#int L dt = 29.4 (Run-1) + 4.2 (2015) + 20.1 (2016) fb^{-1}",""); 
  //  leg4->AddEntry(h[0],"#int L dt = 53.8 fb^{-1}",""); 
  for(int i =0; i < comparisons; i++)
    {
      if(drawing[i]) leg4->AddEntry(h[i],labels[i],"ep"); // +", "+additional[i], "ep");
    }
  leg4->Draw();



  TPaveText pt(0.6,0.82,0.8,0.89,"NDC");
  pt.SetTextSize(0.06);
  pt.AddText(layer);
  pt.SetBorderSize(0);
  pt.SetFillColor(kWhite);
  pt.Draw();  

   TDR2(c);
  c->cd();          // Go back to the main canvas before defining pad2                                                                                                                                                                                                  
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.3);//0.2);                                                                                                                                                                                                                                        
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad                                                                                                                                                                                                                         

  // Define the ratio plot                                                                                                                                                                                                                                                  

  for(int i =1; i<comparisons; i++)
    {
      TH1F *h_ratio_dataMC = (TH1F*)h[0]->Clone("h_ratio_dataMC");

      h_ratio_dataMC->SetMarkerColor(kRed-1+comparisons-i);
      h_ratio_dataMC->SetLineColor(kRed-1+comparisons-i);
      h_ratio_dataMC->SetMinimum(0.5);  // Define Y ..                                                                                                                                                                                                                    
      h_ratio_dataMC->SetMaximum(1.5); // .. range                                                                                                                                                                                                                        
      h_ratio_dataMC->Sumw2();
      h_ratio_dataMC->SetStats(0);      // No statistics on lower plot                                                                                                                                                                                                    
      h_ratio_dataMC->Divide(h[i]);
      h_ratio_dataMC->SetMarkerStyle(20+i);
      if(i==1) h_ratio_dataMC->SetMarkerStyle(25);
      h_ratio_dataMC->Draw("epsame");       // Draw the ratio plot                                                                                                                                                                                                            
      // Y axis ratio plot settings                                                                                                                                                                                                                                             
      h_ratio_dataMC->GetYaxis()->SetTitle("data/MC ");
      h_ratio_dataMC->GetYaxis()->SetNdivisions(505);
      h_ratio_dataMC->GetYaxis()->SetTitleSize(0.1);
      h_ratio_dataMC->GetYaxis()->SetTitleFont(42);
      h_ratio_dataMC->GetYaxis()->SetTitleOffset(0.5);
      h_ratio_dataMC->GetYaxis()->SetLabelFont(42); // Absolute font size in pixel (precision 3)                                                                                                                                                                          
      h_ratio_dataMC->GetYaxis()->SetLabelSize(0.1);
      
      // X axis ratio plot settings                                                                                                                                                                                                                                             
      h_ratio_dataMC->GetXaxis()->SetTitleSize(0.11);
      h_ratio_dataMC->GetXaxis()->SetTitleFont(42);
      h_ratio_dataMC->GetXaxis()->SetTitleOffset(0.9);
      h_ratio_dataMC->GetXaxis()->SetLabelFont(42); // Absolute font size in pixel (precision 3)                                                                                                                                                                          
      h_ratio_dataMC->GetXaxis()->SetLabelSize(0.1);
      h_ratio_dataMC->GetXaxis()->SetTitle("Depth [#mum]");
    }
  
  TLine *line = new TLine(xmin,1,xmax,1);
  line->SetLineColor(kRed);
  line->Draw("same");
      


  TString plotname = "chargeprofiles_"+templates+"_"+name;  
  c->SaveAs(output_folder+plotname+"_ratio.eps");
  c->SaveAs(output_folder+plotname+"_ratio.pdf");
  c->SaveAs(output_folder+plotname+"_ratio.png");
  c->SaveAs(output_folder+plotname+"_ratio.root");
  c->SaveAs(output_folder+plotname+"_ratio.C");
  if(print) std::cout << " saved " <<  std::endl;
  

}//drawfournorm55
                                                                                                                                                                                            



void TDR(TString lumi)
{
  setTDRStyle();
  writeExtraText = true;       // if extra text                                                                                                                                                             
  extraText  = "Internal";  // default extra text is "Preliminary"                                                                                                                                        
  lumi_sqrtS = lumi;//"Integ. Lumi. =  ";
  bool drawLogo = true;
  int H_ref = 600;
  int W_ref = 600;
  float T = 0.07*H_ref;//0.07                                                                                                                                                                                
  float B = 0.11*H_ref;//0.12                                                                                                                                                                                
  float L = 0.12*W_ref;
  float R = 0.01*W_ref;


}

void TDR2(TCanvas * c_all)
{
  //    CMS_lumi( c_all, 4, 11);
  CMS_lumi( c_all, 4, 0);
  //  CMS_lumi( c_all, 0, 11);
  c_all->Update();
  c_all->RedrawAxis();
}


void norm55(TH1F *h)
{
  int bin = h->FindFixBin(55);
  float scale = h->GetBinContent(bin);
  h->Scale(1/scale);
  if(print) cout << scale << endl;

}

void norm1(TH1F *h)
{
  //int bin = h->FindFixBin(55);
  float scale = h->GetEntries();
  h->Scale(1/scale);
  if(print) cout << scale << endl;

}
