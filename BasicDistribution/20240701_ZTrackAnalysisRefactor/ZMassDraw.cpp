#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TLine.h>
#include <cmath>
#include <iostream>
#include <vector>

//#include "CommandLine.h"

using namespace std;

void style(){

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLabelFont(43,"xyz");
  gStyle->SetTitleFont(43);
  gStyle->SetTitleFont(43,"xyz");
  gStyle->SetStatFont(43);

  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetOptStat(0); /*don't show statistics box*/
  gStyle->SetOptTitle(0); /*don't show histogram titles*/
  gStyle->SetTitleSize(24, "xyz");
  gStyle->SetTitleOffset(1, "xz");
  gStyle->SetTitleOffset(1.5, "y");
  gStyle->SetLabelSize(18, "xyz");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFillColor(kWhite);

  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.05);

  gStyle->SetLineScalePS(1.5);

  gROOT->ForceStyle();
}

std::vector<int> GetPrimaryColors();
std::vector<int> GetCVDColors6();
std::vector<int> GetCVDColors8();
std::vector<int> GetCVDColors10();

void DrawRatioPlot(TH1D *h_deno, std::vector<TH1D *> h_nums, std::vector<string> LegendText, std::vector<TLatex> pts, 
   string XTitle, string YTitle, string RTitle, float XMin, float XMax, float YMin, float YMax, float RMax, float RMin, 
   bool isLog, bool isRatio, string OutputBase, string OutputName );

int main(int argc, char *argv[]){

   style();

   // Prepare for histograms from skim.

   //TChain *n_Zee_MC   = new TChain("n_Zee");
   //TChain *n_Zee_Data = new TChain("n_Zee");

   string filebase = "/eos/cms/store/group/phys_heavyions/pchou/SkimZHadron2024/"; 

   //n_Zee_MC->Add((filebase + "OutputMC_v1b_ee/*.root" ).c_str());
   //n_Zee_Data->Add((filebase + "OutputData_v1b_ee/*.root" ).c_str());


   TChain *Tree_MC   = new TChain("Tree");
   TChain *Tree_Data = new TChain("Tree");

   Tree_MC->Add((filebase + "OutputMC_v1b_ee_eidcorrected/*.root" ).c_str());
   Tree_Data->Add((filebase + "OutputData_v1b_ee_eidcorrected/*.root" ).c_str());

   TH1D* hMC = new TH1D("hMC", "Z candidate mass", 30, 60, 120);
   TH1D* hData = new TH1D("hData", "Z candidate mass", 30, 60, 120);

   Tree_MC->Draw("zMass>>hMC","NCollWeight*(zPt>30)","goff");
   Tree_Data->Draw("zMass>>hData","NCollWeight*(zPt>30)","goff");

   TH1F HNSig("HNSig","Normalization", 1, 0, 1);
   TH1F HNBkg("HNBkg","Normalization", 1, 0, 1);

   Tree_MC->Draw("0>>HNSig", "NCollWeight*(zPt>30)","goff");
   Tree_Data->Draw("0>>HNBkg", "NCollWeight*(zPt>30)","goff");

   //float MC_ents = hMC->GetEntries();
   //float Data_ents = hData->GetEntries();

   float MC_ents = HNSig.GetBinContent(1);
   float Data_ents = HNBkg.GetBinContent(1);

   hMC->Scale(1./MC_ents*Data_ents);
   //hData->Scale(1.);

/*
   // Prepare for histograms.

   string HistName = "HZMass";
   string filebase = "/eos/home-p/pchou/BasicPlots/ZHadron2024/";

   TFile *file_sigMC   = TFile::Open((filebase + "GraphMCSig_v1_ee_noZw.root" ).c_str(), "read");
   TFile *file_sigData = TFile::Open((filebase + "GraphDataSig_v1_ee_noZw.root"   ).c_str(), "read");

   float ptL=30, ptH=200, centL=0, centH=100, TptL=1, TptH=1000;

   std::string FolderName = Form("Plot_ZPT_%.0f_%.0f_Cent_%.0f_%.0f_TrackPT_%.2f_%.2f",ptL,ptH,centL,centH,TptL,TptH);
   std::replace(FolderName.begin(), FolderName.end(), '.', 'p');

   TH1D* hMC = (TH1D*) file_sigMC->Get(Form("%s/%s", FolderName.c_str(), HistName.c_str()));
   TH1D* hData = (TH1D*) file_sigData->Get(Form("%s/%s", FolderName.c_str(), HistName.c_str()));

   TNamed *NMC   = (TNamed *) file_sigMC->Get(Form("%s/EntryCount", FolderName.c_str()));
   TNamed *NData = (TNamed *) file_sigData->Get(Form("%s/EntryCount", FolderName.c_str()));

   // These are "Ncoll*Zw*VZw" entries.
   float IMC   = atof(NMC->GetTitle());
   float IData = atof(NData->GetTitle());

   //HZMass: 100, 0, 150

   int rebinnum = 1;
   int bin_width = hMC->GetBinWidth(1); // GeV.

   hMC->Rebin(rebinnum);
   hData->Rebin(rebinnum);
   
   bin_width*=rebinnum;
   bin_width = 1;
   
   float MC_ents = hMC->GetEntries();
   float Data_ents = hData->GetEntries();

   std::cout << "MC_ents = "<<MC_ents<<", Data_ents = "<<Data_ents<<std::endl;
   std::cout << "IMC = "<<IMC<<", IData = "<<IData<<std::endl;

   hMC->Scale(1./bin_width/IMC*IData);
   hData->Scale(1./bin_width);
   */




   std::cout<<"max = "<<hMC->GetMaximum()<<", "<<hData->GetMaximum()<<std::endl;

   std::vector<TH1D *> h_nums;
   h_nums.push_back(hData);

   // Prepare for texts.

   std::vector<string> LegendText;
   LegendText.push_back("Pythia+Hydjet MC");
   LegendText.push_back("PbPb Run2 data");

   std::vector<TLatex> pts;

   TLatex pt0(0.25,0.88,"Z #rightarrow e^{+}e^{-}");
   pt0.SetTextFont(42);
   pt0.SetTextSize(0.04);
   pt0.SetNDC(kTRUE);

   TLatex pt(0.25,0.82,"Cent: 0-100 %");
   pt.SetTextFont(42);
   pt.SetTextSize(0.04);
   pt.SetNDC(kTRUE);

   TLatex pt2(0.25,0.76,"p^{e}_{T} > 20 GeV/c");
   pt2.SetTextFont(42);
   pt2.SetTextSize(0.04);
   pt2.SetNDC(kTRUE);

   TLatex pt3(0.25,0.70,"p^{Z}_{T} > 30 GeV/c");
   pt3.SetTextFont(42);
   pt3.SetTextSize(0.04);
   pt3.SetNDC(kTRUE);

   pts.push_back(pt0);
   pts.push_back(pt);
   pts.push_back(pt2);
   pts.push_back(pt3);

   TLatex pt1(0.99,0.99,"Run2 PbPb, #sqrt{s_{NN}} = 5.02 TeV, 1618 #mub^{-1}");
   pt1.SetTextFont(42);
   pt1.SetTextSize(0.03);
   pt1.SetTextAlign(33);//3 = right/top, 2 = center, 1 = left/bottom
   pt1.SetNDC(kTRUE);

   pts.push_back(pt1);

   string XTitle = "M^{ee} [GeV/c^{2}]";
   string YTitle = "Entries / [2 GeV/c^{2}]";
   string RTitle = "Data / MC";
   float XMin = 60,  XMax = 120,  YMin = 0,  YMax = 1.5* hMC->GetMaximum(),  RMin = 0.5,  RMax = 1.5;

   bool isLog = false, isRatio = true;

   string OutputBase = "/eos/user/p/pchou/figs/ZHadron2024/ZMass/ov1_v1b_Reco_noZw_eidcorrected/20240705/";
   string OutputName = "ZMass_ee_PbPb_run2_fromskim.png";

   DrawRatioPlot(hMC, h_nums, LegendText, pts, XTitle, YTitle, RTitle, XMin, XMax, YMin, YMax, RMin, RMax,
      isLog, isRatio, OutputBase, OutputName);

}

void DrawRatioPlot(TH1D *h_deno, std::vector<TH1D *> h_nums, std::vector<string> LegendText, std::vector<TLatex> pts, 
   string XTitle, string YTitle, string RTitle, float XMin, float XMax, float YMin, float YMax, float RMin, float RMax, 
   bool isLog, bool isRatio, string OutputBase, string OutputName ){

   vector<int> Colors = GetCVDColors8();
   int marks[10] = {20,21,33,24,25,27,22,23,26,32};

   TCanvas *c = new TCanvas("c", "", 500, 600);
   TPad *Pad  = new TPad("Pad" , "", 0, 0.3, 1,   1);
   TPad *RPad = new TPad("RPad", "", 0, 0. , 1, 0.3);

   //c->SetRightMargin(0);
   //Pad->SetRightMargin(0);
   //RPad->SetRightMargin(0);

   //RPad->SetTopMargin(0);
   //Pad->SetBottomMargin(0);

   Pad->Draw();
   RPad->Draw();

   if(isLog)
      Pad->SetLogy();

   Pad->cd();
   h_deno->SetFillColor(kOrange);
   h_deno->SetLineColor(kRed);
   
   h_deno->SetMaximum(YMax);
   h_deno->SetMinimum(YMin);
   h_deno->SetYTitle(YTitle.c_str());
   h_deno->SetXTitle(XTitle.c_str());
   //Pad->SetPad(XMin,XMax,YMin,YMax);
   h_deno->GetXaxis()->SetRangeUser(XMin,XMax);
   h_deno->Draw("hist");

   for(int iH = 0; iH < h_nums.size(); iH++){
      h_nums[iH]->SetLineColor(Colors[iH]);
      h_nums[iH]->SetMarkerColor(Colors[iH]);
      h_nums[iH]->SetMarkerStyle(marks[iH]);
      h_nums[iH]->SetMarkerSize(0.5);
      h_nums[iH]->GetXaxis()->SetRangeUser(XMin,XMax);
      h_nums[iH]->Draw("ep same");
   }

   TLegend Legend(0.58, 0.78, 0.98, 0.90);
   Legend.AddEntry(h_deno, LegendText[0].c_str(), "f");
   for(int iH = 0; iH < h_nums.size(); iH++){
      Legend.AddEntry(h_nums[iH], LegendText[iH+1].c_str(), "lep");
   }

   Legend.SetFillStyle(0);
   Legend.SetLineColor(kBlack);
   Legend.SetLineWidth(1);
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.035);
   Legend.Draw();

   for(int iP = 0; iP < pts.size(); iP++){
      pts[iP].Draw();
   }

   RPad->cd();

   int y_line = isRatio ? 1 : 0;

   TLine horiz_line(XMin,y_line,XMax,y_line);
   horiz_line.SetNDC(kFALSE);
   horiz_line.SetLineColor(kBlack);
   horiz_line.SetLineStyle(kDashed);
   
/*
   TH1D *horiz_line = (TH1D*) h_deno->Clone("horiz_line");

   if(isRatio)
      horiz_line->Divide(h_deno);
   else
      horiz_line->Add(h_deno,-1);
   
   horiz_line->SetLineColor(kBlack);
   horiz_line->SetLineStyle(kDashed);
   
   horiz_line->SetMaximum(RMax);
   horiz_line->SetMinimum(RMin);
   horiz_line->SetYTitle(RTitle.c_str());
   horiz_line->SetXTitle("");
   horiz_line->Draw("lX0");
   */

   std::vector<TH1D *> h_rats;

   for(int iH = 0; iH < h_nums.size(); iH++){
      h_rats.push_back( (TH1D*) h_nums[iH]->Clone(Form("h_rats%d",iH)) );
      if(isRatio)
         h_rats[iH]->Divide(h_deno);
      else
         h_rats[iH]->Add(h_deno,-1);

      h_rats[iH]->SetMaximum(RMax);
      h_rats[iH]->SetMinimum(RMin);
      h_rats[iH]->SetYTitle(RTitle.c_str());
      h_rats[iH]->SetMarkerColor(Colors[iH]);
      h_rats[iH]->SetMarkerStyle(marks[iH]);
      h_rats[iH]->SetMarkerSize(0.5);
      if(iH==0)
         h_rats[iH]->Draw("ep");
      else
         h_rats[iH]->Draw("ep same");
   }

   horiz_line.Draw("same");
   gSystem->Exec(("mkdir -p " + OutputBase).c_str());
   c->SaveAs((OutputBase + OutputName).c_str());
   c->Clear();
}



std::vector<int> GetPrimaryColors()
{
   static std::vector<int> Colors;
   if(Colors.size() > 0)
      return Colors;
   //Colors.push_back(TColor::GetColor("#E74C3C")); // Alizarin
   //Colors.push_back(TColor::GetColor("#3498DB")); // Peter River
   //Colors.push_back(TColor::GetColor("#F1C40F")); // Sum Flower
   //Colors.push_back(TColor::GetColor("#2ECC71")); // Emerald
   //Colors.push_back(TColor::GetColor("#7F8C8D")); // Asbestos
   //Colors.push_back(TColor::GetColor("#8E44AD")); // Wisteria
   //Colors.push_back(TColor::GetColor("#2C3E50")); // Green Sea (dark)
   //Colors.push_back(TColor::GetColor("#16A085")); // Green Sea
   //Colors.push_back(TColor::GetColor("#E67E22")); // Carrot

   Colors.push_back(TColor::GetColor("#377eb8")); // (Blue)
   Colors.push_back(TColor::GetColor("#ff7f00")); // (Orange)
   Colors.push_back(TColor::GetColor("#4daf4a")); // (Green)
   Colors.push_back(TColor::GetColor("#f781bf")); // (Pink)
   Colors.push_back(TColor::GetColor("#a65628")); // (Brown)
   Colors.push_back(TColor::GetColor("#984ea3")); // (Purple)
   Colors.push_back(TColor::GetColor("#e41a1c")); // (Red)

   return Colors;
}

std::vector<int> GetCVDColors6()
{
   static std::vector<int> Colors;
   if(Colors.size() > 0)
      return Colors;

   std::string ColorStrings[6] = {"#5790fc", "#f89c20", "#e42536", "#964a8b", "#9c9ca1", "#7a21dd"};
   for(int i = 0; i < 6; i++)
      Colors.push_back(TColor::GetColor(ColorStrings[i].c_str()));

   return Colors;
}

std::vector<int> GetCVDColors8()
{
   static std::vector<int> Colors;
   if(Colors.size() > 0)
      return Colors;

   std::string ColorStrings[8] = {"#1845fb", "#ff5e02", "#c91f16", "#c849a9", "#adad7d", "#86c8dd", "#578dff", "#656364"};
   for(int i = 0; i < 8; i++)
      Colors.push_back(TColor::GetColor(ColorStrings[i].c_str()));

   return Colors;
}

std::vector<int> GetCVDColors10()
{
   static std::vector<int> Colors;
   if(Colors.size() > 0)
      return Colors;

   std::string ColorStrings[10] = {"#3f90da", "#ffa90e", "#bd1f01", "#94a4a2", "#832db6", "#a96b59", "#e76300", "#b9ac70", "#717581", "#92dadd"};
   for(int i = 0; i < 10; i++)
      Colors.push_back(TColor::GetColor(ColorStrings[i].c_str()));

   return Colors;
}