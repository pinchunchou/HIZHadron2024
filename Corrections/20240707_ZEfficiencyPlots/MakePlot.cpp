#include <iostream>
using namespace std;

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include <TSystem.h>

#include "SetStyle.h"

int main();
void MakePlot(double XMin, double XMax, double YMin, double YMax, string Title,
   TProfile &P1, TProfile &P2, bool LogX, string Output);

string OutputBase = "/eos/home-p/pchou/figs/ZHadron2024/ZeeEfficiency/20240905/vF/";

int main()
{
   SetThumbStyle();

   TFile FPbPb("../20240706_ZeeEfficiency/ZEfficiencyPbPb.root");
   TFile FPP("../20240706_ZeeEfficiency/ZEfficiencyPP.root");

   //TTree *TPbPb = (TTree *)FPbPb.Get("Tree");
   //TTree *TPP = (TTree *)FPP.Get("Tree");

   TChain* TPbPb = new TChain("Tree"); TPbPb->Add("../20240706_ZeeEfficiency/PbPb_0905_vF/*.root");
   TChain* TPP = new TChain("Tree"); TPP->Add("/eos/cms/store/group/phys_heavyions/pchou/HIZHadron2024/Corrections/20240706_ZeeEfficiency/pp_0723_vF/*.root");

   const int binnum = 12;

   double PTBins[binnum+1] = {0};
   //for(int i = 0; i <= 30; i++)
   //   PTBins[i] = exp(log(1) + (log(200) - log(1)) / 30 * i);

   for(int i = 0; i <= binnum; i++)
      PTBins[i] = exp(log(40) + (log(200) - log(40)) / binnum * i);

   TProfile HPbPbMCYRaw            ("HPbPbMCYRaw",             ";;", binnum, -2.1, 2.1);
   TProfile HPbPbMCYCorrected      ("HPbPbMCYCorrected",       ";;", binnum, -2.1, 2.1);
   TProfile HPbPbMCPTRaw           ("HPbPbMCPTRaw",            ";;", binnum, PTBins);
   TProfile HPbPbMCPTCorrected     ("HPbPbMCPTCorrected",      ";;", binnum, PTBins);
   TProfile HPbPbMCHiBinRaw        ("HPbPbMCHiBinRaw",         ";;", binnum, 0, 100);
   TProfile HPbPbMCHiBinCorrected  ("HPbPbMCHiBinCorrected",   ";;", binnum, 0, 100);
   TProfile HPbPbDataYRaw          ("HPbPbDataYRaw",           ";;", binnum, -2.1, 2.1);
   TProfile HPbPbDataYCorrected    ("HPbPbDataYCorrected",     ";;", binnum, -2.1, 2.1);
   TProfile HPbPbDataPTRaw         ("HPbPbDataPTRaw",          ";;", binnum, PTBins);
   TProfile HPbPbDataPTCorrected   ("HPbPbDataPTCorrected",    ";;", binnum, PTBins);

   TProfile HPbPbMCYRaw_cent1            ("HPbPbMCYRaw_cent1",             ";;", binnum, -2.1, 2.1);
   TProfile HPbPbMCYCorrected_cent1      ("HPbPbMCYCorrected_cent1",       ";;", binnum, -2.1, 2.1);
   TProfile HPbPbMCPTRaw_cent1           ("HPbPbMCPTRaw_cent1",            ";;", binnum, PTBins);
   TProfile HPbPbMCPTCorrected_cent1     ("HPbPbMCPTCorrected_cent1",      ";;", binnum, PTBins);

   TProfile HPbPbMCYRaw_pher1            ("HPbPbMCYRaw_pher1",             ";;", binnum, -2.1, 2.1);
   TProfile HPbPbMCYCorrected_pher1      ("HPbPbMCYCorrected_pher1",       ";;", binnum, -2.1, 2.1);
   TProfile HPbPbMCPTRaw_pher1           ("HPbPbMCPTRaw_pher1",            ";;", binnum, PTBins);
   TProfile HPbPbMCPTCorrected_pher1     ("HPbPbMCPTCorrected_pher1",      ";;", binnum, PTBins);

   TProfile HPbPbMCYRaw_cent2            ("HPbPbMCYRaw_cent2",             ";;", binnum, -2.1, 2.1);
   TProfile HPbPbMCYCorrected_cent2      ("HPbPbMCYCorrected_cent2",       ";;", binnum, -2.1, 2.1);
   TProfile HPbPbMCPTRaw_cent2           ("HPbPbMCPTRaw_cent2",            ";;", binnum, PTBins);
   TProfile HPbPbMCPTCorrected_cent2     ("HPbPbMCPTCorrected_cent2",      ";;", binnum, PTBins);

   TProfile HPbPbMCYRaw_pher2            ("HPbPbMCYRaw_pher2",             ";;", binnum, -2.1, 2.1);
   TProfile HPbPbMCYCorrected_pher2      ("HPbPbMCYCorrected_pher2",       ";;", binnum, -2.1, 2.1);
   TProfile HPbPbMCPTRaw_pher2           ("HPbPbMCPTRaw_pher2",            ";;", binnum, PTBins);
   TProfile HPbPbMCPTCorrected_pher2     ("HPbPbMCPTCorrected_pher2",      ";;", binnum, PTBins);


   TProfile HPbPbDataHiBinRaw      ("HPbPbDataHiBinRaw",       ";;", binnum, 0, 100);
   TProfile HPbPbDataHiBinCorrected("HPbPbDataHiBinCorrected", ";;", binnum, 0, 100);
   TProfile HPPMCYRaw              ("HPPMCYRaw",               ";;", binnum, -2.1, 2.1);
   TProfile HPPMCYCorrected        ("HPPMCYCorrected",         ";;", binnum, -2.1, 2.1);
   TProfile HPPMCPTRaw             ("HPPMCPTRaw",              ";;", binnum, PTBins);
   TProfile HPPMCPTCorrected       ("HPPMCPTCorrected",        ";;", binnum, PTBins);
   TProfile HPPDataYRaw            ("HPPDataYRaw",             ";;", binnum, -2.1, 2.1);
   TProfile HPPDataYCorrected      ("HPPDataYCorrected",       ";;", binnum, -2.1, 2.1);
   TProfile HPPDataPTRaw           ("HPPDataPTRaw",            ";;", binnum, PTBins);
   TProfile HPPDataPTCorrected     ("HPPDataPTCorrected",      ";;", binnum, PTBins);

   TPbPb->Draw("HasReco:Y>>HPbPbMCYRaw_cent1",                             "PT > 40 && PT < 200 && HiBin <= 20", "prof");
   TPbPb->Draw("HasReco*TreeWPbPbMC:Y>>HPbPbMCYCorrected_cent1",                     "PT > 40 && PT < 200 && HiBin <= 20", "prof");
   TPbPb->Draw("HasReco:PT>>HPbPbMCPTRaw_cent1",                           "PT > 40 && PT < 200 && HiBin <= 20", "prof");
   TPbPb->Draw("HasReco*TreeWPbPbMC:PT>>HPbPbMCPTCorrected_cent1",                   "PT > 40 && PT < 200 && HiBin <= 20", "prof");

   TPbPb->Draw("HasReco:Y>>HPbPbMCYRaw_cent2",                             "PT > 40 && PT < 200 && HiBin > 20 && HiBin <= 60", "prof");
   TPbPb->Draw("HasReco*TreeWPbPbMC:Y>>HPbPbMCYCorrected_cent2",                     "PT > 40 && PT < 200 && HiBin > 20 && HiBin <= 60", "prof");
   TPbPb->Draw("HasReco:PT>>HPbPbMCPTRaw_cent2",                           "PT > 40 && PT < 200 && HiBin > 20 && HiBin <= 60", "prof");
   TPbPb->Draw("HasReco*TreeWPbPbMC:PT>>HPbPbMCPTCorrected_cent2",                   "PT > 40 && PT < 200 && HiBin > 20 && HiBin <= 60", "prof");

   TPbPb->Draw("HasReco:Y>>HPbPbMCYRaw_pher1",                             "PT > 40 && PT < 200 && HiBin > 60 && HiBin <= 100", "prof");
   TPbPb->Draw("HasReco*TreeWPbPbMC:Y>>HPbPbMCYCorrected_pher1",                     "PT > 40 && PT < 200 && HiBin > 60 && HiBin <= 100", "prof");
   TPbPb->Draw("HasReco:PT>>HPbPbMCPTRaw_pher1",                           "PT > 40 && PT < 200 && HiBin > 60 && HiBin <= 100", "prof");
   TPbPb->Draw("HasReco*TreeWPbPbMC:PT>>HPbPbMCPTCorrected_pher1",                   "PT > 40 && PT < 200 && HiBin > 60 && HiBin <= 100", "prof");

   TPbPb->Draw("HasReco:Y>>HPbPbMCYRaw_pher2",                             "PT > 40 && PT < 200 && HiBin > 100 && HiBin <= 180", "prof");
   TPbPb->Draw("HasReco*TreeWPbPbMC:Y>>HPbPbMCYCorrected_pher2",                     "PT > 40 && PT < 200 && HiBin > 100 && HiBin <= 180", "prof");
   TPbPb->Draw("HasReco:PT>>HPbPbMCPTRaw_pher2",                           "PT > 40 && PT < 200 && HiBin > 100 && HiBin <= 180", "prof");
   TPbPb->Draw("HasReco*TreeWPbPbMC:PT>>HPbPbMCPTCorrected_pher2",                   "PT > 40 && PT < 200 && HiBin > 100 && HiBin <= 180", "prof");

   TPbPb->Draw("HasReco:Y>>HPbPbMCYRaw",                             "PT > 40 && PT < 200", "prof");
   TPbPb->Draw("HasReco*TreeWPbPbMC:Y>>HPbPbMCYCorrected",                     "PT > 40 && PT < 200", "prof");
   TPbPb->Draw("HasReco:PT>>HPbPbMCPTRaw",                           "PT > 40 && PT < 200", "prof");
   TPbPb->Draw("HasReco*TreeWPbPbMC:PT>>HPbPbMCPTCorrected",                   "PT > 40 && PT < 200", "prof");

   TPbPb->Draw("HasReco:HiBin/2>>HPbPbMCHiBinRaw",                   "PT > 40 && PT < 200", "prof");
   TPbPb->Draw("HasReco*TreeWPbPbMC:HiBin/2>>HPbPbMCHiBinCorrected",           "PT > 40 && PT < 200", "prof");
   
   TPP  ->Draw("HasReco:Y>>HPPMCYRaw",                               "PT > 40 && PT < 200", "prof");
   TPP  ->Draw("HasReco*TreeWPPMC:Y>>HPPMCYCorrected",                   "PT > 40 && PT < 200", "prof");
   TPP  ->Draw("HasReco:PT>>HPPMCPTRaw",                             "PT > 40 && PT < 200", "prof");
   TPP  ->Draw("HasReco*TreeWPPMC:PT>>HPPMCPTCorrected",                 "PT > 40 && PT < 200", "prof");

   MakePlot(-2.1, 2.1, 0., 1.2, ";y^{Z} (0 #leq Centrality #leq 10%);Efficiency", HPbPbMCYRaw_cent1, HPbPbMCYCorrected_cent1, false, "PbPbMCY_cent1.png");
   MakePlot(40, 200, 0., 1.2, ";p_{T}^{Z} (0 #leq Centrality #leq 10%);Efficiency", HPbPbMCPTRaw_cent1, HPbPbMCPTCorrected_cent1, true, "PbPbMCPT_cent1.png");
   MakePlot(-2.1, 2.1, 0., 1.2, ";y^{Z} (30% < Centrality #leq 50%);Efficiency", HPbPbMCYRaw_pher1, HPbPbMCYCorrected_pher1, false, "PbPbMCY_pher1.png");
   MakePlot(40, 200, 0., 1.2, ";p_{T}^{Z} (30% < Centrality #leq 50%);Efficiency", HPbPbMCPTRaw_pher1, HPbPbMCPTCorrected_pher1, true, "PbPbMCPT_pher1.png");

   MakePlot(-2.1, 2.1, 0., 1.2, ";y^{Z} (10% < Centrality #leq 30%);Efficiency", HPbPbMCYRaw_cent2, HPbPbMCYCorrected_cent2, false, "PbPbMCY_cent2.png");
   MakePlot(40, 200, 0., 1.2, ";p_{T}^{Z} (10% < Centrality #leq 30%);Efficiency", HPbPbMCPTRaw_cent2, HPbPbMCPTCorrected_cent2, true, "PbPbMCPT_cent2.png");
   MakePlot(-2.1, 2.1, 0., 1.2, ";y^{Z} (50% < Centrality #leq 90%);Efficiency", HPbPbMCYRaw_pher2, HPbPbMCYCorrected_pher2, false, "PbPbMCY_pher2.png");
   MakePlot(40, 200, 0., 1.2, ";p_{T}^{Z} (50% < Centrality #leq 90%);Efficiency", HPbPbMCPTRaw_pher2, HPbPbMCPTCorrected_pher2, true, "PbPbMCPT_pher2.png");


   MakePlot(-2.1, 2.1, 0., 1.2, ";y^{Z};Efficiency", HPbPbMCYRaw, HPbPbMCYCorrected, false, "PbPbMCY.png");
   MakePlot(40, 200, 0., 1.2, ";p_{T}^{Z};Efficiency", HPbPbMCPTRaw, HPbPbMCPTCorrected, true, "PbPbMCPT.png");
   MakePlot(0, 100, 0., 1.2, ";Centrality (%);Efficiency", HPbPbMCHiBinRaw, HPbPbMCHiBinCorrected, false, "PbPbMCHiBin.png");
   //MakePlot(-2.1, 2.1, 0.4, 1.2, ";y^{Z};Efficiency", HPbPbDataYRaw, HPbPbDataYCorrected, false, "PbPbDataY_pt40.png");
   //MakePlot(1, 200, 0.4, 1.2, ";p_{T}^{Z};Efficiency", HPbPbDataPTRaw, HPbPbDataPTCorrected, true, "PbPbDataPT_pt40.png");
   //MakePlot(0, 100, 0.4, 1.2, ";Centrality (%);Efficiency", HPbPbDataHiBinRaw, HPbPbDataHiBinCorrected, false, "PbPbDataHiBin_pt40.png");
   MakePlot(-2.1, 2.1, 0., 1.2, ";y^{Z};Efficiency", HPPMCYRaw, HPPMCYCorrected, false, "PPMCY.png");
   MakePlot(40, 200, 0., 1.2, ";p_{T}^{Z};Efficiency", HPPMCPTRaw, HPPMCPTCorrected, true, "PPMCPT.png");
   //MakePlot(-2.1, 2.1, 0.4, 1.2, ";y^{Z};Efficiency", HPPDataYRaw, HPPDataYCorrected, false, "PPDataY.png");
   //MakePlot(40, 200, 0.4, 1.2, ";p_{T}^{Z};Efficiency", HPPDataPTRaw, HPPDataPTCorrected, true, "PPDataPT.png");

   FPP.Close();
   FPbPb.Close();

   return 0;
}

void MakePlot(double XMin, double XMax, double YMin, double YMax, string Title,
   TProfile &P1, TProfile &P2, bool LogX, string Output)
{
   static vector<int> Colors = GetPrimaryColors();

   TH2D HWorld("HWorld", Title.c_str(), 100, XMin, XMax, 100, YMin, YMax);
   HWorld.SetStats(0);

   P1.SetMarkerStyle(20);
   P2.SetMarkerStyle(20);
   P1.SetMarkerColor(Colors[0]);
   P2.SetMarkerColor(Colors[1]);
   P1.SetMarkerSize(2);
   P2.SetMarkerSize(2);
   P1.SetLineWidth(2);
   P2.SetLineWidth(2);
   P1.SetLineColor(Colors[0]);
   P2.SetLineColor(Colors[1]);

   TLegend Legend(0.5, 0.2, 0.8, 0.4);
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.06);
   Legend.SetFillStyle(0);
   Legend.SetBorderSize(0);
   Legend.AddEntry(&P1, "Efficiency", "lp");
   Legend.AddEntry(&P2, "Corrected", "lp");

   TCanvas Canvas;

   TGraph G;
   G.SetPoint(0, XMin, 1);
   G.SetPoint(1, XMax, 1);

   HWorld.Draw("axis");
   G.Draw("l");
   P1.Draw("same");
   P2.Draw("same");
   Legend.Draw();

   if(LogX == true)
      Canvas.SetLogx();

   gSystem->Exec(Form("mkdir -p %s",OutputBase.c_str()));

   Canvas.SaveAs((OutputBase + Output).c_str());
}


