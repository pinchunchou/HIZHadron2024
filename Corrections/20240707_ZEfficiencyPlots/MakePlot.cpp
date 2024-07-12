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

string OutputBase = "/eos/home-p/pchou/figs/ZHadron2024/ZeeEfficiency/20240712/";

int main()
{
   SetThumbStyle();

   TFile FPbPb("../20240706_ZeeEfficiency/ZEfficiencyPbPb.root");
   TFile FPP("../20240706_ZeeEfficiency/ZEfficiencyPP.root");

   //TTree *TPbPb = (TTree *)FPbPb.Get("Tree");
   //TTree *TPP = (TTree *)FPP.Get("Tree");

   TChain* TPbPb = new TChain("Tree"); TPbPb->Add("../20240706_ZeeEfficiency/PbPb/*.root");
   TChain* TPP = new TChain("Tree"); TPP->Add("../20240706_ZeeEfficiency/pp/*.root");

   double PTBins[21] = {0};
   //for(int i = 0; i <= 30; i++)
   //   PTBins[i] = exp(log(1) + (log(200) - log(1)) / 30 * i);

   for(int i = 0; i <= 20; i++)
      PTBins[i] = exp(log(40) + (log(200) - log(40)) / 20 * i);

   TProfile HPbPbMCYRaw            ("HPbPbMCYRaw",             ";;", 20, -2.1, 2.1);
   TProfile HPbPbMCYCorrected      ("HPbPbMCYCorrected",       ";;", 20, -2.1, 2.1);
   TProfile HPbPbMCPTRaw           ("HPbPbMCPTRaw",            ";;", 20, PTBins);
   TProfile HPbPbMCPTCorrected     ("HPbPbMCPTCorrected",      ";;", 20, PTBins);
   TProfile HPbPbMCHiBinRaw        ("HPbPbMCHiBinRaw",         ";;", 20, 0, 100);
   TProfile HPbPbMCHiBinCorrected  ("HPbPbMCHiBinCorrected",   ";;", 20, 0, 100);
   TProfile HPbPbDataYRaw          ("HPbPbDataYRaw",           ";;", 20, -2.1, 2.1);
   TProfile HPbPbDataYCorrected    ("HPbPbDataYCorrected",     ";;", 20, -2.1, 2.1);
   TProfile HPbPbDataPTRaw         ("HPbPbDataPTRaw",          ";;", 20, PTBins);
   TProfile HPbPbDataPTCorrected   ("HPbPbDataPTCorrected",    ";;", 20, PTBins);

   TProfile HPbPbMCYRaw_cent            ("HPbPbMCYRaw_cent",             ";;", 20, -2.1, 2.1);
   TProfile HPbPbMCYCorrected_cent      ("HPbPbMCYCorrected_cent",       ";;", 20, -2.1, 2.1);
   TProfile HPbPbMCPTRaw_cent           ("HPbPbMCPTRaw_cent",            ";;", 20, PTBins);
   TProfile HPbPbMCPTCorrected_cent     ("HPbPbMCPTCorrected_cent",      ";;", 20, PTBins);

   TProfile HPbPbMCYRaw_pher            ("HPbPbMCYRaw_pher",             ";;", 20, -2.1, 2.1);
   TProfile HPbPbMCYCorrected_pher      ("HPbPbMCYCorrected_pher",       ";;", 20, -2.1, 2.1);
   TProfile HPbPbMCPTRaw_pher           ("HPbPbMCPTRaw_pher",            ";;", 20, PTBins);
   TProfile HPbPbMCPTCorrected_pher     ("HPbPbMCPTCorrected_pher",      ";;", 20, PTBins);


   TProfile HPbPbDataHiBinRaw      ("HPbPbDataHiBinRaw",       ";;", 20, 0, 100);
   TProfile HPbPbDataHiBinCorrected("HPbPbDataHiBinCorrected", ";;", 20, 0, 100);
   TProfile HPPMCYRaw              ("HPPMCYRaw",               ";;", 20, -2.1, 2.1);
   TProfile HPPMCYCorrected        ("HPPMCYCorrected",         ";;", 20, -2.1, 2.1);
   TProfile HPPMCPTRaw             ("HPPMCPTRaw",              ";;", 20, PTBins);
   TProfile HPPMCPTCorrected       ("HPPMCPTCorrected",        ";;", 20, PTBins);
   TProfile HPPDataYRaw            ("HPPDataYRaw",             ";;", 20, -2.1, 2.1);
   TProfile HPPDataYCorrected      ("HPPDataYCorrected",       ";;", 20, -2.1, 2.1);
   TProfile HPPDataPTRaw           ("HPPDataPTRaw",            ";;", 20, PTBins);
   TProfile HPPDataPTCorrected     ("HPPDataPTCorrected",      ";;", 20, PTBins);

   TPbPb->Draw("HasReco:Y>>HPbPbMCYRaw_cent",                             "PT > 40 && HiBin <= 60", "prof");
   TPbPb->Draw("HasReco*TreeWPbPbMC:Y>>HPbPbMCYCorrected_cent",                     "PT > 40 && HiBin <= 60", "prof");
   TPbPb->Draw("HasReco:PT>>HPbPbMCPTRaw_cent",                           "PT > 40 && HiBin <= 60", "prof");
   TPbPb->Draw("HasReco*TreeWPbPbMC:PT>>HPbPbMCPTCorrected_cent",                   "PT > 1 && HiBin <= 60", "prof");

   TPbPb->Draw("HasReco:Y>>HPbPbMCYRaw_pher",                             "PT > 40 && HiBin > 60", "prof");
   TPbPb->Draw("HasReco*TreeWPbPbMC:Y>>HPbPbMCYCorrected_pher",                     "PT > 1 && HiBin > 60", "prof");
   TPbPb->Draw("HasReco:PT>>HPbPbMCPTRaw_pher",                           "PT > 40 && HiBin > 60", "prof");
   TPbPb->Draw("HasReco*TreeWPbPbMC:PT>>HPbPbMCPTCorrected_pher",                   "PT > 40 && HiBin > 60", "prof");

   TPbPb->Draw("HasReco:Y>>HPbPbMCYRaw",                             "PT > 40", "prof");
   TPbPb->Draw("HasReco*TreeWPbPbMC:Y>>HPbPbMCYCorrected",                     "PT > 40", "prof");
   TPbPb->Draw("HasReco:PT>>HPbPbMCPTRaw",                           "PT > 40", "prof");
   TPbPb->Draw("HasReco*TreeWPbPbMC:PT>>HPbPbMCPTCorrected",                   "PT > 40", "prof");

   TPbPb->Draw("HasReco:HiBin/2>>HPbPbMCHiBinRaw",                   "PT > 40", "prof");
   TPbPb->Draw("HasReco*TreeWPbPbMC:HiBin/2>>HPbPbMCHiBinCorrected",           "PT > 40", "prof");
   
   TPP  ->Draw("HasReco:Y>>HPPMCYRaw",                               "PT > 40", "prof");
   TPP  ->Draw("HasReco*TreeWPPMC:Y>>HPPMCYCorrected",                   "PT > 40", "prof");
   TPP  ->Draw("HasReco:PT>>HPPMCPTRaw",                             "PT > 40", "prof");
   TPP  ->Draw("HasReco*TreeWPPMC:PT>>HPPMCPTCorrected",                 "PT > 40", "prof");

   MakePlot(-2.1, 2.1, 0.4, 1.2, ";y^{Z} (0 #leq Centrality #leq 30%);Efficiency", HPbPbMCYRaw_cent, HPbPbMCYCorrected_cent, false, "PbPbMCY_cent.png");
   MakePlot(40, 200, 0.4, 1.2, ";p_{T}^{Z} (0 #leq Centrality #leq 30%);Efficiency", HPbPbMCPTRaw_cent, HPbPbMCPTCorrected_cent, true, "PbPbMCPT_cent.png");
   MakePlot(-2.1, 2.1, 0.4, 1.2, ";y^{Z} (30% < Centrality #leq 100%);Efficiency", HPbPbMCYRaw_pher, HPbPbMCYCorrected_pher, false, "PbPbMCY_pher.png");
   MakePlot(40, 200, 0.4, 1.2, ";p_{T}^{Z} (30% < Centrality #leq 100%);Efficiency", HPbPbMCPTRaw_pher, HPbPbMCPTCorrected_pher, true, "PbPbMCPT_pher.png");

   MakePlot(-2.4, 2.4, 0.4, 1.2, ";y^{Z};Efficiency", HPbPbMCYRaw, HPbPbMCYCorrected, false, "PbPbMCY.png");
   MakePlot(40, 200, 0.4, 1.2, ";p_{T}^{Z};Efficiency", HPbPbMCPTRaw, HPbPbMCPTCorrected, true, "PbPbMCPT.png");
   MakePlot(0, 100, 0.4, 1.2, ";Centrality (%);Efficiency", HPbPbMCHiBinRaw, HPbPbMCHiBinCorrected, false, "PbPbMCHiBin.png");
   //MakePlot(-2.4, 2.4, 0.4, 1.2, ";y^{Z};Efficiency", HPbPbDataYRaw, HPbPbDataYCorrected, false, "PbPbDataY_pt40.png");
   //MakePlot(1, 200, 0.4, 1.2, ";p_{T}^{Z};Efficiency", HPbPbDataPTRaw, HPbPbDataPTCorrected, true, "PbPbDataPT_pt40.png");
   //MakePlot(0, 100, 0.4, 1.2, ";Centrality (%);Efficiency", HPbPbDataHiBinRaw, HPbPbDataHiBinCorrected, false, "PbPbDataHiBin_pt40.png");
   MakePlot(-2.1, 2.1, 0.4, 1.2, ";y^{Z};Efficiency", HPPMCYRaw, HPPMCYCorrected, false, "PPMCY.png");
   MakePlot(40, 200, 0.4, 1.2, ";p_{T}^{Z};Efficiency", HPPMCPTRaw, HPPMCPTCorrected, true, "PPMCPT.png");
   MakePlot(-2.1, 2.1, 0.4, 1.2, ";y^{Z};Efficiency", HPPDataYRaw, HPPDataYCorrected, false, "PPDataY.png");
   MakePlot(40, 200, 0.4, 1.2, ";p_{T}^{Z};Efficiency", HPPDataPTRaw, HPPDataPTCorrected, true, "PPDataPT.png");

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


