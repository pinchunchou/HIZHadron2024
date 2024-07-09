#include <iostream>
using namespace std;

#include "TFile.h"
#include "TTree.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"

#include "SetStyle.h"

int main();
void MakePlot(double XMin, double XMax, double YMin, double YMax, string Title,
   TProfile &P1, TProfile &P2, bool LogX, string Output);

int main()
{
   SetThumbStyle();

   TFile FPbPb("../20240706_ZeeEfficiency/ZEfficiencyAll.root");
   TFile FPP("../20240706_ZeeEfficiency/ZEfficiencyAllPP.root");

   TTree *TPbPb = (TTree *)FPbPb.Get("Tree");
   TTree *TPP = (TTree *)FPP.Get("Tree");

   double PTBins[31] = {0};
   for(int i = 0; i <= 30; i++)
      PTBins[i] = exp(log(1) + (log(200) - log(1)) / 30 * i);

   TProfile HPbPbMCYRaw            ("HPbPbMCYRaw",             ";;", 30, -2.4, 2.4);
   TProfile HPbPbMCYCorrected      ("HPbPbMCYCorrected",       ";;", 30, -2.4, 2.4);
   TProfile HPbPbMCPTRaw           ("HPbPbMCPTRaw",            ";;", 30, PTBins);
   TProfile HPbPbMCPTCorrected     ("HPbPbMCPTCorrected",      ";;", 30, PTBins);
   TProfile HPbPbMCHiBinRaw        ("HPbPbMCHiBinRaw",         ";;", 30, 0, 100);
   TProfile HPbPbMCHiBinCorrected  ("HPbPbMCHiBinCorrected",   ";;", 30, 0, 100);
   TProfile HPbPbDataYRaw          ("HPbPbDataYRaw",           ";;", 30, -2.4, 2.4);
   TProfile HPbPbDataYCorrected    ("HPbPbDataYCorrected",     ";;", 30, -2.4, 2.4);
   TProfile HPbPbDataPTRaw         ("HPbPbDataPTRaw",          ";;", 30, PTBins);
   TProfile HPbPbDataPTCorrected   ("HPbPbDataPTCorrected",    ";;", 30, PTBins);
   TProfile HPbPbDataHiBinRaw      ("HPbPbDataHiBinRaw",       ";;", 30, 0, 100);
   TProfile HPbPbDataHiBinCorrected("HPbPbDataHiBinCorrected", ";;", 30, 0, 100);
   TProfile HPPMCYRaw              ("HPPMCYRaw",               ";;", 30, -2.4, 2.4);
   TProfile HPPMCYCorrected        ("HPPMCYCorrected",         ";;", 30, -2.4, 2.4);
   TProfile HPPMCPTRaw             ("HPPMCPTRaw",              ";;", 30, PTBins);
   TProfile HPPMCPTCorrected       ("HPPMCPTCorrected",        ";;", 30, PTBins);
   TProfile HPPDataYRaw            ("HPPDataYRaw",             ";;", 30, -2.4, 2.4);
   TProfile HPPDataYCorrected      ("HPPDataYCorrected",       ";;", 30, -2.4, 2.4);
   TProfile HPPDataPTRaw           ("HPPDataPTRaw",            ";;", 30, PTBins);
   TProfile HPPDataPTCorrected     ("HPPDataPTCorrected",      ";;", 30, PTBins);

   TPbPb->Draw("HasReco:Y>>HPbPbMCYRaw",                             "PT > 1", "prof");
   TPbPb->Draw("HasReco*WPbPbMC:Y>>HPbPbMCYCorrected",                     "PT > 1", "prof");
   TPbPb->Draw("HasReco:PT>>HPbPbMCPTRaw",                           "PT > 1", "prof");
   TPbPb->Draw("HasReco*WPbPbMC:PT>>HPbPbMCPTCorrected",                   "PT > 1", "prof");
   TPbPb->Draw("HasReco:HiBin/2>>HPbPbMCHiBinRaw",                   "PT > 1", "prof");
   TPbPb->Draw("HasReco*WPbPbMC:HiBin/2>>HPbPbMCHiBinCorrected",           "PT > 1", "prof");
   TPP  ->Draw("HasReco:Y>>HPPMCYRaw",                               "PT > 1", "prof");
   TPP  ->Draw("HasReco*WPPMC:Y>>HPPMCYCorrected",                   "PT > 1", "prof");
   TPP  ->Draw("HasReco:PT>>HPPMCPTRaw",                             "PT > 1", "prof");
   TPP  ->Draw("HasReco*WPPMC:PT>>HPPMCPTCorrected",                 "PT > 1", "prof");
   
   MakePlot(-2.4, 2.4, 0.6, 1.1, ";y^{Z};Efficiency", HPbPbMCYRaw, HPbPbMCYCorrected, false, "PbPbMCY.pdf");
   MakePlot(1, 200, 0.6, 1.1, ";p_{T}^{Z};Efficiency", HPbPbMCPTRaw, HPbPbMCPTCorrected, true, "PbPbMCPT.pdf");
   MakePlot(0, 100, 0.6, 1.1, ";Centrality (%);Efficiency", HPbPbMCHiBinRaw, HPbPbMCHiBinCorrected, false, "PbPbMCHiBin.pdf");
   MakePlot(-2.4, 2.4, 0.6, 1.1, ";y^{Z};Efficiency", HPbPbDataYRaw, HPbPbDataYCorrected, false, "PbPbDataY.pdf");
   MakePlot(1, 200, 0.6, 1.1, ";p_{T}^{Z};Efficiency", HPbPbDataPTRaw, HPbPbDataPTCorrected, true, "PbPbDataPT.pdf");
   MakePlot(0, 100, 0.6, 1.1, ";Centrality (%);Efficiency", HPbPbDataHiBinRaw, HPbPbDataHiBinCorrected, false, "PbPbDataHiBin.pdf");
   MakePlot(-2.4, 2.4, 0.6, 1.1, ";y^{Z};Efficiency", HPPMCYRaw, HPPMCYCorrected, false, "PPMCY.pdf");
   MakePlot(1, 200, 0.6, 1.1, ";p_{T}^{Z};Efficiency", HPPMCPTRaw, HPPMCPTCorrected, true, "PPMCPT.pdf");
   MakePlot(-2.4, 2.4, 0.6, 1.1, ";y^{Z};Efficiency", HPPDataYRaw, HPPDataYCorrected, false, "PPDataY.pdf");
   MakePlot(1, 200, 0.6, 1.1, ";p_{T}^{Z};Efficiency", HPPDataPTRaw, HPPDataPTCorrected, true, "PPDataPT.pdf");

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

   Canvas.SaveAs(Output.c_str());
}


