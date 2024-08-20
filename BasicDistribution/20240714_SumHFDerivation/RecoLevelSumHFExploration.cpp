#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TLine.h>
#include <cmath>
#include <iostream>
#include <TF1.h> 

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
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);

  gStyle->SetLineScalePS(1.5);

  gROOT->ForceStyle();
}

void PVCheck(TTree* Tree, string OutputBase, string OutputName);
void ZPTCheck(TTree* Tree, string OutputBase, string OutputName);



int main(int argc, char *argv[]){
   
   style();
   string OutputBase = "/eos/user/p/pchou/figs/ZHadron2024/SumHFCheck/v3c_PPMCReco/20240817/";
   string filebase = "/eos/cms/store/group/phys_heavyions/pchou/SkimZHadron2024/"; 

   TChain *TreePPGen  = new TChain("Tree");
   TChain *TreePPReco = new TChain("Tree");
   //TChain *TreePPData = new TChain("Tree");

   TreePPGen->Add((filebase + "OutputPPMCGen_v3c_ee/*.root" ).c_str());
   TreePPReco->Add((filebase + "OutputPPMC_v3c_ee/*.root" ).c_str());
   //TreePPData->Add((filebase + "OutputPPData_v1d_ee/*.root" ).c_str());

   PVCheck(TreePPGen, OutputBase, "SumHFCheck_v3c_PPMCGen_NPU.png");
   PVCheck(TreePPReco, OutputBase, "SumHFCheck_v3c_PPMCReco_NPU.png");

   ZPTCheck(TreePPGen, OutputBase, "SumHFCheck_v3c_PPMCGen_ZPT.png");
   ZPTCheck(TreePPReco, OutputBase, "SumHFCheck_v3c_PPMCReco_ZPT.png");


   return 0;

}

void PVCheck(TTree* Tree, string OutputBase, string OutputName)
{
   TCanvas c1("c1", "", 800, 600);

   TProfile PSignalHFVsNPV = TProfile("PSignalHFVsNPV", ";NPV or NPU;<SignalHF>", 14, 0, 14);
   TProfile PSignalHFVsNPU = TProfile("PSignalHFVsNPU", ";NPV or NPU;<SignalHF>", 14, 0, 14);

   Tree->Draw("SignalHF:NVertex>>PSignalHFVsNPV", "bestZidx>=0 && zPt[bestZidx] > 40", "prof text0");
   Tree->Draw("SignalHF:NPU>>PSignalHFVsNPU", "bestZidx>=0 && zPt[bestZidx] > 40", "prof text0 same");

   PSignalHFVsNPV.SetStats(0);
   PSignalHFVsNPV.SetMarkerStyle(20);
   PSignalHFVsNPV.SetMarkerColor(kRed);
   PSignalHFVsNPV.SetLineColor(kRed);

   PSignalHFVsNPU.SetMarkerStyle(20);
   PSignalHFVsNPU.SetMarkerColor(kBlue);
   PSignalHFVsNPU.SetLineColor(kBlue);

   PSignalHFVsNPV.SetMinimum(500);
   PSignalHFVsNPU.SetMinimum(500);
   
   TLegend Legend(0.18, 0.78, 0.48, 0.90);

   Legend.SetFillStyle(0);
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.035);
   Legend.SetBorderSize(0);
   Legend.AddEntry(&PSignalHFVsNPV, "NPV", "lp");
   Legend.AddEntry(&PSignalHFVsNPU, "NPU", "lp");
   Legend.Draw();

   //c1.SetLogy(0);

   gSystem->mkdir(OutputBase.c_str(), kTRUE);
   c1.SaveAs((OutputBase + OutputName).c_str());
}


void ZPTCheck(TTree* Tree, string OutputBase, string OutputName){
	
	TCanvas c1("c1", "", 800, 600);

	TProfile PSignalHFVsZPT_NPV1 = TProfile("PSignalHFVsZPT_NPV1", ";Z p_{T} (GeV);SignalHF", 100, 0, 100);
	TProfile PSignalHFVsZPT_NPU0 = TProfile("PSignalHFVsZPT_NPU0", "", 100, 0, 100);

	PSignalHFVsZPT_NPV1.SetStats(0);
	PSignalHFVsZPT_NPU0.SetStats(0);
	
	PSignalHFVsZPT_NPV1.SetMarkerStyle(20);
	PSignalHFVsZPT_NPV1.SetMarkerColor(kRed);
	PSignalHFVsZPT_NPV1.SetLineColor(kRed);
	PSignalHFVsZPT_NPU0.SetMarkerStyle(20);
	PSignalHFVsZPT_NPU0.SetMarkerColor(kBlue);
	PSignalHFVsZPT_NPU0.SetLineColor(kBlue);

	Tree->Draw("SignalHF:zPt>>PSignalHFVsZPT_NPV1", "NVertex == 1 && zPt < 100", "prof");
	Tree->Draw("SignalHF:zPt>>PSignalHFVsZPT_NPU0", "NPU == 0 && zPt < 100", "prof same");

	TF1 FSignalHFVsZPT_NPV1("FSignalHFVsZPT_NPV1", "[0]-[1]*exp(-[2]*x)", 0, 100);
	TF1 FSignalHFVsZPT_NPU0("FSignalHFVsZPT_NPU0", "[0]-[1]*exp(-[2]*x)", 0, 100);

	FSignalHFVsZPT_NPV1.SetLineColor(kRed + 1);
	FSignalHFVsZPT_NPU0.SetLineColor(kBlue + 1);

	PSignalHFVsZPT_NPV1.Fit(&FSignalHFVsZPT_NPV1,"E");
	PSignalHFVsZPT_NPU0.Fit(&FSignalHFVsZPT_NPU0,"E");

	PSignalHFVsZPT_NPV1.SetMinimum(600);
	PSignalHFVsZPT_NPU0.SetMinimum(600);

	TLegend Legend(0.18, 0.78, 0.48, 0.90);
	Legend.SetFillStyle(0);
	Legend.SetTextFont(42);
	Legend.SetTextSize(0.035);
	Legend.SetBorderSize(0);
	Legend.AddEntry(&PSignalHFVsZPT_NPV1, "NPV = 1", "lp");
	Legend.AddEntry(&PSignalHFVsZPT_NPU0, "NPU = 0", "lp");
	Legend.Draw();

	gSystem->mkdir(OutputBase.c_str(), kTRUE);
    c1.SaveAs((OutputBase + OutputName).c_str());
}
