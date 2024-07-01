#include <map>
#include <cmath>
#include <iostream>
using namespace std;

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TMath.h"
#include "TROOT.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPad.h"

#include "SetStyle.h"

#define TptL_min 0.5
#define typeofdata "20230518"
#define typeofdatatext "single muon"

struct Files;
struct Plots;
struct Setting;
int main(int argc, char *argv[]);
void ZtrackDraw_single(Files &File, Setting S, string OutputBase);
void DivideByBinWidth(TH1D *H = nullptr);
void DivideByBinWidth(TH2D *H = nullptr);
void Draw1DPlot(TH1D *H1, TH1D *H2, string XTitle, string YTitle, bool WithMin,
   Setting &S, string Base, string Tag);
void Draw2DPlot2Panel(TH2D *H1, TH2D *H2, string XTitle, string YTitle, string ZTitle,
   Setting &S, string Base, string Tag, string Option, vector<string> Identifier = {"MC", "Data"});
void Draw2DPlot3Panel(TH2D *H1, TH2D *H2, TH2D *H3, string XTitle, string YTitle, string ZTitle,
   Setting &S, string Base, string Tag, string Option, vector<string> Identifier = {"PbPB MC", "Data", "pp MC"});

struct Files
{
public:
   TFile *SignalMC;
   TFile *BackgroundMC;
   TFile *SignalData;
   TFile *BackgroundData;
   TFile *ppMC;
   TFile *ppData;
   TFile *SignalMCGen;
   TFile *BackgroundMCGen;
public:
   Files() : SignalMC(nullptr), BackgroundMC(nullptr), SignalData(nullptr), BackgroundData(nullptr),
      ppMC(nullptr), ppData(nullptr), SignalMCGen(nullptr), BackgroundMCGen(nullptr)
   {
   }
   ~Files()
   {
      if(SignalMC != nullptr)          {SignalMC->Close(); delete SignalMC;}
      if(BackgroundMC != nullptr)      {BackgroundMC->Close(); delete BackgroundMC;}
      if(SignalData != nullptr)        {SignalData->Close(); delete SignalData;}
      if(BackgroundData != nullptr)    {BackgroundData->Close(); delete BackgroundData;}
      if(ppMC != nullptr)              {ppMC->Close(); delete ppMC;}
      if(ppData != nullptr)            {ppData->Close(); delete ppData;}
      if(SignalMCGen != nullptr)       {SignalMCGen->Close(); delete SignalMCGen;}
      if(BackgroundMCGen != nullptr)   {BackgroundMCGen->Close(); delete BackgroundMCGen;}
   }
};

struct Plots
{
public:
   vector<string> N1;
   vector<string> N2;
   map<string, TH1D *> H1;
   map<string, TH2D *> H2;
public:
   Plots() {}
   Plots(TFile *File, string FolderName, int Rebin1D = 1, int Rebin2D = 1) :
      N1{"HEta", "HPhi", "HTrackMuonDEta", "HTrackMuonDPhi"},
      N2{"HEtaPhi", "HTrackMuonDEtaDPhi", "HMaxHadronEtaPhi", "HMaxOppositeHadronEtaPhi",
         "HWTAEtaPhi", "HWTAMoreEtaPhi", "HZMaxHadronEtaPhi", "HZMaxOppositeHadronEtaPhi",
         "HZWTAEtaPhi", "HZWTAMoreEtaPhi"}
   {
      for(string N : N1) H1[N] = Prepare1DHistogram(File, FolderName, N, Rebin1D);
      for(string N : N2) H2[N] = Prepare2DHistogram(File, FolderName, N, Rebin2D);
   }
   Plots(Plots &other) : N1(other.N1), N2(other.N2)
   {
      for(string N : N1)
         if(other.H1[N] != nullptr)
            H1[N] = (TH1D *)other.H1[N]->Clone();
      for(string N : N2)
         if(other.H2[N] != nullptr)
            H2[N] = (TH2D *)other.H2[N]->Clone();
   }
   void HistogramStyle(int Color, int Marker)
   {
      for(string N : N1)
      {
         if(H1[N] != nullptr)
         {
            H1[N]->SetLineColor(Color);
            H1[N]->SetMarkerColor(Color);
            H1[N]->SetLineWidth(2);
            H1[N]->SetMarkerStyle(Marker);
         }
      }
      for(string N : N2)
      {
         if(H2[N] != nullptr)
         {
            H2[N]->SetMarkerColor(Color);
            H2[N]->SetMarkerStyle(Marker);
         }
      }
   }
   void Subtract(Plots &other)
   {
      for(string N : N1)
         if(other.H1.find(N) != other.H1.end() && other.H1[N] != nullptr && H1[N] != nullptr)
            H1[N]->Add(other.H1[N], -1);
      for(string N : N2)
         if(other.H2.find(N) != other.H2.end() && other.H2[N] != nullptr && H2[N] != nullptr)
            H2[N]->Add(other.H2[N], -1);
   }
   void Divide(Plots &other)
   {
      for(string N : N1)
         if(other.H1.find(N) != other.H1.end() && other.H1[N] != nullptr && H1[N] != nullptr)
            H1[N]->Divide(other.H1[N]);
      for(string N : N2)
         if(other.H2.find(N) != other.H2.end() && other.H2[N] != nullptr && H2[N] != nullptr)
            H2[N]->Divide(other.H2[N]);
   }
   TH1D *Prepare1DHistogram(TFile *File, string FolderName, string HistogramName, int Rebin)
   {
      if(File == nullptr)
         return nullptr;

      TH1D *H = (TH1D *)File->Get(Form("%s/%s", FolderName.c_str(), HistogramName.c_str()));

      TNamed *N  = (TNamed *)File->Get(Form("%s/EntryCount", FolderName.c_str()));
      float I = atof(N->GetTitle());
      H->Scale(1 / I);

      if(Rebin > 1)
         H->Rebin(Rebin);

      DivideByBinWidth(H);

      return H;
   }
   TH2D *Prepare2DHistogram(TFile *File, string FolderName, string HistogramName, int Rebin)
   {
      if(File == nullptr)
         return nullptr;

      TH2D *H = (TH2D *)File->Get(Form("%s/%s", FolderName.c_str(), HistogramName.c_str()));

      TNamed *N  = (TNamed *)File->Get(Form("%s/EntryCount", FolderName.c_str()));
      float I = atof(N->GetTitle());
      H->Scale(1 / I);

      if(Rebin > 1)
         H->Rebin2D(Rebin, Rebin);

      DivideByBinWidth(H);

      return H;
   }
};

struct Setting
{
public:
   int binnum;
   float ptL;
   float ptH;
   float centL;
   float centH;
   float TptL;
   float TptH;
public:
   Setting() {}
   Setting(int b, float p1, float p2, float c1, float c2, float T1, float T2)
      : binnum(b), ptL(p1), ptH(p2), centL(c1), centH(c2), TptL(T1), TptH(T2)
   {
   }
};

int main(int argc, char *argv[])
{
   SetCorrelationStyle();

   Files File;

   string filebase = "/eos/home-p/pchou/BasicPlots/ZHadron2024/"
   
   File.SignalData      = TFile::Open((filebase + "GraphDataSig_v1_ee.root" ).c_str(), "read");
   File.SignalMC        = TFile::Open((filebase + "GraphMCSig_v1_ee.root"   ).c_str(), "read");
   File.BackgroundData  = TFile::Open((filebase + "GraphDataBkg_v1_ee.root" ).c_str(), "read");
   File.BackgroundMC    = TFile::Open((filebase + "GraphMCBkg_v1_ee.root"   ).c_str(), "read");
   File.ppData          = TFile::Open((filebase + "GraphPPData_v1_ee.root"  ).c_str(), "read");
   File.ppMC            = TFile::Open((filebase + "GraphPPMC_v1_ee.root"    ).c_str(), "read");
   File.SignalMCGen     = TFile::Open((filebase + "GraphMCSigGen_v1_ee.root").c_str(), "read");
   File.BackgroundMCGen = TFile::Open((filebase + "GraphMCBkgGen_v1_ee.root").c_str(), "read");

   string OutputBase = "/eos/user/p/pchou/figs/ZHadron2024/";
   // string OutputBase = ".";

   ZtrackDraw_single(File, Setting(40, 40,  200,  0,  10, 1, 2), OutputBase);
   ZtrackDraw_single(File, Setting(40, 30,  200,  0,  10, 1, 2), OutputBase);
   ZtrackDraw_single(File, Setting(40, 40,  200,  0,  10, 1, 1000), OutputBase);
   ZtrackDraw_single(File, Setting(40, 30,  200,  0,  10, 1, 1000), OutputBase);

   ZtrackDraw_single(File, Setting(40, 40,  200,  0, 100, 1, 2), OutputBase);
   ZtrackDraw_single(File, Setting(40, 30,  200,  0, 100, 1, 2), OutputBase);
   ZtrackDraw_single(File, Setting(40, 40,  200,  0, 100, 1, 1000), OutputBase);
   ZtrackDraw_single(File, Setting(40, 30,  200,  0, 100, 1, 1000), OutputBase);

   ZtrackDraw_single(File, Setting(40, 40,  200, 10,  30, 1, 2), OutputBase);
   ZtrackDraw_single(File, Setting(40, 40,  200, 30,  50, 1, 2), OutputBase);
   ZtrackDraw_single(File, Setting(40, 40,  200, 50,  90, 1, 2), OutputBase);

   ZtrackDraw_single(File, Setting(40, 40,  200, 10,  30, 1, 1000), OutputBase);
   ZtrackDraw_single(File, Setting(40, 40,  200, 30,  50, 1, 1000), OutputBase);
   ZtrackDraw_single(File, Setting(40, 40,  200, 50,  90, 1, 1000), OutputBase);

   ZtrackDraw_single(File, Setting(40, 40,  200,  0,  10, 2, 4), OutputBase);
   ZtrackDraw_single(File, Setting(40, 40,  200, 10,  30, 2, 4), OutputBase);
   ZtrackDraw_single(File, Setting(40, 40,  200, 30,  50, 2, 4), OutputBase);
   ZtrackDraw_single(File, Setting(40, 40,  200, 50,  90, 2, 4), OutputBase);
   ZtrackDraw_single(File, Setting(40, 40,  200,  0, 100, 2, 4), OutputBase);

   ZtrackDraw_single(File, Setting(40, 40,  200,  0,  10, 4, 10), OutputBase);
   ZtrackDraw_single(File, Setting(40, 40,  200, 10,  30, 4, 10), OutputBase);
   ZtrackDraw_single(File, Setting(40, 40,  200, 30,  50, 4, 10), OutputBase);
   ZtrackDraw_single(File, Setting(40, 40,  200, 50,  90, 4, 10), OutputBase);
   ZtrackDraw_single(File, Setting(40, 40,  200,  0, 100, 4, 10), OutputBase);

   return 0;

}

void ZtrackDraw_single(Files &File, Setting S, string OutputBase)
{
   cout << "ptL = " << S.ptL << ", ptH = " << S.ptH
      << ", centL = " << S.centL << ", centH = " << S.centH
      << ", TptL = " << S.TptL << ", TptH = " << S.TptH << endl;
   
   TCanvas *c = new TCanvas("Canvas", "", 800, 800);

   cout << "Getting histograms..." << endl;

   string FolderName = Form("Plot_ZPT_%.0f_%.0f_Cent_%.0f_%.0f_TrackPT_%.2f_%.2f",S.ptL,S.ptH,S.centL,S.centH,S.TptL,S.TptH);
   replace(FolderName.begin(), FolderName.end(), '.', 'p');

   Plots HSignalData(File.SignalData, FolderName, 1, 5);
   Plots HSignalMC(File.SignalMC, FolderName, 1, 4);
   Plots HppData(File.ppData, FolderName, 1, 5);
   Plots HppMC(File.ppMC, FolderName, 1, 4);
   Plots HBackgroundData(File.BackgroundData, FolderName, 1, 5);
   Plots HBackgroundMC(File.BackgroundMC, FolderName, 1, 4);
   Plots HSignalMCGen(File.SignalMCGen, FolderName, 1, 4);
   Plots HBackgroundMCGen(File.BackgroundMCGen, FolderName, 1, 4);

   HSignalData.HistogramStyle(kBlack, 24);
   HSignalMC.HistogramStyle(kRed, 24);
   HppData.HistogramStyle(kBlack, 24);
   HppMC.HistogramStyle(kRed, 24);
   HBackgroundData.HistogramStyle(kBlack, 24);
   HBackgroundMC.HistogramStyle(kRed, 24);

   Plots HSignalDataSB(HSignalData);
   Plots HSignalMCSB(HSignalMC);
   Plots HSignalMCGenSB(HSignalMCGen);
   
   HSignalDataSB.Subtract(HBackgroundData);
   HSignalMCSB.Subtract(HBackgroundMC);
   HSignalMCGenSB.Subtract(HBackgroundMCGen);

   Plots HSignalDataSBR(HSignalData);
   Plots HSignalMCSBR(HSignalMC);
   Plots HSignalMCGenSBR(HSignalMCGen);
   
   HSignalDataSBR.Divide(HBackgroundData);
   HSignalMCSBR.Divide(HBackgroundMC);
   HSignalMCGenSBR.Divide(HBackgroundMCGen);
   
   int countD0 = HSignalData.H1["HEta"]->GetEntries();
   cout<<"Data 0 = "<<countD0<<endl;
   int countM0 = HSignalMC.H1["HEta"]->GetEntries();
   cout<<"MC 0 = "<<countM0<<endl;

   int countDb = HBackgroundData.H1["HEta"]->GetEntries();
   cout<<"Data Bkg = "<<countDb<<endl;
   int countMb = HBackgroundMC.H1["HEta"]->GetEntries();
   cout<<"MC Bkg = "<<countMb<<endl;

   if(S.TptL==0) S.TptL=TptL_min;

   TLegend Legend(0.58, 0.78, 0.98, 0.90);
   Legend.AddEntry(HSignalMC.H1["HEta"], "Monte Carlo: DYLL","lep");
   Legend.AddEntry(HSignalData.H1["HEta"], Form("Data: %s",typeofdatatext),"lep");
   Legend.SetFillStyle(0);
   Legend.SetLineColor(kBlack);
   Legend.SetLineWidth(1);
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.035);

   // == Start drawing == //

   //// // Draw eta 

   string Identifier = Form("%.0f_%.0f_%.0f_%.0f_%.0f_%.0f", S.ptL, S.ptH, S.centL, S.centH, S.TptL, S.TptH);

   Draw1DPlot(HSignalMC.H1["HEta"], HSignalData.H1["HEta"],
      "Signal |#Delta#eta_{Z,track}|", "dN/d#Delta#eta", false, S,
      OutputBase + "/" + typeofdata + "/Deta",
      Form("Ztrack_%s_sig_%s_Deta", typeofdata, Identifier.c_str()));
   
   Draw1DPlot(HBackgroundMC.H1["HEta"], HBackgroundData.H1["HEta"],
      "Background |#Delta#eta_{Z,track}|", "dN/d#Delta#eta", false, S,
      OutputBase + "/" + typeofdata + "/Deta",
      Form("Ztrack_%s_bkg_%s_Deta", typeofdata, Identifier.c_str()));
   
   Draw1DPlot(HSignalMCSB.H1["HEta"], HSignalDataSB.H1["HEta"],
      "Signal - Background |#Delta#eta_{Z,track}|", "dN/d#Delta#eta", false, S,
      OutputBase + "/" + typeofdata + "/Deta",
      Form("Ztrack_%s_sb_%s_Deta", typeofdata, Identifier.c_str()));
   
   Draw1DPlot(HSignalMCSBR.H1["HEta"], HSignalDataSBR.H1["HEta"],
      "Signal - Background |#Delta#eta_{Z,track}|", "dN/d#Delta#eta", false, S,
      OutputBase + "/" + typeofdata + "/Deta",
      Form("Ztrack_%s_sb_%s_Deta", typeofdata, Identifier.c_str()));
   
   Draw1DPlot(HSignalMC.H1["HTrackMuonEta"], HSignalData.H1["HTrackMuonEta"],
      "Signal |#Delta#eta_{#mu#mu}|", "dN/d#Delta#eta", false, S,
      OutputBase + "/" + typeofdata + "/muD",
      Form("Ztrack_%s_sig_%s_muDeta", typeofdata, Identifier.c_str()));
   
   Draw1DPlot(HBackgroundMC.H1["HTrackMuonEta"], HBackgroundData.H1["HTrackMuonEta"],
      "Background |#Delta#eta_{#mu#mu}|", "dN/d#Delta#eta", false, S,
      OutputBase + "/" + typeofdata + "/muD",
      Form("Ztrack_%s_bkg_%s_muDeta", typeofdata, Identifier.c_str()));
   
   Draw1DPlot(HSignalMC.H1["HPhi"], HSignalData.H1["HPhi"],
      "Signal |#Delta#phi_{Z,track}|", "dN/d#Delta#phi", false, S,
      OutputBase + "/" + typeofdata + "/Dphi",
      Form("Ztrack_%s_sig_%s_Dphi", typeofdata, Identifier.c_str()));
   
   Draw1DPlot(HBackgroundMC.H1["HPhi"], HBackgroundData.H1["HPhi"],
      "Background |#Delta#phi_{Z,track}|", "dN/d#Delta#phi", false, S,
      OutputBase + "/" + typeofdata + "/Dphi",
      Form("Ztrack_%s_bkg_%s_Dphi", typeofdata, Identifier.c_str()));
   
   Draw1DPlot(HSignalMCSB.H1["HPhi"], HSignalDataSB.H1["HPhi"],
      "Signal - Background |#Delta#phi_{Z,track}|", "dN/d#Delta#phi", false, S,
      OutputBase + "/" + typeofdata + "/Dphi",
      Form("Ztrack_%s_sb_%s_Dphi", typeofdata, Identifier.c_str()));

   
   // TODO: do something about this 4-curve plot
   //// TH1D *hMC_phi_com = (TH1D*) hMC_phi->Clone("hMC_phi_com");
   //// TH1D *hMC_bkg_phi_com = (TH1D*) hMC_bkg_phi->Clone("hMC_bkg_phi_com");
   //// TH1D *hMC_sb_phi_com = (TH1D*) hMC_sb_phi->Clone("hMC_sb_phi_com");
   //// TH1D *hpp_phi_com = (TH1D*) hpp_phi->Clone("hpp_phi_com");
   

   Draw1DPlot(HSignalMCSBR.H1["HPhi"], HSignalDataSBR.H1["HPhi"],
         "Signal / Background |#Delta#phi_{Z,track}|", "dN/d#Delta#phi", true, S,
         OutputBase + "/" + typeofdata + "/Dphi",
         Form("Ztrack_%s_sbr_%s_Dphi", typeofdata, Identifier.c_str()));

   // TODO Check this?
   //// hMC_MuDphi->SetMaximum(3.0/binnum); 
   //// hData_MuDphi->SetMaximum(3.0/binnum); 
   Draw1DPlot(HSignalMC.H1["HTrackMuonPhi"], HSignalData.H1["HTrackMuonPhi"],
         "Signal |#Delta#phi_{#mu#mu}|", "dN/d#Delta#phi", false, S,
         OutputBase + "/" + typeofdata + "/muD",
         Form("Ztrack_%s_sig_%s_muDphi", typeofdata, Identifier.c_str()));
   
   //// hMC_bkg_MuDphi->SetMaximum(3.0/binnum); 
   //// hData_bkg_MuDphi->SetMaximum(3.0/binnum); 
   Draw1DPlot(HSignalMC.H1["HTrackMuonPhi"], HSignalData.H1["HTrackMuonPhi"],
      "Background |#Delta#phi_{#mu#mu}|", "dN/d#Delta#phi", false, S,
      OutputBase + "/" + typeofdata + "/muD",
      Form("Ztrack_%s_bkg_%s_muDphi", typeofdata, Identifier.c_str()));


   // Now we move to 2D plots

   Draw2DPlot2Panel(HSignalMC.H2["HEtaPhi"], HSignalData.H2["HEtaPhi"],
      "Signal #Delta#eta_{Z,track}", "Signal #Delta#phi_{Z,track}", "dN/d#Delta#etad#Delta#phi", S,
      OutputBase + "/" + typeofdata + "/3D",
      Form("Ztrack_%s_sig_%s_Detaphi_3D", typeofdata, Identifier.c_str()), "lego20");
   
   Draw2DPlot2Panel(HBackgroundMC.H2["HEtaPhi"], HBackgroundData.H2["HEtaPhi"],
      "Background #Delta#eta_{Z,track}", "Background #Delta#phi_{Z,track}", "dN/d#Delta#etad#Delta#phi", S,
      OutputBase + "/" + typeofdata + "/3D",
      Form("Ztrack_%s_bkg_%s_Detaphi_3D", typeofdata, Identifier.c_str()), "lego20");

   Draw2DPlot2Panel(HSignalMCSB.H2["HEtaPhi"], HSignalDataSB.H2["HEtaPhi"],
      "Signal - Background #Delta#eta_{Z,track}", "Signal - Background #Delta#phi_{Z,track}", "dN/d#Delta#etad#Delta#phi", S,
      OutputBase + "/" + typeofdata + "/3D",
      Form("Ztrack_%s_sb_%s_Detaphi_3D", typeofdata, Identifier.c_str()), "lego20");
   
   Draw2DPlot2Panel(HSignalMCSBR.H2["HEtaPhi"], HSignalDataSBR.H2["HEtaPhi"],
      "Signal/Background #Delta#eta_{Z,track}", "Signal/Background #Delta#phi_{Z,track}", "dN/d#Delta#etad#Delta#phi", S,
      OutputBase + "/" + typeofdata + "/3D",
      Form("Ztrack_%s_sbr_%s_Detaphi_3D", typeofdata, Identifier.c_str()), "lego20");
   
   Draw2DPlot2Panel(HSignalMCSBR.H2["HEtaPhi"], HSignalDataSBR.H2["HEtaPhi"],
      "Signal/Background #Delta#eta_{Z,track}", "Signal/Background #Delta#phi_{Z,track}", "dN/d#Delta#etad#Delta#phi", S,
      OutputBase + "/" + typeofdata + "/3D",
      Form("Ztrack_%s_sbr_%s_Detaphi_COLZ", typeofdata, Identifier.c_str()), "colz");
   
   Draw2DPlot2Panel(HSignalMCGen.H2["HEtaPhi"], HSignalMCGen.H2["HEtaPhi"],
      "MC Gen #Delta#eta_{Z,track}", "MC Gen #Delta#phi_{Z,track}", "dN/d#Delta#etad#Delta#phi", S,
      OutputBase + "/" + typeofdata + "/gen",
      Form("Ztrack_%s_sig_%s_Detaphi_gen", typeofdata, Identifier.c_str()), "lego20",
      {"Signal", "Background"});
   
   Draw2DPlot2Panel(HSignalMCGenSB.H2["HEtaPhi"], HSignalMCGenSBR.H2["HEtaPhi"],
      "MC Gen #Delta#eta_{Z,track}", "MC Gen #Delta#phi_{Z,track}", "dN/d#Delta#etad#Delta#phi", S,
      OutputBase + "/" + typeofdata + "/gen",
      Form("Ztrack_%s_sig_%s_Detaphi_gensub", typeofdata, Identifier.c_str()), "lego20",
      {"Signal - Background", "Signal/Background"});
   
   Draw2DPlot2Panel(HSignalMCGen.H2["HEtaPhi"], HSignalMCGen.H2["HEtaPhi"],
      "MC Gen #Delta#eta_{Z,track}", "MC Gen #Delta#phi_{Z,track}", "dN/d#Delta#etad#Delta#phi", S,
      OutputBase + "/" + typeofdata + "/gen",
      Form("Ztrack_%s_sig_%s_Detaphi_gen_COLZ", typeofdata, Identifier.c_str()), "colz",
      {"Signal", "Background"});

   Draw2DPlot2Panel(HSignalMCGenSB.H2["HEtaPhi"], HSignalMCGenSBR.H2["HEtaPhi"],
      "MC Gen #Delta#eta_{Z,track}", "MC Gen #Delta#phi_{Z,track}", "dN/d#Delta#etad#Delta#phi", S,
      OutputBase + "/" + typeofdata + "/gen",
      Form("Ztrack_%s_sig_%s_Detaphi_gensub_COLZ", typeofdata, Identifier.c_str()), "colz",
      {"Signal - Background", "Signal/Background"});
  
   // Now we move on to some 3-panel plots
   
   Draw2DPlot3Panel(HSignalMCSB.H2["HEtaPhi"], HSignalDataSB.H2["HEtaPhi"], HppMC.H2["HEtaPhi"],
      "#Delta#eta_{Z,track}", "#Delta#phi_{Z,track}", "dN/d#Delta#eta#Delta#phi", S,
      OutputBase + "/" + typeofdata + "/pp",
      Form("Ztrack_%s_sb_%s_Detaphi_pp", typeofdata, Identifier.c_str()), "lego20",
      {"Signal - Background MC", "Signal - Background Data", "pp MC (NPU = 0)"});
   
   Draw2DPlot3Panel(HSignalMCSB.H2["HEtaPhi"], HSignalDataSB.H2["HEtaPhi"], HppMC.H2["HEtaPhi"],
      "#Delta#eta_{Z,track}", "#Delta#phi_{Z,track}", "dN/d#Delta#eta#Delta#phi", S,
      OutputBase + "/" + typeofdata + "/pp",
      Form("Ztrack_%s_sb_%s_Detaphi_pp_COLZ", typeofdata, Identifier.c_str()), "colz",
      {"Signal - Background MC", "Signal - Background Data", "pp MC (NPU = 0)"});

   Draw2DPlot3Panel(HSignalMCSB.H2["HEtaPhi"], HSignalMCSBR.H2["HEtaPhi"], HppMC.H2["HEtaPhi"],
      "#Delta#eta_{Z,track}", "#Delta#phi_{Z,track}", "dN/d#Delta#eta#Delta#phi", S,
      OutputBase + "/" + typeofdata + "/pp",
      Form("Ztrack_%s_sbr_%s_Detaphi_pp", typeofdata, Identifier.c_str()), "lego20",
      {"Signal - Background MC", "Signal/Background Data", "pp MC (NPU = 0)"});
   
   Draw2DPlot3Panel(HSignalMCSB.H2["HEtaPhi"], HSignalMCSBR.H2["HEtaPhi"], HppMC.H2["HEtaPhi"],
      "#Delta#eta_{Z,track}", "#Delta#phi_{Z,track}", "dN/d#Delta#eta#Delta#phi", S,
      OutputBase + "/" + typeofdata + "/pp",
      Form("Ztrack_%s_sbr_%s_Detaphi_pp_COLZ", typeofdata, Identifier.c_str()), "colz",
      {"Signal - Background MC", "Signal/Background Data", "pp MC (NPU = 0)"});
  
   // Back to 2 panel?
   
   Draw2DPlot2Panel(HSignalMC.H2["HTrackMuonDEtaDPhi"], HSignalData.H2["HTrackMuonDEtaDPhi"],
      "#Delta#eta_{#mu#mu}", "#Delta#phi_{#mu#mu}", "dN/d#Delta#etad#Delta#phi", S,
      OutputBase + "/" + typeofdata + "/muD",
      Form("Ztrack_%s_sig_%s_muDetaphi", typeofdata, Identifier.c_str()), "lego20",
      {"Signal MC", "Signal Data"});

   Draw2DPlot2Panel(HBackgroundMC.H2["HTrackMuonDEtaDPhi"], HBackgroundData.H2["HTrackMuonDEtaDPhi"],
      "#Delta#eta_{#mu#mu}", "#Delta#phi_{#mu#mu}", "dN/d#Delta#etad#Delta#phi", S,
      OutputBase + "/" + typeofdata + "/muD",
      Form("Ztrack_%s_bkg_%s_muDetaphi", typeofdata, Identifier.c_str()), "lego20",
      {"Background MC", "Background Data"});

   // TODO: do something about these projected 2-panel plots


}


void DivideByBinWidth(TH1D *H)
{
   if(H == nullptr)
      return;

   int NBin = H->GetNbinsX();
   for(int i = 1; i <= NBin; i++)
   {
      double XL = H->GetXaxis()->GetBinLowEdge(i);
      double XH = H->GetXaxis()->GetBinUpEdge(i);
      H->SetBinContent(i, H->GetBinContent(i) / (XH - XL));
      H->SetBinError(i, H->GetBinError(i) / (XH - XL));
   }
}

void DivideByBinWidth(TH2D *H)
{
   if(H == nullptr)
      return;

   int NBinX = H->GetNbinsX();
   int NBinY = H->GetNbinsY();
   for(int iX = 1; iX <= NBinX; iX++)
   {
      for(int iY = 1; iY <= NBinY; iY++)
      {
         double XL = H->GetXaxis()->GetBinLowEdge(iY);
         double XH = H->GetXaxis()->GetBinUpEdge(iY);
         double YL = H->GetYaxis()->GetBinLowEdge(iY);
         double YH = H->GetYaxis()->GetBinUpEdge(iY);
         H->SetBinContent(iX, iY, H->GetBinContent(iX, iY) / (XH - XL) / (YH - YL));
         H->SetBinError(iX, iY, H->GetBinError(iX, iY) / (XH - XL) / (YH - YL));
      }
   }
}

void Draw1DPlot(TH1D *H1, TH1D *H2, string XTitle, string YTitle, bool WithMin,
   Setting &S, string Base, string Tag)
{
   if(H1 == nullptr || H2 == nullptr)
      return;
   
   TCanvas Canvas("Canvas", "", 800, 800);

   double Max1 = H1->GetMaximum();
   double Max2 = H2->GetMaximum();
   double Min1 = H1->GetMinimum();
   double Min2 = H2->GetMinimum();
   
   if(Max1 < Max2)
   {
      H1->Draw();
      H2->Draw("same");
   }
   else
   {
      H2->Draw();
      H1->Draw("same");
   }

   double Max = (Max1 < Max2) ? Max2 : Max1;
   double Min = (Min1 < Min2) ? Min1 : Min2;

   H1->SetXTitle(XTitle.c_str());
   H1->SetYTitle(YTitle.c_str());
   H2->SetXTitle(XTitle.c_str());
   H2->SetYTitle(YTitle.c_str());

   TLegend Legend(0.58, 0.78, 0.98, 0.90);
   Legend.AddEntry(H1, "Monte Carlo: DYLL","lep");
   Legend.AddEntry(H2, Form("Data: %s",typeofdatatext),"lep");
   Legend.SetFillStyle(0);
   Legend.SetLineColor(kBlack);
   Legend.SetLineWidth(1);
   Legend.SetTextFont(42);
   Legend.SetTextSize(0.035);
   
   Legend.Draw();

   TLatex Latex;
   Latex.SetTextFont(42);
   Latex.SetTextSize(0.03);
   Latex.SetNDC();

   Latex.DrawLatex(0.18, 0.88, Form("%.0f %%< Centrality < %.0f %%", S.centL, S.centH));
   Latex.DrawLatex(0.18, 0.82, Form("%.1f < Z p_{T} < %.1f", S.ptL, S.ptH));
   Latex.DrawLatex(0.18, 0.76, Form("%.1f < Track p_{T} < %.1f", S.TptL, S.TptH));
   // Latex.DrawLatex(0.10, 0.97, Form("Signal N_{MC}^{Z} = %.1f, N_{Data}^{Z} = %.1f", tM_tN, tD_tN));

   H1->SetMaximum(1.6 * Max); 
   H2->SetMaximum(1.6 * Max); 

   if(WithMin == false)
   {
      H1->SetMinimum(0);
      H2->SetMinimum(0);
   }
   else
   {
      H1->SetMinimum(Min);
      H2->SetMinimum(Min);
   }
   
   gSystem->Exec(Form("mkdir -p %s", Base.c_str()));
   Canvas.SaveAs(Form("%s/%s.png", Base.c_str(), Tag.c_str())); 
   // Canvas.SaveAs(Form("%s/pdf/%s.pdf", Base.c_str(), Tag.c_str())); 
   // Canvas.SaveAs(Form("%s/C/%s.C", Base.c_str(), Tag.c_str())); 
}

void Draw2DPlot2Panel(TH2D *H1, TH2D *H2, string XTitle, string YTitle, string ZTitle,
   Setting &S, string Base, string Tag, string Option, vector<string> Identifier)
{
   if(H1 == nullptr || H2 == nullptr)
      return;

   TCanvas Canvas("Canvas", "", 1400, 800);

   Canvas.Divide(2);
   TPad *P1 = (TPad *)Canvas.cd(1);
   P1->SetTheta(60.839);
   P1->SetPhi(38.0172);

   if(Identifier.size() >= 1)
      H1->SetTitle(Identifier[0].c_str());
   else
      H1->SetTitle("MC");
   H1->Draw(Option.c_str());
   H1->GetXaxis()->SetTitle(XTitle.c_str());
   H1->GetYaxis()->SetTitle(YTitle.c_str());
   H1->GetXaxis()->SetTitleSize(30);
   H1->GetYaxis()->SetTitleSize(30);
   H1->GetXaxis()->SetTitleOffset(3.0);
   H1->GetYaxis()->SetTitleOffset(2.5);
   H1->GetXaxis()->SetNdivisions(50205, kFALSE);
   H1->GetZaxis()->SetTitle(ZTitle.c_str());

   TLatex Latex;
   Latex.SetTextFont(42);
   Latex.SetTextSize(0.03);
   Latex.SetNDC();

   Latex.DrawLatex(0.03, 0.94, Form("%.0f %% < Centrality < %.0f %%", S.centL, S.centH));
   Latex.DrawLatex(0.03, 0.88, Form("%.1f < Z p_{T} < %.1f", S.ptL, S.ptH));
   Latex.DrawLatex(0.03, 0.82, Form("%.1f < Track p_{T} < %.1f", S.TptL, S.TptH));

   Latex.SetTextSize(0.05);
   Latex.SetTextAlign(32);
   Latex.DrawLatex(0.97, 0.94, H1->GetTitle());
   
   TPad *P2 = (TPad *)Canvas.cd(2);
   P2->SetTheta(60.839);
   P2->SetPhi(38.0172);

   if(Identifier.size() >= 2)
      H2->SetTitle(Identifier[1].c_str());
   else
      H2->SetTitle("Data");
   H2->Draw(Option.c_str());
   H2->GetXaxis()->SetTitle(XTitle.c_str());
   H2->GetYaxis()->SetTitle(YTitle.c_str());
   H2->GetXaxis()->SetTitleSize(30);
   H2->GetYaxis()->SetTitleSize(30);
   H2->GetXaxis()->SetTitleOffset(3.0);
   H2->GetYaxis()->SetTitleOffset(2.5);
   H2->GetXaxis()->SetNdivisions(50205, kFALSE);
   H2->GetZaxis()->SetTitle(ZTitle.c_str());

   Latex.SetTextSize(0.05);
   Latex.SetTextAlign(32);
   Latex.DrawLatex(0.97, 0.94, H2->GetTitle());
   
   // ptN0->Draw();

   gSystem->Exec(Form("mkdir -p %s", Base.c_str()));
   Canvas.SaveAs(Form("%s/%s.png", Base.c_str(), Tag.c_str())); 
   // Canvas.SaveAs(Form("%s/pdf/%s.pdf", Base.c_str(), Tag.c_str())); 
   // Canvas.SaveAs(Form("%s/C/%s.C", Base.c_str(), Tag.c_str())); 
}

void Draw2DPlot3Panel(TH2D *H1, TH2D *H2, TH2D *H3, string XTitle, string YTitle, string ZTitle,
   Setting &S, string Base, string Tag, string Option, vector<string> Identifier)
{
   if(H1 == nullptr || H2 == nullptr || H3 == nullptr)
      return;

   TCanvas Canvas("Canvas", "", 2000, 800);

   Canvas.Divide(3);
   TPad *P1 = (TPad *)Canvas.cd(1);
   P1->SetTheta(60.839);
   P1->SetPhi(38.0172);

   if(Identifier.size() >= 1)
      H1->SetTitle(Identifier[0].c_str());
   else
      H1->SetTitle("PbPb MC");
   H1->Draw(Option.c_str());
   H1->GetXaxis()->SetTitle(XTitle.c_str());
   H1->GetYaxis()->SetTitle(YTitle.c_str());
   H1->GetXaxis()->SetTitleSize(24);
   H1->GetYaxis()->SetTitleSize(24);
   H1->GetXaxis()->SetTitleOffset(3.0);
   H1->GetYaxis()->SetTitleOffset(2.5);
   H1->GetXaxis()->SetNdivisions(50205, kFALSE);
   H1->GetZaxis()->SetTitle(ZTitle.c_str());

   TLatex Latex;
   Latex.SetTextFont(42);
   Latex.SetTextSize(0.03);
   Latex.SetNDC();

   Latex.DrawLatex(0.03, 0.94, Form("%.0f %% < Centrality < %.0f %%", S.centL, S.centH));
   Latex.DrawLatex(0.03, 0.88, Form("%.1f < Z p_{T} < %.1f", S.ptL, S.ptH));
   Latex.DrawLatex(0.03, 0.82, Form("%.1f < Track p_{T} < %.1f", S.TptL, S.TptH));

   Latex.SetTextSize(0.05);
   Latex.SetTextAlign(32);
   Latex.DrawLatex(0.97, 0.94, H1->GetTitle());

   TPad *P2 = (TPad *)Canvas.cd(2);
   P2->SetTheta(60.839);
   P2->SetPhi(38.0172);

   if(Identifier.size() >= 2)
      H2->SetTitle(Identifier[1].c_str());
   else
      H2->SetTitle("Data");
   H2->Draw(Option.c_str());
   H2->GetXaxis()->SetTitle(XTitle.c_str());
   H2->GetYaxis()->SetTitle(YTitle.c_str());
   H2->GetXaxis()->SetTitleSize(24);
   H2->GetYaxis()->SetTitleSize(24);
   H2->GetXaxis()->SetTitleOffset(3.0);
   H2->GetYaxis()->SetTitleOffset(2.5);
   H2->GetXaxis()->SetNdivisions(50205, kFALSE);
   H2->GetZaxis()->SetTitle(ZTitle.c_str());
   
   Latex.SetTextSize(0.05);
   Latex.SetTextAlign(32);
   Latex.DrawLatex(0.97, 0.94, H2->GetTitle());

   TPad *P3 = (TPad *)Canvas.cd(3);
   P3->SetTheta(60.839);
   P3->SetPhi(38.0172);

   if(Identifier.size() >= 3)
      H3->SetTitle(Identifier[1].c_str());
   else
      H3->SetTitle("pp MC");
   H3->Draw(Option.c_str());
   H3->GetXaxis()->SetTitle(XTitle.c_str());
   H3->GetYaxis()->SetTitle(YTitle.c_str());
   H3->GetXaxis()->SetTitleSize(24);
   H3->GetYaxis()->SetTitleSize(24);
   H3->GetXaxis()->SetTitleOffset(3.0);
   H3->GetYaxis()->SetTitleOffset(2.5);
   H3->GetXaxis()->SetNdivisions(50205, kFALSE);
   H3->GetZaxis()->SetTitle(ZTitle.c_str());

   Latex.SetTextSize(0.05);
   Latex.SetTextAlign(32);
   Latex.DrawLatex(0.97, 0.94, H3->GetTitle());
   
   // ptN0->Draw();

   gSystem->Exec(Form("mkdir -p %s", Base.c_str()));
   Canvas.SaveAs(Form("%s/%s.png", Base.c_str(), Tag.c_str())); 
   // Canvas.SaveAs(Form("%s/pdf/%s.pdf", Base.c_str(), Tag.c_str())); 
   // Canvas.SaveAs(Form("%s/C/%s.C", Base.c_str(), Tag.c_str())); 
}
