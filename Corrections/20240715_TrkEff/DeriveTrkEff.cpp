#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TTree.h"

#include "CommandLine.h"
#include "ProgressBar.h"
#include "PlotHelper4.h"
#include "Messenger.h"
#include "CommonFunctions.h"
#include "SetStyle.h"

// Include necessary headers for weight calculations and corrections
#include "trackingEfficiency2017pp.h"
#include "trackingEfficiency2018PbPb.h"
#include "TrackResidualCorrector.h"

#include "CustomAssert.h"

#define MAX 10000
#define E_MASS 0.0005111

TH1D* ProjectX(TH3* h3, const char* name = "_px", Option_t* option = "")
{
    return h3->ProjectionX(name, 1, h3->GetNbinsY(), 1, h3->GetNbinsZ(), option);
}

TH1D* ProjectY(TH3* h3, const char* name = "_py", Option_t* option = "")
{
    return h3->ProjectionY(name, 1, h3->GetNbinsX(), 1, h3->GetNbinsZ(), option);
}

TH1D* ProjectZ(TH3* h3, const char* name = "_pz", Option_t* option = "")
{
    return h3->ProjectionZ(name, 1, h3->GetNbinsX(), 1, h3->GetNbinsY(), option);
}

int main(int argc, char *argv[])
{
   SetThesisStyle();

   CommandLine CL(argc, argv);

   std::vector<double> PTs(49);
   std::generate_n(PTs.begin(), 49, [n = 0.2]() mutable { return n += 0.2; });

   std::vector<double> Etas(51);
   std::generate_n(Etas.begin(), 51, [n = -2.4-0.096]() mutable { return n += 0.096; });

   std::vector<double> Phis(51);
   std::generate_n(Phis.begin(), 51, [n = -0.04*M_PI]() mutable { return n += 0.04*M_PI; });

   vector<string> InputFileNames = CL.GetStringVector("Input");
   string OutputFileName         = CL.Get("Output", "TrackEfficiencyPlots.pdf");
   string RootOutputFileName     = CL.Get("RootOutput", "TrkEffPbPb.root");
   double Fraction               = CL.GetDouble("Fraction", 1.00);
   bool IsPP                     = CL.GetBool("IsPP", false);
   bool DoTrackEfficiency        = CL.GetBool("DoTrackEfficiency", true);
   bool DoTrackResidual          = CL.GetBool("DoTrackResidual", true);
   string TrackEfficiencyPath         = (DoTrackEfficiency == true) ? CL.Get("TrackEfficiencyPath") : "";
   vector<string> TrackResidualPath   = (DoTrackResidual == true) ? CL.GetStringVector("TrackResidualPath") : vector<string>{"", "", "", ""};
   bool DoMCHiBinShift                = IsPP ? CL.GetBool("DoMCHiBinShift", true): false;
   double MCHiBinShift                = DoMCHiBinShift ? CL.GetDouble("MCHiBinShift", 3) : 0;
   double MinTrackPT                   = CL.GetDouble("MinTrackPT", 0.500);
   double MaxTrackPT                   = CL.GetDouble("MaxTrackPT", 1000.00);
   int MinHiBin                        = CL.GetInt("MinHiBin", 0);
   int MaxHiBin                        = CL.GetInt("MaxHiBin", 200);
   bool DoIteration                    = (DoTrackResidual == true) ? CL.GetBool("DoIteration", false): false;

   string GGTreeName                  = IsPP ? "ggHiNtuplizerGED/EventTree" : "ggHiNtuplizer/EventTree";
   GGTreeName                         = CL.Get("GGTree", GGTreeName);

   bool DoZSelection                  = CL.GetBool("DoZSelection", false);
   bool DoElectron                    = CL.GetBool("DoElectron", true);


   if(DoTrackResidual == true)
      Assert(TrackResidualPath.size() == 1 || TrackResidualPath.size() == 4, "You need 1 file for residual correction or 4 files for centrality-dependence");

   sort(Etas.begin(), Etas.end());
   sort(Phis.begin(), Phis.end());
   sort(PTs.begin(), PTs.end());

   TFile *OutputFile = new TFile(RootOutputFileName.c_str(), "RECREATE");

   PdfFileHelper PdfFile(OutputFileName);
   PdfFile.AddTextPage("Track efficiency derivation");

   // Initialize track efficiency and residual correction objects
   TrkEff2017pp *TrackEfficiencyPP = nullptr;
   TrkEff2018PbPb *TrackEfficiencyPbPb = nullptr;
   if(DoTrackEfficiency)
   {
      if(IsPP)
         TrackEfficiencyPP = new TrkEff2017pp(false, TrackEfficiencyPath);
      else
         TrackEfficiencyPbPb = new TrkEff2018PbPb("general", "", false, TrackEfficiencyPath);
   }

   TrackResidualCentralityCorrector TrackResidual(TrackResidualPath);

   TFile* ResFile = nullptr;
   TH1D *hPtCorrTotal_old = nullptr;
   TH1D *hEtaCorrTotal_old = nullptr;
   TH1D *hPhiCorrTotal_old = nullptr;

   int ResNum = 0;
   if(MinHiBin == 0) ResNum = 0;
   else if(MinHiBin == 20) ResNum = 1; 
   else if(MinHiBin == 60)  ResNum = 2;
   else if(MinHiBin == 100) ResNum = 3;
   else ResNum = 0;       

   if(DoTrackResidual == true){
      ResFile = TFile::Open(TrackResidualPath[ResNum].c_str(), "READ");
      hPtCorrTotal_old = (TH1D*)ResFile->Get("hPtCorrTotal");
      hEtaCorrTotal_old = (TH1D*)ResFile->Get("hEtaCorrTotal");
      hPhiCorrTotal_old = (TH1D*)ResFile->Get("hPhiCorrTotal");
      hPtCorrTotal_old->SetName("hPtCorrTotal_old");
      hEtaCorrTotal_old->SetName("hEtaCorrTotal_old");
      hPhiCorrTotal_old->SetName("hPhiCorrTotal_old");
   }

   int NEta = Etas.size() - 1;
   int NPhi = Phis.size() - 1;
   int NPT = PTs.size() - 1;

   double EtaBins[MAX], PhiBins[MAX], PTBins[MAX];
   for(int i = 0; i <= NEta; i++)
      EtaBins[i] = Etas[i];
   for(int i = 0; i <= NPhi; i++)
      PhiBins[i] = Phis[i];
   for(int i = 0; i <= NPT; i++)
      PTBins[i] = PTs[i];

   //PTBins[0] = PTBins[1] * 0.5;
   //PTBins[NPT] = PTBins[NPT-1] * 1.25;

   TH3D HGenTrack("HGenTrack", "Gen Track;#eta;#phi;p_{T}", NEta, EtaBins, NPhi, PhiBins, NPT, PTBins);
   TH3D HRecoTrack("HRecoTrack", "Reco Track;#eta;#phi;p_{T}", NEta, EtaBins, NPhi, PhiBins, NPT, PTBins);
   TH3D HRecoTrackCorrected("HRecoTrackCorrected", "Reco Track Corrected;#eta;#phi;p_{T}", NEta, EtaBins, NPhi, PhiBins, NPT, PTBins);
   TH3D HRecoTrackResCorrected("HRecoTrackResCorrected", "Reco Track Residual Corrected;#eta;#phi;p_{T}", NEta, EtaBins, NPhi, PhiBins, NPT, PTBins);
   TH3D HEfficiency("HEfficiency", "Track Efficiency;#eta;#phi;p_{T}", NEta, EtaBins, NPhi, PhiBins, NPT, PTBins);
   TH3D HEfficiencyCorrected("HEfficiencyCorrected", "Corrected Track Efficiency;#eta;#phi;p_{T}", NEta, EtaBins, NPhi, PhiBins, NPT, PTBins);
   TH3D HEfficiencyResCorrected("HEfficiencyResCorrected", "Residual Corrected Track Efficiency;#eta;#phi;p_{T}", NEta, EtaBins, NPhi, PhiBins, NPT, PTBins);

   TH1D *hPtCorrTotal = new TH1D("hPtCorrTotal", "p_{T} Correction Factor", NPT, PTBins);
   TH1D *hEtaCorrTotal = new TH1D("hEtaCorrTotal", "#eta Correction Factor", NEta, EtaBins);
   TH1D *hPhiCorrTotal = new TH1D("hPhiCorrTotal", "#phi Correction Factor", NPhi, PhiBins);

   HGenTrack.SetStats(0);
   HRecoTrack.SetStats(0);
   HRecoTrackCorrected.SetStats(0);
   HRecoTrackResCorrected.SetStats(0);

   HEfficiency.SetStats(0);
   HEfficiencyCorrected.SetStats(0);
   HEfficiencyResCorrected.SetStats(0);

   HGenTrack.Sumw2();
   HRecoTrack.Sumw2();

   for(string InputFileName : InputFileNames)
   {
      cout << "Processing file " << InputFileName << endl;

      TFile InputFile(InputFileName.c_str());

      HiEventTreeMessenger MEvent(InputFile);
      GenParticleTreeMessenger MGen(InputFile);
      TrackTreeMessenger MSignalTrackPP(InputFile);
      PbPbTrackTreeMessenger MSignalTrack(InputFile);
      GGTreeMessenger          MSignalGG(InputFile, GGTreeName);

      int EntryCount = MEvent.GetEntries() * Fraction;
      ProgressBar Bar(cout, EntryCount);
      Bar.SetStyle(-1);
      for(int iE = 0; iE < EntryCount; iE++)
      {
         if(EntryCount < 200 || (iE % (EntryCount / 200)) == 0)
         {
            Bar.Update(iE);
            Bar.Print();
         }

         MEvent.GetEntry(iE);
         MGen.GetEntry(iE);
         MSignalGG.GetEntry(iE);

         if(IsPP == true)
            MSignalTrackPP.GetEntry(iE);
         else
            MSignalTrack.GetEntry(iE);

         if(IsPP == false && DoMCHiBinShift == true)
         {
            MEvent.hiBin = MEvent.hiBin - MCHiBinShift;   // MC shift
            if((MEvent.hiBin < MinHiBin) || (MEvent.hiBin >= MaxHiBin))   // out of range after shifting.  Skip!
               continue;
         }


         //Do Z selection

         bool isRecoZ = true, isGenZ = true;

         TLorentzVector VGenZ, VGenEle1, VGenEle2;

         if( DoZSelection == true && DoElectron == true && MSignalGG.NMC > 1 )
         {
            for(int igen1 = 0; igen1 < MSignalGG.NMC; igen1++){

               // We only want muon from Z's
               if(fabs(MSignalGG.MCPID->at(igen1)) != 11)
                  isGenZ = false;
               if(MSignalGG.MCMomPID->at(igen1) != 23)
                  isGenZ = false;
               if(MSignalGG.MCPt->at(igen1) < 20)
                  isGenZ = false;
               if(fabs(MSignalGG.MCEta->at(igen1)) > 2.1)
                  isGenZ = false;
   
               VGenEle1.SetPtEtaPhiM(MSignalGG.MCPt->at(igen1),
                     MSignalGG.MCEta->at(igen1),
                     MSignalGG.MCPhi->at(igen1),
                     E_MASS);

               for(int igen2 = igen1 + 1; igen2 < MSignalGG.NMC; igen2++){

                  // We only want electron from Z's
                  if(MSignalGG.MCPID->at(igen2) != -MSignalGG.MCPID->at(igen1))
                     isGenZ = false;
                  if(MSignalGG.MCMomPID->at(igen2) != 23)
                     isGenZ = false;
                  if(MSignalGG.MCPt->at(igen2) < 20)
                     isGenZ = false;
                  if(fabs(MSignalGG.MCEta->at(igen2)) > 2.1)
                     isGenZ = false;
   
                  VGenEle2.SetPtEtaPhiM(MSignalGG.MCPt->at(igen2),
                        MSignalGG.MCEta->at(igen2),
                        MSignalGG.MCPhi->at(igen2),
                        E_MASS);
   
                  VGenZ = VGenEle1 + VGenEle2;

                  if(VGenZ.M() < 60 || VGenZ.M() > 120)
                     isGenZ = false;

               }
            }// GenZ loop end
            
         }

         int N_eles;
            
         if(DoZSelection == true && DoElectron == true)
            N_eles = MSignalGG.NEle;
         else
            N_eles = 0;

         for(int iele1 = 0; iele1 < N_eles; iele1++)
         {

            // Some basic electron kinematic cuts
            if(fabs(MSignalGG.EleSCEta->at(iele1)) > 2.5)  isRecoZ = false;
            if(fabs(MSignalGG.EleEta->at(iele1)) > 2.1)    isRecoZ = false;
            if(fabs(MSignalGG.ElePt->at(iele1)) < 20)      isRecoZ = false;
            if(IsPP == false && MSignalGG.DielectronPassVetoCut(iele1, MZHadron.hiBin) == false) isRecoZ = false;
            if(IsPP == true  && MSignalGG.DielectronPassVetoCutPP(iele1) == false) isRecoZ = false;

            if(IsPP == false){ // per Kaya, HCAL failure gives rise to misidentified electrons.
               if(MSignalGG.EleSCEta->at(iele1) < -1.39 && MSignalGG.EleSCPhi->at(iele1) > -1.6 &&  MSignalGG.EleSCPhi->at(iele1) < -0.9 ) isRecoZ = false;
            }

            TLorentzVector Ele1;  
            Ele1.SetPtEtaPhiM(MSignalGG.ElePt->at(iele1), MSignalGG.EleEta->at(iele1), MSignalGG.ElePhi->at(iele1), E_MASS);

            // Loop over 2nd reco electron
            for(int iele2 = iele1+1; iele2 < N_eles; iele2++)
            {
               // We want opposite-charge electrons with some basic kinematic cuts
               if(MSignalGG.EleCharge->at(iele1) == MSignalGG.EleCharge->at(iele2))  isRecoZ = false;
               if(fabs(MSignalGG.EleSCEta->at(iele2)) > 2.5)                         isRecoZ = false;
               if(fabs(MSignalGG.EleEta->at(iele2)) > 2.1)                           isRecoZ = false;
               if(fabs(MSignalGG.ElePt->at(iele2)) < 20)                             isRecoZ = false;
               if(IsPP == false && MSignalGG.DielectronPassVetoCut(iele2, MZHadron.hiBin) == false)   isRecoZ = false;
               if(IsPP == true  && MSignalGG.DielectronPassVetoCutPP(iele2) == false)   isRecoZ = false;

               if(IsPP == false){ // per Kaya, HCAL failure gives rise to misidentified electrons.
                  if(MSignalGG.EleSCEta->at(iele2) < -1.39 && MSignalGG.EleSCPhi->at(iele2) > -1.6 &&  MSignalGG.EleSCPhi->at(iele2) < -0.9 ) isRecoZ = false;
               }

               TLorentzVector Ele2;  
               Ele2.SetPtEtaPhiM(MSignalGG.ElePt->at(iele2), MSignalGG.EleEta->at(iele2), MSignalGG.ElePhi->at(iele2), E_MASS);

               TLorentzVector Z = Ele1+Ele2;
               double Zmass = Z.M();

               if(Zmass < 60 || Zmass > 120) isRecoZ = false;
               //if(fabs(Z.Rapidity()) > 2.4) continue;

            } // Loop over 2nd reco electron ended

         } // Loop over 1st reco electron ended

         if(DoZSelection == true && DoElectron == true)
         {
            if(isGenZ == false || isRecoZ == false)
               continue;
         }


         // Loop over gen particles
         for(int iG = 0; iG < MGen.Mult; iG++)
         {
            // Apply gen-level cuts (e.g., charged particles only, stable particles only)
            if(MGen.Charge->at(iG) == 0)
               continue;
            if(MGen.DaughterCount->at(iG) > 0)
               continue;

            if(MGen.PT->at(iG) < MinTrackPT || MGen.PT->at(iG) > MaxTrackPT)
               continue;
            if(MGen.Eta->at(iG) < -2.4 || MGen.Eta->at(iG) > 2.4)
               continue;
     
            HGenTrack.Fill(MGen.Eta->at(iG), MGen.Phi->at(iG) < 0 ? MGen.Phi->at(iG)+2*M_PI : MGen.Phi->at(iG), MGen.PT->at(iG));
         }

         // Loop over reco tracks

         if(IsPP)
         {
            for(int iT = 0; iT < MSignalTrackPP.nTrk; iT++)
            {
               
               if(MSignalTrackPP.trkPt[iT] < MinTrackPT || MSignalTrackPP.trkPt[iT] > MaxTrackPT)
                  continue;

               if(!MSignalTrackPP.PassZHadron2022Cut(iT))
                  continue;

               double RecoEta = MSignalTrackPP.trkEta[iT];
               double RecoPhi = MSignalTrackPP.trkPhi[iT];
               double RecoPT = MSignalTrackPP.trkPt[iT];

               HRecoTrack.Fill(RecoEta, RecoPhi < 0 ? RecoPhi+2*M_PI : RecoPhi, RecoPT);

               // Calculate track efficiency correction
               double TreeTrackCorrection = 1;
               if(DoTrackEfficiency)
               {
                  TreeTrackCorrection = TrackEfficiencyPP->getCorrection(RecoPT, RecoEta);
               }

               // Calculate track residual correction
               double TreeTrackResidualCorrection = 1;
               if(DoTrackResidual)
               {
                  TreeTrackResidualCorrection = TrackResidual.GetCorrectionFactor(RecoPT, RecoEta, RecoPhi, MEvent.hiBin);
               }
               HRecoTrackCorrected.Fill(RecoEta, RecoPhi < 0 ? RecoPhi+2*M_PI : RecoPhi, RecoPT, TreeTrackCorrection );
               HRecoTrackResCorrected.Fill(RecoEta, RecoPhi < 0 ? RecoPhi+2*M_PI : RecoPhi, RecoPT, TreeTrackCorrection*TreeTrackResidualCorrection);
            }
         }
         else
         {
            for(int iT = 0; iT < MSignalTrack.TrackPT->size(); iT++)
            {

               if(MSignalTrack.TrackPT->at(iT) < MinTrackPT || MSignalTrack.TrackPT->at(iT) > MaxTrackPT)
                  continue;

               if(!MSignalTrack.PassZHadron2022Cut(iT))
                  continue;

               double RecoEta = MSignalTrack.TrackEta->at(iT);
               double RecoPhi = MSignalTrack.TrackPhi->at(iT);
               double RecoPT = MSignalTrack.TrackPT->at(iT);

               HRecoTrack.Fill(RecoEta, RecoPhi < 0 ? RecoPhi+2*M_PI : RecoPhi, RecoPT);

               // Calculate track efficiency correction
               double TreeTrackCorrection = 1;
               if(DoTrackEfficiency)
               {
                  TreeTrackCorrection = TrackEfficiencyPbPb->getCorrection(RecoPT, RecoEta, MEvent.hiBin + MCHiBinShift );
               }

               // Calculate track residual correction
               double TreeTrackResidualCorrection = 1;
               if(DoTrackResidual)
               {
                  TreeTrackResidualCorrection = TrackResidual.GetCorrectionFactor(RecoPT, RecoEta, RecoPhi, MEvent.hiBin );
               }

               HRecoTrackCorrected.Fill(RecoEta, RecoPhi < 0 ? RecoPhi+2*M_PI : RecoPhi, RecoPT, TreeTrackCorrection);
               HRecoTrackResCorrected.Fill(RecoEta, RecoPhi < 0 ? RecoPhi+2*M_PI : RecoPhi, RecoPT, TreeTrackCorrection*TreeTrackResidualCorrection);

            }
         }
      }
      Bar.Update(EntryCount);
      Bar.Print();
      Bar.PrintLine();

      InputFile.Close();
   }

   // Calculate efficiency
   HEfficiency.Divide(&HRecoTrack, &HGenTrack);
   HEfficiencyCorrected.Divide(&HRecoTrackCorrected, &HGenTrack);
   HEfficiencyResCorrected.Divide(&HRecoTrackResCorrected, &HGenTrack);

   TH1D *HEfficiency_eta = ProjectX(&HRecoTrack, "HEfficiency_eta"); 
   HEfficiency_eta->Divide(ProjectX(&HGenTrack));
   HEfficiency_eta->SetTitle("#eta efficiency");
   HEfficiency_eta->SetStats(0);
   TH1D *HEfficiency_phi = ProjectY(&HRecoTrack, "HEfficiency_phi");
   HEfficiency_phi->Divide(ProjectY(&HGenTrack));
   HEfficiency_phi->SetTitle("#phi efficiency");
   HEfficiency_phi->SetStats(0);
   TH1D *HEfficiency_pt  = ProjectZ(&HRecoTrack, "HEfficiency_pt");
   HEfficiency_pt->Divide(ProjectZ(&HGenTrack));
   HEfficiency_pt->SetTitle("p_{T} efficiency");
   HEfficiency_pt->SetStats(0);

   TH1D *HEfficiencyCorrected_eta = ProjectX(&HRecoTrackCorrected, "HEfficiencyCorrected_eta");
   HEfficiencyCorrected_eta->Divide(ProjectX(&HGenTrack));
   HEfficiencyCorrected_eta->SetTitle("#eta corrected efficiency");
   HEfficiencyCorrected_eta->SetStats(0);
   TH1D *HEfficiencyCorrected_phi = ProjectY(&HRecoTrackCorrected, "HEfficiencyCorrected_phi");
   HEfficiencyCorrected_phi->Divide(ProjectY(&HGenTrack));
   HEfficiencyCorrected_phi->SetTitle("#phi corrected efficiency");
   HEfficiencyCorrected_phi->SetStats(0);
   TH1D *HEfficiencyCorrected_pt = ProjectZ(&HRecoTrackCorrected, "HEfficiencyCorrected_pt");
   HEfficiencyCorrected_pt->Divide(ProjectZ(&HGenTrack));
   HEfficiencyCorrected_pt->SetTitle("p_{T} corrected efficiency");
   HEfficiencyCorrected_pt->SetStats(0);
   
   
   TH1D *HEfficiencyResCorrected_eta = ProjectX(&HRecoTrackResCorrected,"HEfficiencyResCorrected_eta");
   HEfficiencyResCorrected_eta->Divide(ProjectX(&HGenTrack));
   HEfficiencyResCorrected_eta->SetTitle("#eta residual corrected efficiency");
   HEfficiencyResCorrected_eta->SetStats(0);
   TH1D *HEfficiencyResCorrected_phi = ProjectY(&HRecoTrackResCorrected,"HEfficiencyResCorrected_phi");
   HEfficiencyResCorrected_phi->Divide(ProjectY(&HGenTrack));
   HEfficiencyResCorrected_phi->SetTitle("#phi residual corrected efficiency");
   HEfficiencyResCorrected_phi->SetStats(0);
   TH1D *HEfficiencyResCorrected_pt = ProjectZ(&HRecoTrackResCorrected,"HEfficiencyResCorrected_pt");
   HEfficiencyResCorrected_pt->Divide(ProjectZ(&HGenTrack));
   HEfficiencyResCorrected_pt->SetTitle("p_{T} residual corrected efficiency");
   HEfficiencyResCorrected_pt->SetStats(0);

   // Calculate correction factors  

   for(int i = 1; i <= HEfficiencyResCorrected_eta->GetNbinsX(); i++)
      if(DoIteration && HEfficiencyResCorrected_eta->GetNbinsX() == hEtaCorrTotal_old->GetNbinsX())
         hEtaCorrTotal->SetBinContent(i, HEfficiencyResCorrected_eta->GetBinContent(i) > 0 ? hEtaCorrTotal_old->GetBinContent(i) / HEfficiencyResCorrected_eta->GetBinContent(i) : 1.0);
      else
         hEtaCorrTotal->SetBinContent(i, HEfficiencyCorrected_eta->GetBinContent(i) > 0 ? 1.0 / HEfficiencyCorrected_eta->GetBinContent(i) : 1.0);

   for(int i = 1; i <= HEfficiencyResCorrected_phi->GetNbinsX(); i++)
      if(DoIteration && HEfficiencyResCorrected_phi->GetNbinsX() == hPhiCorrTotal_old->GetNbinsX())
         hPhiCorrTotal->SetBinContent(i, HEfficiencyResCorrected_phi->GetBinContent(i) > 0 ? hPhiCorrTotal_old->GetBinContent(i) / HEfficiencyResCorrected_phi->GetBinContent(i) : 1.0);
      else
         hPhiCorrTotal->SetBinContent(i, HEfficiencyCorrected_phi->GetBinContent(i) > 0 ? 1.0 / HEfficiencyCorrected_phi->GetBinContent(i) : 1.0);

   for(int i = 1; i <= HEfficiencyResCorrected_pt->GetNbinsX(); i++)
      if(DoIteration && HEfficiencyResCorrected_pt->GetNbinsX() == hPtCorrTotal_old->GetNbinsX())
         hPtCorrTotal->SetBinContent(i, HEfficiencyResCorrected_pt->GetBinContent(i) > 0 ? hPtCorrTotal_old->GetBinContent(i) / HEfficiencyResCorrected_pt->GetBinContent(i) : 1.0);
      else
         hPtCorrTotal->SetBinContent(i, HEfficiencyCorrected_pt->GetBinContent(i) > 0 ? 1.0 / HEfficiencyCorrected_pt->GetBinContent(i) : 1.0);


   // Normalize correction factors

   //double ptIntegral_i = HEfficiencyResCorrected_pt->Integral();
   //double etaIntegral_i = HEfficiencyResCorrected_eta->Integral();
   //double phiIntegral_i = HEfficiencyResCorrected_phi->Integral();

   double ptIntegral = hPtCorrTotal->Integral();
   double etaIntegral = hEtaCorrTotal->Integral();
   double phiIntegral = hPhiCorrTotal->Integral();

   double recoEffIntegral = HRecoTrackCorrected.Integral();
   double recoResIntegral = HRecoTrackResCorrected.Integral();
   double genIntegral = HGenTrack.Integral();

   double CorrTotal = (hPtCorrTotal->GetNbinsX() / ptIntegral) * (hEtaCorrTotal->GetNbinsX() / etaIntegral) * (hPhiCorrTotal->GetNbinsX() / phiIntegral);
   std::cout<<"CorrTotal = "<<1/CorrTotal<<std::endl;

   //hPtCorrTotal ->Scale(pow(CorrTotal,2/9.));
   //hEtaCorrTotal->Scale(pow(CorrTotal,2/9.));
   //hPhiCorrTotal->Scale(pow(CorrTotal,2/9.));

   //hPtCorrTotal ->Scale(pow((hEtaCorrTotal->GetNbinsX() / etaIntegral) * (hPhiCorrTotal->GetNbinsX() / phiIntegral),1/3.));
   //hEtaCorrTotal->Scale(pow((hPtCorrTotal->GetNbinsX() / ptIntegral) * (hPhiCorrTotal->GetNbinsX() / phiIntegral),1/3.));
   //hPhiCorrTotal->Scale(pow((hPtCorrTotal->GetNbinsX() / ptIntegral) * (hEtaCorrTotal->GetNbinsX() / etaIntegral),1/3.));

   hEtaCorrTotal->Scale(hEtaCorrTotal->GetNbinsX() / etaIntegral);
   hPhiCorrTotal->Scale(hPhiCorrTotal->GetNbinsX() / phiIntegral);

   // Add plots to PDF
   PdfFile.AddTextPage("Track efficiency plots: #eta");

   PdfFile.AddPlot(HGenTrack.ProjectionX(), "hist text00", true);
   PdfFile.AddPlot(HRecoTrack.ProjectionX(), "hist text00", true);
   PdfFile.AddPlot(HEfficiency_eta, "hist e text00", false);
   PdfFile.AddPlot(HEfficiencyCorrected_eta, "hist e text00", false);
   PdfFile.AddPlot(HEfficiencyResCorrected_eta, "hist e text00", false);

   PdfFile.AddTextPage("Track efficiency plots: #phi");

   PdfFile.AddPlot(HGenTrack.ProjectionY(), "hist text00", true);
   PdfFile.AddPlot(HRecoTrack.ProjectionY(), "hist text00", true);
   PdfFile.AddPlot(HEfficiency_phi, "hist e text00", false);
   PdfFile.AddPlot(HEfficiencyCorrected_phi, "hist e text00", false);
   PdfFile.AddPlot(HEfficiencyResCorrected_phi, "hist e text00", false);

   PdfFile.AddTextPage("Track efficiency plots: p_{T}");

   PdfFile.AddPlot(HGenTrack.ProjectionZ(), "hist text00", true);
   PdfFile.AddPlot(HRecoTrack.ProjectionZ(), "hist text00", true);
   PdfFile.AddPlot(HEfficiency_pt, "hist e text00", false);
   PdfFile.AddPlot(HEfficiencyCorrected_pt, "hist e text00", false);
   PdfFile.AddPlot(HEfficiencyResCorrected_pt, "hist e text00", false);

   PdfFile.AddTimeStampPage();
   PdfFile.Close();

   OutputFile->cd();
   HGenTrack.Write();
   HRecoTrack.Write();
   HRecoTrackCorrected.Write();
   HRecoTrackResCorrected.Write();
   
   HEfficiency.Write();
   HEfficiencyCorrected.Write();
   HEfficiencyResCorrected.Write();

   hPtCorrTotal->Write();
   hEtaCorrTotal->Write();
   hPhiCorrTotal->Write();

   hEtaCorrTotal_old->Write();
   hPhiCorrTotal_old->Write();
   hPtCorrTotal_old->Write();
   
   OutputFile->Close();

   delete OutputFile;  // Clean up

   return 0;
}