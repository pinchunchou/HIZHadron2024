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
#include "TLorentzVector.h"

// Include necessary headers for weight calculations and corrections
#include "trackingEfficiency2017pp.h"
#include "trackingEfficiency2018PbPb.h"
#include "TrackResidualCorrector.h"

#include "CustomAssert.h"

#define MAX 10000
#define E_MASS 0.0005111
#define Z_MASS 91.1880

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

// Updated function for error calculation
TH1D* CalculateEfficiencyAndError(const char* name, TH1D* hNum, TH1D* hDen) {

   TH1D* hResult = (TH1D*)hNum->Clone(name);
    for (int i = 1; i <= hNum->GetNbinsX(); ++i) {
        double num = hNum->GetBinContent(i);
        double den = hDen->GetBinContent(i);
        
        // Efficiency
        double eff = (den > 0) ? num / den : 0;
        
        // Error calculation using binomial error propagation
        double err = (den > 0 && eff > 0 && eff < 1) ? sqrt(eff * (1 - eff) / den) : 0;
        
        hResult->SetBinContent(i, eff);
        hResult->SetBinError(i, err);
    }

    hResult->SetStats(0);
    return hResult;
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

   double MinZPT                   = CL.GetDouble("MinZPT", 40.00);
   double MaxZPT                   = CL.GetDouble("MaxZPT", 200.00);

   double LeptonVeto                    = CL.GetDouble("LeptonVeto", 0.01);


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

         bool isRecoZ = false, isGenZ = false;

         TLorentzVector VGenZ, VGenEle1, VGenEle2;

         //std::cout << "MSignalGG.NMC = "<<MSignalGG.NMC<<std::endl;

         std::vector<double> ZEta_gen, ZPhi_gen, ZMass_gen, Zpt_gen;
         std::vector<double> ZEta_reco, ZPhi_reco, ZMass_reco, Zpt_reco;
         std::vector<double> eleEta1_gen, elePhi1_gen, elePt1_gen;
         std::vector<double> eleEta2_gen, elePhi2_gen, elePt2_gen;
         std::vector<double> eleEta1_reco, elePhi1_reco, elePt1_reco;
         std::vector<double> eleEta2_reco, elePhi2_reco, elePt2_reco;

         if( DoZSelection == true && DoElectron == true && MSignalGG.NMC > 1 )
         {
            for(int igen1 = 0; igen1 < MSignalGG.NMC; igen1++){


               //std::cout << "MSignalGG.MCPID->at(igen1) = "<<MSignalGG.MCPID->at(igen1)<<std::endl;
               //std::cout << "MSignalGG.MCMomPID->at(igen1) = "<<MSignalGG.MCMomPID->at(igen1)<<std::endl;
               //std::cout << "MSignalGG.MCPt->at(igen1) = "<<MSignalGG.MCPt->at(igen1)<<std::endl;
               //std::cout << "MSignalGG.MCEta->at(igen1) = "<<MSignalGG.MCEta->at(igen1)<<std::endl;

               //if(fabs(MSignalGG.MCPID->at(igen1)) == 11 && MSignalGG.MCMomPID->at(igen1) == 23)
               //   std::cout<<"hi1!"<<std::endl;

               if(igen1 > MSignalGG.MCPID->size() - 1){
                  cerr<<"Warning: MSignalGG.NMC and igen1 > MSignalGG.MCPID->size(): "<< MSignalGG.NMC <<" or "<<igen1<<" > "<< MSignalGG.MCPID->size()<<endl;
                  break;
               }

               // We only want electron from Z's
               if(fabs(MSignalGG.MCPID->at(igen1)) != 11)
                  continue;
               if(MSignalGG.MCMomPID->at(igen1) != 23 && ((MSignalGG.MCMomPID->at(igen1) != MSignalGG.MCPID->at(igen1)) || MSignalGG.MCGMomPID->at(igen1) != 23) )
                  continue;
               if(MSignalGG.MCPt->at(igen1) < 20)
                  continue;
               if(fabs(MSignalGG.MCEta->at(igen1)) > 2.1)
                  continue;
   
               VGenEle1.SetPtEtaPhiM(MSignalGG.MCPt->at(igen1),
                     MSignalGG.MCEta->at(igen1),
                     MSignalGG.MCPhi->at(igen1),
                     E_MASS);

               for(int igen2 = igen1 + 1; igen2 < MSignalGG.NMC; igen2++){

                  //std::cout << "MSignalGG.MCPID->at(igen2) = "<<MSignalGG.MCPID->at(igen2)<<std::endl;
                  //std::cout << "MSignalGG.MCMomPID->at(igen2) = "<<MSignalGG.MCMomPID->at(igen2)<<std::endl;
                  //std::cout << "MSignalGG.MCPt->at(igen2) = "<<MSignalGG.MCPt->at(igen2)<<std::endl;
                  //std::cout << "MSignalGG.MCEta->at(igen2) = "<<MSignalGG.MCEta->at(igen2)<<std::endl;

                  //if(fabs(MSignalGG.MCPID->at(igen2)) == 11 && MSignalGG.MCMomPID->at(igen2) == 23)
                  //   std::cout<<"hi2!"<<std::endl;

                  if(igen2 > MSignalGG.MCPID->size() - 1){
                     cout<<"Warning: MSignalGG.NMC and igen1 > MSignalGG.MCPID->size():"<< MSignalGG.NMC <<" or "<<igen2<<" > "<< MSignalGG.MCPID->size()<<endl;
                     break;
                  }

                  // We only want electron from Z's
                  if(MSignalGG.MCPID->at(igen2) != -MSignalGG.MCPID->at(igen1))
                     continue;
                  if(MSignalGG.MCMomPID->at(igen2) != 23 && ((MSignalGG.MCMomPID->at(igen2) != MSignalGG.MCPID->at(igen2)) || MSignalGG.MCGMomPID->at(igen2) != 23) )
                     continue;
                  if(MSignalGG.MCPt->at(igen2) < 20)
                     continue;
                  if(fabs(MSignalGG.MCEta->at(igen2)) > 2.1)
                     continue;
   
                  VGenEle2.SetPtEtaPhiM(MSignalGG.MCPt->at(igen2),
                        MSignalGG.MCEta->at(igen2),
                        MSignalGG.MCPhi->at(igen2),
                        E_MASS);
   
                  VGenZ = VGenEle1 + VGenEle2;

                  //std::cout << "VGenZ.M() = "<<VGenZ.M()<<std::endl;

                  if(VGenZ.M() < 60 || VGenZ.M() > 120)
                     continue;

                  if(VGenZ.Pt() < MinZPT || VGenZ.Pt() > MaxZPT) continue;

                  ZEta_gen.push_back(VGenZ.Eta());
                  ZPhi_gen.push_back(VGenZ.Phi());
                  ZMass_gen.push_back(VGenZ.M());
                  Zpt_gen.push_back(VGenZ.Pt());

                  eleEta1_gen.push_back(MSignalGG.MCEta->at(igen1));
                  elePhi1_gen.push_back(MSignalGG.MCPhi->at(igen1));
                  elePt1_gen.push_back(MSignalGG.MCPt->at(igen1));

                  eleEta2_gen.push_back(MSignalGG.MCEta->at(igen2));
                  elePhi2_gen.push_back(MSignalGG.MCPhi->at(igen2));
                  elePt2_gen.push_back(MSignalGG.MCPt->at(igen2));


                  isGenZ = true;
               }
            }// GenZ loop end
            
         }

         int N_eles;
            
         if(DoZSelection == true && DoElectron == true)
            N_eles = MSignalGG.NEle;
         else
            N_eles = 0;

         if(DoZSelection == true && DoElectron == true && N_eles > MSignalGG.ElePt->size() ){
            cerr<<"Warning: MSignalGG.NEle and N_eles > MSignalGG.ElePt->size(): "<< MSignalGG.NEle <<" or "<<N_eles<<" > "<< MSignalGG.ElePt->size()<<endl;
            N_eles = MSignalGG.ElePt->size();
         }

         for(int iele1 = 0; iele1 < N_eles; iele1++)
         {

            // Some basic electron kinematic cuts
            if(fabs(MSignalGG.EleSCEta->at(iele1)) > 2.5)  continue;
            if(fabs(MSignalGG.EleEta->at(iele1)) > 2.1)    continue;
            if(fabs(MSignalGG.ElePt->at(iele1)) < 20)      continue;
            if(IsPP == false && MSignalGG.DielectronPassVetoCut(iele1, MEvent.hiBin) == false) continue;
            if(IsPP == true  && MSignalGG.DielectronPassVetoCutPP(iele1) == false) continue;

            if(IsPP == false){ // per Kaya, HCAL failure gives rise to misidentified electrons.
               if(MSignalGG.EleSCEta->at(iele1) < -1.39 && MSignalGG.EleSCPhi->at(iele1) > -1.6 &&  MSignalGG.EleSCPhi->at(iele1) < -0.9 ) continue;
            }

            TLorentzVector Ele1;  
            Ele1.SetPtEtaPhiM(MSignalGG.ElePt->at(iele1), MSignalGG.EleEta->at(iele1), MSignalGG.ElePhi->at(iele1), E_MASS);

            // Loop over 2nd reco electron
            for(int iele2 = iele1+1; iele2 < N_eles; iele2++)
            {
               // We want opposite-charge electrons with some basic kinematic cuts
               if(MSignalGG.EleCharge->at(iele1) == MSignalGG.EleCharge->at(iele2))  continue;
               if(fabs(MSignalGG.EleSCEta->at(iele2)) > 2.5)                         continue;
               if(fabs(MSignalGG.EleEta->at(iele2)) > 2.1)                           continue;
               if(fabs(MSignalGG.ElePt->at(iele2)) < 20)                             continue;
               if(IsPP == false && MSignalGG.DielectronPassVetoCut(iele2, MEvent.hiBin) == false)   continue;
               if(IsPP == true  && MSignalGG.DielectronPassVetoCutPP(iele2) == false)   continue;

               if(IsPP == false){ // per Kaya, HCAL failure gives rise to misidentified electrons.
                  if(MSignalGG.EleSCEta->at(iele2) < -1.39 && MSignalGG.EleSCPhi->at(iele2) > -1.6 &&  MSignalGG.EleSCPhi->at(iele2) < -0.9 ) continue;
               }

               TLorentzVector Ele2;  
               Ele2.SetPtEtaPhiM(MSignalGG.ElePt->at(iele2), MSignalGG.EleEta->at(iele2), MSignalGG.ElePhi->at(iele2), E_MASS);

               TLorentzVector Z = Ele1+Ele2;
               double Zmass = Z.M();

               //std::cout << "Z.M() = "<<Z.M()<<std::endl;

               if(Zmass < 60 || Zmass > 120) continue;
               //if(fabs(Z.Rapidity()) > 2.4) continue;

               if(Z.Pt() < MinZPT || Z.Pt() > MaxZPT) continue;

               ZEta_reco.push_back(Z.Eta());
               ZPhi_reco.push_back(Z.Phi());
               ZMass_reco.push_back(Z.M());
               Zpt_reco.push_back(Z.Pt());

               eleEta1_reco.push_back(MSignalGG.EleEta->at(iele1));
               elePhi1_reco.push_back(MSignalGG.ElePhi->at(iele1));
               elePt1_reco.push_back(MSignalGG.ElePt->at(iele1));

               eleEta2_reco.push_back(MSignalGG.EleEta->at(iele2));
               elePhi2_reco.push_back(MSignalGG.ElePhi->at(iele2));
               elePt2_reco.push_back(MSignalGG.ElePt->at(iele2));

               isRecoZ = true;

            } // Loop over 2nd reco electron ended

         } // Loop over 1st reco electron ended


         

         if(DoZSelection == true && DoElectron == true)
         {
            if(isGenZ == false || isRecoZ == false)
               continue;
         }

         float Zmass_temp = -1;

         float bestEtaGen, bestPhiGen, bestEtaReco, bestPhiReco;

         int bestZidx = -1, bestZgenIdx = -1;

         float bestEta1Gen , bestPhi1Gen , bestPt1Gen ;
         float bestEta2Gen , bestPhi2Gen , bestPt2Gen ;
         float bestEta1Reco, bestPhi1Reco, bestPt1Reco;
         float bestEta2Reco, bestPhi2Reco, bestPt2Reco;

         // Loop to select the best Z mass
         for(int idx=0; idx < size(ZMass_reco) ; idx++){
            if((Zpt_reco[idx] < MinZPT) || (Zpt_reco[idx] > MaxZPT)){
                cout<<"ZpT might be wrong!"<<endl;
                continue;
            }

            if(abs(Zmass_temp - Z_MASS) > abs(ZMass_reco[idx] - Z_MASS)){
               Zmass_temp = ZMass_reco[idx];
               bestZidx = idx;
               bestEtaReco = ZEta_reco[idx];
               bestPhiReco = ZPhi_reco[idx];

               bestEta1Reco = eleEta1_reco[idx];
               bestPhi1Reco = elePhi1_reco[idx];
               bestPt1Reco = elePt1_reco[idx];

               bestEta2Reco = eleEta2_reco[idx];
               bestPhi2Reco = elePhi2_reco[idx];
               bestPt2Reco = elePt2_reco[idx];


            }
         }


         Zmass_temp = -1;

         for(int idx=0; idx < size(ZMass_gen) ; idx++){
            if((Zpt_gen[idx] < MinZPT) || (Zpt_gen[idx] > MaxZPT)){
               cout<<"ZpT might be wrong!"<<endl;  
               continue;
            }

            if(abs(Zmass_temp - Z_MASS) > abs(ZMass_gen[idx] - Z_MASS)){
               Zmass_temp = ZMass_gen[idx];
               bestZgenIdx = idx;
               bestEtaGen = ZEta_gen[idx];
               bestPhiGen = ZPhi_gen[idx];

               bestEta1Gen = eleEta1_gen[idx];
               bestPhi1Gen = elePhi1_gen[idx];
               bestPt1Gen = elePt1_gen[idx];

               bestEta2Gen = eleEta2_gen[idx];
               bestPhi2Gen = elePhi2_gen[idx];
               bestPt2Gen = elePt2_gen[idx];

            }
         } 


         // Loop over gen particles

         int NGenTrack = (MGen.PT->size() > MGen.Mult) ? MGen.Mult : MGen.PT->size(); // To prevent some weird out_of_range errors.
         if(MGen.PT->size() < MGen.Mult){
            cerr << "Warning: Less Gen tracks than Mult: " << MGen.PT->size() << " < " << MGen.Mult << endl;
         }

         for(int iG = 0; iG < NGenTrack; iG++)
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

            if(DoZSelection == true){

               if(bestZgenIdx == -1) continue;

               double DeltaEtaEle1 = bestEta1Gen - MGen.Eta->at(iG);
               double DeltaPhiEle1 = DeltaPhi(bestPhi1Gen, MGen.Phi->at(iG));
               double DeltaR1 = sqrt(DeltaEtaEle1*DeltaEtaEle1 + DeltaPhiEle1*DeltaPhiEle1);
            
               double DeltaEtaEle2 = bestEta2Gen - MGen.Eta->at(iG);
               double DeltaPhiEle2 = DeltaPhi(bestPhi2Gen, MGen.Phi->at(iG));
               double DeltaR2 = sqrt(DeltaEtaEle2*DeltaEtaEle2 + DeltaPhiEle2*DeltaPhiEle2);
            
               if(DeltaR1 <  LeptonVeto || DeltaR2 < LeptonVeto)
                  continue;
            }

            HGenTrack.Fill(MGen.Eta->at(iG), MGen.Phi->at(iG) < 0 ? MGen.Phi->at(iG)+2*M_PI : MGen.Phi->at(iG), MGen.PT->at(iG));

         }

         // Loop over reco tracks

         //cout<<"Loop over reco tracks"<<endl;

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

               if(DoZSelection == true){

                  if(bestZidx == -1) continue;

                  double DeltaEtaEle1 = bestEta1Reco - RecoEta;
                  double DeltaPhiEle1 = DeltaPhi(bestPhi1Reco, RecoPhi);
                  double DeltaR1 = sqrt(DeltaEtaEle1*DeltaEtaEle1 + DeltaPhiEle1*DeltaPhiEle1);
               
                  double DeltaEtaEle2 = bestEta2Reco - RecoEta;
                  double DeltaPhiEle2 = DeltaPhi(bestPhi2Reco, RecoPhi);
                  double DeltaR2 = sqrt(DeltaEtaEle2*DeltaEtaEle2 + DeltaPhiEle2*DeltaPhiEle2);
               
                  if(DeltaR1 <  LeptonVeto || DeltaR2 < LeptonVeto)
                     continue;
               }

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
                  TreeTrackResidualCorrection = TrackResidual.GetCorrectionFactor(RecoPT, RecoEta, RecoPhi, MEvent.hiBin, 0);
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

               if(DoZSelection == true){

                  if(bestZidx == -1) continue;

                  double DeltaEtaEle1 = bestEta1Reco - RecoEta;
                  double DeltaPhiEle1 = DeltaPhi(bestPhi1Reco, RecoPhi);
                  double DeltaR1 = sqrt(DeltaEtaEle1*DeltaEtaEle1 + DeltaPhiEle1*DeltaPhiEle1);
   
                  double DeltaEtaEle2 = bestEta2Reco - RecoEta;
                  double DeltaPhiEle2 = DeltaPhi(bestPhi2Reco, RecoPhi);
                  double DeltaR2 = sqrt(DeltaEtaEle2*DeltaEtaEle2 + DeltaPhiEle2*DeltaPhiEle2);
   
                  if(DeltaR1 <  LeptonVeto || DeltaR2 < LeptonVeto)
                     continue;
               }

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
                  TreeTrackResidualCorrection = TrackResidual.GetCorrectionFactor(RecoPT, RecoEta, RecoPhi, MEvent.hiBin, 0);
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

   TH1D *HGenTrack_px = ProjectX(&HGenTrack, "HGenTrack_px"); 
   TH1D *HGenTrack_py = ProjectY(&HGenTrack, "HGenTrack_py"); 
   TH1D *HGenTrack_pz = ProjectZ(&HGenTrack, "HGenTrack_pz"); 


   // Calculate efficiencies and errors
   //TH1D* HEfficiency_eta = new TH1D("HEfficiency_eta", "#eta Efficiency", NEta, EtaBins);
   //TH1D* HEfficiency_phi = new TH1D("HEfficiency_phi", "#phi Efficiency", NPhi, PhiBins);
   //TH1D* HEfficiency_pt  = new TH1D("HEfficiency_pt", "p_{T} Efficiency", NPT, PTBins);
   
   TH1D* HEfficiency_eta = CalculateEfficiencyAndError("HEfficiency_eta", ProjectX(&HRecoTrack,"hr_px"), HGenTrack_px);
   TH1D* HEfficiency_phi = CalculateEfficiencyAndError("HEfficiency_phi", ProjectY(&HRecoTrack,"hr_py"), HGenTrack_py);
   TH1D* HEfficiency_pt  = CalculateEfficiencyAndError("HEfficiency_pt" , ProjectZ(&HRecoTrack,"hr_pz"), HGenTrack_pz);

   //TH1D *HEfficiencyCorrected_eta = ProjectX(&HRecoTrackCorrected, "HEfficiencyCorrected_eta");
   //TH1D *HEfficiencyCorrected_phi = ProjectY(&HRecoTrackCorrected, "HEfficiencyCorrected_phi");
   //TH1D *HEfficiencyCorrected_pt = ProjectZ(&HRecoTrackCorrected, "HEfficiencyCorrected_pt");

   //TH1D *HEfficiencyCorrected_eta = new TH1D("HEfficiencyCorrected_eta", "#eta Corrected Efficiency", NEta, EtaBins);
   //TH1D *HEfficiencyCorrected_phi = new TH1D("HEfficiencyCorrected_phi", "#phi Corrected Efficiency", NPhi, PhiBins);
   //TH1D *HEfficiencyCorrected_pt  = new TH1D("HEfficiencyCorrected_pt", "p_{T} Corrected Efficiency", NPT, PTBins);

   TH1D *HEfficiencyCorrected_eta = CalculateEfficiencyAndError("HEfficiencyCorrected_eta", ProjectX(&HRecoTrackCorrected,"hrc_px"), HGenTrack_px);
   TH1D *HEfficiencyCorrected_phi = CalculateEfficiencyAndError("HEfficiencyCorrected_phi", ProjectY(&HRecoTrackCorrected,"hrc_py"), HGenTrack_py);
   TH1D *HEfficiencyCorrected_pt  = CalculateEfficiencyAndError("HEfficiencyCorrected_pt" , ProjectZ(&HRecoTrackCorrected,"hrc_pz"), HGenTrack_pz);

   //TH1D *HEfficiencyResCorrected_eta = ProjectX(&HRecoTrackResCorrected,"HEfficiencyResCorrected_eta");
   //TH1D *HEfficiencyResCorrected_phi = ProjectY(&HRecoTrackResCorrected,"HEfficiencyResCorrected_phi");
   //TH1D *HEfficiencyResCorrected_pt = ProjectZ(&HRecoTrackResCorrected,"HEfficiencyResCorrected_pt");


   //TH1D *HEfficiencyResCorrected_eta = new TH1D("HEfficiencyResCorrected_eta", "#eta Residual Corrected Efficiency", NEta, EtaBins);
   //TH1D *HEfficiencyResCorrected_phi = new TH1D("HEfficiencyResCorrected_phi", "#phi Residual Corrected Efficiency", NPhi, PhiBins);
   //TH1D *HEfficiencyResCorrected_pt  = new TH1D("HEfficiencyResCorrected_pt", "p_{T} Residual Corrected Efficiency", NPT, PTBins);

   TH1D *HEfficiencyResCorrected_eta = CalculateEfficiencyAndError("HEfficiencyResCorrected_eta", ProjectX(&HRecoTrackResCorrected,"hrr_px"), HGenTrack_px);
   TH1D *HEfficiencyResCorrected_phi = CalculateEfficiencyAndError("HEfficiencyResCorrected_phi", ProjectY(&HRecoTrackResCorrected,"hrr_py"), HGenTrack_py);
   TH1D *HEfficiencyResCorrected_pt  = CalculateEfficiencyAndError("HEfficiencyResCorrected_pt" , ProjectZ(&HRecoTrackResCorrected,"hrr_pz"), HGenTrack_pz);

   cout<<"HEfficiency_eta->GetEntries() = "<<HEfficiency_eta->GetEntries()<<endl;
   cout<<"HEfficiencyCorrected_eta->GetEntries() = "<<HEfficiencyCorrected_eta->GetEntries()<<endl;
   cout<<"HEfficiencyResCorrected_eta->GetEntries() = "<<HEfficiencyResCorrected_eta->GetEntries()<<endl;

   cout<<"HEfficiency_eta->Integral() = "<<HEfficiency_eta->Integral()<<endl;
   cout<<"HEfficiencyCorrected_eta->Integral() = "<<HEfficiencyCorrected_eta->Integral()<<endl;
   cout<<"HEfficiencyResCorrected_eta->Integral() = "<<HEfficiencyResCorrected_eta->Integral()<<endl;

   cout<<"HEfficiency_phi->GetEntries() = "<<HEfficiency_phi->GetEntries()<<endl;
   cout<<"HEfficiencyCorrected_phi->GetEntries() = "<<HEfficiencyCorrected_phi->GetEntries()<<endl;
   cout<<"HEfficiencyResCorrected_phi->GetEntries() = "<<HEfficiencyResCorrected_phi->GetEntries()<<endl;

   cout<<"HEfficiency_phi->Integral() = "<<HEfficiency_phi->Integral()<<endl;
   cout<<"HEfficiencyCorrected_phi->Integral() = "<<HEfficiencyCorrected_phi->Integral()<<endl;
   cout<<"HEfficiencyResCorrected_phi->Integral() = "<<HEfficiencyResCorrected_phi->Integral()<<endl;
   

   cout<<"HEfficiency_pt->GetEntries() = "<<HEfficiency_pt->GetEntries()<<endl;
   cout<<"HEfficiencyCorrected_pt->GetEntries() = "<<HEfficiencyCorrected_pt->GetEntries()<<endl;
   cout<<"HEfficiencyResCorrected_pt->GetEntries() = "<<HEfficiencyResCorrected_pt->GetEntries()<<endl;

   cout<<"HEfficiency_pt->Integral() = "<<HEfficiency_pt->Integral()<<endl;
   cout<<"HEfficiencyCorrected_pt->Integral() = "<<HEfficiencyCorrected_pt->Integral()<<endl;
   cout<<"HEfficiencyResCorrected_pt->Integral() = "<<HEfficiencyResCorrected_pt->Integral()<<endl;




   // Calculate correction factors  

   for (int i = 1; i <= HEfficiencyResCorrected_eta->GetNbinsX(); ++i) {
      bool DoEtaIteration = DoIteration && HEfficiencyResCorrected_eta->GetNbinsX() == hEtaCorrTotal_old->GetNbinsX();
      double eff = DoEtaIteration ? HEfficiencyResCorrected_eta->GetBinContent(i)/hEtaCorrTotal_old->GetBinContent(i) : HEfficiencyCorrected_eta->GetBinContent(i);
      double err = (HGenTrack_px->GetBinContent(i) > 0 && eff>0 && eff<1) ? sqrt(eff * (1-eff) / (HGenTrack_px->GetBinContent(i))) : 0;
      hEtaCorrTotal->SetBinContent(i, eff > 0 ? 1.0 / eff : 1.0);
      hEtaCorrTotal->SetBinError(i, eff > 0 ? err / (eff * eff) : 0.0);
   }

      
   for(int i = 1; i <= HEfficiencyResCorrected_phi->GetNbinsX(); i++){
      bool DoPhiIteration = DoIteration && HEfficiencyResCorrected_phi->GetNbinsX() == hPhiCorrTotal_old->GetNbinsX();
      double eff = DoPhiIteration ? HEfficiencyResCorrected_phi->GetBinContent(i)/hPhiCorrTotal_old->GetBinContent(i) : HEfficiencyCorrected_phi->GetBinContent(i);
      double err = (HGenTrack_py->GetBinContent(i) > 0 && eff>0 && eff<1) ? sqrt(eff * (1-eff) / (HGenTrack_py->GetBinContent(i))) : 0;
      hPhiCorrTotal->SetBinContent(i, eff > 0 ? 1.0 / eff : 1.0);
      hPhiCorrTotal->SetBinError(i, eff > 0 ? err / (eff * eff) : 0.0);
   }


   for(int i = 1; i <= HEfficiencyResCorrected_pt->GetNbinsX(); i++){
      bool DoPtIteration = DoIteration && HEfficiencyResCorrected_pt->GetNbinsX() == hPtCorrTotal_old->GetNbinsX();
      double eff = DoPtIteration ? HEfficiencyResCorrected_pt->GetBinContent(i)/hPtCorrTotal_old->GetBinContent(i) : HEfficiencyCorrected_pt->GetBinContent(i);
      double err = (HGenTrack_pz->GetBinContent(i) > 0 && eff>0 && eff<1) ? sqrt(eff * (1-eff) / (HGenTrack_pz->GetBinContent(i))) : 0;
      hPtCorrTotal->SetBinContent(i, eff > 0 ? 1.0 / eff : 1.0);
      hPtCorrTotal->SetBinError(i, eff > 0 ? err / (eff * eff) : 0.0);
   }

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


   cout<<"==================="<<endl;

   cout<<"HEfficiency_eta->GetEntries() = "<<HEfficiency_eta->GetEntries()<<endl;
   cout<<"HEfficiencyCorrected_eta->GetEntries() = "<<HEfficiencyCorrected_eta->GetEntries()<<endl;
   cout<<"HEfficiencyResCorrected_eta->GetEntries() = "<<HEfficiencyResCorrected_eta->GetEntries()<<endl;

   cout<<"HEfficiency_eta->Integral() = "<<HEfficiency_eta->Integral()<<endl;
   cout<<"HEfficiencyCorrected_eta->Integral() = "<<HEfficiencyCorrected_eta->Integral()<<endl;
   cout<<"HEfficiencyResCorrected_eta->Integral() = "<<HEfficiencyResCorrected_eta->Integral()<<endl;

   cout<<"HEfficiency_phi->GetEntries() = "<<HEfficiency_phi->GetEntries()<<endl;
   cout<<"HEfficiencyCorrected_phi->GetEntries() = "<<HEfficiencyCorrected_phi->GetEntries()<<endl;
   cout<<"HEfficiencyResCorrected_phi->GetEntries() = "<<HEfficiencyResCorrected_phi->GetEntries()<<endl;

   cout<<"HEfficiency_phi->Integral() = "<<HEfficiency_phi->Integral()<<endl;
   cout<<"HEfficiencyCorrected_phi->Integral() = "<<HEfficiencyCorrected_phi->Integral()<<endl;
   cout<<"HEfficiencyResCorrected_phi->Integral() = "<<HEfficiencyResCorrected_phi->Integral()<<endl;
   

   cout<<"HEfficiency_pt->GetEntries() = "<<HEfficiency_pt->GetEntries()<<endl;
   cout<<"HEfficiencyCorrected_pt->GetEntries() = "<<HEfficiencyCorrected_pt->GetEntries()<<endl;
   cout<<"HEfficiencyResCorrected_pt->GetEntries() = "<<HEfficiencyResCorrected_pt->GetEntries()<<endl;

   cout<<"HEfficiency_pt->Integral() = "<<HEfficiency_pt->Integral()<<endl;
   cout<<"HEfficiencyCorrected_pt->Integral() = "<<HEfficiencyCorrected_pt->Integral()<<endl;
   cout<<"HEfficiencyResCorrected_pt->Integral() = "<<HEfficiencyResCorrected_pt->Integral()<<endl;

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

   HEfficiency_eta->Write();
   HEfficiency_phi->Write();
   HEfficiency_pt->Write();

   HEfficiencyCorrected_eta->Write();
   HEfficiencyCorrected_phi->Write();
   HEfficiencyCorrected_pt->Write();

   HEfficiencyResCorrected_eta->Write();
   HEfficiencyResCorrected_phi->Write();
   HEfficiencyResCorrected_pt->Write();


   
   OutputFile->Close();

   delete OutputFile;  // Clean up

   return 0;
}