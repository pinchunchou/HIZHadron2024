#include <vector>
#include <algorithm>
#include <iostream>
using namespace std;

#include "TH1D.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TTree.h"

#include "CommandLine.h"
#include "ProgressBar.h"
#include "PlotHelper4.h"
#include "Messenger.h"
#include "CommonFunctions.h"
#include "SetStyle.h"

//#include "tnp_weight.h"

#define MAX 10000
#define E_MASS 0.0005111

int main(int argc, char *argv[]);
int FindBin(double Value, vector<double> &Bins);

int main(int argc, char *argv[])
{
   SetThesisStyle();

   CommandLine CL(argc, argv);

   vector<string> InputFileNames = CL.GetStringVector("Input");
   string OutputFileName         = CL.Get("Output", "EfficiencyPlots.pdf");
   string RootOutputFileName     = CL.Get("RootOutput", "ZEfficiency.root");
   vector<double> Ys             = CL.GetDoubleVector("Y", vector<double>{-2.1, -1.9, -1.7, -1.5, -1.3, -1.1, -0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1});
   vector<double> PTs            = CL.GetDoubleVector("PT", vector<double>{0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 125, 150, 1000});
   double Fraction               = CL.GetDouble("Fraction", 1.00);
   bool IsPP                     = CL.GetBool("IsPP", false);
   string GGTreeName             = IsPP ? "ggHiNtuplizerGED/EventTree" : "ggHiNtuplizer/EventTree";
   GGTreeName                    = CL.Get("GGTree", GGTreeName);
   bool DoMCHiBinShift           = CL.GetBool("DoMCHiBinShift", false);
   double MCHiBinShift           = DoMCHiBinShift ? CL.GetDouble("MCHiBinShift", 3) : 0;


   sort(Ys.begin(), Ys.end());
   sort(PTs.begin(), PTs.end());

   TFile OutputFile(RootOutputFileName.c_str(), "RECREATE");

   PdfFileHelper PdfFile(OutputFileName);
   PdfFile.AddTextPage("Z efficiency derivation");

   TTree Tree("Tree", "Z efficiency tree");
   double TreeZPT, TreeZY, TreeZPhi, TreeZMass;   bool TreeZHasReco;   double TreeHiBin;
   double TreeWPbPbMC, TreeWPbPbData, TreeWPbPbDataTrigger, TreeWPPMC, TreeWPPData, TreeWPPDataTrigger;
   
   Tree.Branch("PT", &TreeZPT, "PT/D");
   Tree.Branch("Y", &TreeZY, "Y/D");
   Tree.Branch("Phi", &TreeZPhi, "Phi/D");
   Tree.Branch("Mass", &TreeZMass, "Mass/D");
   Tree.Branch("HasReco", &TreeZHasReco, "Y/O");
   Tree.Branch("HiBin", &TreeHiBin, "HiBin/D");
   Tree.Branch("WPbPbMC", &TreeWPbPbMC, "TreeWPbPbMC/D");
   Tree.Branch("WPbPbData", &TreeWPbPbData, "TreeWPbPbData/D");
   Tree.Branch("WPbPbDataTrigger", &TreeWPbPbDataTrigger, "TreeWPbPbDataTrigger/D");
   Tree.Branch("WPPMC", &TreeWPPMC, "TreeWPPMC/D");
   Tree.Branch("WPPData", &TreeWPPData, "TreeWPPData/D");
   Tree.Branch("WPPDataTrigger", &TreeWPPDataTrigger, "TreeWPPDataTrigger/D");

   int NY = Ys.size() - 1;
   int NPT = PTs.size() - 1;
   double YBins[MAX], PTBins[MAX];
   for(int i = 0; i <= NY; i++)
      YBins[i] = Ys[i];
   for(int i = 0; i <= NPT; i++)
      PTBins[i] = PTs[i];
   PTBins[0] = PTBins[1] * 0.5;
   PTBins[NPT] = PTBins[NPT-1] * 1.25;

   TH1D HNGen("HNGen", ";Number of Gen Z found;", 25, 0, 25);
   TH1D HNReco("HNReco", ";Number of Reco Z found;", 25, 0, 25);
   TH2D HGenZ("HGenZ", ";y;p_{T}", NY, YBins, NPT, PTBins);
   TH2D HRecoZ("HRecoZ", ";y;p_{T}", NY, YBins, NPT, PTBins);
   TH1D HMatchDR("HMatchDR", ";#DeltaR(y-#phi);", 100, 0, 1.5);
   TH2D HMatchZ("HMatchZ", ";y;p_{T}", NY, YBins, NPT, PTBins);
   TH2D HGenZHasReco("HGenZHasReco", ";y;p_{T}", NY, YBins, NPT, PTBins);
   TH2D HEfficiency("HEfficiency", ";y;p_{T}", NY, YBins, NPT, PTBins);
   TH2D HGenZHasRecoWeighted("HGenZHasRecoWeighted", ";y;p_{T}", NY, YBins, NPT, PTBins);

   HNGen.SetStats(0);
   HNReco.SetStats(0);
   HGenZ.SetStats(0);
   HRecoZ.SetStats(0);
   HMatchDR.SetStats(0);
   HMatchZ.SetStats(0);
   HGenZHasReco.SetStats(0);
   HEfficiency.SetStats(0);
   HGenZHasRecoWeighted.SetStats(0);

   for(string InputFileName : InputFileNames)
   {
      cout << "Processing file " << InputFileName << endl;

      TFile InputFile(InputFileName.c_str());

      HiEventTreeMessenger MEvent(InputFile);
      GGTreeMessenger MSignalGG(InputFile, GGTreeName);

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
         MSignalGG.GetEntry(iE);

         if(IsPP == false)
         {
            MEvent.hiBin = MEvent.hiBin -  MCHiBinShift;   // MC shift
            if(MEvent.hiBin < 0)   // out of range after shifting.  Skip!
               continue;
         }

         if(MSignalGG.NMC == 0)
            continue;

         // First find all the Zs in the event
         vector<TLorentzVector> PGenZ, PGenEle1, PGenEle2;
         for(int iEle1 = 0; iEle1 < MSignalGG.NMC; iEle1++)
         {
            if(iEle1 > MSignalGG.MCPID->size() - 1){
               cerr<<"Warning: MSignalGG.NMC and iEle1 > MSignalGG.MCPID->size(): "<< MSignalGG.NMC <<" or "<<iEle1<<" > "<< MSignalGG.MCPID->size()<<endl;
               break;
            }

            // We only want electron from Z's
            if(fabs(MSignalGG.MCPID->at(iEle1)) != 11)
                continue;
            if(MSignalGG.MCMomPID->at(iEle1) != 23 && ((MSignalGG.MCMomPID->at(iEle1) != MSignalGG.MCPID->at(iEle1)) || MSignalGG.MCGMomPID->at(iEle1) != 23) )
               continue;
            if(MSignalGG.MCPt->at(iEle1) < 20)
               continue;
            if(fabs(MSignalGG.MCEta->at(iEle1)) > 2.1)
               continue;

            TLorentzVector Ele1;
            Ele1.SetPtEtaPhiM(MSignalGG.MCPt->at(iEle1), MSignalGG.MCEta->at(iEle1), MSignalGG.MCPhi->at(iEle1), E_MASS);

            for(int iEle2 = iEle1 + 1; iEle2 < MSignalGG.NMC; iEle2++)
            {
               if(iEle2 > MSignalGG.MCPID->size() - 1){
                  cerr<<"Warning: MSignalGG.NMC and iEle2 > MSignalGG.MCPID->size(): "<< MSignalGG.NMC <<" or "<<iEle2<<" > "<< MSignalGG.MCPID->size()<<endl;
                  break;
               }
               // We only want electron from Z's
               if(MSignalGG.MCPID->at(iEle2) != -MSignalGG.MCPID->at(iEle1))
                   continue;
               if(MSignalGG.MCMomPID->at(iEle2) != 23 && ((MSignalGG.MCMomPID->at(iEle2) != MSignalGG.MCPID->at(iEle2)) || MSignalGG.MCGMomPID->at(iEle2) != 23) )
                  continue;
               if(MSignalGG.MCPt->at(iEle2) < 20)
                  continue;
               if(fabs(MSignalGG.MCEta->at(iEle2)) > 2.1)
                  continue;

               TLorentzVector Ele2;
               Ele2.SetPtEtaPhiM(MSignalGG.MCPt->at(iEle2), MSignalGG.MCEta->at(iEle2), MSignalGG.MCPhi->at(iEle2), E_MASS);

               TLorentzVector Z = Ele1 + Ele2;

               if(Z.M() < 60 || Z.M() > 120)
                  continue;

               if(fabs(Z.Rapidity()) > 2.1) 
                  continue;

               PGenZ.push_back(Ele1 + Ele2);
               PGenEle1.push_back(Ele1);
               PGenEle2.push_back(Ele2);
            }
         }

         int NGen = PGenZ.size();
         HNGen.Fill(NGen);
         for(int i = 0; i < NGen; i++)
         {
            int iY = FindBin(PGenZ[i].Rapidity(), Ys);
            int iPT = FindBin(PGenZ[i].Pt(), PTs);
            HGenZ.SetBinContent(iY + 1, iPT + 1, HGenZ.GetBinContent(iY + 1, iPT + 1) + 1);
         }

         if(NGen == 0)   // no need to continue without a gen Z
            continue;

			// Then find all the reco Zs			
         vector<TLorentzVector> PRecoZ, PRecoEle1, PRecoEle2;

         int N_eles = MSignalGG.NEle ;

         if( N_eles > MSignalGG.ElePt->size() ){
            cerr<<"Warning: MSignalGG.NEle and N_eles > MSignalGG.ElePt->size(): "<< MSignalGG.NEle <<" or "<<N_eles<<" > "<< MSignalGG.ElePt->size()<<endl;
            N_eles = MSignalGG.ElePt->size();
         }

         for(int iele1 = 0; iele1 < N_eles; iele1++)
         {
            // Some basic electron kinematic cuts
            if(fabs(MSignalGG.EleSCEta->at(iele1)) > 2.5)                       continue;
            if(fabs(MSignalGG.EleEta->at(iele1)) > 2.1)                         continue;
            if(fabs(MSignalGG.ElePt->at(iele1)) < 20)                           continue;
            if(IsPP == false && MSignalGG.DielectronPassVetoCut(iele1, MEvent.hiBin + MCHiBinShift) == false) continue;
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
               if(IsPP == false && MSignalGG.DielectronPassVetoCut(iele2, MEvent.hiBin + MCHiBinShift) == false)   continue;
               if(IsPP == true  && MSignalGG.DielectronPassVetoCutPP(iele2) == false)   continue;

               if(IsPP == false){ // per Kaya, HCAL failure gives rise to misidentified electrons.
                  if(MSignalGG.EleSCEta->at(iele2) < -1.39 && MSignalGG.EleSCPhi->at(iele2) > -1.6 &&  MSignalGG.EleSCPhi->at(iele2) < -0.9 ) continue;
               }

               TLorentzVector Ele2;  
               Ele2.SetPtEtaPhiM(MSignalGG.ElePt->at(iele2), MSignalGG.EleEta->at(iele2), MSignalGG.ElePhi->at(iele2), E_MASS);

               TLorentzVector Z = Ele1+Ele2;
               double Zmass = Z.M();

               if(Zmass < 60 || Zmass > 120) continue;
               

               PRecoZ.push_back(Z);
               PRecoEle1.push_back(Ele1);
               PRecoEle2.push_back(Ele2);
            }
         }

         
         int NReco = PRecoZ.size();
         HNReco.Fill(NReco);
         for(int i = 0; i < NReco; i++)
         {
            int iY = FindBin(PRecoZ[i].Rapidity(), Ys);
            int iPT = FindBin(PRecoZ[i].Pt(), PTs);
            HRecoZ.SetBinContent(iY + 1, iPT + 1, HRecoZ.GetBinContent(iY + 1, iPT + 1) + 1);
         }
         
         if(NReco > 0)
         {
            for(int i = 0; i < NGen; i++)
            {
               int iY = FindBin(PGenZ[i].Rapidity(), Ys);
               int iPT = FindBin(PGenZ[i].Pt(), PTs);
               HGenZHasReco.SetBinContent(iY + 1, iPT + 1, HGenZHasReco.GetBinContent(iY + 1, iPT + 1) + 1);

               // double W = GetZWeight(PGenZ[i].Pt(), PGenZ[i].Rapidity(), MEvent.hiBin);
            }
         }
 
         // Now we match
         for(int i = 0; i < NGen; i++)
         {
            double BestDR = -1;

            for(int j = 0; j < NReco; j++)
            {
               double DY = PGenZ[i].Rapidity() - PRecoZ[j].Rapidity();
               double DPhi = DeltaPhi(PGenZ[i].Phi(), PRecoZ[j].Phi());
               double DR = sqrt(DY * DY + DPhi * DPhi);

               if(BestDR < 0 || DR < BestDR)
                  BestDR = DR;
            }

            HMatchDR.Fill(BestDR);
            if(BestDR < 5)   // Matched!
            {
               int iY = FindBin(PGenZ[i].Rapidity(), Ys);
               int iPT = FindBin(PGenZ[i].Pt(), PTs);
               HMatchZ.SetBinContent(iY + 1, iPT + 1, HMatchZ.GetBinContent(iY + 1, iPT + 1) + 1);
            }
         }

         // Finally fill trees
         for(int i = 0; i < NGen; i++)
         {
            TreeZPT = PGenZ[i].Pt();
            TreeZY = PGenZ[i].Rapidity();
            TreeZPhi = PGenZ[i].Phi();
            TreeZMass = PGenZ[i].M();
            TreeZHasReco = (NReco > 0);
            TreeHiBin = MEvent.hiBin;
            TreeWPbPbMC = GetZeeWeightPbPbMC(TreeZPT, TreeZY, TreeHiBin);
            TreeWPbPbData = GetZeeWeightPbPbData(TreeZPT, TreeZY, TreeHiBin);
            TreeWPbPbDataTrigger = GetZeeWeightPbPbDataTrigger(TreeZPT, TreeZY, TreeHiBin);
            TreeWPPMC = GetZeeWeightPPMC(TreeZPT, TreeZY);
            TreeWPPData = GetZeeWeightPPData(TreeZPT, TreeZY);
            TreeWPPDataTrigger = GetZeeWeightPPDataTrigger(TreeZPT, TreeZY);
            
            Tree.Fill();
         }
      }
      Bar.Update(EntryCount);
      Bar.Print();
      Bar.PrintLine();

      InputFile.Close();
   }

   for(int iY = 1; iY <= NY; iY++)
   {
      for(int iPT = 1; iPT <= NPT; iPT++)
      {
         double All = HGenZ.GetBinContent(iY, iPT);
         double Pass = HGenZHasReco.GetBinContent(iY, iPT);
         if(All == 0)
            continue;
         HEfficiency.SetBinContent(iY, iPT, Pass / All);
      }
   }

   PdfFile.AddPlot(HNGen, "hist text00", true);
   PdfFile.AddPlot(HNReco, "hist text00", true);
   PdfFile.AddTextPage("Gen vs reco vs matched");
   PdfFile.AddPlot(HGenZ, "colz text00", true);
   PdfFile.AddPlot(HRecoZ, "colz text00", true);
   PdfFile.AddPlot(HMatchDR, "", true);
   PdfFile.AddPlot(HMatchZ, "colz text00", true);
   PdfFile.AddTextPage("Gen vs has reco");
   PdfFile.AddPlot(HGenZ, "colz text00", true);
   PdfFile.AddPlot(HGenZHasReco, "colz text00", true);
   PdfFile.AddPlot(HEfficiency, "colz, text00", true);

   PdfFile.AddTimeStampPage();
   PdfFile.Close();

   OutputFile.cd();
   HNGen.Write();
   HNReco.Write();
   HGenZ.Write();
   HRecoZ.Write();
   HMatchDR.Write();
   HMatchZ.Write();
   HGenZHasReco.Write();
   HEfficiency.Write();
   Tree.Write();
   OutputFile.Close();

   return 0;
}

int FindBin(double Value, vector<double> &Bins)
{
   for(int i = 0; i < (int)Bins.size() - 1; i++)
      if(Value >= Bins[i] && Value < Bins[i+1])
         return i;
   return -1;
}


