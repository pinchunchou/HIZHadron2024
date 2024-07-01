#include <cstdlib>
#include <vector>
#include <iostream>
#include <algorithm>
#include <ctime>
using namespace std;

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TMath.h"
#include "TSystem.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TLorentzVector.h"

#include "TRandom.h"

#include "Messenger.h"
#include "CommandLine.h"
#include "CommonFunctions.h"
#include "ProgressBar.h"
#include "CustomAssert.h"

#include "tnp_weight.h"
#include "trackingEfficiency2017pp.h"
#include "trackingEfficiency2018PbPb.h"
#include "TrackResidualCorrector.h"

struct EventIndex;
int main(int argc, char *argv[]);
int FindFirstAbove(vector<EventIndex> &Indices, double X);
double GetHFSum(PFTreeMessenger *M, double MinPFPT);
double GetGenHFSum(GenParticleTreeMessenger *M, double MinGenTrackPT);

struct EventIndex
{
public:
   double HF;
   double VZ;
   int File;
   int Event;
public:
   bool operator <(const EventIndex &other) const
   {
      if(HF < other.HF)         return true;
      if(HF > other.HF)         return false;
      if(VZ < other.VZ)         return true;
      if(VZ > other.VZ)         return false;
      if(File < other.File)     return true;
      if(File > other.File)     return false;
      if(Event < other.Event)   return true;
      if(Event > other.Event)   return false;
      return false;
   }
};


#define E_MASS 0.0005111
#define Z_MASS 91.1880

int main(int argc, char *argv[]){

   string Version = "V1";
   CommandLine CL(argc, argv);

   vector<string> InputFileNames      = CL.GetStringVector("Input");
   string OutputFileName              = CL.Get("Output");

   bool DoMuon                    	  = CL.GetBool("DoMuon", true);
   bool DoElectron                    = CL.GetBool("DoElectron", true);

   bool DoGenLevel                    = CL.GetBool("DoGenLevel", true);
   double Fraction                    = CL.GetDouble("Fraction", 1.00);
   double MinZPT                      = CL.GetDouble("MinZPT", 20.00);
   double MaxZPT                      = CL.GetDouble("MaxZPT", 2000000.0);
   double MinTrackPT                  = CL.GetDouble("MinTrackPT", 1.00);
   double MaxTrackPT                  = CL.GetDouble("MaxTrackPT", 10000.00);
   double MinGenTrackPT               = CL.GetDouble("MinGenTrackPT", 0.40);
   double MinPFPT                     = CL.GetDouble("MinPFPT", 0);

   bool IsData                        = CL.GetBool("IsData", false);
   bool IsPP                          = CL.GetBool("IsPP", false);
   bool DoGenCorrelation              = CL.GetBool("DoGenCorrelation", false);
   bool GenCorrelationCharged         = CL.GetBool("GenCorrelationCharged", false);
   bool DoBackground                  = CL.GetBool("DoBackground", false);
   bool DoSumET                       = CL.GetBool("DoSumET", true);
   double MuonVeto                    = CL.GetDouble("MuonVeto", 0.01);

   bool DoMCHiBinShift                = CL.GetBool("DoMCHiBinShift", true);
   double MCHiBinShift                = DoMCHiBinShift ? CL.GetDouble("MCHiBinShift", 3) : 0;

   bool DoAlternateTrackSelection     = CL.GetBool("DoAlternateTrackSelection", false);
   int AlternateTrackSelection        = DoAlternateTrackSelection ? CL.GetInt("AlternateTrackSelection") : 0;

   bool DoTrackEfficiency             = CL.GetBool("DoTrackEfficiency", true);
   string TrackEfficiencyPath         = (DoTrackEfficiency == true) ? CL.Get("TrackEfficiencyPath") : "";
   bool DoTrackResidual               = CL.GetBool("DoTrackResidual", false);
   vector<string> TrackResidualPath   = (DoTrackResidual == true) ? CL.GetStringVector("TrackResidualPath") : vector<string>{"", "", "", ""};

   string PFTreeName                  = IsPP ? "pfcandAnalyzer/pfTree" : "particleFlowAnalyser/pftree";
   PFTreeName                         = CL.Get("PFTree", PFTreeName);

   string GGTreeName                  = IsPP ? "ggHiNtuplizerGED/EventTree" : "ggHiNtuplizer/EventTree";
   GGTreeName                         = CL.Get("GGTree", GGTreeName);

   bool WithProgressBar               = CL.GetBool("WithProgressBar", false);

   Assert(!(DoGenCorrelation == true && DoGenLevel == false), "You need to turn on gen level to do gen correlation!");
   if(DoTrackResidual == true)
      Assert(TrackResidualPath.size() == 1 || TrackResidualPath.size() == 4, "You need 1 file for residual correction or 4 files for centrality-dependence");

   if(DoBackground == true)
   {
      cerr << "=== WARNING ===" << endl;
      cerr << "Background mixing mode not yet supported with this code." << endl;
      cerr << "=== WARNING ===" << endl;
   }

   vector<string> BackgroundFileNames;
   int NBackground = 0;
   double HFShift = 0;
   double HFTolerance = 0;
   double HFToleranceFraction = 0;
   double HFCeiling = -1;
   double VZTolerance = 0;
   int Oversample = 1;
   bool ReuseBackground = false;
   bool ForceGenMatch = false;

   if(DoBackground == true)
   {
      BackgroundFileNames         = CL.GetStringVector("Background");
      VZTolerance                 = CL.GetDouble("VZTolerance", 2);
      HFShift                     = CL.GetDouble("HFShift");
      HFTolerance                 = CL.GetDouble("Tolerance");
      HFToleranceFraction         = CL.GetDouble("ToleranceFraction");
      HFCeiling                   = CL.GetDouble("HFCeiling", -1);
      NBackground                 = BackgroundFileNames.size();
      Oversample                  = CL.GetInteger("Oversample", 1);
      ReuseBackground             = CL.GetBool("ReuseBackground", false);
      ForceGenMatch               = CL.GetBool("ForceGenMatch", false);

   }

   TrkEff2017pp *TrackEfficiencyPP = nullptr;
   TrkEff2018PbPb *TrackEfficiencyPbPb = nullptr;
   if(DoTrackEfficiency == true)
   {
      if(IsPP == true)
         TrackEfficiencyPP = new TrkEff2017pp(false, TrackEfficiencyPath);
      else
      {
         if(DoAlternateTrackSelection == false)
            TrackEfficiencyPbPb = new TrkEff2018PbPb("general", "", false, TrackEfficiencyPath);
         if(DoAlternateTrackSelection == true && AlternateTrackSelection == 0)
            TrackEfficiencyPbPb = new TrkEff2018PbPb("general", "", false, TrackEfficiencyPath);
         if(DoAlternateTrackSelection == true && AlternateTrackSelection == 1)
            TrackEfficiencyPbPb = new TrkEff2018PbPb("general", "Loose", false, TrackEfficiencyPath);
         if(DoAlternateTrackSelection == true && AlternateTrackSelection == 2)
            TrackEfficiencyPbPb = new TrkEff2018PbPb("general", "Tight", false, TrackEfficiencyPath);
      }
   }

   TrackResidualCentralityCorrector TrackResidual(TrackResidualPath);

   // Do some pre-caching if we read background files.
   // Later on if speed is an issue we can do some optimizations
   vector<TFile *>                    BackgroundFiles;
   vector<HiEventTreeMessenger *>     MBackgroundEvent;
   vector<GenParticleTreeMessenger *> MBackgroundGen;
   vector<TrackTreeMessenger *>       MBackgroundTrackPP;
   vector<PbPbTrackTreeMessenger *>   MBackgroundTrack;
   vector<PFTreeMessenger *>          MBackgroundPF;
   vector<EventIndex>                 BackgroundIndices;

   if(DoBackground == true)
   {
      for(int iB = 0; iB < NBackground; iB++)
      {
         BackgroundFiles.push_back(new TFile(BackgroundFileNames[iB].c_str()));
         MBackgroundEvent.push_back(new HiEventTreeMessenger(BackgroundFiles[iB]));
         MBackgroundGen.push_back(new GenParticleTreeMessenger(BackgroundFiles[iB]));
         MBackgroundTrackPP.push_back(new TrackTreeMessenger(BackgroundFiles[iB]));
         MBackgroundTrack.push_back(new PbPbTrackTreeMessenger(BackgroundFiles[iB]));
         MBackgroundPF.push_back(new PFTreeMessenger(BackgroundFiles[iB], PFTreeName));
         //MBackgroundRho.push_back(new RhoTreeMessenger(BackgroundFiles[iB], RhoTreeName));

         int EntryCount = MBackgroundEvent[iB]->GetEntries();
         for(int iE = 0; iE < EntryCount; iE++)
         {
            MBackgroundEvent[iB]->GetEntry(iE);
            MBackgroundGen[iB]->GetEntry(iE);
            MBackgroundPF[iB]->GetEntry(iE);
            EventIndex E;
            E.HF = DoGenCorrelation ? GetGenHFSum(MBackgroundGen[iB], MinGenTrackPT) : (DoSumET ? MBackgroundEvent[iB]->hiHF : GetHFSum(MBackgroundPF[iB], MinPFPT));
            E.VZ = MBackgroundEvent[iB]->vz;
            E.File = iB;
            E.Event = iE;
            BackgroundIndices.push_back(E);
         }
      }
   }
   sort(BackgroundIndices.begin(), BackgroundIndices.end());

   // Declare output files
   TFile OutputFile(OutputFileName.c_str(), "RECREATE");

   TNtuple   n_Zee( "n_Zee"  , "Invariant Mass of e^{+}e{-}",   "mass:pt:eta:phi");
   TNtuple n_Zmumu( "n_Zmumu", "Invariant Mass of mu^{+}mu{-}", "mass:pt:eta:phi");
   TTree Tree("Tree", Form("Tree for ZHadron analysis, %s", Version.c_str()));
   TTree InfoTree("InfoTree", "Information");

   string Key, Value;
   InfoTree.Branch("Key", &Key);
   InfoTree.Branch("Value", &Value);

   time_t CurrentTime = time(NULL);
   string StringTime = ctime(&CurrentTime);
   replace(StringTime.begin(), StringTime.end(), '\n', ' ');

   Key = "Version";                 Value = InfoString(Version);                 InfoTree.Fill();
   Key = "CurrentTime";             Value = InfoString(StringTime);              InfoTree.Fill();
   Key = "Input";                   Value = InfoString(InputFileNames);          InfoTree.Fill();
   Key = "Output";                  Value = InfoString(OutputFileName);          InfoTree.Fill();
   Key = "DoMuon";                  Value = InfoString(DoMuon);               	InfoTree.Fill();
   Key = "DoElectron";              Value = InfoString(DoElectron);     	      InfoTree.Fill();
   Key = "DoGenLevel";              Value = InfoString(DoGenLevel);              InfoTree.Fill();
   Key = "Fraction";                Value = InfoString(Fraction);                InfoTree.Fill();
   Key = "MinZPT";                  Value = InfoString(MinZPT);                  InfoTree.Fill();
   Key = "MaxZPT";                  Value = InfoString(MaxZPT);                  InfoTree.Fill();
   Key = "MinTrackPT";              Value = InfoString(MinTrackPT);              InfoTree.Fill();
   Key = "MaxTrackPT";              Value = InfoString(MaxTrackPT);              InfoTree.Fill();
   Key = "MinGenTrackPT";           Value = InfoString(MinGenTrackPT);           InfoTree.Fill();
   Key = "MinPFPT";                 Value = InfoString(MinPFPT);                 InfoTree.Fill();
   Key = "IsData";                  Value = InfoString(IsData);                  InfoTree.Fill();
   Key = "IsPP";                    Value = InfoString(IsPP);                    InfoTree.Fill();
   Key = "DoGenCorrelation";        Value = InfoString(DoGenCorrelation);        InfoTree.Fill();
   Key = "GenCorrelationCharged";   Value = InfoString(GenCorrelationCharged);   InfoTree.Fill();
   Key = "DoBackground";            Value = InfoString(DoBackground);            InfoTree.Fill();
   Key = "DoSumET";                 Value = InfoString(DoSumET);                 InfoTree.Fill();
   Key = "MuonVeto";                Value = InfoString(MuonVeto);                InfoTree.Fill();
   Key = "DoTrackEfficiency";       Value = InfoString(DoTrackEfficiency);       InfoTree.Fill();
   Key = "TrackEfficiencyPath";     Value = InfoString(TrackEfficiencyPath);     InfoTree.Fill();
   Key = "DoTrackResidual";         Value = InfoString(DoTrackResidual);         InfoTree.Fill();
   Key = "TrackResidualPath";       Value = InfoString(TrackResidualPath);       InfoTree.Fill();
   Key = "PFTreeName";              Value = InfoString(PFTreeName);              InfoTree.Fill();
   Key = "GGTreeName";              Value = InfoString(GGTreeName);              InfoTree.Fill();
   Key = "Background";              Value = InfoString(BackgroundFileNames);     InfoTree.Fill();
   Key = "VZTolerance";             Value = InfoString(VZTolerance);             InfoTree.Fill();
   Key = "HFShift";                 Value = InfoString(HFShift);                 InfoTree.Fill();
   Key = "Tolerance";               Value = InfoString(HFTolerance);             InfoTree.Fill();
   Key = "ToleranceFraction";       Value = InfoString(HFToleranceFraction);     InfoTree.Fill();
   Key = "HFCeiling";               Value = InfoString(HFCeiling);               InfoTree.Fill();
   Key = "Oversample";              Value = InfoString(Oversample);              InfoTree.Fill();
   Key = "ReuseBackground";         Value = InfoString(ReuseBackground);         InfoTree.Fill();
   Key = "ForceGenMatch";           Value = InfoString(ForceGenMatch);           InfoTree.Fill();
   Key = "MCHiBinShift";            Value = InfoString(MCHiBinShift);            InfoTree.Fill();

   TH2D H2D("H2D", "", 100, -6, 6, 100, -M_PI, M_PI);

   ZHadronMessenger MZHadron;
   MZHadron.SetBranch(&Tree);
   Tree.SetAlias("zP", "(zPt*cosh(zEta))");
   Tree.SetAlias("zPz", "(zPt*sinh(zEta))");
   Tree.SetAlias("zE", "sqrt(zP*zP+zMass*zMass)");
   Tree.SetAlias("zY", "(0.5*log((zE+zPz)/(zE-zPz)))");

   // Loop over signal files
   for(string InputFileName : InputFileNames)
   {
      MZHadron.Clear();

      // Get the input file
      TFile InputFile(InputFileName.c_str());

      // Setup all the messengers.  In the future we'll add more for triggers etc.
      HiEventTreeMessenger     MSignalEvent(InputFile);
      TrackTreeMessenger       MSignalTrackPP(InputFile);
      PbPbTrackTreeMessenger   MSignalTrack(InputFile);
      GenParticleTreeMessenger MSignalGen(InputFile);
      PFTreeMessenger          MSignalPF(InputFile, PFTreeName);
      MuTreeMessenger          MSignalMu(InputFile);
      SkimTreeMessenger        MSignalSkim(InputFile);
      TriggerTreeMessenger     MSignalTrigger(InputFile);
      GGTreeMessenger          MSignalGG(InputFile, GGTreeName);

      // Start looping over events
      int EntryCount = MSignalEvent.GetEntries() * Fraction;
      ProgressBar Bar(cout, EntryCount);
      Bar.SetStyle(-1);
      // Bar.SetStyle(6);

      for(int iE = 0; iE < EntryCount; iE++)
      {
         // Progress bar stuff
         if(WithProgressBar && (EntryCount < 300 || (iE % (EntryCount / 250)) == 0))
         {
            Bar.Update(iE);
            Bar.Print();
         }

         TLorentzVector VGenZ, VGenMu1, VGenMu2, VGenEle1, VGenEle2;

         MSignalEvent.GetEntry(iE);
         MSignalGen.GetEntry(iE);
         if(IsPP == true)
            MSignalTrackPP.GetEntry(iE);
         else
            MSignalTrack.GetEntry(iE);
         MSignalMu.GetEntry(iE);
         MSignalSkim.GetEntry(iE);
         MSignalTrigger.GetEntry(iE);
         MSignalPF.GetEntry(iE);
         MSignalGG.GetEntry(iE);

         if(IsPP == false && IsData == false && DoMCHiBinShift == true)   // PbPb MC, we shift 1.5% as per Kaya
         {
            MSignalEvent.hiBin = MSignalEvent.hiBin - MCHiBinShift;
            if(MSignalEvent.hiBin < 0)   // too central, skip
               continue;
         }

         MZHadron.Run   = MSignalEvent.Run;
         MZHadron.Lumi  = MSignalEvent.Lumi;
         MZHadron.Event = MSignalEvent.Event;
         MZHadron.hiBin = MSignalEvent.hiBin;

         if(IsPP == false && IsData == true)   // need hiBin shifts!
         {
            MZHadron.hiBinUp   = GetHiBin(MSignalEvent.hiHF, 1);
            MZHadron.hiBinDown = GetHiBin(MSignalEvent.hiHF, -1);
         }

         MZHadron.hiHF  = MSignalEvent.hiHF;
         MZHadron.NPU   = 0;
         if(MSignalEvent.npus->size() == 9)
            MZHadron.NPU = MSignalEvent.npus->at(5);
         else if(MSignalEvent.npus->size() > 1)
            MZHadron.NPU = MSignalEvent.npus->at(0);
         else
            MZHadron.NPU = 0;

         // Fill vertex information
         MZHadron.NVertex = 0;
         int BestVertex = -1;
         for(int i = 0; i < (IsPP ? MSignalTrackPP.nVtx : MSignalTrack.VX->size()); i++)
         {
            // TODO: Add vertex selections here

            if(IsPP == true && (BestVertex < 0 || MSignalTrackPP.sumPtVtx[i] > MSignalTrackPP.sumPtVtx[BestVertex]))
               BestVertex = i;
            if(IsPP == false && (BestVertex < 0 || MSignalTrack.VPTSum->at(i) > MSignalTrack.VPTSum->at(BestVertex)))
               BestVertex = i;

            MZHadron.NVertex = MZHadron.NVertex + 1;
         }

         if(BestVertex >= 0)
         {
            MZHadron.VX      = IsPP ? MSignalTrackPP.xVtx[BestVertex] : MSignalTrack.VX->at(BestVertex);
            MZHadron.VY      = IsPP ? MSignalTrackPP.yVtx[BestVertex] : MSignalTrack.VY->at(BestVertex);
            MZHadron.VZ      = IsPP ? MSignalTrackPP.zVtx[BestVertex] : MSignalTrack.VZ->at(BestVertex);
            MZHadron.VXError = IsPP ? MSignalTrackPP.xVtxErr[BestVertex] : MSignalTrack.VXError->at(BestVertex);
            MZHadron.VYError = IsPP ? MSignalTrackPP.yVtxErr[BestVertex] : MSignalTrack.VYError->at(BestVertex);
            MZHadron.VZError = IsPP ? MSignalTrackPP.zVtxErr[BestVertex] : MSignalTrack.VZError->at(BestVertex);

            if(IsData == false)
               MZHadron.VZWeight = IsPP ? GetVZWeightPP(MZHadron.VZ) : GetVZWeightPbPb(MZHadron.VZ);
            else
               MZHadron.VZWeight = 1;
         }

         // Do event selection and triggers
         
         if(IsPP == true)
         {
            if(IsData == true)
            {
               int pprimaryVertexFilter = MSignalSkim.PVFilter;
               int beamScrapingFilter = MSignalSkim.BeamScrapingFilter;
      
               // Event selection criteria
               //    see https://twiki.cern.ch/twiki/bin/viewauth/CMS/HIPhotonJe5TeVpp2017PbPb2018
               if(pprimaryVertexFilter == 0 || beamScrapingFilter == 0)
                  continue;
      
               //HLT trigger to select dimuon events, see Kaya's note: AN2019_143_v12, p.5
               int HLT_HIL2Mu12 = MSignalTrigger.CheckTriggerStartWith("HLT_HIL2Mu12");
               int HLT_HIL3Mu12 = MSignalTrigger.CheckTriggerStartWith("HLT_HIL3Mu12");
               int HLT_HIL3SingleMu12 = MSignalTrigger.CheckTriggerStartWith("HLT_HIL3SingleMu12");// 2023 trigger.

               //e+e-
               int HLT_HIEle20_WPLoose_Gsf = MSignalTrigger.CheckTriggerStartWith("HLT_HIEle20_WPLoose_Gsf");// run2 ppref
               int HLT_HIEle30_WPLoose_Gsf = MSignalTrigger.CheckTriggerStartWith("HLT_HIEle30_WPLoose_Gsf");
               int HLT_HIEle20Gsf = 0; // MSignalTrigger.CheckTriggerStartWith("HLT_HIEle20Gsf"); // Check this after we have run3 ppref data.

               bool muonTrigger = (HLT_HIL3Mu12 != 0 || HLT_HIL2Mu12 != 0 || HLT_HIL3SingleMu12 != 0);
               bool electronTrigger = (HLT_HIEle20_WPLoose_Gsf != 0 || HLT_HIEle30_WPLoose_Gsf != 0 || HLT_HIEle20Gsf != 0 );

               if ((DoMuon && !DoElectron && !muonTrigger) || (!DoMuon && DoElectron && !electronTrigger) || (DoMuon && DoElectron && !muonTrigger && !electronTrigger))
                   continue;
               

               MZHadron.NCollWeight = 1;
            }
            else
               MZHadron.NCollWeight = 1;
         }
         else
         {
            if(IsData == true)
            {
               int pprimaryVertexFilter = MSignalSkim.PVFilter;
               int phfCoincFilter2Th4 = MSignalSkim.HFCoincidenceFilter2Th4;
               int pclusterCompatibilityFilter = MSignalSkim.ClusterCompatibilityFilter;
      
               // Event selection criteria
               //    see https://twiki.cern.ch/twiki/bin/viewauth/CMS/HIPhotonJe5TeVpp2017PbPb2018
               if(pprimaryVertexFilter == 0 || phfCoincFilter2Th4 == 0 || pclusterCompatibilityFilter == 0)
                  continue;
      
               //HLT trigger to select dimuon events, see Kaya's note: AN2019_143_v12, p.5
               int HLT_HIL3Mu12 = MSignalTrigger.CheckTriggerStartWith("HLT_HIL3Mu12");
               int HLT_HIL3SingleMu12 = MSignalTrigger.CheckTriggerStartWith("HLT_HIL3SingleMu12");// 2023 trigger.

               int HLT_HIEle20Gsf = MSignalTrigger.CheckTriggerStartWith("HLT_HIEle20Gsf");

               bool muonTrigger = (HLT_HIL3Mu12 != 0 || HLT_HIL3SingleMu12 != 0);
               bool electronTrigger = (HLT_HIEle20Gsf != 0 );

               if ((DoMuon && !DoElectron && !muonTrigger) || (!DoMuon && DoElectron && !electronTrigger) || (DoMuon && DoElectron && !muonTrigger && !electronTrigger))
                   continue;

               MZHadron.NCollWeight = 1;
            }
            else
               MZHadron.NCollWeight = FindNColl(MSignalEvent.hiBin);
           	   //MZHadron.NCollWeight = 1;
         } // End event selection and triggers

         // Oversample if needed

         

         for(int iS = 0; iS < Oversample; iS++)
         {
            // Loop over gen muons

            if(DoGenLevel == true && DoMuon == true && MSignalMu.NGen > 1)
            {
               // Loop over 1st gen muon
               for(int igen1 = 0; igen1 < MSignalMu.NGen; igen1++)
               {
                  // We only want muon from Z's
                  if(MSignalMu.GenMom[igen1] != 23)
                     continue;
                  if(MSignalMu.GenPT[igen1] < 20)
                     continue;
                  if(fabs(MSignalMu.GenEta[igen1]) > 2.4)
                     continue;
   
                  VGenMu1.SetPtEtaPhiM(MSignalMu.GenPT[igen1],
                        MSignalMu.GenEta[igen1],
                        MSignalMu.GenPhi[igen1],
                        M_MU);

                  // Loop over 2nd gen muon
                  for(int igen2 = igen1 + 1; igen2 < MSignalMu.NGen; igen2++)
                  {
                     // We only want muon from Z's
                     if(MSignalMu.GenMom[igen2] != 23)
                        continue;
                     if(MSignalMu.GenPT[igen2] < 20)
                        continue;
                     if(fabs(MSignalMu.GenEta[igen2]) > 2.4)
                        continue;
   
                     VGenMu2.SetPtEtaPhiM(MSignalMu.GenPT[igen2],
                           MSignalMu.GenEta[igen2],
                           MSignalMu.GenPhi[igen2],
                           M_MU);
   
                     VGenZ = VGenMu1 + VGenMu2;

                     if(VGenZ.M() < 60 || VGenZ.M() > 120)
                        continue;
                     if(fabs(VGenZ.Rapidity()) > 2.4)
                        continue;

                     MZHadron.genZMass->push_back(VGenZ.M());
                     MZHadron.genZPt->push_back  (VGenZ.Pt());
                     MZHadron.genZPhi->push_back (VGenZ.Phi());
                     MZHadron.genZEta->push_back (VGenZ.Eta());
   
                     MZHadron.genMuPt1->push_back(MSignalMu.GenPT[igen1]);
                     MZHadron.genMuPt2->push_back(MSignalMu.GenPT[igen2]);
                     MZHadron.genMuEta1->push_back(MSignalMu.GenEta[igen1]);
                     MZHadron.genMuEta2->push_back(MSignalMu.GenEta[igen2]);
                     MZHadron.genMuPhi1->push_back(MSignalMu.GenPhi[igen1]);
                     MZHadron.genMuPhi2->push_back(MSignalMu.GenPhi[igen2]);
   
                     double genDeltaMuEta = MSignalMu.GenEta[igen1] - MSignalMu.GenEta[igen2];
                     double genDeltaMuPhi = PhiRangePositive(DeltaPhi(MSignalMu.GenPhi[igen1], MSignalMu.GenPhi[igen2]));
   
                     MZHadron.genMuDeta->push_back(genDeltaMuEta);
                     MZHadron.genMuDphi->push_back(genDeltaMuPhi);
                     MZHadron.genMuDR->push_back(sqrt(genDeltaMuEta * genDeltaMuEta + genDeltaMuPhi * genDeltaMuPhi));

                     MZHadron.isGenMuon->push_back(true);
   
                     double genDeltaPhiStar = tan((M_PI - genDeltaMuPhi) / 2)
                        * sqrt(1 - tanh(genDeltaMuEta / 2) * tanh(genDeltaMuEta / 2));
                     MZHadron.genMuDphiS->push_back(genDeltaPhiStar);

                  } // Loop over 2nd gen muon ended

               } // Loop over 1st gen muon ended

            } // Loop over gen muons ended

            // Loop over gen electrons

            if(DoGenLevel == true && DoElectron == true && MSignalGG.NMC > 1)
            {
               // Loop over 1st gen electron
               for(int igen1 = 0; igen1 < MSignalGG.NMC; igen1++)
               {
                  // We only want electron from Z's
               	  if(fabs(MSignalGG.MCPID->at(igen1)) != 11)
                     continue;
                  if(MSignalGG.MCMomPID->at(igen1) != 23)
                     continue;
                  if(MSignalGG.MCPt->at(igen1) < 20)
                     continue;
                  if(fabs(MSignalGG.MCEta->at(igen1)) > 2.1)
                     continue;
   
                  VGenEle1.SetPtEtaPhiM(MSignalGG.MCPt->at(igen1),
                        MSignalGG.MCEta->at(igen1),
                        MSignalGG.MCPhi->at(igen1),
                        E_MASS);

                  // Loop over 2nd gen electron
                  for(int igen2 = igen1 + 1; igen2 < MSignalGG.NMC; igen2++)
                  {
                     // We only want electron from Z's
                     if(MSignalGG.MCPID->at(igen2) != -MSignalGG.MCPID->at(igen1))
                     	continue;
                  	 if(MSignalGG.MCMomPID->at(igen2) != 23)
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

                     if(VGenZ.M() < 60 || VGenZ.M() > 120)
                        continue;
                     if(fabs(VGenZ.Rapidity()) > 2.4)
                        continue;

                     MZHadron.genZMass->push_back(VGenZ.M());
                     MZHadron.genZPt->push_back  (VGenZ.Pt());
                     MZHadron.genZPhi->push_back (VGenZ.Phi());
                     MZHadron.genZEta->push_back (VGenZ.Eta());
   
                     MZHadron.genMuPt1->push_back(MSignalGG.MCPt->at(igen1));
                     MZHadron.genMuPt2->push_back(MSignalGG.MCPt->at(igen2));
                     MZHadron.genMuEta1->push_back(MSignalGG.MCEta->at(igen1));
                     MZHadron.genMuEta2->push_back(MSignalGG.MCEta->at(igen2));
                     MZHadron.genMuPhi1->push_back(MSignalGG.MCPhi->at(igen1));
                     MZHadron.genMuPhi2->push_back(MSignalGG.MCPhi->at(igen2));
   
                     double genDeltaEleEta = MSignalGG.MCEta->at(igen1) - MSignalGG.MCEta->at(igen2);
                     double genDeltaElePhi = PhiRangePositive(DeltaPhi(MSignalGG.MCPhi->at(igen1), MSignalGG.MCPhi->at(igen2)));
   
                     MZHadron.genMuDeta->push_back(genDeltaEleEta);
                     MZHadron.genMuDphi->push_back(genDeltaElePhi);
                     MZHadron.genMuDR->push_back(sqrt(genDeltaEleEta * genDeltaEleEta + genDeltaElePhi * genDeltaElePhi));

                     MZHadron.isGenMuon->push_back(false);
   
                     double genDeltaPhiStar = tan((M_PI - genDeltaElePhi) / 2)
                        * sqrt(1 - tanh(genDeltaEleEta / 2) * tanh(genDeltaEleEta / 2));
                     MZHadron.genMuDphiS->push_back(genDeltaPhiStar);

                  } // Loop over 2nd gen electron ended

               } // Loop over 1st gen electron ended

            } // Loop over gen electrons ended

            // Loop over reco dimuon pairs
            
            int N_pairs;
            
            if(DoMuon == true)
            	N_pairs = MSignalMu.NDi ;
            else
            	N_pairs = 0;

            for(int ipair = 0; ipair < N_pairs; ipair++)
            {
               // We want opposite-charge muons with some basic kinematic cuts
               if(MSignalMu.DiCharge1[ipair] == MSignalMu.DiCharge2[ipair])        continue;
               if(fabs(MSignalMu.DiEta1[ipair]) > 2.4)                             continue;
               if(fabs(MSignalMu.DiEta2[ipair]) > 2.4)                             continue;
               if(fabs(MSignalMu.DiPT1[ipair]) < 20)                               continue;
               if(fabs(MSignalMu.DiPT2[ipair]) < 20)                               continue;
               if(MSignalMu.DimuonPassTightCut(ipair) == false)                    continue;
               if(MSignalMu.DiMass[ipair] < 60 || MSignalMu.DiMass[ipair] > 120)   continue;

               TLorentzVector Mu1, Mu2;
               Mu1.SetPtEtaPhiM(MSignalMu.DiPT1[ipair], MSignalMu.DiEta1[ipair], MSignalMu.DiPhi1[ipair], M_MU);
               Mu2.SetPtEtaPhiM(MSignalMu.DiPT2[ipair], MSignalMu.DiEta2[ipair], MSignalMu.DiPhi2[ipair], M_MU);
               TLorentzVector Z = Mu1 + Mu2;
               if(fabs(Z.Rapidity()) > 2.4)
                  continue;

               MZHadron.zMass->push_back(MSignalMu.DiMass[ipair]);
               MZHadron.zEta->push_back(MSignalMu.DiEta[ipair]);
               MZHadron.zPhi->push_back(MSignalMu.DiPhi[ipair]);
               MZHadron.zPt->push_back(MSignalMu.DiPT[ipair]);
   
               MZHadron.muEta1->push_back(MSignalMu.DiEta1[ipair]);
               MZHadron.muEta2->push_back(MSignalMu.DiEta2[ipair]);
               MZHadron.muPhi1->push_back(MSignalMu.DiPhi1[ipair]);
               MZHadron.muPhi2->push_back(MSignalMu.DiPhi2[ipair]);
   
               MZHadron.muPt1->push_back(MSignalMu.DiPT1[ipair]);
               MZHadron.muPt2->push_back(MSignalMu.DiPT2[ipair]);
   
               double deltaMuEta = MSignalMu.DiEta1[ipair] - MSignalMu.DiEta2[ipair];
               double deltaMuPhi = PhiRangePositive(DeltaPhi(MSignalMu.DiPhi1[ipair], MSignalMu.DiPhi2[ipair]));
   
               MZHadron.muDeta->push_back(deltaMuEta);
               MZHadron.muDphi->push_back(deltaMuPhi);
               MZHadron.muDR->push_back(sqrt(deltaMuEta * deltaMuEta + deltaMuPhi * deltaMuPhi));
   
               double deltaPhiStar = tan((M_PI - deltaMuPhi) / 2) * sqrt(1 - tanh(deltaMuEta / 2) * tanh(deltaMuEta / 2));
   
               MZHadron.muDphiS->push_back(deltaPhiStar);

               MZHadron.isMuon->push_back(true);
   
               n_Zmumu.Fill(MSignalMu.DiMass[ipair], MSignalMu.DiPT[ipair], MSignalMu.DiEta[ipair], MSignalMu.DiPhi[ipair]);

            } // Loop over reco dimuon pairs ended

            // Loop over 1st reco electron
             
            int N_eles;
            
            if(DoElectron == true)
            	N_eles = MSignalGG.NEle;
            else
            	N_eles = 0;

            for(int iele1 = 0; iele1 < N_eles; iele1++)
            {
            	if(DoElectron == false) break;

            	// Some basic electron kinematic cuts
            	if(fabs(MSignalGG.EleEta->at(iele1)) > 2.1)               continue;
            	if(fabs(MSignalGG.ElePt->at(iele1)) < 20)                 continue;
            	if(MSignalGG.DielectronPassVetoCut(iele1) == false)   continue;

            	if(IsPP == false){ // per Kaya, HCAL failure gives rise to misidentified electrons.
            		if(MSignalGG.EleEta->at(iele1) < -1.39 && MSignalGG.ElePhi->at(iele1) > -1.6 &&  MSignalGG.ElePhi->at(iele1) < -0.9 ) continue;
            	}

            	TLorentzVector Ele1;  
            	Ele1.SetPtEtaPhiM(MSignalGG.ElePt->at(iele1), MSignalGG.EleEta->at(iele1), MSignalGG.ElePhi->at(iele1), E_MASS);

            	// Loop over 2nd reco electron
            	for(int iele2 = iele1+1; iele2 < N_eles; iele2++)
            	{
            		// We want opposite-charge electrons with some basic kinematic cuts
            		if(MSignalGG.EleCharge->at(iele1) == MSignalGG.EleCharge->at(iele2))  continue;
            		if(fabs(MSignalGG.EleEta->at(iele2)) > 2.1)               		  continue;
            		if(fabs(MSignalGG.ElePt->at(iele2)) < 20)                 		  continue;
            		if(MSignalGG.DielectronPassVetoCut(iele2) == false)  		  continue;

            		if(IsPP == false){ // per Kaya, HCAL failure gives rise to misidentified electrons.
            			if(MSignalGG.EleEta->at(iele2) < -1.39 && MSignalGG.ElePhi->at(iele2) > -1.6 &&  MSignalGG.ElePhi->at(iele2) < -0.9 ) continue;
            		}

            		TLorentzVector Ele2;  
            		Ele2.SetPtEtaPhiM(MSignalGG.ElePt->at(iele2), MSignalGG.EleEta->at(iele2), MSignalGG.ElePhi->at(iele2), E_MASS);

            		TLorentzVector Z = Ele1+Ele2;
            		double Zmass = Z.M();

            		if(Zmass < 60 || Zmass > 120) continue;
            		if(fabs(Z.Rapidity()) > 2.4) continue;

            		MZHadron.zMass->push_back(Zmass);
               		MZHadron.zEta->push_back(Z.Eta());
               		MZHadron.zPhi->push_back(Z.Phi());
               		MZHadron.zPt->push_back(Z.Pt());
   		
               		MZHadron.muEta1->push_back(MSignalGG.EleEta->at(iele1));
               		MZHadron.muEta2->push_back(MSignalGG.EleEta->at(iele2));
               		MZHadron.muPhi1->push_back(MSignalGG.ElePhi->at(iele1));
               		MZHadron.muPhi2->push_back(MSignalGG.ElePhi->at(iele2));
   		
               		MZHadron.muPt1->push_back(MSignalGG.ElePt->at(iele1));
               		MZHadron.muPt2->push_back(MSignalGG.ElePt->at(iele2));

               		double deltaEleEta = MSignalGG.EleEta->at(iele1) - MSignalGG.EleEta->at(iele2);
               		double deltaElePhi = PhiRangePositive(DeltaPhi(MSignalGG.ElePhi->at(iele1), MSignalGG.ElePhi->at(iele2)));
   		
               		MZHadron.muDeta->push_back(deltaEleEta);
               		MZHadron.muDphi->push_back(deltaElePhi);
               		MZHadron.muDR->push_back(sqrt(deltaEleEta * deltaEleEta + deltaElePhi * deltaElePhi));
   		
               		double deltaPhiStar = tan((M_PI - deltaElePhi) / 2) * sqrt(1 - tanh(deltaEleEta / 2) * tanh(deltaEleEta / 2));
   		
               		MZHadron.muDphiS->push_back(deltaPhiStar);

               		MZHadron.isMuon->push_back(false);
   		
               		n_Zee.Fill(Zmass, Z.Pt(), Z.Eta(), Z.Phi());

            	} // Loop over 2nd reco electron ended

            } // Loop over 1st reco electron ended

            MZHadron.SignalHF = DoGenCorrelation ? GetGenHFSum(&MSignalGen, MinGenTrackPT) : (DoSumET ? MSignalEvent.hiHF : GetHFSum(&MSignalPF, MinPFPT));

            if(DoGenLevel == true){
               MZHadron.SignalGenHF = GetGenHFSum(&MSignalGen, MinGenTrackPT);
            }

            MZHadron.SignalVZ = MSignalEvent.vz;

            // Select the best Z mass. 
            // MZHadron.zMass->at(i) the best one, not 0.
            float Zmass_temp = -1;

            // Loop to select the best Z mass
            for(int idx=0; idx < MZHadron.zPt->size() ; idx++){
            	if((MZHadron.zPt->at(idx) < MinZPT) || (MZHadron.zPt->at(idx) > MaxZPT))
            		continue;

            	if(abs(Zmass_temp - Z_MASS) > abs(MZHadron.zMass->at(idx) - Z_MASS)){
            		Zmass_temp = MZHadron.zMass->at(idx);
            		MZHadron.bestZidx = idx;
            	}
            }

            Zmass_temp = -1;

            for(int idx=0; idx < MZHadron.genZPt->size() ; idx++){
            	if((MZHadron.genZPt->at(idx) < MinZPT) || (MZHadron.genZPt->at(idx) > MaxZPT))
            		continue;

            	if(abs(Zmass_temp - Z_MASS) > abs(MZHadron.genZMass->at(idx) - Z_MASS)){
            		Zmass_temp = MZHadron.genZMass->at(idx);
            		MZHadron.bestZgenIdx = idx;
            	}
            } // or should we pick the reco, or both? or just save the index?

            bool GoodGenZ = MZHadron.bestZgenIdx >= 0 && MZHadron.genZPt->size() > 0 && (MZHadron.genZPt->at(MZHadron.bestZgenIdx) > MinZPT) && (MZHadron.genZPt->at(MZHadron.bestZgenIdx) < MaxZPT);
            bool GoodRecoZ = MZHadron.bestZidx >= 0 && MZHadron.zPt->size() > 0 && (MZHadron.zPt->at(MZHadron.bestZidx) > MinZPT) && (MZHadron.zPt->at(MZHadron.bestZidx) < MaxZPT);

            // Z-track correlation
            // TODO: We can skip bkg first, and then come back later before skimming.
            // Remember to check the latest bkg skim code.

            if((DoGenCorrelation == true && GoodGenZ == true) || (DoGenCorrelation == false && GoodRecoZ == true))
            {
               EventIndex Location;

               PbPbTrackTreeMessenger *MTrack = DoBackground ? MBackgroundTrack[Location.File] : &MSignalTrack;
               TrackTreeMessenger *MTrackPP   = DoBackground ? MBackgroundTrackPP[Location.File] : &MSignalTrackPP;
               GenParticleTreeMessenger *MGen = DoBackground ? MBackgroundGen[Location.File] : &MSignalGen;
               PFTreeMessenger *MPF           = &MSignalPF;

               // Loop over reco tracks and build the correlation function
               int NTrack = DoGenCorrelation ? MGen->Mult : (IsPP ? MTrackPP->nTrk : MTrack->TrackPT->size());
               for(int itrack = 0; itrack < NTrack; itrack++)
               {
                  if(DoGenCorrelation == false)   // track selection on reco
                  {
                     if(IsPP == true)
                     {
                        if(DoAlternateTrackSelection == false && MTrackPP->PassZHadron2022Cut(itrack) == false)
                           continue;
                        if(DoAlternateTrackSelection == true && AlternateTrackSelection == 0 && MTrackPP->PassZHadron2022Cut(itrack) == false)
                           continue;
                        if(DoAlternateTrackSelection == true && AlternateTrackSelection == 1 && MTrackPP->PassZHadron2022CutLoose(itrack) == false)
                           continue;
                        if(DoAlternateTrackSelection == true && AlternateTrackSelection == 2 && MTrackPP->PassZHadron2022CutTight(itrack) == false)
                           continue;
                     }
                     if(IsPP == false)
                     {
                        if(DoAlternateTrackSelection == false && MTrack->PassZHadron2022Cut(itrack) == false)
                           continue;
                        if(DoAlternateTrackSelection == true && AlternateTrackSelection == 0 && MTrack->PassZHadron2022Cut(itrack) == false)
                           continue;
                        if(DoAlternateTrackSelection == true && AlternateTrackSelection == 1 && MTrack->PassZHadron2022CutLoose(itrack) == false)
                           continue;
                        if(DoAlternateTrackSelection == true && AlternateTrackSelection == 2 && MTrack->PassZHadron2022CutTight(itrack) == false)
                           continue;
                     }

                     if((IsPP ? MTrackPP->trkPt[itrack] : MTrack->TrackPT->at(itrack)) < MinTrackPT)
                        continue;

                     if((IsPP ? MTrackPP->trkPt[itrack] : MTrack->TrackPT->at(itrack)) > MaxTrackPT)
                        continue;
                  }// track selection on reco ended

                  if(DoGenCorrelation == true)
                  {
                     if(MGen->PT->at(itrack) < MinTrackPT || MGen->PT->at(itrack) > MaxTrackPT)
                        continue;
                     if(MGen->Eta->at(itrack) < -2.4)
                        continue;
                     if(MGen->Eta->at(itrack) > +2.4)
                        continue;
                     if(MGen->DaughterCount->at(itrack) > 0)
                        continue;
                     if(GenCorrelationCharged == true && MGen->Charge->at(itrack) == 0)
                        continue;
                  }

                  double TrackEta = DoGenCorrelation ? MGen->Eta->at(itrack) : (IsPP ? MTrackPP->trkEta[itrack] : MTrack->TrackEta->at(itrack));
                  double TrackPhi = DoGenCorrelation ? MGen->Phi->at(itrack) : (IsPP ? MTrackPP->trkPhi[itrack] : MTrack->TrackPhi->at(itrack));
                  double TrackPT  = DoGenCorrelation ? MGen->PT->at(itrack) : (IsPP ? MTrackPP->trkPt[itrack] : MTrack->TrackPT->at(itrack));
                  int TrackCharge = DoGenCorrelation ? MGen->Charge->at(itrack) : (IsPP ? MTrackPP->trkCharge[itrack] : MTrack->TrackCharge->at(itrack));
                  int SubEvent    = DoGenCorrelation ? (MGen->SubEvent->at(itrack) + DoBackground) : (IsPP ? 0 : DoBackground);

                  double Mu1Eta = DoGenCorrelation ? MZHadron.genMuEta1->at(MZHadron.bestZgenIdx) : MZHadron.muEta1->at(MZHadron.bestZidx);
                  double Mu1Phi = DoGenCorrelation ? MZHadron.genMuPhi1->at(MZHadron.bestZgenIdx) : MZHadron.muPhi1->at(MZHadron.bestZidx);
                  double Mu2Eta = DoGenCorrelation ? MZHadron.genMuEta2->at(MZHadron.bestZgenIdx) : MZHadron.muEta2->at(MZHadron.bestZidx);
                  double Mu2Phi = DoGenCorrelation ? MZHadron.genMuPhi2->at(MZHadron.bestZgenIdx) : MZHadron.muPhi2->at(MZHadron.bestZidx);

                  double DeltaEtaMu1 = TrackEta - Mu1Eta;
                  double DeltaEtaMu2 = TrackEta - Mu2Eta;
                  double DeltaPhiMu1 = DeltaPhi(TrackPhi, Mu1Phi);
                  double DeltaPhiMu2 = DeltaPhi(TrackPhi, Mu2Phi);

                  double DeltaRMu1 = sqrt(DeltaEtaMu1 * DeltaEtaMu1 + DeltaPhiMu1 * DeltaPhiMu1);
                  double DeltaRMu2 = sqrt(DeltaEtaMu2 * DeltaEtaMu2 + DeltaPhiMu2 * DeltaPhiMu2);

                  bool MuTagged = false;
                  if(DeltaRMu1 < MuonVeto)   MuTagged = true;
                  if(DeltaRMu2 < MuonVeto)   MuTagged = true;

                  double ZEta = DoGenCorrelation ? MZHadron.genZEta->at(MZHadron.bestZgenIdx) : MZHadron.zEta->at(MZHadron.bestZidx);
                  double ZPhi = DoGenCorrelation ? MZHadron.genZPhi->at(MZHadron.bestZgenIdx) : MZHadron.zPhi->at(MZHadron.bestZidx);

                  double deltaEta = TrackEta - ZEta;
                  double deltaPhi = DeltaPhi(TrackPhi, ZPhi);

                  H2D.Fill(+deltaEta, +deltaPhi, 0.25);
                  H2D.Fill(-deltaEta, +deltaPhi, 0.25);
                  H2D.Fill(-deltaEta, -deltaPhi, 0.25);
                  H2D.Fill(+deltaEta, -deltaPhi, 0.25);

                  MZHadron.trackDphi->push_back(deltaPhi);
                  MZHadron.trackDeta->push_back(deltaEta);
                  MZHadron.trackPt->push_back(TrackPT);
                  MZHadron.trackMuTagged->push_back(MuTagged);
                  MZHadron.trackMuDR->push_back(min(DeltaRMu1, DeltaRMu2));
                  MZHadron.subevent->push_back(SubEvent);

                  MZHadron.trackEta->push_back(TrackEta);
                  MZHadron.trackPhi->push_back(TrackPhi);
                  MZHadron.trackCharge->push_back(TrackCharge);

                  double TrackCorrection = 1;
                  if(DoTrackEfficiency == true && DoGenCorrelation == false)
                  {
                     if(IsPP == true)
                        TrackCorrection = TrackEfficiencyPP->getCorrection(TrackPT, TrackEta);
                     else
                        TrackCorrection = TrackEfficiencyPbPb->getCorrection(TrackPT, TrackEta, MZHadron.hiBin + MCHiBinShift);
                  }
                  double TrackResidualCorrection = 1;
                  if(DoTrackResidual == true && DoGenCorrelation == false)
                  {
                     TrackResidualCorrection = TrackResidual.GetCorrectionFactor(TrackPT, TrackEta, TrackPhi, MZHadron.hiBin + MCHiBinShift );
                  }
                  MZHadron.trackWeight->push_back(TrackCorrection);
                  MZHadron.trackResidualWeight->push_back(TrackResidualCorrection);

               } // Loop over reco tracks ended

            } // Z-track correlation loop ended


            // Loop over gen tracks
            if(DoGenLevel == true && GoodGenZ == true)
            {
               EventIndex Location;
               GenParticleTreeMessenger *MGen = DoBackground ? MBackgroundGen[Location.File] : &MSignalGen;

               //std::cout<<"MGen->Mult = "<<MGen->Mult<<std::endl;
               for(int itrack = 0; itrack < MGen->Mult; itrack++)
               {
                  //std::cout<<"MGen->PT->at(itrack) = "<<MGen->PT->at(itrack)<<std::endl;
                  if(MGen->PT->at(itrack) < MinGenTrackPT )
                     continue;

                  if(MGen->PT->at(itrack) < MinTrackPT || MGen->PT->at(itrack) > MaxTrackPT)
                     continue;
                  if(MGen->Eta->at(itrack) < -2.4)
                     continue;
                  if(MGen->Eta->at(itrack) > +2.4)
                     continue;
                  if(MGen->DaughterCount->at(itrack) > 0)
                     continue;
                  if(MGen->Charge->at(itrack) == 0)
                     continue;
                  
                  //std::cout<<"passed gen selection"<<std::endl;
                  double TrackEta = MGen->Eta->at(itrack) ;
                  double TrackPhi = MGen->Phi->at(itrack) ;
                  double TrackPT  = MGen->PT->at(itrack);
                  int TrackCharge = MGen->Charge->at(itrack) ;
                  int SubEvent    = MGen->SubEvent->at(itrack) + DoBackground;

                  //std::cout<<"MZHadron.genMuEta1->size()= "<<MZHadron.genMuEta1->size()<<std::endl;
                  //std::cout<<"MZHadron.genZPt->size()= "<<MZHadron.genZPt->size()<<std::endl;

                  //int idx = DoGenCorrelation ? MZHadron.bestZidx  : MZHadron.bestZgenIdx;
                  double Mu1Eta = MZHadron.genMuEta1->at(MZHadron.bestZgenIdx);
                  double Mu1Phi = MZHadron.genMuPhi1->at(MZHadron.bestZgenIdx);
                  double Mu2Eta = MZHadron.genMuEta2->at(MZHadron.bestZgenIdx);
                  double Mu2Phi = MZHadron.genMuPhi2->at(MZHadron.bestZgenIdx);

                  double DeltaEtaMu1 = TrackEta - Mu1Eta;
                  double DeltaEtaMu2 = TrackEta - Mu2Eta;
                  double DeltaPhiMu1 = DeltaPhi(TrackPhi, Mu1Phi);
                  double DeltaPhiMu2 = DeltaPhi(TrackPhi, Mu2Phi);

                  double DeltaRMu1 = sqrt(DeltaEtaMu1 * DeltaEtaMu1 + DeltaPhiMu1 * DeltaPhiMu1);
                  double DeltaRMu2 = sqrt(DeltaEtaMu2 * DeltaEtaMu2 + DeltaPhiMu2 * DeltaPhiMu2);

                  bool MuTagged = false;
                  if(DeltaRMu1 < MuonVeto)   MuTagged = true;
                  if(DeltaRMu2 < MuonVeto)   MuTagged = true;

                  //std::cout<<"MZHadron.genZEta->size()= "<<MZHadron.genZEta->size()<<std::endl;

                  double ZEta = MZHadron.genZEta->at(MZHadron.bestZgenIdx);
                  double ZPhi = MZHadron.genZPhi->at(MZHadron.bestZgenIdx);

                  double deltaEta = TrackEta - ZEta;
                  double deltaPhi = DeltaPhi(TrackPhi, ZPhi);

                  //std::cout<<"push back"<<std::endl;
                  MZHadron.GenTrackDphi->push_back(deltaPhi);
                  MZHadron.GenTrackDeta->push_back(deltaEta);
                  MZHadron.GenTrackPt->push_back(TrackPT);
                  MZHadron.GenTrackMuTagged->push_back(MuTagged);
                  MZHadron.GenTrackMuDR->push_back(min(DeltaRMu1, DeltaRMu2));
                  MZHadron.GenSubevent->push_back(SubEvent);

                  MZHadron.GenTrackEta->push_back(TrackEta);
                  MZHadron.GenTrackPhi->push_back(TrackPhi);
                  MZHadron.GenTrackCharge->push_back(TrackCharge);
                  //std::cout<<"pushed back"<<std::endl;
               }
               
            } // Loop over gen tracks ended


            //Do Z weights
            MZHadron.ZWeight = 1;
            if(DoGenCorrelation == false && MZHadron.zPt->size() > 0 && MZHadron.bestZidx >= 0 )
            {
               TLorentzVector Z;
               Z.SetPtEtaPhiM(MZHadron.zPt->at(MZHadron.bestZidx), MZHadron.zEta->at(MZHadron.bestZidx), MZHadron.zPhi->at(MZHadron.bestZidx), MZHadron.zMass->at(MZHadron.bestZidx));
               if(IsPP == false)
               {
                  if(IsData == false)
                     MZHadron.ZWeight = GetZWeightPbPbMC(Z.Pt(), Z.Rapidity(), MZHadron.hiBin);
                  else
                  {
                     MZHadron.ZWeight = GetZWeightPbPbDataTrigger(Z.Pt(), Z.Rapidity(), MZHadron.hiBin);
                     
                     double Mu1Eta = MZHadron.muEta1->at(MZHadron.bestZidx);
                     double Mu1PT = MZHadron.muPt1->at(MZHadron.bestZidx);
                     double Mu2Eta = MZHadron.muEta1->at(MZHadron.bestZidx);
                     double Mu2PT = MZHadron.muPt1->at(MZHadron.bestZidx);
                     double Centrality = MZHadron.hiBin * 0.5;

                     MZHadron.ExtraZWeight[0] =
                        tnp_weight_glbPFtrk_pbpb(Mu1Eta, Centrality, -1)
                        / tnp_weight_glbPFtrk_pbpb(Mu1Eta, Centrality, 0)
                        * tnp_weight_glbPFtrk_pbpb(Mu2Eta, Centrality, -1)
                        / tnp_weight_glbPFtrk_pbpb(Mu2Eta, Centrality, 0);
                     MZHadron.ExtraZWeight[1] =
                        tnp_weight_glbPFtrk_pbpb(Mu1Eta, Centrality, -2)
                        / tnp_weight_glbPFtrk_pbpb(Mu1Eta, Centrality, 0)
                        * tnp_weight_glbPFtrk_pbpb(Mu2Eta, Centrality, -2)
                        / tnp_weight_glbPFtrk_pbpb(Mu2Eta, Centrality, 0);
                     MZHadron.ExtraZWeight[2] =
                        tnp_weight_muid_pbpb(Mu1Eta, -1)
                        / tnp_weight_muid_pbpb(Mu1Eta, 0)
                        * tnp_weight_muid_pbpb(Mu2Eta, -1)
                        / tnp_weight_muid_pbpb(Mu2Eta, 0);
                     MZHadron.ExtraZWeight[3] =
                        tnp_weight_muid_pbpb(Mu1Eta, -2)
                        / tnp_weight_muid_pbpb(Mu1Eta, 0)
                        * tnp_weight_muid_pbpb(Mu2Eta, -2)
                        / tnp_weight_muid_pbpb(Mu2Eta, 0);
                     MZHadron.ExtraZWeight[4] =
                        tnp_weight_trig_double_pbpb(Mu1PT, Mu1Eta, Centrality, Mu2PT, Mu2Eta, Centrality, -1)
                        / tnp_weight_trig_double_pbpb(Mu1PT, Mu1Eta, Centrality, Mu2PT, Mu2Eta, Centrality, 0);
                     MZHadron.ExtraZWeight[5] =
                        tnp_weight_trig_double_pbpb(Mu1PT, Mu1Eta, Centrality, Mu2PT, Mu2Eta, Centrality, -2)
                        / tnp_weight_trig_double_pbpb(Mu1PT, Mu1Eta, Centrality, Mu2PT, Mu2Eta, Centrality, 0);

                     MZHadron.ExtraZWeight[6] =
                        tnp_weight_glbPFtrk_pbpb(Mu1Eta, Centrality, 1)
                        / tnp_weight_glbPFtrk_pbpb(Mu1Eta, Centrality, 0)
                        * tnp_weight_glbPFtrk_pbpb(Mu2Eta, Centrality, 1)
                        / tnp_weight_glbPFtrk_pbpb(Mu2Eta, Centrality, 0);
                     MZHadron.ExtraZWeight[7] =
                        tnp_weight_glbPFtrk_pbpb(Mu1Eta, Centrality, 2)
                        / tnp_weight_glbPFtrk_pbpb(Mu1Eta, Centrality, 0)
                        * tnp_weight_glbPFtrk_pbpb(Mu2Eta, Centrality, 2)
                        / tnp_weight_glbPFtrk_pbpb(Mu2Eta, Centrality, 0);
                     MZHadron.ExtraZWeight[8] =
                        tnp_weight_muid_pbpb(Mu1Eta, 1)
                        / tnp_weight_muid_pbpb(Mu1Eta, 0)
                        * tnp_weight_muid_pbpb(Mu2Eta, 1)
                        / tnp_weight_muid_pbpb(Mu2Eta, 0);
                     MZHadron.ExtraZWeight[9] =
                        tnp_weight_muid_pbpb(Mu1Eta, 2)
                        / tnp_weight_muid_pbpb(Mu1Eta, 0)
                        * tnp_weight_muid_pbpb(Mu2Eta, 2)
                        / tnp_weight_muid_pbpb(Mu2Eta, 0);
                     MZHadron.ExtraZWeight[10] =
                        tnp_weight_trig_double_pbpb(Mu1PT, Mu1Eta, Centrality, Mu2PT, Mu2Eta, Centrality, 1)
                        / tnp_weight_trig_double_pbpb(Mu1PT, Mu1Eta, Centrality, Mu2PT, Mu2Eta, Centrality, 0);
                     MZHadron.ExtraZWeight[11] =
                        tnp_weight_trig_double_pbpb(Mu1PT, Mu1Eta, Centrality, Mu2PT, Mu2Eta, Centrality, 2)
                        / tnp_weight_trig_double_pbpb(Mu1PT, Mu1Eta, Centrality, Mu2PT, Mu2Eta, Centrality, 0);
                  }
               }

               else
               {
                  if(IsData == false)
                     MZHadron.ZWeight = GetZWeightPPMC(Z.Pt(), Z.Rapidity());
                  else
                  {
                     MZHadron.ZWeight = GetZWeightPPDataTrigger(Z.Pt(), Z.Rapidity());
                     // Extra Z weight for systematics

                     double Mu1Eta = MZHadron.muEta1->at(MZHadron.bestZidx);
                     double Mu1PT = MZHadron.muPt1->at(MZHadron.bestZidx);
                     double Mu2Eta = MZHadron.muEta1->at(MZHadron.bestZidx);
                     double Mu2PT = MZHadron.muPt1->at(MZHadron.bestZidx);

                     MZHadron.ExtraZWeight[0] =
                        tnp_weight_TightID_pp(Mu1Eta, 1)
                        / tnp_weight_TightID_pp(Mu1Eta, 0)
                        * tnp_weight_TightID_pp(Mu2Eta, 1)
                        / tnp_weight_TightID_pp(Mu2Eta, 0);
                     MZHadron.ExtraZWeight[1] =
                        tnp_weight_TightID_pp(Mu1Eta, -1)
                        / tnp_weight_TightID_pp(Mu1Eta, 0)
                        * tnp_weight_TightID_pp(Mu2Eta, -1)
                        / tnp_weight_TightID_pp(Mu2Eta, 0);
                     MZHadron.ExtraZWeight[2] =
                        tnp_weight_L3Mu12_double_pp(Mu1Eta, Mu2Eta, 1)
                        / tnp_weight_L3Mu12_double_pp(Mu1Eta, Mu2Eta, 0);
                     MZHadron.ExtraZWeight[3] =
                        tnp_weight_L3Mu12_double_pp(Mu1Eta, Mu2Eta, -1)
                        / tnp_weight_L3Mu12_double_pp(Mu1Eta, Mu2Eta, 0);

                     for(int i = 4; i < 12; i++)
                        MZHadron.ExtraZWeight[i] = 1;
                  }
               }
            } // Z weights ended
            MZHadron.FillEntry();
         } //Oversample loop ended
      } // End looping over events

      if(WithProgressBar){
         Bar.Update(EntryCount);
         Bar.Print();
         Bar.PrintLine();
      }

      InputFile.Close();

   } // Loop over signal files ended


   OutputFile.cd();
   H2D.Write();
   n_Zee.Write();
   n_Zmumu.Write();
   Tree.Write();
   InfoTree.Write();

   OutputFile.Close();

   if(DoBackground == true)
   {
      for(TFile *F : BackgroundFiles)
      {
         if(F == nullptr)
            continue;
         F->Close();
         delete F;
      }

      for(HiEventTreeMessenger *M : MBackgroundEvent)
         if(M != nullptr)
            delete M;

      for(PbPbTrackTreeMessenger *M : MBackgroundTrack)
         if(M != nullptr)
            delete M;

      for(PFTreeMessenger *M : MBackgroundPF)
         if(M != nullptr)
            delete M;
      
      //for(RhoTreeMessenger *M : MBackgroundRho)
      //   if(M != nullptr)
      //      delete M;
   }

   if(DoTrackEfficiency == true)
   {
      if(TrackEfficiencyPP != nullptr)
         delete TrackEfficiencyPP;
      if(TrackEfficiencyPbPb != nullptr)
         delete TrackEfficiencyPbPb;
   }

   return 0;
}


int FindFirstAbove(vector<EventIndex> &Indices, double X)
{
   if(X < Indices[0].HF)
      return 0;

   if(X >= Indices[Indices.size()-1].HF)
      return Indices.size();

   int Low = 0;
   int High = Indices.size();

   while(High - Low > 1)
   {
      int Middle = (High + Low) / 2;
      if(X < Indices[Middle].HF)
         High = Middle;
      else
         Low = Middle;
   }

   return Low;
}

double GetHFSum(PFTreeMessenger *M, double MinPFPT)
{
   if(M == nullptr)
      return -1;
   if(M->Tree == nullptr)
      return -1;

   double Sum = 0;
   for(int iPF = 0; iPF < M->ID->size(); iPF++)
      if(fabs(M->Eta->at(iPF)) > 3 && fabs(M->Eta->at(iPF)) < 5 && M->PT->at(iPF) > MinPFPT )
         Sum = Sum + M->E->at(iPF);

   //cout << Sum << endl;

   return Sum;
}

double GetGenHFSum(GenParticleTreeMessenger *M, double MinGenTrackPT)
{
   if(M == nullptr)
      return -1;
   if(M->Tree == nullptr)
      return -1;

   double Sum = 0;
   for(int iGen = 0; iGen < M->Mult; iGen++)
   {
      if(fabs(M->Eta->at(iGen)) < 3)
         continue;
      if(fabs(M->Eta->at(iGen)) > 5)
         continue;
      if(fabs(M->PT->at(iGen)) < MinGenTrackPT)
         continue;
      if(M->DaughterCount->at(iGen) > 0)
         continue;
      Sum = Sum + M->PT->at(iGen) * cosh(M->Eta->at(iGen));
   }

   return Sum;
}