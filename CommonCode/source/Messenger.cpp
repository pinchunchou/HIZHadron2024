#include <algorithm>
#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"

#include "Messenger.h"

namespace EmptyVectors
{
   std::vector<double> EmptyVectorDouble;
   std::vector<float> EmptyVectorFloat;
   std::vector<int> EmptyVectorInt;
   std::vector<std::vector<float> > EmptyVectorVectorFloat;
}

HiEventTreeMessenger::HiEventTreeMessenger(TFile &File)
{
   Tree = (TTree *)File.Get("hiEvtAnalyzer/HiTree");
   Initialize();
}

HiEventTreeMessenger::HiEventTreeMessenger(TFile *File)
{
   if(File != nullptr)
      Tree = (TTree *)File->Get("hiEvtAnalyzer/HiTree");
   else
      Tree = nullptr;
   Initialize();
}

HiEventTreeMessenger::HiEventTreeMessenger(TTree *HiEventTree)
{
   Initialize(HiEventTree);
}
   
bool HiEventTreeMessenger::Initialize(TTree *HiEventTree)
{
   Tree = HiEventTree;
   return Initialize();
}

bool HiEventTreeMessenger::Initialize()
{
   if(Tree == nullptr)
      return false;

   npus = nullptr;
   tnpus = nullptr;

   if(Tree->GetBranch("hiHF"))          Tree->SetBranchAddress("hiHF", &hiHF);
   else                                 hiHF = 0;
   if(Tree->GetBranch("hiBin"))         Tree->SetBranchAddress("hiBin", &hiBin);
   else                                 hiBin = 0;
   if(Tree->GetBranch("run"))           Tree->SetBranchAddress("run", &Run);
   else                                 Run = 1;
   if(Tree->GetBranch("evt"))           Tree->SetBranchAddress("evt", &Event);
   else                                 Event = 1;
   if(Tree->GetBranch("lumi"))          Tree->SetBranchAddress("lumi", &Lumi);
   else                                 Lumi = 1;
   if(Tree->GetBranch("weight"))        Tree->SetBranchAddress("weight", &weight);
   else                                 weight = 1;
   if(Tree->GetBranch("vx"))            Tree->SetBranchAddress("vx", &vx);
   else                                 vx = -999;
   if(Tree->GetBranch("vy"))            Tree->SetBranchAddress("vy", &vy);
   else                                 vy = -999;
   if(Tree->GetBranch("vz"))            Tree->SetBranchAddress("vz", &vz);
   else                                 vz = -999;
   if(Tree->GetBranch("npus"))          Tree->SetBranchAddress("npus", &npus);
   else                                 npus = &EmptyVectors::EmptyVectorInt;
   if(Tree->GetBranch("tnpus"))         Tree->SetBranchAddress("tnpus", &tnpus);
   else                                 tnpus = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("hiHFplus"))      Tree->SetBranchAddress("hiHFplus", &hiHFplus);
   else                                 hiHFplus = 0;
   if(Tree->GetBranch("hiHFminus"))     Tree->SetBranchAddress("hiHFminus", &hiHFminus);
   else                                 hiHFminus = 0;
   if(Tree->GetBranch("hiHFplusEta4"))  Tree->SetBranchAddress("hiHFplusEta4", &hiHFplusEta4);
   else                                 hiHFplusEta4 = 0;
   if(Tree->GetBranch("hiHFminusEta4")) Tree->SetBranchAddress("hiHFminusEta4", &hiHFminusEta4);
   else                                 hiHFminusEta4 = 0;
   if(Tree->GetBranch("hiNevtPlane"))   Tree->SetBranchAddress("hiNevtPlane", &hiNevtPlane);
   else                                 hiNevtPlane = 0;
   if(Tree->GetBranch("hiEvtPlanes"))   Tree->SetBranchAddress("hiEvtPlanes", &hiEvtPlanes);
   
   return true;
}

bool HiEventTreeMessenger::GetEntry(int iEntry)
{
   if(Tree == nullptr)
      return false;

   Tree->GetEntry(iEntry);
   return true;
}
   
int HiEventTreeMessenger::GetEntries()
{
   if(Tree == nullptr)
      return 0;

   return Tree->GetEntries();
}

GGTreeMessenger::GGTreeMessenger(TFile &File, std::string TreeName)
{
   Tree = (TTree *)File.Get(TreeName.c_str());
   Initialize();
}

GGTreeMessenger::GGTreeMessenger(TFile *File, std::string TreeName)
{
   if(File != nullptr)
      Tree = (TTree *)File->Get(TreeName.c_str());
   else
      Tree = nullptr;
   Initialize();
}

GGTreeMessenger::GGTreeMessenger(TTree *EventTree)
{
   Tree = EventTree;
   Initialize();
}

bool GGTreeMessenger::Initialize(TTree *EventTree)
{
   Tree = EventTree;
   return Initialize();
}

bool GGTreeMessenger::Initialize()
{
   if(Tree == nullptr)
      return false;

   PUCount = nullptr;
   PUBX = nullptr;
   PUTrue = nullptr;

   if(Tree->GetBranch("nPUInfo"))   Tree->SetBranchAddress("nPUInfo", &NPUInfo);
   else                             NPUInfo = 0;
   if(Tree->GetBranch("nPU"))       Tree->SetBranchAddress("nPU", &PUCount);
   else                             PUCount = &EmptyVectors::EmptyVectorInt;
   if(Tree->GetBranch("puBX"))      Tree->SetBranchAddress("puBX", &PUBX);
   else                             PUBX = &EmptyVectors::EmptyVectorInt;
   if(Tree->GetBranch("puTrue"))    Tree->SetBranchAddress("puTrue", &PUTrue);
   else                             PUTrue = &EmptyVectors::EmptyVectorFloat;
   
   PFJetPT = nullptr;
   PFJetEta = nullptr;
   PFJetPhi = nullptr;
   CaloJetPT = nullptr;
   CaloJetEta = nullptr;
   CaloJetPhi = nullptr;
   GenJetPT = nullptr;
   GenJetEta = nullptr;
   GenJetPhi = nullptr;

   if(Tree->GetBranch("nPFJt"))         Tree->SetBranchAddress("nPFJt", &PFJetCount);
   else                                 PFJetCount = 0;
   if(Tree->GetBranch("pfJtPt"))        Tree->SetBranchAddress("pfJtPt", &PFJetPT);
   else                                 PFJetPT = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("pfJtEta"))       Tree->SetBranchAddress("pfJtEta", &PFJetEta);
   else                                 PFJetEta = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("pfJtPhi"))       Tree->SetBranchAddress("pfJtPhi", &PFJetPhi);
   else                                 PFJetPhi = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("nCaloJt"))       Tree->SetBranchAddress("nCaloJt", &CaloJetCount);
   else                                 CaloJetCount = 0;
   if(Tree->GetBranch("caloJtPt"))      Tree->SetBranchAddress("caloJtPt", &CaloJetPT);
   else                                 CaloJetPT = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("caloJtEta"))     Tree->SetBranchAddress("caloJtEta", &CaloJetEta);
   else                                 CaloJetEta = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("caloJtPhi"))     Tree->SetBranchAddress("caloJtPhi", &CaloJetPhi);
   else                                 CaloJetPhi = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("nGenJt"))        Tree->SetBranchAddress("nGenJt", &GenJetCount);
   else                                 GenJetCount = 0;
   if(Tree->GetBranch("genJtPt"))       Tree->SetBranchAddress("genJtPt", &GenJetPT);
   else                                 GenJetPT = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("genJtEta"))      Tree->SetBranchAddress("genJtEta", &GenJetEta);
   else                                 GenJetEta = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("genJtPhi"))      Tree->SetBranchAddress("genJtPhi", &GenJetPhi);
   else                                 GenJetPhi = &EmptyVectors::EmptyVectorFloat;

   MCPID = nullptr;
   MCStatus = nullptr;
   MCPt = nullptr;
   MCEta = nullptr;
   MCPhi = nullptr;
   MCE = nullptr;
   MCEt = nullptr;
   MCMass = nullptr;
   MCParentage = nullptr;
   MCMomPID = nullptr;
   MCMomPt = nullptr;
   MCMomEta = nullptr;
   MCMomPhi = nullptr;
   MCMomMass = nullptr;
   MCGMomPID = nullptr;
   MCIndex = nullptr;
   MCSube = nullptr;
   MCCalIsoDR03 = nullptr;
   MCCalIsoDR04 = nullptr;
   MCTrkIsoDR03 = nullptr;
   MCTrkIsoDR04 = nullptr;

   if(Tree->GetBranch("nMC"))           Tree->SetBranchAddress("nMC", &NMC);
   else                                 NMC = 0;
   if(Tree->GetBranch("mcPID"))         Tree->SetBranchAddress("mcPID", &MCPID);
   else                                 MCPID = &EmptyVectors::EmptyVectorInt;
   if(Tree->GetBranch("mcStatus"))      Tree->SetBranchAddress("mcStatus", &MCStatus);
   else                                 MCStatus = &EmptyVectors::EmptyVectorInt;
   if(Tree->GetBranch("mcPt"))          Tree->SetBranchAddress("mcPt", &MCPt);
   else                                 MCPt = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("mcEta"))         Tree->SetBranchAddress("mcEta", &MCEta);
   else                                 MCEta = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("mcPhi"))         Tree->SetBranchAddress("mcPhi", &MCPhi);
   else                                 MCPhi = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("mcE"))           Tree->SetBranchAddress("mcE", &MCE);
   else                                 MCE = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("mcEt"))          Tree->SetBranchAddress("mcEt", &MCEt);
   else                                 MCEt = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("mcMass"))        Tree->SetBranchAddress("mcMass", &MCMass);
   else                                 MCMass = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("mcParentage"))   Tree->SetBranchAddress("mcParentage", &MCParentage);
   else                                 MCParentage = &EmptyVectors::EmptyVectorInt;
   if(Tree->GetBranch("mcMomPID"))      Tree->SetBranchAddress("mcMomPID", &MCMomPID);
   else                                 MCMomPID = &EmptyVectors::EmptyVectorInt;
   if(Tree->GetBranch("mcMomPt"))       Tree->SetBranchAddress("mcMomPt", &MCMomPt);
   else                                 MCMomPt = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("mcMomEta"))      Tree->SetBranchAddress("mcMomEta", &MCMomEta);
   else                                 MCMomEta = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("mcMomPhi"))      Tree->SetBranchAddress("mcMomPhi", &MCMomPhi);
   else                                 MCMomPhi = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("mcMomMass"))     Tree->SetBranchAddress("mcMomMass", &MCMomMass);
   else                                 MCMomMass = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("mcGMomPID"))     Tree->SetBranchAddress("mcGMomPID", &MCGMomPID);
   else                                 MCGMomPID = &EmptyVectors::EmptyVectorInt;
   if(Tree->GetBranch("mcIndex"))       Tree->SetBranchAddress("mcIndex", &MCIndex);
   else                                 MCIndex = &EmptyVectors::EmptyVectorInt;
   if(Tree->GetBranch("mcSube"))        Tree->SetBranchAddress("mcSube", &MCSube);
   else                                 MCSube = &EmptyVectors::EmptyVectorInt;
   if(Tree->GetBranch("mcCalIsoDR03"))  Tree->SetBranchAddress("mcCalIsoDR03", &MCCalIsoDR03);
   else                                 MCCalIsoDR03 = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("mcCalIsoDR04"))  Tree->SetBranchAddress("mcCalIsoDR04", &MCCalIsoDR04);
   else                                 MCCalIsoDR04 = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("mcTrkIsoDR03"))  Tree->SetBranchAddress("mcTrkIsoDR03", &MCTrkIsoDR03);
   else                                 MCTrkIsoDR03 = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("mcTrkIsoDR04"))  Tree->SetBranchAddress("mcTrkIsoDR04", &MCTrkIsoDR04);
   else                                 MCTrkIsoDR04 = &EmptyVectors::EmptyVectorFloat;


   ElePt = nullptr;
   EleEta = nullptr;
   EleHoverE = nullptr;
   ElePhi = nullptr;
   EleCharge = nullptr;
   EleSigmaIEtaIEta_2012 = nullptr;
   EledEtaSeedAtVtx = nullptr;
   EledPhiAtVtx = nullptr;
   EleHoverEBc = nullptr;
   EleEoverPInv = nullptr;
   EleMissHits = nullptr;
   EleIP3D = nullptr;
   EleSCEta = nullptr;
   EleSCPhi = nullptr;
   EleSCRawEn = nullptr;

   if(Tree->GetBranch("nEle"))                  Tree->SetBranchAddress("nEle", &NEle);
   else                                         NEle = 0;
   if(Tree->GetBranch("elePt"))                 Tree->SetBranchAddress("elePt", &ElePt);
   else                                         ElePt = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("eleEta"))                Tree->SetBranchAddress("eleEta", &EleEta);
   else                                         EleEta = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("eleHoverE"))             Tree->SetBranchAddress("eleHoverE", &EleHoverE);
   else                                         EleHoverE = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("elePhi"))                Tree->SetBranchAddress("elePhi", &ElePhi);
   else                                         ElePhi = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("eleCharge"))             Tree->SetBranchAddress("eleCharge", &EleCharge);
   else                                         EleCharge = &EmptyVectors::EmptyVectorInt;
   if(Tree->GetBranch("eleSigmaIEtaIEta_2012")) Tree->SetBranchAddress("eleSigmaIEtaIEta_2012", &EleSigmaIEtaIEta_2012);
   else                                         EleSigmaIEtaIEta_2012 = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("eledEtaSeedAtVtx"))      Tree->SetBranchAddress("eledEtaSeedAtVtx", &EledEtaSeedAtVtx);
   else                                         EledEtaSeedAtVtx = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("eledPhiAtVtx"))          Tree->SetBranchAddress("eledPhiAtVtx", &EledPhiAtVtx);
   else                                         EledPhiAtVtx = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("eleHoverEBc"))           Tree->SetBranchAddress("eleHoverEBc", &EleHoverEBc);
   else                                         EleHoverEBc = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("eleEoverPInv"))          Tree->SetBranchAddress("eleEoverPInv", &EleEoverPInv);
   else                                         EleEoverPInv = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("eleMissHits"))           Tree->SetBranchAddress("eleMissHits", &EleMissHits);
   else                                         EleMissHits = &EmptyVectors::EmptyVectorInt;
   if(Tree->GetBranch("eleIP3D"))               Tree->SetBranchAddress("eleIP3D", &EleIP3D);
   else                                         EleIP3D = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("eleSCEta"))              Tree->SetBranchAddress("eleSCEta", &EleSCEta);
   else                                         EleSCEta = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("eleSCPhi"))              Tree->SetBranchAddress("eleSCPhi", &EleSCPhi);
   else                                         EleSCPhi = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("eleSCRawEn"))            Tree->SetBranchAddress("eleSCRawEn", &EleSCRawEn);
   else                                         EleSCRawEn = &EmptyVectors::EmptyVectorFloat;

   int NEle;

   return true;
}

bool GGTreeMessenger::GetEntry(int iEntry)
{
   if(Tree == nullptr)
      return false;

   Tree->GetEntry(iEntry);
   return true;
}

bool GGTreeMessenger::DielectronPassVetoCut(int index, int hiBin)
{
   // If somebody requested a pair that does not exist, we return false
   if(index >= NEle)
      return false;

   if(EleMissHits->at(index) > 3 || EleIP3D->at(index) >= 0.03)
      return false;

   bool VetoCut = true;

   if(fabs(EleEta->at(index)) < 1.44 && hiBin <= 60 ){
      if(EleSigmaIEtaIEta_2012->at(index) >= 0.0147)
         VetoCut = false;
      if(EledEtaSeedAtVtx->at(index) >= 0.0041)
         VetoCut = false;
      if(EledPhiAtVtx->at(index) >= 0.0853)
         VetoCut = false;
      if(EleHoverEBc->at(index) >= 0.2733)
         VetoCut = false;
      if(EleEoverPInv->at(index) >= 0.0367)
         VetoCut = false;

   }else if(fabs(EleEta->at(index)) > 1.57 && fabs(EleEta->at(index)) < 2.1 && hiBin <= 60 ){
      if(EleSigmaIEtaIEta_2012->at(index) >= 0.048)
         VetoCut = false;
      if(EledEtaSeedAtVtx->at(index) >= 0.0097)
         VetoCut = false;
      if(EledPhiAtVtx->at(index) >= 0.2348)
         VetoCut = false;
      if(EleHoverEBc->at(index) >= 0.1898)
         VetoCut = false;
      if(EleEoverPInv->at(index) >= 0.0300)
         VetoCut = false;

   }else if(fabs(EleEta->at(index)) < 1.44 && hiBin > 60 ){
      if(EleSigmaIEtaIEta_2012->at(index) >= 0.0113)
         VetoCut = false;
      if(EledEtaSeedAtVtx->at(index) >= 0.0037)
         VetoCut = false;
      if(EledPhiAtVtx->at(index) >= 0.1280)
         VetoCut = false;
      if(EleHoverEBc->at(index) >= 0.1814)
         VetoCut = false;
      if(EleEoverPInv->at(index) >= 0.1065)
         VetoCut = false;

   }else if(fabs(EleEta->at(index)) > 1.57 && fabs(EleEta->at(index)) < 2.1 && hiBin > 60 ){
      if(EleSigmaIEtaIEta_2012->at(index) >= 0.0376)
         VetoCut = false;
      if(EledEtaSeedAtVtx->at(index) >= 0.0074)
         VetoCut = false;
      if(EledPhiAtVtx->at(index) >= 0.2085)
         VetoCut = false;
      if(EleHoverEBc->at(index) >= 0.1138)
         VetoCut = false;
      if(EleEoverPInv->at(index) >= 0.0237)
         VetoCut = false;

   }else{
      VetoCut = false;
   }

   // Note that the PF electron is not checked right now
   
   return VetoCut;
}

bool GGTreeMessenger::DielectronPassVetoCutPP(int index)
{
   // If somebody requested a pair that does not exist, we return false
   if(index >= NEle)
      return false;

   if(EleMissHits->at(index) > 3 || EleIP3D->at(index) >= 0.03)
      return false;

   bool VetoCut = true;

   if(fabs(EleEta->at(index)) < 1.44 ){
      if(EleSigmaIEtaIEta_2012->at(index) >= 0.0113)
         VetoCut = false;
      if(EledEtaSeedAtVtx->at(index) >= 0.0037)
         VetoCut = false;
      if(EledPhiAtVtx->at(index) >= 0.1280)
         VetoCut = false;
      if(EleHoverEBc->at(index) >= 0.1814)
         VetoCut = false;
      if(EleEoverPInv->at(index) >= 0.1065)
         VetoCut = false;

   }else if(fabs(EleEta->at(index)) > 1.57 && fabs(EleEta->at(index)) < 2.1 ){
      if(EleSigmaIEtaIEta_2012->at(index) >= 0.0376)
         VetoCut = false;
      if(EledEtaSeedAtVtx->at(index) >= 0.0074)
         VetoCut = false;
      if(EledPhiAtVtx->at(index) >= 0.2085)
         VetoCut = false;
      if(EleHoverEBc->at(index) >= 0.1138)
         VetoCut = false;
      if(EleEoverPInv->at(index) >= 0.0237)
         VetoCut = false;

   }else{
      VetoCut = false;
   }

   // Note that the PF electron is not checked right now
   
   return VetoCut;
}

RhoTreeMessenger::RhoTreeMessenger(TFile &File, std::string TreeName)
{
   Tree = (TTree *)File.Get(TreeName.c_str());
   Initialize();
}

RhoTreeMessenger::RhoTreeMessenger(TFile *File, std::string TreeName)
{
   if(File != nullptr)
      Tree = (TTree *)File->Get(TreeName.c_str());
   else
      Tree = nullptr;
   Initialize();
}

RhoTreeMessenger::RhoTreeMessenger(TTree *RhoTree)
{
   Initialize(RhoTree);
}
   
bool RhoTreeMessenger::Initialize(TTree *RhoTree)
{
   Tree = RhoTree;
   return Initialize();
}

bool RhoTreeMessenger::Initialize()
{
   if(Tree == nullptr)
      return false;

   EtaMin = nullptr;
   EtaMax = nullptr;
   Rho = nullptr;
   RhoM = nullptr;

   if(Tree->GetBranch("etaMin")) Tree->SetBranchAddress("etaMin", &EtaMin);
   else                          EtaMin = &EmptyVectors::EmptyVectorDouble;
   if(Tree->GetBranch("etaMax")) Tree->SetBranchAddress("etaMax", &EtaMax);
   else                          EtaMax = &EmptyVectors::EmptyVectorDouble;
   if(Tree->GetBranch("rho"))    Tree->SetBranchAddress("rho", &Rho);
   else                          Rho = &EmptyVectors::EmptyVectorDouble;
   if(Tree->GetBranch("rhom"))   Tree->SetBranchAddress("rhom", &RhoM);
   else                          RhoM = &EmptyVectors::EmptyVectorDouble;

   return true;
}

bool RhoTreeMessenger::GetEntry(int iEntry)
{
   if(Tree == nullptr)
      return false;

   Tree->GetEntry(iEntry);
   return true;
}

SkimTreeMessenger::SkimTreeMessenger(TFile &File)
{
   Tree = (TTree *)File.Get("skimanalysis/HltTree");
   Initialize();
}

SkimTreeMessenger::SkimTreeMessenger(TFile *File)
{
   if(File != nullptr)
      Tree = (TTree *)File->Get("skimanalysis/HltTree");
   else
      Tree = nullptr;
   Initialize();
}

SkimTreeMessenger::SkimTreeMessenger(TTree *SkimTree)
{
   Initialize(SkimTree);
}

bool SkimTreeMessenger::Initialize(TTree *SkimTree)
{
   Tree = SkimTree;
   return Initialize();
}

bool SkimTreeMessenger::Initialize()
{
   if(Tree == nullptr)
      return false;

   if(Tree->GetBranch("HBHENoiseFilterResult"))
      Tree->SetBranchAddress("HBHENoiseFilterResult", &HBHENoise);
   else
      HBHENoise = 1;
   if(Tree->GetBranch("HBHENoiseFilterResultRun2Loose"))
      Tree->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseRun2Loose);
   else
      HBHENoiseRun2Loose = 1;
   if(Tree->GetBranch("pprimaryVertexFilter"))
      Tree->SetBranchAddress("pprimaryVertexFilter", &PVFilter);
   else if(Tree->GetBranch("pPAprimaryVertexFilter"))
      Tree->SetBranchAddress("pPAprimaryVertexFilter", &PVFilter);
   else
      PVFilter = 1;
   if(Tree->GetBranch("pBeamScrapingFilter"))
      Tree->SetBranchAddress("pBeamScrapingFilter", &BeamScrapingFilter);
   else
      BeamScrapingFilter = 1;
   if(Tree->GetBranch("pclusterCompatibilityFilter"))
      Tree->SetBranchAddress("pclusterCompatibilityFilter", &ClusterCompatibilityFilter);
   else
      ClusterCompatibilityFilter = 1;
   if(Tree->GetBranch("phfCoincFilter3"))
      Tree->SetBranchAddress("phfCoincFilter3", &HFCoincidenceFilter);
   else
      HFCoincidenceFilter = 1;
   if(Tree->GetBranch("phfCoincFilter2Th4"))
      Tree->SetBranchAddress("phfCoincFilter2Th4", &HFCoincidenceFilter2Th4);
   else
      HFCoincidenceFilter2Th4 = 1;
   if(Tree->GetBranch("pcollisionEventSelection"))
      Tree->SetBranchAddress("pcollisionEventSelection", &CollisionEventSelection);
   else
      CollisionEventSelection = 1;

   return true;
}

bool SkimTreeMessenger::PassBasicFilter()
{
   if(HBHENoise == 0)
      return false;
   if(HBHENoiseRun2Loose == 0)
      return false;
   if(PVFilter == 0)
      return false;
   if(ClusterCompatibilityFilter == 0)
      return false;
   if(HFCoincidenceFilter2Th4 == 0)
      return false;
   if(CollisionEventSelection == 0)
      return false;

   return true;
}

bool SkimTreeMessenger::PassBasicFilterLoose()
{
   if(HBHENoiseRun2Loose == 0)
      return false;
   if(PVFilter == 0)
      return false;
   if(ClusterCompatibilityFilter == 0)
      return false;
   if(HFCoincidenceFilter2Th4 == 0)
      return false;
   if(CollisionEventSelection == 0)
      return false;

   return true;
}

bool SkimTreeMessenger::GetEntry(int iEntry)
{
   if(Tree == nullptr)
      return false;

   Tree->GetEntry(iEntry);
   return true;
}

JetTreeMessenger::JetTreeMessenger(TFile &File, std::string TreeName)
{
   Tree = (TTree *)File.Get(TreeName.c_str());
   Initialize();
}

JetTreeMessenger::JetTreeMessenger(TFile *File, std::string TreeName)
{
   if(File != nullptr)
      Tree = (TTree *)File->Get(TreeName.c_str());
   else
      Tree = nullptr;
   Initialize();
}

JetTreeMessenger::JetTreeMessenger(TTree *JetTree)
{
   Initialize(JetTree);
}

bool JetTreeMessenger::Initialize(TTree *JetTree)
{
   Tree = JetTree;
   return Initialize();
}

bool JetTreeMessenger::Initialize()
{
   if(Tree == nullptr)
      return false;

   JetSubJetPT = nullptr;
   JetSubJetEta = nullptr;
   JetSubJetPhi = nullptr;
   JetSubJetM = nullptr;
   RefGSubJetPT = nullptr;
   RefGSubJetEta = nullptr;
   RefGSubJetPhi = nullptr;
   RefGSubJetM = nullptr;

   if(Tree->GetBranch("pthat"))   Tree->SetBranchAddress("pthat", &PTHat);

   if(Tree->GetBranch("nref"))    Tree->SetBranchAddress("nref", &JetCount);
   else                           JetCount = 0;
   if(Tree->GetBranch("rawpt"))   Tree->SetBranchAddress("rawpt", &JetRawPT);
   if(Tree->GetBranch("jtpt"))    Tree->SetBranchAddress("jtpt", &JetPT);
   if(Tree->GetBranch("jteta"))   Tree->SetBranchAddress("jteta", &JetEta);
   if(Tree->GetBranch("jty"))     Tree->SetBranchAddress("jty", &JetY);
   if(Tree->GetBranch("jtphi"))   Tree->SetBranchAddress("jtphi", &JetPhi);
   if(Tree->GetBranch("jtpu"))    Tree->SetBranchAddress("jtpu", &JetPU);
   if(Tree->GetBranch("jtm"))     Tree->SetBranchAddress("jtm", &JetM);
   if(Tree->GetBranch("jttau1"))  Tree->SetBranchAddress("jttau1", &JetTau1);
   if(Tree->GetBranch("jttau2"))  Tree->SetBranchAddress("jttau2", &JetTau2);
   if(Tree->GetBranch("jttau3"))  Tree->SetBranchAddress("jttau3", &JetTau3);
   if(Tree->GetBranch("jtarea"))  Tree->SetBranchAddress("jtarea", &JetArea);
   if(Tree->GetBranch("discr_csvV1"))  Tree->SetBranchAddress("discr_csvV1", &JetCSVV1);
   if(Tree->GetBranch("discr_csvV2"))  Tree->SetBranchAddress("discr_csvV2", &JetCSVV2);
   if(Tree->GetBranch("ndiscr_csvV1")) Tree->SetBranchAddress("ndiscr_csvV1", &JetCSVV1N);
   if(Tree->GetBranch("ndiscr_csvV2")) Tree->SetBranchAddress("ndiscr_csvV2", &JetCSVV2N);
   if(Tree->GetBranch("pdiscr_csvV1")) Tree->SetBranchAddress("pdiscr_csvV1", &JetCSVV1P);
   if(Tree->GetBranch("pdiscr_csvV2")) Tree->SetBranchAddress("pdiscr_csvV2", &JetCSVV2P);

   if(Tree->GetBranch("jtSubJetPt"))
      Tree->SetBranchAddress("jtSubJetPt", &JetSubJetPT);
   else
      JetSubJetPT = &EmptyVectors::EmptyVectorVectorFloat;
   if(Tree->GetBranch("jtSubJetEta"))
      Tree->SetBranchAddress("jtSubJetEta", &JetSubJetEta);
   else
      JetSubJetEta = &EmptyVectors::EmptyVectorVectorFloat;
   if(Tree->GetBranch("jtSubJetPhi"))
      Tree->SetBranchAddress("jtSubJetPhi", &JetSubJetPhi);
   else
      JetSubJetPhi = &EmptyVectors::EmptyVectorVectorFloat;
   if(Tree->GetBranch("jtSubJetM"))
      Tree->SetBranchAddress("jtSubJetM", &JetSubJetM);
   else
      JetSubJetM = &EmptyVectors::EmptyVectorVectorFloat;
   if(Tree->GetBranch("jtsym"))
      Tree->SetBranchAddress("jtsym", &JetSym);
   else
      for(int i = 0; i < JETCOUNTMAX; i++)
         JetSym[i] = -1;
   if(Tree->GetBranch("jtdroppedBranches"))
      Tree->SetBranchAddress("jtdroppedBranches", &JetDroppedBranches);
   else
      for(int i = 0; i < JETCOUNTMAX; i++)
         JetDroppedBranches[i] = -1;
   
   if(Tree->GetBranch("refpt"))   Tree->SetBranchAddress("refpt", &RefPT);
   if(Tree->GetBranch("refeta"))  Tree->SetBranchAddress("refeta", &RefEta);
   if(Tree->GetBranch("refy"))    Tree->SetBranchAddress("refy", &RefY);
   if(Tree->GetBranch("refphi"))  Tree->SetBranchAddress("refphi", &RefPhi);
   if(Tree->GetBranch("refm"))    Tree->SetBranchAddress("refm", &RefM);
   if(Tree->GetBranch("refarea")) Tree->SetBranchAddress("refarea", &RefArea);
   
   if(Tree->GetBranch("refptG"))  Tree->SetBranchAddress("refptG", &RefGPT);
   if(Tree->GetBranch("refetaG")) Tree->SetBranchAddress("refetaG", &RefGEta);
   if(Tree->GetBranch("refphiG")) Tree->SetBranchAddress("refphiG", &RefGPhi);
   if(Tree->GetBranch("refmG"))   Tree->SetBranchAddress("refmG", &RefGM);

   if(Tree->GetBranch("refparton_pt"))
      Tree->SetBranchAddress("refparton_pt", &RefPartonPT);
   if(Tree->GetBranch("refparton_flavor"))
      Tree->SetBranchAddress("refparton_flavor", &RefPartonFlavor);
   if(Tree->GetBranch("refparton_flavorForB"))
      Tree->SetBranchAddress("refparton_flavorForB", &RefPartonFlavorForB);
   
   if(Tree->GetBranch("refSubJetPt"))
      Tree->SetBranchAddress("refSubJetPt", &RefGSubJetPT);
   else
      RefGSubJetPT = &EmptyVectors::EmptyVectorVectorFloat;
   if(Tree->GetBranch("refSubJetEta"))
      Tree->SetBranchAddress("refSubJetEta", &RefGSubJetEta);
   else
      RefGSubJetEta = &EmptyVectors::EmptyVectorVectorFloat;
   if(Tree->GetBranch("refSubJetPhi"))
      Tree->SetBranchAddress("refSubJetPhi", &RefGSubJetPhi);
   else
      RefGSubJetPhi = &EmptyVectors::EmptyVectorVectorFloat;
   if(Tree->GetBranch("refSubJetM"))
      Tree->SetBranchAddress("refSubJetM", &RefGSubJetM);
   else
      RefGSubJetM = &EmptyVectors::EmptyVectorVectorFloat;
   
   if(Tree->GetBranch("ngen"))    Tree->SetBranchAddress("ngen", &GenCount);
   else                           GenCount = 0;
   
   if(Tree->GetBranch("genmatchindex"))
      Tree->SetBranchAddress("genmatchindex", &GenMatchIndex);
   else
      for(int i = 0; i < JETCOUNTMAX; i++)
         GenMatchIndex[i] = -1;
   if(Tree->GetBranch("genpt"))   Tree->SetBranchAddress("genpt", &GenPT);
   if(Tree->GetBranch("geneta"))  Tree->SetBranchAddress("geneta", &GenEta);
   if(Tree->GetBranch("geny"))    Tree->SetBranchAddress("geny", &GenY);
   if(Tree->GetBranch("genphi"))  Tree->SetBranchAddress("genphi", &GenPhi);
   if(Tree->GetBranch("genm"))    Tree->SetBranchAddress("genm", &GenM);

   if(Tree->GetBranch("hcalSum")) Tree->SetBranchAddress("hcalSum", &HcalSum);
   if(Tree->GetBranch("ecalSum")) Tree->SetBranchAddress("ecalSum", &EcalSum);

   if(Tree->GetBranch("jtPfCHF")) Tree->SetBranchAddress("jtPfCHF", &JetPFCHF);
   if(Tree->GetBranch("jtPfNHF")) Tree->SetBranchAddress("jtPfNHF", &JetPFNHF);
   if(Tree->GetBranch("jtPfCEF")) Tree->SetBranchAddress("jtPfCEF", &JetPFCEF);
   if(Tree->GetBranch("jtPfNEF")) Tree->SetBranchAddress("jtPfNEF", &JetPFNEF);
   if(Tree->GetBranch("jtPfMUF")) Tree->SetBranchAddress("jtPfMUF", &JetPFMUF);
   if(Tree->GetBranch("jtPfCHM")) Tree->SetBranchAddress("jtPfCHM", &JetPFCHM);
   if(Tree->GetBranch("jtPfNHM")) Tree->SetBranchAddress("jtPfNHM", &JetPFNHM);
   if(Tree->GetBranch("jtPfCEM")) Tree->SetBranchAddress("jtPfCEM", &JetPFCEM);
   if(Tree->GetBranch("jtPfNEM")) Tree->SetBranchAddress("jtPfNEM", &JetPFNEM);
   if(Tree->GetBranch("jtPfMUM")) Tree->SetBranchAddress("jtPfMUM", &JetPFMUM);

   return true;
}

bool JetTreeMessenger::GetEntry(int iEntry)
{
   if(Tree == nullptr)
      return false;

   Tree->GetEntry(iEntry);
   return true;
}

GenParticleTreeMessenger::GenParticleTreeMessenger(TFile &File)
{
   Tree = (TTree *)File.Get("HiGenParticleAna/hi");
   Initialize();
}

GenParticleTreeMessenger::GenParticleTreeMessenger(TFile *File)
{
   if(File != nullptr)
      Tree = (TTree *)File->Get("HiGenParticleAna/hi");
   else
      Tree = nullptr;
   Initialize();
}

GenParticleTreeMessenger::GenParticleTreeMessenger(TTree *GenParticleTree)
{
   Initialize(GenParticleTree);
}

bool GenParticleTreeMessenger::Initialize(TTree *GenParticleTree)
{
   Tree = GenParticleTree;
   return Initialize();
}

bool GenParticleTreeMessenger::Initialize()
{
   if(Tree == nullptr)
      return false;

   PT = nullptr;
   Eta = nullptr;
   Phi = nullptr;
   ID = nullptr;
   Charge = nullptr;
   DaughterCount = nullptr;
   SubEvent = nullptr;

   if(Tree->GetBranch("npart")) Tree->SetBranchAddress("npart", &NPart);
   else                         NPart = 0;
   if(Tree->GetBranch("ncoll")) Tree->SetBranchAddress("ncoll", &NColl);
   else                         NColl = 0;
   if(Tree->GetBranch("nhard")) Tree->SetBranchAddress("nhard", &NHard);
   else                         NHard = 0;
   if(Tree->GetBranch("phi0"))  Tree->SetBranchAddress("phi0", &Phi0);
   else                         Phi0 = 0;

   if(Tree->GetBranch("mult"))  Tree->SetBranchAddress("mult", &Mult);
   else                         Mult = 0;
   if(Tree->GetBranch("pt"))    Tree->SetBranchAddress("pt", &PT);
   else                         PT = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("eta"))   Tree->SetBranchAddress("eta", &Eta);
   else                         Eta = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("phi"))   Tree->SetBranchAddress("phi", &Phi);
   else                         Phi = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("pdg"))   Tree->SetBranchAddress("pdg", &ID);
   else                         ID = &EmptyVectors::EmptyVectorInt;
   if(Tree->GetBranch("chg")) Tree->SetBranchAddress("chg", &Charge);
   else                         Charge = &EmptyVectors::EmptyVectorInt;
   if(Tree->GetBranch("nDaughters")) Tree->SetBranchAddress("nDaughters", &DaughterCount);
   else                         DaughterCount = &EmptyVectors::EmptyVectorInt;
   if(Tree->GetBranch("sube"))  Tree->SetBranchAddress("sube", &SubEvent);
   else                         SubEvent = &EmptyVectors::EmptyVectorInt;

   return true;
}

bool GenParticleTreeMessenger::GetEntry(int iEntry)
{
   if(Tree == nullptr)
      return false;

   Tree->GetEntry(iEntry);
   return true;
}

PFTreeMessenger::PFTreeMessenger(TFile &File, std::string TreeName)
{
   Initialize(File, TreeName);
}

PFTreeMessenger::PFTreeMessenger(TFile *File, std::string TreeName)
{
   Initialize(File, TreeName);
}

PFTreeMessenger::PFTreeMessenger(TTree *PFTree)
{
   Initialize(PFTree);
}

bool PFTreeMessenger::Initialize(TFile &File, std::string TreeName)
{
   Tree = (TTree *)File.Get(TreeName.c_str());
   return Initialize();
}

bool PFTreeMessenger::Initialize(TFile *File, std::string TreeName)
{
   if(File != nullptr)
      Tree = (TTree *)File->Get(TreeName.c_str());
   else
      Tree = nullptr;
   return Initialize();
}

bool PFTreeMessenger::Initialize(TTree *PFTree)
{
   Tree = PFTree;
   return Initialize();
}

bool PFTreeMessenger::Initialize()
{
   if(Tree == nullptr)
      return false;

   ID = nullptr;
   PT = nullptr;
   E = nullptr;
   Eta = nullptr;
   Phi = nullptr;
   M = nullptr;

   if(Tree->GetBranch("pfId"))     Tree->SetBranchAddress("pfId", &ID);
   else                            ID = &EmptyVectors::EmptyVectorInt;
   if(Tree->GetBranch("pfPt"))     Tree->SetBranchAddress("pfPt", &PT);
   else                            PT = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("pfEnergy")) Tree->SetBranchAddress("pfEnergy", &E);
   else if(Tree->GetBranch("pfE")) Tree->SetBranchAddress("pfE", &E);
   else                            E = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("pfEta"))    Tree->SetBranchAddress("pfEta", &Eta);
   else                            Eta = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("pfPhi"))    Tree->SetBranchAddress("pfPhi", &Phi);
   else                            Phi = &EmptyVectors::EmptyVectorFloat;
   if(Tree->GetBranch("pfM"))      Tree->SetBranchAddress("pfM", &M);
   else                            M = &EmptyVectors::EmptyVectorFloat;

   return true;
}

bool PFTreeMessenger::GetEntry(int iEntry)
{
   if(Tree == nullptr)
      return false;

   Tree->GetEntry(iEntry);
   return true;
}

TriggerTreeMessenger::TriggerTreeMessenger(TFile &File, std::string TreeName)
{
   Tree = (TTree *)File.Get(TreeName.c_str());
   Initialize();
}

TriggerTreeMessenger::TriggerTreeMessenger(TFile *File, std::string TreeName)
{
   if(File != nullptr)
      Tree = (TTree *)File->Get(TreeName.c_str());
   else
      Tree = nullptr;
   Initialize();
}

TriggerTreeMessenger::TriggerTreeMessenger(TTree *TriggerTree)
{
   Initialize(TriggerTree);
}

bool TriggerTreeMessenger::Initialize(TTree *TriggerTree)
{
   Tree = TriggerTree;
   return Initialize();
}

bool TriggerTreeMessenger::Initialize()
{
   if(Tree == nullptr)
      return false;

   if(Tree->GetBranch("Run"))   Tree->SetBranchAddress("Run", &Run);
   else                         Run = 0;
   if(Tree->GetBranch("LumiBlock")) Tree->SetBranchAddress("LumiBlock", &Lumi);
   else                         Lumi = 0;
   if(Tree->GetBranch("Event")) Tree->SetBranchAddress("Event", &Event);
   else                         Event = 0;

   FillTriggerNames();


   for(int i = 0; i < (int)Name.size(); i++)
   {
      if(Tree->GetBranch(Name[i].c_str()))
      {
         Tree->SetBranchAddress(Name[i].c_str(), &Decision[i]);
         Exist[i] = true;
         
         if(Tree->GetBranch((Name[i] + "_Prescl").c_str()))
         {
            Tree->SetBranchAddress((Name[i] + "_Prescl").c_str(), &Prescale[i]);
            PrescaleExist[i] = true;
         }
         else if(( Tree->GetBranch((Name[i] + "_PrescaleNumerator").c_str()) ) && ( Tree->GetBranch((Name[i] + "_PrescaleDenominator").c_str()) ) )
         {
            Tree->SetBranchAddress((Name[i] + "_PrescaleNumerator").c_str(), &PrescaleNumerator[i]);
            Tree->SetBranchAddress((Name[i] + "_PrescaleDenominator").c_str(), &PrescaleDenominator[i]);
            // TODO: We should change it to float at some time? But it seems that we are not using this info yet.

            if(PrescaleDenominator[i] > 0 && PrescaleNumerator[i] > 0)
               Prescale[i] = (int) (PrescaleNumerator[i]/PrescaleDenominator[i]); 
            else
               Prescale[i] = 0;
            PrescaleExist[i] = true;
         }
         else
         {
            Prescale[i] = -1;
            PrescaleExist[i] = false;
         }
      }
      else
      {
         Decision[i] = false;
         Exist[i] = false;
      }
   }

   return true;
}


bool TriggerTreeMessenger::GetEntry(int iEntry)
{
   if(Tree == nullptr)
      return false;

   Tree->GetEntry(iEntry);
   return true;
}
   
void TriggerTreeMessenger::FillTriggerNames()
{
   Name.clear();
   Decision.clear();
   
   // 2018 triggers
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part0_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part1_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part2_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part3_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part4_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part5_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part6_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part7_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part8_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part9_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part10_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part11_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part12_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part13_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part14_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part15_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part16_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part17_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part18_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixBypass_part19_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part0_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part1_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part2_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part3_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part4_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part5_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part6_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part7_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part8_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part9_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part10_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part11_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part12_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part13_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part14_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part15_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part16_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part17_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part18_v1");
   Name.push_back("HLT_HIMinimumBias_SinglePixelTrack_NpixGated_part19_v1");

   // 2018 AA trigger
   Name.push_back("HLT_HIPuAK4CaloJet100Eta5p1_v1");
   Name.push_back("HLT_HIPuAK4CaloJet80Eta5p1_v1");
   Name.push_back("HLT_HIPuAK4CaloJet60Eta5p1_v1");
   Name.push_back("HLT_HIPuAK4CaloJet60Eta5p1_v1");
   Name.push_back("HLT_HIPuAK4CaloJet40Eta5p1_v1");

   // AA trigger
   Name.push_back("HLT_HIPuAK4CaloJet100_Eta5p1_v1");
   Name.push_back("HLT_HIPuAK4CaloJet80_Eta5p1_v1");
   Name.push_back("HLT_HIPuAK4CaloJet60_Eta5p1_v1");
   Name.push_back("HLT_HIPuAK4CaloJet60_Eta5p1_v1");
   Name.push_back("HLT_HIPuAK4CaloJet40_Eta5p1_v1");

   Name.push_back("HLT_PAAK4PFJet40_Eta5p1v2");
   Name.push_back("HLT_PAAK4PFJet60_Eta5p1_v2");
   Name.push_back("HLT_PAAK4PFJet80_Eta5p1_v2");

   // pp trigger
   Name.push_back("HLT_AK4PFJet80_Eta5p1_v1");
   Name.push_back("HLT_HIAK4PFJet40_v1");
   Name.push_back("HLT_HIAK4PFJet60_v1");
   Name.push_back("HLT_HIAK4PFJet80_v1");
   Name.push_back("HLT_HIAK4PFJet100_v1");
   Name.push_back("HLT_HIAK4PFJet120_v1");

   // L1 pass through
   Name.push_back("HLT_L1SingleJet8_v1_BptxAND_v1");
   Name.push_back("HLT_L1SingleJet12_v1_BptxAND_v1");
   Name.push_back("HLT_L1SingleJet16_v1");
   Name.push_back("HLT_L1SingleJet20_v1");
   Name.push_back("HLT_L1SingleJet35_v1");
   Name.push_back("HLT_L1SingleJet40_v1");
   Name.push_back("HLT_L1SingleJet44_v1");
   Name.push_back("HLT_L1SingleJet60_v1");

   // pA temporary triggers - v16
   // Name.push_back("HLT_PAAK4CaloBJetCSVv160_Eta2p1_v1");   // same in v18
   // Name.push_back("HLT_PAAK4CaloBJetCSVv260_Eta2p1_v1");   // same in v18
   Name.push_back("HLT_PAAK4PFJetCSV60_Eta2p1_v2");

   // pA temporary triggers - v22
   Name.push_back("HLT_PAAK4CaloJet60_Eta2p1_v1");
   Name.push_back("HLT_PAAK4CaloJet80_Eta2p1_v1");
   Name.push_back("HLT_PAAK4CaloJet100_Eta2p1_v1");
   Name.push_back("HLT_PAAK4CaloBJetCSVv160_Eta2p1_v1");
   Name.push_back("HLT_PAAK4CaloBJetCSVv180_Eta2p1_v1");
   Name.push_back("HLT_PAAK4CaloBJetCSVv1100_Eta2p1_v1");
   Name.push_back("HLT_PAAK4CaloBJetCSVv260_Eta2p1_v1");
   Name.push_back("HLT_PAAK4CaloBJetCSVv280_Eta2p1_v1");
   Name.push_back("HLT_PAAK4CaloBJetCSVv2100_Eta2p1_v1");
   Name.push_back("HLT_PAAK4PFJet60_Eta2p1_v1");
   Name.push_back("HLT_PAAK4PFJet80_Eta2p1_v1");
   Name.push_back("HLT_PAAK4PFJet100_Eta2p1_v1");
   Name.push_back("HLT_PAAK4PFBJetCSV60_Eta2p1_v1");
   Name.push_back("HLT_PAAK4PFBJetCSV80_Eta2p1_v1");
   Name.push_back("HLT_PAAK4PFBJetCSV100_Eta2p1_v1");

   // pA temporary triggers - v24
   Name.push_back("HLT_L1ZeroBias_BptxAND_v1");
   
   // pA temporary triggers - v27
   Name.push_back("HLT_PAAK4CaloJet40_Eta2p1_v1");
   Name.push_back("HLT_PAAK4CaloBJetCSV40_Eta2p1_v1");
   Name.push_back("HLT_PAAK4CaloBJetCSV60_Eta2p1_v1");
   Name.push_back("HLT_PAAK4CaloBJetCSV80_Eta2p1_v1");
   Name.push_back("HLT_PAAK4PFJet40_Eta2p1_v1");
   Name.push_back("HLT_PAAK4PFBJetCSV40_Eta2p1_v1");
   
   // pA temporary triggers - v31 / JCv1
   Name.push_back("HLT_PAAK4CaloBJetCSV40_CommonTracking_Eta2p1_v1");
   Name.push_back("HLT_PAAK4CaloBJetCSV60_CommonTracking_Eta2p1_v1");
   Name.push_back("HLT_PAAK4CaloBJetCSV80_CommonTracking_Eta2p1_v1");

   // Chris V14
   Name.push_back("HLT_PAAK4CaloJet40_Eta5p1_v2");
   Name.push_back("HLT_PAAK4CaloJet60_Eta5p1_v2");
   Name.push_back("HLT_PAAK4CaloJet80_Eta5p1_v2");
   Name.push_back("HLT_PAAK4CaloJet40_Eta1p9toEta5p1_v2");
   Name.push_back("HLT_PAAK4CaloJet60_Eta1p9toEta5p1_v2");
   Name.push_back("HLT_PAAK4CaloJet40_Eta2p9toEta5p1_v2");
   Name.push_back("HLT_PAAK4CaloJet40_Eta5p1_PAL3Mu3_v2");
   Name.push_back("HLT_PAAK4CaloJet40_Eta5p1_PAL3Mu5_v2");
   Name.push_back("HLT_PAAK4CaloJet60_Eta5p1_PAL3Mu3_v2");
   Name.push_back("HLT_PAAK4CaloJet60_Eta5p1_PAL3Mu5_v2");
   Name.push_back("HLT_PAAK4PFJet40_Eta5p1_v2");
   Name.push_back("HLT_PAAK4PFJet60_Eta5p1_v2");
   Name.push_back("HLT_PAAK4PFJet80_Eta5p1_v2");
   Name.push_back("HLT_PAAK4PFJet60_Eta1p9toEta5p1_v2");
   Name.push_back("HLT_PAAK4PFJet40_Eta1p9toEta5p1_v2");
   Name.push_back("HLT_PAAK4PFJet40_Eta2p9toEta5p1_v2");
   Name.push_back("HLT_PADiAK4CaloJetAve40_Eta5p1_v2");
   Name.push_back("HLT_PADiAK4CaloJetAve60_Eta5p1_v2");
   Name.push_back("HLT_PADiAK4CaloJetAve80_Eta5p1_v2");
   Name.push_back("HLT_PADiAK4PFJetAve40_Eta5p1_v2");
   Name.push_back("HLT_PADiAK4PFJetAve60_Eta5p1_v2");
   Name.push_back("HLT_PADiAK4PFJetAve80_Eta5p1_v2");
   Name.push_back("HLT_PASinglePhoton10_Eta3p1_v2");
   Name.push_back("HLT_PASinglePhoton15_Eta3p1_v2");
   Name.push_back("HLT_PASinglePhoton20_Eta3p1_v2");
   Name.push_back("HLT_PASinglePhoton30_Eta3p1_v2");
   Name.push_back("HLT_PADoublePhoton15_Eta3p1_Mass50_1000_v2");
   Name.push_back("HLT_PASinglePhoton10_Eta3p1_PAL3Mu3_v2");
   Name.push_back("HLT_PASinglePhoton10_Eta3p1_PAL3Mu5_v2");
   Name.push_back("HLT_PASinglePhoton15_Eta3p1_PAL3Mu3_v2");
   Name.push_back("HLT_PASinglePhoton15_Eta3p1_PAL3Mu5_v2");
   Name.push_back("HLT_PASinglePhoton20_Eta3p1_PAL3Mu3_v2");
   Name.push_back("HLT_PASinglePhoton20_Eta3p1_PAL3Mu5_v2");

   Name.push_back("L1_SingleJet12_BptxAND");
   
   Name.push_back("HLT_PAAK4CaloJet40_Eta5p1_v1");
   Name.push_back("HLT_PAAK4CaloJet60_Eta5p1_v1");
   Name.push_back("HLT_PAAK4CaloJet80_Eta5p1_v1");
   Name.push_back("HLT_PAAK4CaloJet100_Eta5p1_v1");
   Name.push_back("HLT_PAAK4PFJet40_Eta5p1_v1");
   Name.push_back("HLT_PAAK4PFJet60_Eta5p1_v1");
   Name.push_back("HLT_PAAK4PFJet60_Eta5p1_v3");
   Name.push_back("HLT_PAAK4PFJet80_Eta5p1_v1");
   Name.push_back("HLT_PAAK4PFJet100_Eta5p1_v1");

   Name.push_back("L1_MinimumBiasHF0_OR_BptxAND_Initial");
   Name.push_back("L1_MinimumBiasHF0_OR_BptxAND_Final");
   Name.push_back("L1_MinimumBiasHF0_OR_BptxAND");
   Name.push_back("HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_v1");
   Name.push_back("HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_v2");
   Name.push_back("HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_v3");
   Name.push_back("HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_v4");
   Name.push_back("HLT_PAL1MinimumBiasHF_AND_SinglePixelTrack_v5");

   Name.push_back("L1_SingleJet12_BptxAND_Initial");
   Name.push_back("L1_SingleJet12_BptxAND_Final");
   Name.push_back("L1_SingleJet16_BptxAND_Initial");
   Name.push_back("L1_SingleJet16_BptxAND_Final");
   Name.push_back("L1_SingleJet24_BptxAND_Initial");
   Name.push_back("L1_SingleJet24_BptxAND_Final");
   Name.push_back("L1_SingleJet36_BptxAND_Initial");
   Name.push_back("L1_SingleJet36_BptxAND_Final");
   Name.push_back("L1_SingleJet52_BptxAND_Initial");
   Name.push_back("L1_SingleJet52_BptxAND_Final");
   
   Name.push_back("L1_SingleJet12_ForEta1p5_BptxAND_Initial");
   Name.push_back("L1_SingleJet12_ForEta1p5_BptxAND_Final");
   Name.push_back("L1_SingleJet12_ForEta2p5_BptxAND_Initial");
   Name.push_back("L1_SingleJet12_ForEta2p5_BptxAND_Final");
   Name.push_back("L1_SingleJet24_ForEta1p5_BptxAND_Initial");
   Name.push_back("L1_SingleJet24_ForEta1p5_BptxAND_Final");
   Name.push_back("L1_SingleJet24_ForEta2p5_BptxAND_Initial");
   Name.push_back("L1_SingleJet24_ForEta2p5_BptxAND_Final");

   Name.push_back("HLT_PAAK4CaloJet40_Eta5p1_v3");
   Name.push_back("HLT_PAAK4CaloJet60_Eta5p1_v3");
   Name.push_back("HLT_PAAK4CaloJet80_Eta5p1_v3");
   Name.push_back("HLT_PAAK4CaloJet40_Eta1p9toEta5p1_v3");
   Name.push_back("HLT_PAAK4CaloJet60_Eta1p9toEta5p1_v3");
   Name.push_back("HLT_PAAK4CaloJet40_Eta2p9toEta5p1_v3");
   Name.push_back("HLT_PAAK4CaloJet40_Eta5p1_PAL3Mu3_v3");
   Name.push_back("HLT_PAAK4CaloJet40_Eta5p1_PAL3Mu5_v3");
   Name.push_back("HLT_PAAK4CaloJet60_Eta5p1_PAL3Mu3_v3");
   Name.push_back("HLT_PAAK4CaloJet60_Eta5p1_PAL3Mu5_v3");
   Name.push_back("HLT_PAAK4PFJet40_Eta5p1_v3");
   Name.push_back("HLT_PAAK4PFJet60_Eta5p1_v3");
   Name.push_back("HLT_PAAK4PFJet60_Eta5p1_v4");
   Name.push_back("HLT_PAAK4PFJet80_Eta5p1_v3");
   Name.push_back("HLT_PAAK4PFJet60_Eta1p9toEta5p1_v3");
   Name.push_back("HLT_PAAK4PFJet40_Eta1p9toEta5p1_v3");
   Name.push_back("HLT_PAAK4PFJet40_Eta2p9toEta5p1_v3");
   Name.push_back("HLT_PADiAK4CaloJetAve40_Eta5p1_v3");
   Name.push_back("HLT_PADiAK4CaloJetAve60_Eta5p1_v3");
   Name.push_back("HLT_PADiAK4CaloJetAve80_Eta5p1_v3");
   Name.push_back("HLT_PADiAK4PFJetAve40_Eta5p1_v3");
   Name.push_back("HLT_PADiAK4PFJetAve60_Eta5p1_v3");
   Name.push_back("HLT_PADiAK4PFJetAve80_Eta5p1_v3");
   Name.push_back("HLT_PASinglePhoton10_Eta3p1_v3");
   Name.push_back("HLT_PASinglePhoton15_Eta3p1_v3");
   Name.push_back("HLT_PASinglePhoton20_Eta3p1_v3");
   Name.push_back("HLT_PASinglePhoton30_Eta3p1_v3");
   Name.push_back("HLT_PADoublePhoton15_Eta3p1_Mass50_1000_v3");
   Name.push_back("HLT_PASinglePhoton10_Eta3p1_PAL3Mu3_v3");
   Name.push_back("HLT_PASinglePhoton10_Eta3p1_PAL3Mu5_v3");
   Name.push_back("HLT_PASinglePhoton15_Eta3p1_PAL3Mu3_v3");
   Name.push_back("HLT_PASinglePhoton15_Eta3p1_PAL3Mu5_v3");
   Name.push_back("HLT_PASinglePhoton20_Eta3p1_PAL3Mu3_v3");
   Name.push_back("HLT_PASinglePhoton20_Eta3p1_PAL3Mu5_v3");

   // PP Reference Runs - v4
   Name.push_back("HLT_AK4PFJet30_v15");
   Name.push_back("HLT_AK4PFJet50_v15");
   Name.push_back("HLT_AK4PFJet80_v15");
   Name.push_back("HLT_AK4PFJet100_v15");
   Name.push_back("HLT_AK4PFJet120_v14");

   Name.push_back("HLT_AK4CaloJet30_v10");
   Name.push_back("HLT_AK4CaloJet40_v9");
   Name.push_back("HLT_AK4CaloJet50_v9");
   Name.push_back("HLT_AK4CaloJet80_v9");
   Name.push_back("HLT_AK4CaloJet100_v9");
   Name.push_back("HLT_AK4CaloJet120_v8");

   Name.push_back("HLT_L1SingleJet16_Eta5p1_v1");
   Name.push_back("HLT_L1SingleJet20_Eta5p1_v1");
   Name.push_back("HLT_L1SingleJet28_Eta5p1_v1");
   Name.push_back("HLT_L1SingleJet40_Eta5p1_v1");
   Name.push_back("HLT_L1SingleJet60_Eta5p1_v1");

   // pp reference runs - v6
   Name.push_back("HLT_L1SingleJet20FWD");
   Name.push_back("HLT_L1SingleJet35FWD");
   Name.push_back("HLT_L1SingleJet40FWD");
   Name.push_back("HLT_L1SingleJet50FWD");
   Name.push_back("HLT_L1SingleJet60FWD");
   
   Name.push_back("HLT_AK4PFJet30FWD_v14");
   Name.push_back("HLT_AK4PFJet50FWD_v14");
   Name.push_back("HLT_AK4PFJet60FWD_v14");
   Name.push_back("HLT_AK4PFJet70FWD_v14");
   Name.push_back("HLT_AK4PFJet80FWD_v14");

   // muon triggers
   Name.push_back("HLT_HIL2Mu12_v1");
   Name.push_back("HLT_HIL2Mu12_v2");
   Name.push_back("HLT_HIL2Mu12_v3");
   Name.push_back("HLT_HIL2Mu12_v4");
   Name.push_back("HLT_HIL2Mu12_v5");
   Name.push_back("HLT_HIL3Mu12_v1");
   Name.push_back("HLT_HIL3Mu12_v2");
   Name.push_back("HLT_HIL3Mu12_v3");
   Name.push_back("HLT_HIL3Mu15_v1");
   Name.push_back("HLT_HIL3Mu15_v2");
   Name.push_back("HLT_HIL3Mu15_v3");
   Name.push_back("HLT_HIL3Mu20_v1");
   Name.push_back("HLT_HIL3Mu20_v2");
   Name.push_back("HLT_HIL3Mu20_v3");

   // electron triggers
   Name.push_back("HLT_HIEle10Gsf_v1");
   Name.push_back("HLT_HIEle15Gsf_v1");
   Name.push_back("HLT_HIEle20Gsf_v1");
   Name.push_back("HLT_HIEle30Gsf_v1");
   Name.push_back("HLT_HIEle40Gsf_v1");
   Name.push_back("HLT_HIEle50Gsf_v1");
   Name.push_back("HLT_HIEle15Ele10Gsf_v1");
   Name.push_back("HLT_HIEle15Ele10GsfMass50_v1");
   Name.push_back("HLT_HIDoubleEle10Gsf_v1");
   Name.push_back("HLT_HIDoubleEle10GsfMass50_v1");
   Name.push_back("HLT_HIDoubleEle15Gsf_v1");
   Name.push_back("HLT_HIDoubleEle15GsfMass50_v1");

   Name.push_back("HLT_HIEle10_WPLoose_Gsf_v1");
   Name.push_back("HLT_HIEle15_WPLoose_Gsf_v1");
   Name.push_back("HLT_HIEle17_WPLoose_Gsf_v1");
   Name.push_back("HLT_HIEle20_WPLoose_Gsf_v1");
   Name.push_back("HLT_HIEle30_WPLoose_Gsf_v1");
   Name.push_back("HLT_HIEle40_WPLoose_Gsf_v1");

   // 2023 triggers
   Name.push_back("HLT_HIL3SingleMu5_v3");
   Name.push_back("HLT_HIL3SingleMu7_v3");
   Name.push_back("HLT_HIL3SingleMu12_v3");
   Name.push_back("HLT_HIL2SingleMu5_v3");
   Name.push_back("HLT_HIL2SingleMu7_v3");
   Name.push_back("HLT_HIL2DoubleMu0_Open_v3");
   Name.push_back("HLT_HIL3DoubleMu0_M0toInf_Open_v3");
   
   Name.push_back("HLT_HIEle10Gsf_v10");
   Name.push_back("HLT_HIEle15Gsf_v10");
   Name.push_back("HLT_HIEle20Gsf_v10");
   Name.push_back("HLT_HIEle30Gsf_v10");
   Name.push_back("HLT_HIEle40Gsf_v10");
   Name.push_back("HLT_HIEle50Gsf_v10");
   Name.push_back("HLT_HIEle15Ele10Gsf_v10");
   Name.push_back("HLT_HIEle15Ele10GsfMass50_v10");
   Name.push_back("HLT_HIDoubleEle10Gsf_v10");
   Name.push_back("HLT_HIDoubleEle10GsfMass50_v10");
   Name.push_back("HLT_HIDoubleEle15Gsf_v10");
   Name.push_back("HLT_HIDoubleEle15GsfMass50_v10");

   std::sort(Name.begin(), Name.end());
   std::vector<std::string>::iterator iter = std::unique(Name.begin(), Name.end());
   Name.erase(iter, Name.end());
   Decision.resize(Name.size());
   Prescale.resize(Name.size());
   PrescaleNumerator.resize(Name.size());
   PrescaleDenominator.resize(Name.size());
   Exist.resize(Name.size());
   PrescaleExist.resize(Name.size());
}

int TriggerTreeMessenger::FindIndex(std::string Trigger)
{
   std::vector<std::string>::iterator iter
      = std::lower_bound(Name.begin(), Name.end(), Trigger);

   if(iter == Name.end())
      return -1;
   if(Trigger == *iter)
      return (iter - Name.begin());
   return -1;
}

int TriggerTreeMessenger::CheckTrigger(std::string Trigger)
{
   return CheckTrigger(FindIndex(Trigger));
}
   
int TriggerTreeMessenger::CheckTriggerStartWith(std::string Trigger)
{
   bool Found = false;
   int Result = 0;

   for(int i = 0; i < (int)Name.size(); i++)
   {
      if(Exist[i] == false || Name[i].find(Trigger) != 0)
         continue;

      Found = true;
      Result = Result + CheckTrigger(i);
   }

   if(Found == true)
      return Result;
   return -1;
}
   
int TriggerTreeMessenger::CheckTrigger(int Index)
{
   if(Index < 0 || Index >= (int)Name.size())
      return -1;
   if(Exist[Index] == false)
      return -1;

   return Decision[Index];
}

int TriggerTreeMessenger::GetPrescale(std::string Trigger)
{
   return GetPrescale(FindIndex(Trigger));
}
   
int TriggerTreeMessenger::GetPrescaleStartWith(std::string Trigger)
{
   for(int i = 0; i < (int)Name.size(); i++)
   {
      if(Exist[i] == false || Name[i].find(Trigger) != 0)
         continue;
      return GetPrescale(i);
   }

   return -1;
}

int TriggerTreeMessenger::GetPrescale(int Index)
{
   if(Index < 0 || Index >= (int)Name.size())
      return -1;
   if(Exist[Index] == false)
      return -1;
   if(PrescaleExist[Index] == false)
      return 1;

   return Prescale[Index];
}

TriggerObjectTreeMessenger::TriggerObjectTreeMessenger(TFile &File, std::string TreeName)
{
   Tree = (TTree *)File.Get(TreeName.c_str());
   Initialize();
}

TriggerObjectTreeMessenger::TriggerObjectTreeMessenger(TFile *File, std::string TreeName)
{
   if(File != nullptr)
      Tree = (TTree *)File->Get(TreeName.c_str());
   else
      Tree = nullptr;
   Initialize();
}

TriggerObjectTreeMessenger::TriggerObjectTreeMessenger(TTree *TriggerTree)
{
   Initialize(TriggerTree);
}

bool TriggerObjectTreeMessenger::Initialize(TTree *TriggerTree)
{
   Tree = TriggerTree;
   return Initialize();
}

bool TriggerObjectTreeMessenger::Initialize()
{
   ID = nullptr;
   PT = nullptr;
   Eta = nullptr;
   Phi = nullptr;
   Mass = nullptr;

   if(Tree == nullptr)
      return false;

   Tree->SetBranchAddress("TriggerObjID", &ID);
   Tree->SetBranchAddress("pt", &PT);
   Tree->SetBranchAddress("eta", &Eta);
   Tree->SetBranchAddress("phi", &Phi);
   Tree->SetBranchAddress("mass", &Mass);

   return true;
}

bool TriggerObjectTreeMessenger::GetEntry(int iEntry)
{
   if(Tree == nullptr)
      return false;

   Tree->GetEntry(iEntry);

   return true;
}

TrackTreeMessenger::TrackTreeMessenger(TFile &File, std::string TreeName)
{
   Tree = (TTree *)File.Get(TreeName.c_str());
   Initialize();
}

TrackTreeMessenger::TrackTreeMessenger(TFile *File, std::string TreeName)
{
   if(File != nullptr)
      Tree = (TTree *)File->Get(TreeName.c_str());
   else
      Tree = nullptr;
   Initialize();
}

TrackTreeMessenger::TrackTreeMessenger(TTree *TrackTree)
{
   Initialize(TrackTree);
}

bool TrackTreeMessenger::Initialize(TTree *TrackTree)
{
   Tree = TrackTree;
   return Initialize();
}

bool TrackTreeMessenger::Initialize()
{
   nVtx = 0;

   if(Tree == nullptr)
      return false;

   Tree->SetBranchAddress("nVtx", &nVtx);
   Tree->SetBranchAddress("nTrkVtx", &nTrkVtx);
   Tree->SetBranchAddress("normChi2Vtx", &normChi2Vtx);
   Tree->SetBranchAddress("sumPtVtx", &sumPtVtx);
   Tree->SetBranchAddress("xVtx", &xVtx);
   Tree->SetBranchAddress("yVtx", &yVtx);
   Tree->SetBranchAddress("zVtx", &zVtx);
   Tree->SetBranchAddress("xVtxErr", &xVtxErr);
   Tree->SetBranchAddress("yVtxErr", &yVtxErr);
   Tree->SetBranchAddress("zVtxErr", &zVtxErr);
   Tree->SetBranchAddress("vtxDist2D", &vtxDist2D);
   Tree->SetBranchAddress("vtxDist2DErr", &vtxDist2DErr);
   Tree->SetBranchAddress("vtxDist2DSig", &vtxDist2DSig);
   Tree->SetBranchAddress("vtxDist3D", &vtxDist3D);
   Tree->SetBranchAddress("vtxDist3DErr", &vtxDist3DErr);
   Tree->SetBranchAddress("vtxDist3DSig", &vtxDist3DSig);
   Tree->SetBranchAddress("nTrk", &nTrk);
   Tree->SetBranchAddress("trkPt", &trkPt);
   Tree->SetBranchAddress("trkPtError", &trkPtError);
   Tree->SetBranchAddress("trkNHit", &trkNHit);
   Tree->SetBranchAddress("trkNlayer", &trkNlayer);
   Tree->SetBranchAddress("trkEta", &trkEta);
   Tree->SetBranchAddress("trkPhi", &trkPhi);
   Tree->SetBranchAddress("trkCharge", &trkCharge);
   Tree->SetBranchAddress("trkNVtx", &trkNVtx);
   Tree->SetBranchAddress("nTrkTimesnVtx", &nTrkTimesnVtx);
   Tree->SetBranchAddress("trkAssocVtx", &trkAssocVtx);
   Tree->SetBranchAddress("trkDxyOverDxyError", &trkDxyOverDxyError);
   Tree->SetBranchAddress("trkDzOverDzError", &trkDzOverDzError);
   Tree->SetBranchAddress("highPurity", &highPurity);
   Tree->SetBranchAddress("tight", &tight);
   Tree->SetBranchAddress("loose", &loose);
   Tree->SetBranchAddress("trkChi2", &trkChi2);
   Tree->SetBranchAddress("trkNdof", &trkNdof);
   Tree->SetBranchAddress("trkDxy1", &trkDxy1);
   Tree->SetBranchAddress("trkDxyError1", &trkDxyError1);
   Tree->SetBranchAddress("trkDz1", &trkDz1);
   Tree->SetBranchAddress("trkDzError1", &trkDzError1);
   Tree->SetBranchAddress("trkFake", &trkFake);
   Tree->SetBranchAddress("trkAlgo", &trkAlgo);
   Tree->SetBranchAddress("trkOriginalAlgo", &trkOriginalAlgo);

   return true;
}

bool TrackTreeMessenger::GetEntry(int iEntry)
{
   if(Tree == nullptr)
      return false;

   nVtx = 0;
   Tree->GetEntry(iEntry);

   return true;
}

int TrackTreeMessenger::GetBestVertexIndex()
{
   if(nVtx == 0)
      return -1;
   if(nVtx == 1)
      return 0;

   int BestIndex = 0;
   double BestSumPT = 0;
   for(int iV = 0; iV < nVtx; iV++)
   {
      if(BestSumPT < sumPtVtx[iV])
      {
         BestIndex = iV;
         BestSumPT = sumPtVtx[iV];
      }
   }

   return BestIndex;
}

bool TrackTreeMessenger::PassZHadron2022Cut(int index)
{
   if(index >= nTrk)
      return false;

   if(highPurity[index] == false)
      return false;

   double RelativeUncertainty = trkPtError[index] / trkPt[index];
   if(RelativeUncertainty > 0.1)
      return false;

   // if(trkDxyOverDxyError[index] > 3)
   //    return false;
   if(fabs(trkDxy1[index]) / trkDxyError1[index] > 3)
      return false;

   // if(trkDzOverDzError[index] > 3)
   //    return false;
   if(fabs(trkDz1[index]) / trkDzError1[index] > 3)
      return false;

   if(trkEta[index] < -2.4 || trkEta[index] > 2.4)
      return false;

   return true;
}

bool TrackTreeMessenger::PassZHadron2022CutLoose(int index)
{
   // TODO

   if(index >= nTrk)
      return false;

   if(highPurity[index] == false)
      return false;

   double RelativeUncertainty = trkPtError[index] / trkPt[index];
   if(RelativeUncertainty > 0.15)
      return false;

   // if(trkDxyOverDxyError[index] > 3)
   //    return false;
   if(fabs(trkDxy1[index]) / trkDxyError1[index] > 5)
      return false;

   // if(trkDzOverDzError[index] > 3)
   //    return false;
   if(fabs(trkDz1[index]) / trkDzError1[index] > 5)
      return false;

   if(trkEta[index] < -2.4 || trkEta[index] > 2.4)
      return false;

   return true;
}

bool TrackTreeMessenger::PassZHadron2022CutTight(int index)
{
   // TODO
   
   if(index >= nTrk)
      return false;

   if(highPurity[index] == false)
      return false;

   double RelativeUncertainty = trkPtError[index] / trkPt[index];
   if(RelativeUncertainty > 0.15)
      return false;

   // if(trkDxyOverDxyError[index] > 3)
   //    return false;
   if(fabs(trkDxy1[index]) / trkDxyError1[index] > 2)
      return false;

   // if(trkDzOverDzError[index] > 3)
   //    return false;
   if(fabs(trkDz1[index]) / trkDzError1[index] > 2)
      return false;

   if(trkEta[index] < -2.4 || trkEta[index] > 2.4)
      return false;

   return true;
}

MuTreeMessenger::MuTreeMessenger(TFile &File, std::string TreeName)
{
   Tree = (TTree *)File.Get(TreeName.c_str());
   Initialize();
}

MuTreeMessenger::MuTreeMessenger(TFile *File, std::string TreeName)
{
   if(File != nullptr)
      Tree = (TTree *)File->Get(TreeName.c_str());
   else
      Tree = nullptr;
   Initialize();
}

MuTreeMessenger::MuTreeMessenger(TTree *MuTree)
{
   Initialize(MuTree);
}
   
bool MuTreeMessenger::Initialize(TTree *MuTree)
{
   Tree = MuTree;
   return Initialize();
}

bool MuTreeMessenger::Initialize()
{
   if(Tree == nullptr)
      return false;

   Tree->SetBranchAddress("Run", &Run);
   Tree->SetBranchAddress("Event", &Event);
   Tree->SetBranchAddress("Lumi", &Lumi);
   Tree->SetBranchAddress("vx", &Vx);
   Tree->SetBranchAddress("vy", &Vy);
   Tree->SetBranchAddress("vz", &Vz);
   Tree->SetBranchAddress("Gen_nptl", &NGen);
   Tree->SetBranchAddress("Gen_pid", &GenPID);
   Tree->SetBranchAddress("Gen_mom", &GenMom);
   Tree->SetBranchAddress("Gen_status", &GenStatus);
   Tree->SetBranchAddress("Gen_p", &GenP);
   Tree->SetBranchAddress("Gen_pt", &GenPT);
   Tree->SetBranchAddress("Gen_eta", &GenEta);
   Tree->SetBranchAddress("Gen_phi", &GenPhi);
   Tree->SetBranchAddress("Glb_nptl", &NGlb);
   Tree->SetBranchAddress("Glb_charge", &GlbCharge);
   Tree->SetBranchAddress("Glb_p", &GlbP);
   Tree->SetBranchAddress("Glb_pt", &GlbPT);
   Tree->SetBranchAddress("Glb_eta", &GlbEta);
   Tree->SetBranchAddress("Glb_phi", &GlbPhi);
   Tree->SetBranchAddress("Glb_dxy", &GlbDxy);
   Tree->SetBranchAddress("Glb_dz", &GlbDz);
   Tree->SetBranchAddress("Glb_nValMuHits", &GlbNValMuHits);
   Tree->SetBranchAddress("Glb_nValTrkHits", &GlbNValTrkHits);
   Tree->SetBranchAddress("Glb_nValPixHits", &GlbNValPixHits);
   Tree->SetBranchAddress("Glb_trkLayerWMeas", &GlbTrkLayerWMeas);
   Tree->SetBranchAddress("Glb_nMatchedStations", &GlbNMatchedStations);
   Tree->SetBranchAddress("Glb_nTrkFound", &GlbNTrkFound);
   Tree->SetBranchAddress("Glb_glbChi2_ndof", &GlbGlbChi2NDof);
   Tree->SetBranchAddress("Glb_trkChi2_ndof", &GlbTrkChi2NDof);
   Tree->SetBranchAddress("Glb_pixLayerWMeas", &GlbPixLayerWMeas);
   Tree->SetBranchAddress("Glb_trkDxy", &GlbTrkDxy);
   Tree->SetBranchAddress("Glb_trkDz", &GlbTrkDz);
   Tree->SetBranchAddress("Sta_nptl", &NSta);
   Tree->SetBranchAddress("Sta_charge", &StaCharge);
   Tree->SetBranchAddress("Sta_p", &StaP);
   Tree->SetBranchAddress("Sta_pt", &StaPT);
   Tree->SetBranchAddress("Sta_eta", &StaEta);
   Tree->SetBranchAddress("Sta_phi", &StaPhi);
   Tree->SetBranchAddress("Sta_dxy", &StaDxy);
   Tree->SetBranchAddress("Sta_dz", &StaDz);
   Tree->SetBranchAddress("Glb_isArbitrated", &GlbIsArbitrated);
   Tree->SetBranchAddress("Di_npair", &NDi);
   Tree->SetBranchAddress("Di_vProb", &DiVProb);
   Tree->SetBranchAddress("Di_mass", &DiMass);
   Tree->SetBranchAddress("Di_e", &DiE);
   Tree->SetBranchAddress("Di_pt", &DiPT);
   Tree->SetBranchAddress("Di_pt1", &DiPT1);
   Tree->SetBranchAddress("Di_pt2", &DiPT2);
   Tree->SetBranchAddress("Di_eta", &DiEta);
   Tree->SetBranchAddress("Di_eta1", &DiEta1);
   Tree->SetBranchAddress("Di_eta2", &DiEta2);
   Tree->SetBranchAddress("Di_rapidity", &DiRapidity);
   Tree->SetBranchAddress("Di_phi", &DiPhi);
   Tree->SetBranchAddress("Di_phi1", &DiPhi1);
   Tree->SetBranchAddress("Di_phi2", &DiPhi2);
   Tree->SetBranchAddress("Di_charge", &DiCharge);
   Tree->SetBranchAddress("Di_charge1", &DiCharge1);
   Tree->SetBranchAddress("Di_charge2", &DiCharge2);
   Tree->SetBranchAddress("Di_isArb1", &DiIsArb1);
   Tree->SetBranchAddress("Di_isArb2", &DiIsArb2);
   Tree->SetBranchAddress("Di_nTrkHit1", &DiNTrkHit1);
   Tree->SetBranchAddress("Di_nTrkHit2", &DiNTrkHit2);
   Tree->SetBranchAddress("Di_nMuHit1", &DiNMuHit1);
   Tree->SetBranchAddress("Di_nMuHit2", &DiNMuHit2);
   Tree->SetBranchAddress("Di_nTrkLayers1", &DiNTrkLayers1);
   Tree->SetBranchAddress("Di_nTrkLayers2", &DiNTrkLayers2);
   Tree->SetBranchAddress("Di_nPixHit1", &DiNPixHit1);
   Tree->SetBranchAddress("Di_nPixHit2", &DiNPixHit2);
   Tree->SetBranchAddress("Di_nMatchedStations1", &DiNMatchedStations1);
   Tree->SetBranchAddress("Di_nMatchedStations2", &DiNMatchedStations2);
   Tree->SetBranchAddress("Di_trkChi2_1", &DiTrkChi21);
   Tree->SetBranchAddress("Di_trkChi2_2", &DiTrkChi22);
   Tree->SetBranchAddress("Di_glbChi2_1", &DiGlbChi21);
   Tree->SetBranchAddress("Di_glbChi2_2", &DiGlbChi22);
   Tree->SetBranchAddress("Di_dxy1", &DiDxy1);
   Tree->SetBranchAddress("Di_dxy2", &DiDxy2);
   Tree->SetBranchAddress("Di_dz1", &DiDz1);
   Tree->SetBranchAddress("Di_dz2", &DiDz2);
   
   return true;
}

bool MuTreeMessenger::GetEntry(int iEntry)
{
   if(Tree == nullptr)
      return false;

   Tree->GetEntry(iEntry);
   return true;
}

bool MuTreeMessenger::DimuonPassTightCut(int index)
{
   // If somebody requested a pair that does not exist, we return false
   if(index >= NDi)
      return false;

   bool TightCut = true;

   if(DiGlbChi21[index] >= 10. || DiGlbChi22[index] >= 10.)
      TightCut = false;
   if(DiNMuHit1[index] <= 0 || DiNMuHit2[index] <= 0)
      TightCut = false;
   if(DiNMatchedStations1[index] <=1 || DiNMatchedStations2[index] <= 1)
      TightCut = false;
   if(fabs(DiDxy1[index]) >= 0.2 || fabs(DiDxy2[index]) >= 0.2)
      TightCut = false;
   if(fabs(DiDz1[index]) >= 0.5 || fabs(DiDz2[index]) >= 0.5)
      TightCut = false;
   if(DiNPixHit1[index] <= 0 || DiNPixHit2[index] <= 0)
      TightCut = false;
   if(DiNTrkLayers1[index] <= 5 || DiNTrkLayers2[index] <= 5)
      TightCut = false;

   // Note that the PF muon is not checked right now due to information not stored in the tree
   
   return TightCut;
}

PbPbTrackTreeMessenger::PbPbTrackTreeMessenger(TFile &File, std::string TreeName)
{
   Tree = (TTree *)File.Get(TreeName.c_str());
   Initialize();
}

PbPbTrackTreeMessenger::PbPbTrackTreeMessenger(TFile *File, std::string TreeName)
{
   if(File != nullptr)
      Tree = (TTree *)File->Get(TreeName.c_str());
   else
      Tree = nullptr;
   Initialize();
}

PbPbTrackTreeMessenger::PbPbTrackTreeMessenger(TTree *PbPbTrackTree)
{
   Initialize(PbPbTrackTree);
}
   
bool PbPbTrackTreeMessenger::Initialize(TTree *PbPbTrackTree)
{
   Tree = PbPbTrackTree;
   return Initialize();
}

bool PbPbTrackTreeMessenger::Initialize()
{
   if(Tree == nullptr)
      return false;

   VX = nullptr;
   VY = nullptr;
   VZ = nullptr;
   VXError = nullptr;
   VYError = nullptr;
   VZError = nullptr;
   VChi2 = nullptr;
   VNDof = nullptr;
   VIsFake = nullptr;
   VNTracks = nullptr;
   VPTSum = nullptr;
   TrackPT = nullptr;
   TrackPTError = nullptr;
   TrackEta = nullptr;
   TrackPhi = nullptr;
   TrackCharge = nullptr;
   TrackPDGID = nullptr;
   TrackNHits = nullptr;
   TrackNPixHits = nullptr;
   TrackNLayers = nullptr;
   TrackNormChi2 = nullptr;
   TrackHighPurity = nullptr;
   PFEnergy = nullptr;
   PFEcal = nullptr;
   PFHcal = nullptr;
   TrackAssociatedVertexIndex = nullptr;
   TrackAssociatedVertexQuality = nullptr;
   TrackAssociatedVertexDz = nullptr;
   TrackAssociatedVertexDzError = nullptr;
   TrackAssociatedVertexDxy = nullptr;
   TrackAssociatedVertexDxyError = nullptr;
   TrackFirstVertexQuality = nullptr;
   TrackFirstVertexDz = nullptr;
   TrackFirstVertexDzError = nullptr;
   TrackFirstVertexDxy = nullptr;
   TrackFirstVertexDxyError = nullptr;   

   Tree->SetBranchAddress("nRun", &Run);
   Tree->SetBranchAddress("nEv", &Event);
   Tree->SetBranchAddress("nLumi", &Lumi);
   Tree->SetBranchAddress("xVtx", &VX);
   Tree->SetBranchAddress("yVtx", &VY);
   Tree->SetBranchAddress("zVtx", &VZ);
   Tree->SetBranchAddress("xErrVtx", &VXError);
   Tree->SetBranchAddress("yErrVtx", &VYError);
   Tree->SetBranchAddress("zErrVtx", &VZError);
   Tree->SetBranchAddress("chi2Vtx", &VChi2);
   Tree->SetBranchAddress("ndofVtx", &VNDof);
   Tree->SetBranchAddress("isFakeVtx", &VIsFake);
   Tree->SetBranchAddress("nTracksVtx", &VNTracks);
   Tree->SetBranchAddress("ptSumVtx", &VPTSum);
   Tree->SetBranchAddress("trkPt", &TrackPT);
   Tree->SetBranchAddress("trkPtError", &TrackPTError);
   Tree->SetBranchAddress("trkEta", &TrackEta);
   Tree->SetBranchAddress("trkPhi", &TrackPhi);
   Tree->SetBranchAddress("trkCharge", &TrackCharge);
   Tree->SetBranchAddress("trkPDGId", &TrackPDGID);
   Tree->SetBranchAddress("trkNHits", &TrackNHits);
   Tree->SetBranchAddress("trkNPixHits", &TrackNPixHits);
   Tree->SetBranchAddress("trkNLayers", &TrackNLayers);
   Tree->SetBranchAddress("trkNormChi2", &TrackNormChi2);
   Tree->SetBranchAddress("highPurity", &TrackHighPurity);
   Tree->SetBranchAddress("pfEnergy", &PFEnergy);
   Tree->SetBranchAddress("pfEcal", &PFEcal);
   Tree->SetBranchAddress("pfHcal", &PFHcal);
   Tree->SetBranchAddress("trkAssociatedVtxIndx", &TrackAssociatedVertexIndex);
   Tree->SetBranchAddress("trkAssociatedVtxQuality", &TrackAssociatedVertexQuality);
   Tree->SetBranchAddress("trkDzAssociatedVtx", &TrackAssociatedVertexDz);
   Tree->SetBranchAddress("trkDzErrAssociatedVtx", &TrackAssociatedVertexDzError);
   Tree->SetBranchAddress("trkDxyAssociatedVtx", &TrackAssociatedVertexDxy);
   Tree->SetBranchAddress("trkDxyErrAssociatedVtx", &TrackAssociatedVertexDxyError);
   Tree->SetBranchAddress("trkFirstVtxQuality", &TrackFirstVertexQuality);
   Tree->SetBranchAddress("trkDzFirstVtx", &TrackFirstVertexDz);
   Tree->SetBranchAddress("trkDzErrFirstVtx", &TrackFirstVertexDzError);
   Tree->SetBranchAddress("trkDxyFirstVtx", &TrackFirstVertexDxy);
   Tree->SetBranchAddress("trkDxyErrFirstVtx", &TrackFirstVertexDxyError);

   return true;
}

bool PbPbTrackTreeMessenger::GetEntry(int iEntry)
{
   if(Tree == nullptr)
      return false;

   Tree->GetEntry(iEntry);
   return true;
}

bool PbPbTrackTreeMessenger::PassZHadron2022Cut(int index)
{
   if(TrackPT == nullptr)
      return false;
   if(index >= TrackPT->size())
      return false;

   if(TrackHighPurity->at(index) == false)
      return false;

   double RelativeUncertainty = TrackPTError->at(index) / TrackPT->at(index);
   if(RelativeUncertainty >= 0.1)
      return false;

   //if(TrackCharge==0)
   //   return false;

   // double XYVertexSignificance = fabs(TrackAssociatedVertexDxy->at(index) / TrackAssociatedVertexDxyError->at(index));
   double XYVertexSignificance = fabs(TrackFirstVertexDxy->at(index) / TrackFirstVertexDxyError->at(index));
   if(XYVertexSignificance >= 3)
      return false;

   // double ZVertexSignificance = fabs(TrackAssociatedVertexDz->at(index) / TrackAssociatedVertexDzError->at(index));
   double ZVertexSignificance = fabs(TrackFirstVertexDz->at(index) / TrackFirstVertexDzError->at(index));
   if(ZVertexSignificance >= 3)
      return false;

   if(TrackNHits->at(index) < 11)
      return false;

   if(TrackNormChi2->at(index) / TrackNLayers->at(index) >= 0.18)
      return false;

   double ECAL = -1, HCAL = -1;
   if(PFEcal != nullptr && PFEcal->size() > index)   ECAL = PFEcal->at(index);
   if(PFHcal != nullptr && PFHcal->size() > index)   HCAL = PFHcal->at(index);
   if(TrackPT->at(index) > 20 && ((ECAL + HCAL) / cosh(TrackEta->at(index)) < 0.5 * TrackPT->at(index)))
      return false;

   if(fabs(TrackEta->at(index)) > 2.4)
      return false;

   return true;
}

bool PbPbTrackTreeMessenger::PassZHadron2022CutNoVZ(int index)
{
   if(TrackPT == nullptr)
      return false;
   if(index >= TrackPT->size())
      return false;

   if(TrackHighPurity->at(index) == false)
      return false;

   double RelativeUncertainty = TrackPTError->at(index) / TrackPT->at(index);
   if(RelativeUncertainty >= 0.1)
      return false;

   // double XYVertexSignificance = fabs(TrackAssociatedVertexDxy->at(index) / TrackAssociatedVertexDxyError->at(index));
   double XYVertexSignificance = fabs(TrackFirstVertexDxy->at(index) / TrackFirstVertexDxyError->at(index));
   if(XYVertexSignificance >= 3)
      return false;

   // double ZVertexSignificance = fabs(TrackAssociatedVertexDz->at(index) / TrackAssociatedVertexDzError->at(index));
   // double ZVertexSignificance = fabs(TrackFirstVertexDz->at(index) / TrackFirstVertexDzError->at(index));
   // if(ZVertexSignificance >= 3)
   //    return false;

   if(TrackNHits->at(index) < 11)
      return false;

   if(TrackNormChi2->at(index) / TrackNLayers->at(index) >= 0.18)
      return false;

   double ECAL = -1, HCAL = -1;
   if(PFEcal != nullptr && PFEcal->size() > index)   ECAL = PFEcal->at(index);
   if(PFHcal != nullptr && PFHcal->size() > index)   HCAL = PFHcal->at(index);
   if(TrackPT->at(index) > 20 && ((ECAL + HCAL) / cosh(TrackEta->at(index)) < 0.5 * TrackPT->at(index)))
       return false;

   if(fabs(TrackEta->at(index)) > 2.4)
      return false;

   return true;
}

bool PbPbTrackTreeMessenger::PassZHadron2022CutLoose(int index)
{
   if(TrackPT == nullptr)
      return false;
   if(index >= TrackPT->size())
      return false;

   if(TrackHighPurity->at(index) == false)
      return false;

   double RelativeUncertainty = TrackPTError->at(index) / TrackPT->at(index);
   if(RelativeUncertainty >= 0.15)
      return false;

   // double XYVertexSignificance = fabs(TrackAssociatedVertexDxy->at(index) / TrackAssociatedVertexDxyError->at(index));
   double XYVertexSignificance = fabs(TrackFirstVertexDxy->at(index) / TrackFirstVertexDxyError->at(index));
   if(XYVertexSignificance >= 5)
      return false;

   // double ZVertexSignificance = fabs(TrackAssociatedVertexDz->at(index) / TrackAssociatedVertexDzError->at(index));
   double ZVertexSignificance = fabs(TrackFirstVertexDz->at(index) / TrackFirstVertexDzError->at(index));
   if(ZVertexSignificance >= 5)
      return false;

   if(TrackNHits->at(index) < 11)
      return false;

   if(TrackNormChi2->at(index) / TrackNLayers->at(index) >= 0.18)
      return false;

   double ECAL = -1, HCAL = -1;
   if(PFEcal != nullptr && PFEcal->size() > index)   ECAL = PFEcal->at(index);
   if(PFHcal != nullptr && PFHcal->size() > index)   HCAL = PFHcal->at(index);
   if(TrackPT->at(index) > 20 && ((ECAL + HCAL) / cosh(TrackEta->at(index)) < 0.5 * TrackPT->at(index)))
       return false;

   if(fabs(TrackEta->at(index)) > 2.4)
      return false;

   return true;
}

bool PbPbTrackTreeMessenger::PassZHadron2022CutTight(int index)
{
   if(TrackPT == nullptr)
      return false;
   if(index >= TrackPT->size())
      return false;

   if(TrackHighPurity->at(index) == false)
      return false;

   double RelativeUncertainty = TrackPTError->at(index) / TrackPT->at(index);
   if(RelativeUncertainty >= 0.05)
      return false;

   // double XYVertexSignificance = fabs(TrackAssociatedVertexDxy->at(index) / TrackAssociatedVertexDxyError->at(index));
   double XYVertexSignificance = fabs(TrackFirstVertexDxy->at(index) / TrackFirstVertexDxyError->at(index));
   if(XYVertexSignificance >= 2)
      return false;

   // double ZVertexSignificance = fabs(TrackAssociatedVertexDz->at(index) / TrackAssociatedVertexDzError->at(index));
   double ZVertexSignificance = fabs(TrackFirstVertexDz->at(index) / TrackFirstVertexDzError->at(index));
   if(ZVertexSignificance >= 2)
      return false;

   if(TrackNHits->at(index) < 11)
      return false;

   if(TrackNormChi2->at(index) / TrackNLayers->at(index) >= 0.15)
      return false;

   double ECAL = -1, HCAL = -1;
   if(PFEcal != nullptr && PFEcal->size() > index)   ECAL = PFEcal->at(index);
   if(PFHcal != nullptr && PFHcal->size() > index)   HCAL = PFHcal->at(index);
   if(TrackPT->at(index) > 20 && ((ECAL + HCAL) / cosh(TrackEta->at(index)) < 0.5 * TrackPT->at(index)))
       return false;

   if(fabs(TrackEta->at(index)) > 2.4)
      return false;

   return true;
}

ZHadronMessenger::ZHadronMessenger(TFile &File, std::string TreeName)
{
   Initialized = false;
   WriteMode = false;
   
   Tree = (TTree *)File.Get(TreeName.c_str());
   Initialize();
}

ZHadronMessenger::ZHadronMessenger(TFile *File, std::string TreeName)
{
   Initialized = false;
   WriteMode = false;
   
   if(File != nullptr)
      Tree = (TTree *)File->Get(TreeName.c_str());
   else
      Tree = nullptr;
   Initialize();
}

ZHadronMessenger::ZHadronMessenger(TTree *ZHadronTree)
{
   Initialized = false;
   WriteMode = false;

   Initialize(ZHadronTree);
}
   
ZHadronMessenger::~ZHadronMessenger()
{
   if(Initialized == true && WriteMode == true)
   {
      delete zMass;
      delete zEta;
      delete zPhi;
      delete zY;
      delete zPt;
      delete genZMass;
      delete genZEta;
      delete genZPhi;
      delete genZY;
      delete genZPt;
      delete trackPt;
      delete trackDeta;
      delete trackDphi;
      delete trackDY;
      delete trackMuTagged;
      delete trackMuDR;
      delete TrackPDGID;
      delete trackEta;
      delete trackPhi;
      delete trackY;
      delete trackCharge;

      delete GenTrackDphi;
      delete GenTrackDeta;
      delete GenTrackDY;
      delete GenTrackPt;
      delete GenTrackMuTagged;
      delete GenTrackMuDR;
      delete GenSubevent;
      delete GenTrackEta;
      delete GenTrackPhi;
      delete GenTrackY;
      delete GenTrackCharge;

      delete jetPt;
      delete jetDeta;
      delete jetDphi;
      delete jetRefPt;
      delete jetRefDeta;
      delete jetRefDphi;
      delete jetMuTagged;
      
      delete genJetPt;
      delete genJetEta;
      delete genJetPhi;

      delete muEta1;
      delete muEta2;
      delete muPhi1;
      delete muPhi2;
      delete muPt1;
      delete muPt2;

      delete muDeta;
      delete muDphi;
      delete muDR;
      delete muDphiS;

      delete genMuPt1;
      delete genMuPt2;
      delete genMuEta1;
      delete genMuEta2;
      delete genMuPhi1;
      delete genMuPhi2;

      delete genMuDeta;
      delete genMuDphi;
      delete genMuDR;
      delete genMuDphiS;

      delete isMuon;
      delete isGenMuon;
   }
}
   
bool ZHadronMessenger::Initialize(TTree *ZHadronTree)
{
   Tree = ZHadronTree;
   return Initialize();
}

bool ZHadronMessenger::Initialize()
{
   if(Tree == nullptr)
      return false;

   Initialized = true;

   zMass = nullptr;
   zEta = nullptr;
   zPhi = nullptr;
   zPt = nullptr;
   zY = nullptr;
   genZMass = nullptr;
   genZEta = nullptr;
   genZPhi = nullptr;
   genZY = nullptr;
   genZPt = nullptr;
   trackDphi = nullptr;
   trackPt = nullptr;
   trackDeta = nullptr;
   trackDY = nullptr;
   TrackPDGID = nullptr;
   trackEta = nullptr;
   trackPhi = nullptr;
   trackY = nullptr;
   trackMuTagged = nullptr;
   trackMuDR = nullptr;
   trackWeight = nullptr;
   trackResidualWeight = nullptr;
   trackCharge = nullptr;
   subevent = nullptr;

   GenTrackDphi = nullptr;
   GenTrackDeta = nullptr;
   GenTrackDY = nullptr;
   GenTrackPt = nullptr;
   GenTrackMuTagged = nullptr;
   GenTrackMuDR = nullptr;
   GenSubevent = nullptr;
   GenTrackEta = nullptr;
   GenTrackPhi = nullptr;
   GenTrackY = nullptr;
   GenTrackCharge = nullptr;

   jetPt = nullptr;
   jetDeta = nullptr;
   jetDphi = nullptr;
   jetRefPt = nullptr;
   jetRefDeta = nullptr;
   jetRefDphi = nullptr;
   jetMuTagged = nullptr;
   
   genJetPt = nullptr;
   genJetEta = nullptr;
   genJetPhi = nullptr;

   muEta1 = nullptr;
   muEta2 = nullptr;
   muPhi1 = nullptr;
   muPhi2 = nullptr;
   muPt1 = nullptr;
   muPt2 = nullptr;

   muDeta = nullptr;
   muDphi = nullptr;
   muDR = nullptr;
   muDphiS = nullptr;

   genMuPt1 = nullptr;
   genMuPt2 = nullptr;
   genMuEta1 = nullptr;
   genMuEta2 = nullptr;
   genMuPhi1 = nullptr;
   genMuPhi2 = nullptr;

   genMuDeta = nullptr;
   genMuDphi = nullptr;
   genMuDR = nullptr;
   genMuDphiS = nullptr;

   isMuon = nullptr;
   isGenMuon = nullptr;

   Tree->SetBranchAddress("run", &Run);
   Tree->SetBranchAddress("event", &Event);
   Tree->SetBranchAddress("lumi", &Lumi);

   Tree->SetBranchAddress("hiBin", &hiBin);
   Tree->SetBranchAddress("hiBinUp", &hiBinUp);
   Tree->SetBranchAddress("hiBinDown", &hiBinDown);
   Tree->SetBranchAddress("hiHF", &hiHF);
   
   Tree->SetBranchAddress("SignalHF", &SignalHF);
   Tree->SetBranchAddress("BackgroundHF", &BackgroundHF);
   Tree->SetBranchAddress("SignalGenHF", &SignalGenHF);
   Tree->SetBranchAddress("BackgroundGenHF", &BackgroundGenHF);
   Tree->SetBranchAddress("SignalVZ", &SignalVZ);

   Tree->SetBranchAddress("NCollWeight", &NCollWeight);
   Tree->SetBranchAddress("ZWeight", &ZWeight);
   Tree->SetBranchAddress("VZWeight", &VZWeight);
   Tree->SetBranchAddress("ExtraZWeight", &ExtraZWeight);

   Tree->SetBranchAddress("NVertex", &NVertex);
   Tree->SetBranchAddress("VX", &VX);
   Tree->SetBranchAddress("VY", &VY);
   Tree->SetBranchAddress("VZ", &VZ);
   Tree->SetBranchAddress("VXError", &VXError);
   Tree->SetBranchAddress("VYError", &VYError);
   Tree->SetBranchAddress("VZError", &VZError);
   Tree->SetBranchAddress("NPU", &NPU);
   
   Tree->SetBranchAddress("zMass", &zMass);
   Tree->SetBranchAddress("zEta", &zEta);
   Tree->SetBranchAddress("zPhi", &zPhi);
   Tree->SetBranchAddress("zY", &zY);
   Tree->SetBranchAddress("zPt", &zPt);
   Tree->SetBranchAddress("genZMass", &genZMass);
   Tree->SetBranchAddress("genZEta", &genZEta);
   Tree->SetBranchAddress("genZPhi", &genZPhi);
   Tree->SetBranchAddress("genZY", &genZY);
   Tree->SetBranchAddress("genZPt", &genZPt);
   Tree->SetBranchAddress("trackDphi", &trackDphi);
   Tree->SetBranchAddress("trackPt", &trackPt);
   Tree->SetBranchAddress("trackDeta", &trackDeta);
   Tree->SetBranchAddress("trackDY", &trackDY);
   Tree->SetBranchAddress("TrackPDGID", &TrackPDGID);
   Tree->SetBranchAddress("trackEta", &trackEta);
   Tree->SetBranchAddress("trackPhi", &trackPhi);
   Tree->SetBranchAddress("trackY", &trackY);
   Tree->SetBranchAddress("trackMuTagged", &trackMuTagged);
   Tree->SetBranchAddress("trackMuDR", &trackMuDR);
   Tree->SetBranchAddress("trackWeight", &trackWeight);
   Tree->SetBranchAddress("trackResidualWeight", &trackResidualWeight);
   Tree->SetBranchAddress("trackCharge", &trackCharge);
   Tree->SetBranchAddress("subevent", &subevent);

   Tree->SetBranchAddress("GenTrackDphi", &GenTrackDphi);
   Tree->SetBranchAddress("GenTrackDeta", &GenTrackDeta);
   Tree->SetBranchAddress("GenTrackDY", &GenTrackDY);
   Tree->SetBranchAddress("GenTrackPt", & GenTrackPt);
   Tree->SetBranchAddress("GenTrackMuTagged", & GenTrackMuTagged);
   Tree->SetBranchAddress("GenTrackMuDR", & GenTrackMuDR);
   Tree->SetBranchAddress("GenSubevent", & GenSubevent);
   Tree->SetBranchAddress("GenTrackEta", & GenTrackEta);
   Tree->SetBranchAddress("GenTrackPhi", & GenTrackPhi);
   Tree->SetBranchAddress("GenTrackY", & GenTrackY);
   Tree->SetBranchAddress("GenTrackCharge", & GenTrackCharge);

   Tree->SetBranchAddress("jetPt", &jetPt);
   Tree->SetBranchAddress("jetDeta", &jetDeta);
   Tree->SetBranchAddress("jetDphi", &jetDphi);
   Tree->SetBranchAddress("jetRefPt", &jetRefPt);
   Tree->SetBranchAddress("jetRefDeta", &jetRefDeta);
   Tree->SetBranchAddress("jetRefDphi", &jetRefDphi);
   Tree->SetBranchAddress("jetMuTagged", &jetMuTagged);
   
   Tree->SetBranchAddress("genJetPt", &genJetPt);
   Tree->SetBranchAddress("genJetEta", &genJetEta);
   Tree->SetBranchAddress("genJetPhi", &genJetPhi);

   Tree->SetBranchAddress("muEta1", &muEta1);
   Tree->SetBranchAddress("muEta2", &muEta2);
   Tree->SetBranchAddress("muPhi1", &muPhi1);
   Tree->SetBranchAddress("muPhi2", &muPhi2);
   Tree->SetBranchAddress("muPt1", &muPt1);
   Tree->SetBranchAddress("muPt2", &muPt2);
 
   Tree->SetBranchAddress("muDeta", &muDeta);
   Tree->SetBranchAddress("muDphi", &muDphi);
   Tree->SetBranchAddress("muDR", &muDR);
   Tree->SetBranchAddress("muDphiS", &muDphiS);
  
   Tree->SetBranchAddress("genMuPt1", &genMuPt1);
   Tree->SetBranchAddress("genMuPt2", &genMuPt2);
   Tree->SetBranchAddress("genMuEta1", &genMuEta1);
   Tree->SetBranchAddress("genMuEta2", &genMuEta2);
   Tree->SetBranchAddress("genMuPhi1", &genMuPhi1);
   Tree->SetBranchAddress("genMuPhi2", &genMuPhi2);
   
   Tree->SetBranchAddress("genMuDeta", &genMuDeta);
   Tree->SetBranchAddress("genMuDphi", &genMuDphi);
   Tree->SetBranchAddress("genMuDR", &genMuDR);
   Tree->SetBranchAddress("genMuDphiS", &genMuDphiS);

   Tree->SetBranchAddress("isMuon", &isMuon);
   Tree->SetBranchAddress("isGenMuon", &isGenMuon);

   Tree->SetBranchAddress("bestZidx", &bestZidx);
   Tree->SetBranchAddress("bestZgenIdx", &bestZgenIdx);
   
   return true;
}

int ZHadronMessenger::GetEntries()
{
   if(Tree == nullptr)
      return 0;
   return Tree->GetEntries();
}

bool ZHadronMessenger::GetEntry(int iEntry)
{
   if(Tree == nullptr)
      return false;

   Tree->GetEntry(iEntry);
   return true;
}

bool ZHadronMessenger::SetBranch(TTree *T)
{
   if(T == nullptr)
      return false;

   Initialized = true;
   WriteMode = true;

   zMass = new std::vector<double>();
   zEta = new std::vector<double>();
   zPhi = new std::vector<double>();
   zY = new std::vector<double>();
   zPt = new std::vector<double>();
   genZMass = new std::vector<double>();
   genZEta = new std::vector<double>();
   genZPhi = new std::vector<double>();
   genZY = new std::vector<double>();
   genZPt = new std::vector<double>();
   trackPt = new std::vector<double>();
   trackDeta = new std::vector<double>();
   trackDphi = new std::vector<double>();
   trackDY = new std::vector<double>();
   TrackPDGID = new std::vector<double>();
   trackEta = new std::vector<double>();
   trackPhi = new std::vector<double>();
   trackY = new std::vector<double>();
   trackMuTagged = new std::vector<bool>();
   trackMuDR = new std::vector<double>();
   trackWeight = new std::vector<double>();
   trackResidualWeight = new std::vector<double>();
   trackCharge = new std::vector<int>();
   subevent = new std::vector<int>();

   GenTrackDphi = new std::vector<double>();
   GenTrackDeta = new std::vector<double>();
   GenTrackDY = new std::vector<double>();
   GenTrackPt = new std::vector<double>();
   GenTrackMuTagged = new std::vector<bool>();
   GenTrackMuDR = new std::vector<double>();
   GenSubevent = new std::vector<int>();
   GenTrackEta = new std::vector<double>();
   GenTrackPhi = new std::vector<double>();
   GenTrackY = new std::vector<double>();
   GenTrackCharge = new std::vector<int>();

   jetPt = new std::vector<double>();
   jetDeta = new std::vector<double>();
   jetDphi = new std::vector<double>();
   jetRefPt = new std::vector<double>();
   jetRefDeta = new std::vector<double>();
   jetRefDphi = new std::vector<double>();
   jetMuTagged = new std::vector<bool>();
   
   genJetPt = new std::vector<double>();
   genJetEta = new std::vector<double>();
   genJetPhi = new std::vector<double>();

   muEta1 = new std::vector<double>();
   muEta2 = new std::vector<double>();
   muPhi1 = new std::vector<double>();
   muPhi2 = new std::vector<double>();
   muPt1 = new std::vector<double>();
   muPt2 = new std::vector<double>();

   muDeta = new std::vector<double>();
   muDphi = new std::vector<double>();
   muDR = new std::vector<double>();
   muDphiS = new std::vector<double>();

   genMuPt1 = new std::vector<double>();
   genMuPt2 = new std::vector<double>();
   genMuEta1 = new std::vector<double>();
   genMuEta2 = new std::vector<double>();
   genMuPhi1 = new std::vector<double>();
   genMuPhi2 = new std::vector<double>();

   genMuDeta = new std::vector<double>();
   genMuDphi = new std::vector<double>();
   genMuDR = new std::vector<double>();
   genMuDphiS = new std::vector<double>();

   isMuon = new std::vector<bool>();
   isGenMuon = new std::vector<bool>();

   Tree = T;

   Tree->Branch("run",                    &Run,          "run/I");
   Tree->Branch("event",                  &Event,        "event/L");
   Tree->Branch("lumi",                   &Lumi,         "lumi/I");
   Tree->Branch("hiBin",                  &hiBin,        "hiBin/I");
   Tree->Branch("hiBinUp",                &hiBinUp,      "hiBinUp/I");
   Tree->Branch("hiBinDown",              &hiBinDown,    "hiBinDown/I");
   Tree->Branch("hiHF",                   &hiHF,         "hiHF/F");

   Tree->Branch("SignalHF",               &SignalHF,     "SignalHF/F");
   Tree->Branch("BackgroundHF",           &BackgroundHF, "BackgroundHF/F");
   Tree->Branch("SignalGenHF",            &SignalGenHF,     "SignalGenHF/F");
   Tree->Branch("BackgroundGenHF",        &BackgroundGenHF, "BackgroundGenHF/F");
   Tree->Branch("SignalVZ",               &SignalVZ,     "SignalVZ/F");

   Tree->Branch("NCollWeight",            &NCollWeight,  "NCollWeight/F");
   Tree->Branch("ZWeight",                &ZWeight,      "ZWeight/F");
   Tree->Branch("VZWeight",               &VZWeight,     "VZWeight/F");
   Tree->Branch("ExtraZWeight",           &ExtraZWeight, "ExtraZWeight[12]/F");
   
   Tree->Branch("NVertex",                &NVertex,      "NVertex/I");
   Tree->Branch("VX",                     &VX,           "VX/F");
   Tree->Branch("VY",                     &VY,           "VY/F");
   Tree->Branch("VZ",                     &VZ,           "VZ/F");
   Tree->Branch("VXError",                &VXError,      "VXError/F");
   Tree->Branch("VYError",                &VYError,      "VYError/F");
   Tree->Branch("VZError",                &VZError,      "VZError/F");
   Tree->Branch("NPU",                    &NPU,          "NPU/I");

   Tree->Branch("zMass",                  &zMass);
   Tree->Branch("zEta",                   &zEta);
   Tree->Branch("zPhi",                   &zPhi);
   Tree->Branch("zY",                     &zY);
   Tree->Branch("zPt",                    &zPt);
   Tree->Branch("genZMass",               &genZMass);
   Tree->Branch("genZEta",                &genZEta);
   Tree->Branch("genZPhi",                &genZPhi);
   Tree->Branch("genZY",                  &genZY);
   Tree->Branch("genZPt",                 &genZPt);
   Tree->Branch("trackPt",                &trackPt);
   Tree->Branch("trackDeta",              &trackDeta);
   Tree->Branch("trackDphi",              &trackDphi);
   Tree->Branch("trackDY",                &trackDY);
   Tree->Branch("TrackPDGID",             &TrackPDGID);
   Tree->Branch("trackPhi",               &trackPhi);
   Tree->Branch("trackEta",               &trackEta);
   Tree->Branch("trackY",                 &trackY);
   Tree->Branch("trackMuTagged",          &trackMuTagged);
   Tree->Branch("trackMuDR",              &trackMuDR);
   Tree->Branch("trackWeight",            &trackWeight);
   Tree->Branch("trackResidualWeight",    &trackResidualWeight);
   Tree->Branch("trackCharge",            &trackCharge);
   Tree->Branch("subevent",               &subevent);

   Tree->Branch("GenTrackDphi",           &GenTrackDphi);
   Tree->Branch("GenTrackDeta",           &GenTrackDeta);
   Tree->Branch("GenTrackDY",             &GenTrackDY);
   Tree->Branch("GenTrackPt",             & GenTrackPt);
   Tree->Branch("GenTrackMuTagged",       & GenTrackMuTagged);
   Tree->Branch("GenTrackMuDR",           & GenTrackMuDR);
   Tree->Branch("GenSubevent",            & GenSubevent);
   Tree->Branch("GenTrackEta",            & GenTrackEta);
   Tree->Branch("GenTrackPhi",            & GenTrackPhi);
   Tree->Branch("GenTrackY",              & GenTrackY);
   Tree->Branch("GenTrackCharge",         & GenTrackCharge);
   
   Tree->Branch("jetPt",                  &jetPt);
   Tree->Branch("jetDeta",                &jetDeta);
   Tree->Branch("jetDphi",                &jetDphi);
   Tree->Branch("jetRefPt",               &jetRefPt);
   Tree->Branch("jetRefDeta",             &jetRefDeta);
   Tree->Branch("jetRefDphi",             &jetRefDphi);
   Tree->Branch("jetMuTagged",            &jetMuTagged);
   
   Tree->Branch("genJetPt",               &genJetPt);
   Tree->Branch("genJetEta",              &genJetEta);
   Tree->Branch("genJetPhi",              &genJetPhi);
 
   Tree->Branch("maxOppositeDEta",        &maxOppositeDEta);
   Tree->Branch("maxOppositeDPhi",        &maxOppositeDPhi);
   Tree->Branch("maxDEta",                &maxDEta);
   Tree->Branch("maxDPhi",                &maxDPhi);
   Tree->Branch("maxOppositeWTADEta",     &maxOppositeWTADEta);
   Tree->Branch("maxOppositeWTADPhi",     &maxOppositeWTADPhi);
   Tree->Branch("maxMoreOppositeWTADEta", &maxMoreOppositeWTADEta);
   Tree->Branch("maxMoreOppositeWTADPhi", &maxMoreOppositeWTADPhi);
   Tree->Branch("maxOppositeCSWTADEta",   &maxOppositeCSWTADEta);
   Tree->Branch("maxOppositeCSWTADPhi",   &maxOppositeCSWTADPhi);
   Tree->Branch("maxOppositeChargedWTADEta",     &maxOppositeChargedWTADEta);
   Tree->Branch("maxOppositeChargedWTADPhi",     &maxOppositeChargedWTADPhi);
   Tree->Branch("maxMoreOppositeChargedWTADEta", &maxMoreOppositeChargedWTADEta);
   Tree->Branch("maxMoreOppositeChargedWTADPhi", &maxMoreOppositeChargedWTADPhi);
   Tree->Branch("maxOppositeHardChargedWTADEta", &maxOppositeHardChargedWTADEta);
   Tree->Branch("maxOppositeHardChargedWTADPhi", &maxOppositeHardChargedWTADPhi);
   
   Tree->Branch("maxOppositeJet12Pt",     &maxOppositeJet12Pt);
   Tree->Branch("maxOppositeJet12DEta",   &maxOppositeJet12DEta);
   Tree->Branch("maxOppositeJet12DPhi",   &maxOppositeJet12DPhi);
   Tree->Branch("maxOppositeJet34Pt",     &maxOppositeJet34Pt);
   Tree->Branch("maxOppositeJet34DEta",   &maxOppositeJet34DEta);
   Tree->Branch("maxOppositeJet34DPhi",   &maxOppositeJet34DPhi);
   Tree->Branch("maxOppositeJet56Pt",     &maxOppositeJet56Pt);
   Tree->Branch("maxOppositeJet56DEta",   &maxOppositeJet56DEta);
   Tree->Branch("maxOppositeJet56DPhi",   &maxOppositeJet56DPhi);
   Tree->Branch("maxOppositeJet78Pt",     &maxOppositeJet78Pt);
   Tree->Branch("maxOppositeJet78DEta",   &maxOppositeJet78DEta);
   Tree->Branch("maxOppositeJet78DPhi",   &maxOppositeJet78DPhi);

   Tree->Branch("muEta1",                 &muEta1);
   Tree->Branch("muEta2",                 &muEta2);
   Tree->Branch("muPhi1",                 &muPhi1);
   Tree->Branch("muPhi2",                 &muPhi2);
   Tree->Branch("muPt1",                  &muPt1);
   Tree->Branch("muPt2",                  &muPt2);

   Tree->Branch("genMuPt1",               &genMuPt1);
   Tree->Branch("genMuPt2",               &genMuPt2);
   Tree->Branch("genMuEta1",              &genMuEta1);
   Tree->Branch("genMuEta2",              &genMuEta2);
   Tree->Branch("genMuPhi1",              &genMuPhi1);
   Tree->Branch("genMuPhi2",              &genMuPhi2);

   Tree->Branch("muDeta",                 &muDeta);
   Tree->Branch("muDphi",                 &muDphi);
   Tree->Branch("muDR",                   &muDR);
   Tree->Branch("muDphiS",                &muDphiS);

   Tree->Branch("genMuDeta",              &genMuDeta);
   Tree->Branch("genMuDphi",              &genMuDphi);
   Tree->Branch("genMuDR",                &genMuDR);
   Tree->Branch("genMuDphiS",             &genMuDphiS);

   Tree->Branch("isMuon",                 &isMuon);
   Tree->Branch("isGenMuon",              &isGenMuon);

   Tree->Branch("bestZidx",               &bestZidx,        "bestZidx/I");
   Tree->Branch("bestZgenIdx",            &bestZgenIdx,     "bestZgenIdx/I");



   return true;
}
   
void ZHadronMessenger::Clear()
{
   if(Initialized == false)
      return;

   SignalHF = -1;
   BackgroundHF = -1;
   SignalGenHF = -1;
   BackgroundGenHF = -1;
   SignalVZ = -999;

   NCollWeight = 1;
   ZWeight = 1;
   VZWeight = 1;
   for(int i = 0; i < 12; i++)
      ExtraZWeight[i] = 1;

   NVertex = 0;
   VX = 0;
   VY = 0;
   VZ = 0;
   VXError = 0;
   VYError = 0;
   VZError = 0;
   NPU = 0;

   zMass->clear();
   zEta->clear();
   zPhi->clear();
   zY->clear();
   zPt->clear();
   genZMass->clear();
   genZEta->clear();
   genZPhi->clear();
   genZY->clear();
   genZPt->clear();
   trackPt->clear();
   trackDeta->clear();
   trackDphi->clear();
   trackDY->clear();
   TrackPDGID->clear();
   trackPhi->clear();
   trackEta->clear();
   trackY->clear();
   trackMuTagged->clear();
   trackMuDR->clear();
   trackWeight->clear();
   trackResidualWeight->clear();
   trackCharge->clear();
   subevent->clear();

   GenTrackDphi->clear();
   GenTrackDeta->clear();
   GenTrackDY->clear();
   GenTrackPt->clear();
   GenTrackMuTagged->clear();
   GenTrackMuDR->clear();
   GenSubevent->clear();
   GenTrackEta->clear();
   GenTrackPhi->clear();
   GenTrackY->clear();
   GenTrackCharge->clear();

   jetPt->clear();
   jetDeta->clear();
   jetDphi->clear();
   jetRefPt->clear();
   jetRefDeta->clear();
   jetRefDphi->clear();
   jetMuTagged->clear();
   
   genJetPt->clear();
   genJetEta->clear();
   genJetPhi->clear();

   maxOppositeDEta = 0;
   maxOppositeDPhi = 0;
   maxDEta = 0;
   maxDPhi = 0;
   maxOppositeWTADEta = 0;
   maxOppositeWTADPhi = 0;
   maxMoreOppositeWTADEta = 0;
   maxMoreOppositeWTADPhi = 0;
   maxOppositeCSWTADEta = 0;
   maxOppositeCSWTADPhi = 0;
   maxOppositeChargedWTADEta = 0;
   maxOppositeChargedWTADPhi = 0;
   maxMoreOppositeChargedWTADEta = 0;
   maxMoreOppositeChargedWTADPhi = 0;
   maxOppositeHardChargedWTADEta = 0;
   maxOppositeHardChargedWTADPhi = 0;
   
   maxOppositeJet12Pt = 0;
   maxOppositeJet12DEta = 0;
   maxOppositeJet12DPhi = 0;
   maxOppositeJet34Pt = 0;
   maxOppositeJet34DEta = 0;
   maxOppositeJet34DPhi = 0;
   maxOppositeJet56Pt = 0;
   maxOppositeJet56DEta = 0;
   maxOppositeJet56DPhi = 0;
   maxOppositeJet78Pt = 0;
   maxOppositeJet78DEta = 0;
   maxOppositeJet78DPhi = 0;
   
   muEta1->clear();
   muEta2->clear();
   muPhi1->clear();
   muPhi2->clear();
   muPt1->clear();
   muPt2->clear();

   genMuPt1->clear();
   genMuPt2->clear();
   genMuEta1->clear();
   genMuEta2->clear();
   genMuPhi1->clear();
   genMuPhi2->clear();

   muDeta->clear();
   muDphi->clear();
   muDR->clear();
   muDphiS->clear();
   genMuDeta->clear();
   genMuDphi->clear();
   genMuDR->clear();
   genMuDphiS->clear();

   isMuon->clear();
   isGenMuon->clear();

   bestZidx = -1;
   bestZgenIdx = -1;

}

void ZHadronMessenger::CopyNonTrack(ZHadronMessenger &M)
{
   Run          = M.Run;
   Event        = M.Event;
   Lumi         = M.Lumi;
   
   hiBin        = M.hiBin;
   hiBinUp      = M.hiBinUp;
   hiBinDown    = M.hiBinDown;
   hiHF         = M.hiHF;

   SignalHF     = M.SignalHF;
   BackgroundHF = M.BackgroundHF;
   SignalGenHF     = M.SignalGenHF;
   BackgroundGenHF = M.BackgroundGenHF;
   SignalVZ     = M.SignalVZ;

   NCollWeight  = M.NCollWeight;
   ZWeight      = M.ZWeight;
   VZWeight     = M.VZWeight;
   for(int i = 0; i < 12; i++)
      ExtraZWeight[i] = M.ExtraZWeight[i];

   NVertex      = M.NVertex;
   VX           = M.VX;
   VY           = M.VY;
   VZ           = M.VZ;
   VXError      = M.VXError;
   VYError      = M.VYError;
   VZError      = M.VZError;
   NPU          = M.NPU;
  
   if(zMass != nullptr && M.zMass != nullptr)   *zMass = *(M.zMass);
   if(zEta != nullptr && M.zEta != nullptr)   *zEta = *(M.zEta);
   if(zPhi != nullptr && M.zPhi != nullptr)   *zPhi = *(M.zPhi);
   if(zY != nullptr && M.zY != nullptr)   *zY = *(M.zY);
   if(zPt != nullptr && M.zPt != nullptr)   *zPt = *(M.zPt);
   if(genZMass != nullptr && M.genZMass != nullptr)   *genZMass = *(M.genZMass);
   if(genZEta != nullptr && M.genZEta != nullptr)   *genZEta = *(M.genZEta);
   if(genZPhi != nullptr && M.genZPhi != nullptr)   *genZPhi = *(M.genZPhi);
   if(genZY != nullptr && M.genZY != nullptr)   *genZY = *(M.genZY);
   if(genZPt != nullptr && M.genZPt != nullptr)   *genZPt = *(M.genZPt);
   
   if(jetPt != nullptr && M.jetPt != nullptr)   *jetPt = *(M.jetPt);
   if(jetDeta != nullptr && M.jetDeta != nullptr)   *jetDeta = *(M.jetDeta);
   if(jetDphi != nullptr && M.jetDphi != nullptr)   *jetDphi = *(M.jetDphi);
   if(jetRefPt != nullptr && M.jetRefPt != nullptr)   *jetRefPt = *(M.jetRefPt);
   if(jetRefDeta != nullptr && M.jetRefDeta != nullptr)   *jetRefDeta = *(M.jetRefDeta);
   if(jetRefDphi != nullptr && M.jetRefDphi != nullptr)   *jetRefDphi = *(M.jetRefDphi);
   if(jetMuTagged != nullptr && M.jetMuTagged != nullptr)   *jetMuTagged = *(M.jetMuTagged);
   
   if(genJetPt != nullptr && M.genJetPt != nullptr)   *genJetPt = *(M.genJetPt);
   if(genJetEta != nullptr && M.genJetEta != nullptr)   *genJetEta = *(M.genJetEta);
   if(genJetPhi != nullptr && M.genJetPhi != nullptr)   *genJetPhi = *(M.genJetPhi);

   maxOppositeDEta = M.maxOppositeDEta;
   maxOppositeDPhi = M.maxOppositeDPhi;
   maxDEta = M.maxDEta;
   maxDPhi = M.maxDPhi;
   maxOppositeWTADEta = M.maxOppositeWTADEta;
   maxOppositeWTADPhi = M.maxOppositeWTADPhi;
   maxMoreOppositeWTADEta = M.maxMoreOppositeWTADEta;
   maxMoreOppositeWTADPhi = M.maxMoreOppositeWTADPhi;
   maxOppositeCSWTADEta = M.maxOppositeCSWTADEta;
   maxOppositeCSWTADPhi = M.maxOppositeCSWTADPhi;
   maxOppositeChargedWTADEta = M.maxOppositeChargedWTADEta;
   maxOppositeChargedWTADPhi = M.maxOppositeChargedWTADPhi;
   maxMoreOppositeChargedWTADEta = M.maxMoreOppositeChargedWTADEta;
   maxMoreOppositeChargedWTADPhi = M.maxMoreOppositeChargedWTADPhi;
   maxOppositeHardChargedWTADEta = M.maxOppositeHardChargedWTADEta;
   maxOppositeHardChargedWTADPhi = M.maxOppositeHardChargedWTADPhi;

   maxOppositeJet12Pt = M.maxOppositeJet12Pt;
   maxOppositeJet12DEta = M.maxOppositeJet12DEta;
   maxOppositeJet12DPhi = M.maxOppositeJet12DPhi;
   maxOppositeJet34Pt = M.maxOppositeJet34Pt;
   maxOppositeJet34DEta = M.maxOppositeJet34DEta;
   maxOppositeJet34DPhi = M.maxOppositeJet34DPhi;
   maxOppositeJet56Pt = M.maxOppositeJet56Pt;
   maxOppositeJet56DEta = M.maxOppositeJet56DEta;
   maxOppositeJet56DPhi = M.maxOppositeJet56DPhi;
   maxOppositeJet78Pt = M.maxOppositeJet78Pt;
   maxOppositeJet78DEta = M.maxOppositeJet78DEta;
   maxOppositeJet78DPhi = M.maxOppositeJet78DPhi;

   if(muEta1 != nullptr && M.muEta1 != nullptr)   *muEta1 = *(M.muEta1);
   if(muEta2 != nullptr && M.muEta2 != nullptr)   *muEta2 = *(M.muEta2);
   if(muPhi1 != nullptr && M.muPhi1 != nullptr)   *muPhi1 = *(M.muPhi1);
   if(muPhi2 != nullptr && M.muPhi2 != nullptr)   *muPhi2 = *(M.muPhi2);
   if(muPt1 != nullptr && M.muPt1 != nullptr)   *muPt1 = *(M.muPt1);
   if(muPt2 != nullptr && M.muPt2 != nullptr)   *muPt2 = *(M.muPt2);

   if(muDeta != nullptr && M.muDeta != nullptr)   *muDeta = *(M.muDeta);
   if(muDphi != nullptr && M.muDphi != nullptr)   *muDphi = *(M.muDphi);
   if(muDR != nullptr && M.muDR != nullptr)   *muDR = *(M.muDR);
   if(muDphiS != nullptr && M.muDphiS != nullptr)   *muDphiS = *(M.muDphiS);

   if(genMuPt1 != nullptr && M.genMuPt1 != nullptr)   *genMuPt1 = *(M.genMuPt1);
   if(genMuPt2 != nullptr && M.genMuPt2 != nullptr)   *genMuPt2 = *(M.genMuPt2);
   if(genMuEta1 != nullptr && M.genMuEta1 != nullptr)   *genMuEta1 = *(M.genMuEta1);
   if(genMuEta2 != nullptr && M.genMuEta2 != nullptr)   *genMuEta2 = *(M.genMuEta2);
   if(genMuPhi1 != nullptr && M.genMuPhi1 != nullptr)   *genMuPhi1 = *(M.genMuPhi1);
   if(genMuPhi2 != nullptr && M.genMuPhi2 != nullptr)   *genMuPhi2 = *(M.genMuPhi2);

   if(genMuDeta != nullptr && M.genMuDeta != nullptr)   *genMuDeta = *(M.genMuDeta);
   if(genMuDphi != nullptr && M.genMuDphi != nullptr)   *genMuDphi = *(M.genMuDphi);
   if(genMuDR != nullptr && M.genMuDR != nullptr)   *genMuDR = *(M.genMuDR);
   if(genMuDphiS != nullptr && M.genMuDphiS != nullptr)   *genMuDphiS = *(M.genMuDphiS);

   if(isMuon != nullptr && M.isMuon != nullptr)   *isMuon = *(M.isMuon);
   if(isGenMuon != nullptr && M.isGenMuon != nullptr)   *isGenMuon = *(M.isGenMuon);

   bestZidx     = M.bestZidx;
   bestZgenIdx  = M.bestZgenIdx;
   
}

bool ZHadronMessenger::FillEntry()
{
   if(Initialized == false)
      return false;
   if(WriteMode == false)
      return false;

   if(Tree == nullptr)
      return false;

   Tree->Fill();
   Clear();

   return true;
}
