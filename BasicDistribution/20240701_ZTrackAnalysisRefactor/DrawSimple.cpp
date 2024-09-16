#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TMath.h>
#include <TROOT.h>
#include <cmath>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TCut.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TLine.h>
#include <iostream>

using namespace std;

#include "CommandLine.h"

int main(int argc, char *argv[]);

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
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.2);
  gStyle->SetPadRightMargin(0.05);

  gStyle->SetLineScalePS(1.5);

  gROOT->ForceStyle();
}

const char *typeofdata = "ZHadron2024/SkimBkgSub/ov1_v4d_sub0_GenMatch/20240907/";
const char *typeofdata1 = "ov1_v4d_sub0_GenMatch";

const char *typeofdataRres = "ZHadron2024/SkimBkgSub/ov1_v4d_sub0_Rres_v3/20240912/";
const char *typeofdataRres1 = "ov1_v4d_sub0_Rres_v3";

const char *typeofdataNoZ = "ZHadron2024/SkimBkgSub/ov1_v4d_sub0_NoZ_GenMatch/20240907/";
const char *typeofdataNoZ1 = "ov1_v4d_sub0_NoZ_GenMatch";


TChain *TreeSig = new TChain("Tree"); 
TChain *TreeBkg = new TChain("Tree"); 
//TChain *TreePP0 = new TChain("Tree"); 
TChain *TreeSgG = new TChain("Tree"); 
TChain *TreeBgG = new TChain("Tree"); 

void DefinePhiRangeCorrelation() {
    gROOT->ProcessLine(
        "double PhiRangeCorrelation(double Phi) {"
        "   if(Phi < -TMath::Pi() * 0.5) Phi = Phi + 2 * TMath::Pi();"
        "   if(Phi > +TMath::Pi() * 1.5) Phi = Phi - 2 * TMath::Pi();"
        "   return Phi;"
        "}"
    );
}

void DrawSimple_single(float ptL=20,float ptH=2000,float centL=0,float centH=90,float TptL=0,float TptH=10000, bool isRres = false, bool noZ = false, double leptonveto = 0.01){

	DefinePhiRangeCorrelation();

	string HistName = "trackDphi";

	TCanvas *c = new TCanvas("c", "", 500, 600);
  TPad *Pad  = new TPad("Pad" , "", 0, 0.3, 1,   1);
  TPad *RPad = new TPad("RPad", "", 0, 0. , 1, 0.3);

  Pad->Draw();
  RPad->Draw();

	TH1D* d10 = new TH1D("d10","",20,-0.5*M_PI,1.5*M_PI);
	TH1D* d20 = new TH1D("d20","",20,-0.5*M_PI,1.5*M_PI);
	TH1D* d3 =  new TH1D( "d3","",20,-0.5*M_PI,1.5*M_PI);
	TH1D* d4 =  new TH1D( "d4","",20,-0.5*M_PI,1.5*M_PI);
	TH1D* d5 =  new TH1D( "d5","",20,-0.5*M_PI,1.5*M_PI);

	//TH1D* d10 = new TH1D("d10","",20,TptL,TptH);
	//TH1D* d20 = new TH1D("d20","",20,TptL,TptH);
	//TH1D* d3  =  new TH1D("d3","",20,TptL,TptH);
	//TH1D* d4  =  new TH1D("d4","",20,TptL,TptH);
	//TH1D* d5  =  new TH1D("d5","",20,TptL,TptH);

	
/*
	TChain *TreeSig = new TChain("Tree"); TreeSig->Add("~/PhysicsHIZHadron2022/Skims/20230310_ZHadronSkims/TestMC_v14.root");
	TChain *TreePP0 = new TChain("Tree"); TreePP0->Add("~/PhysicsHIZHadron2022/Skims/20230310_ZHadronSkims/TestPP_v14.root");
	TChain *TreeSgG = new TChain("Tree"); TreeSgG->Add("~/PhysicsHIZHadron2022/Skims/20230310_ZHadronSkims/TestGenMC_v14.root");
	*/

	//TCut evtCut = "zMass>60&&zPt>5";//

	TCut evtCutPP0 = Form("bestZidx==0&&zMass[0]>60&&zPt[0]>%f&&zPt[0]<%f",ptL,ptH);
	TCut evtCutPP1 = Form("bestZidx==1&&zMass[1]>60&&zPt[1]>%f&&zPt[1]<%f",ptL,ptH);
	TCut evtCutPP2 = Form("bestZidx==2&&zMass[2]>60&&zPt[2]>%f&&zPt[2]<%f",ptL,ptH);
	TCut evtCutPP3 = Form("bestZidx==3&&zMass[3]>60&&zPt[3]>%f&&zPt[3]<%f",ptL,ptH);

	TCut evtCutPP = evtCutPP0;// || evtCutPP1 || evtCutPP2 || evtCutPP3;
	TCut evtCut = evtCutPP && Form("hiBin>%f&&hiBin<%f",2*centL,2*centH);

	TCut evtCutGenPP0 = Form("bestZgenIdx==0&&genZMass[0]>60&&genZPt[0]>%f&&genZPt[0]<%f",ptL,ptH);
	TCut evtCutGenPP1 = Form("bestZgenIdx==1&&genZMass[1]>60&&genZPt[1]>%f&&genZPt[1]<%f",ptL,ptH);
	TCut evtCutGenPP2 = Form("bestZgenIdx==2&&genZMass[2]>60&&genZPt[2]>%f&&genZPt[2]<%f",ptL,ptH);
	TCut evtCutGenPP3 = Form("bestZgenIdx==3&&genZMass[3]>60&&genZPt[3]>%f&&genZPt[3]<%f",ptL,ptH);

	TCut evtCutGenPP = evtCutGenPP0;// || evtCutGenPP1 || evtCutGenPP2 || evtCutGenPP3;
	TCut evtCutGen = evtCutGenPP && Form("hiBin>%f&&hiBin<%f",2*centL,2*centH);

	TCut trkCut = Form("trackMuTagged==0&&trackPt>%f&&trackPt<%f&&trackMuDR>%f",TptL,TptH,leptonveto);
	TCut trkCutSum = Form("Sum$(trackMuTagged==0&&trackPt>%f&&trackPt<%f&&trackMuDR>%f)>0",TptL,TptH,leptonveto);


/*
	TreeSig->SetAlias("TrackEta", "(trackDeta+zEta[0])");
	TreeBkg->SetAlias("TrackEta", "(trackDeta+zEta[0])");
	TreePP0->SetAlias("TrackEta", "(trackDeta+zEta[0])");
	TreeSgG->SetAlias("TrackEta", "(trackDeta+genZEta[0])");
	TreeBgG->SetAlias("TrackEta", "(trackDeta+genZEta[0])");

	TreeSig->SetAlias("RawTrackPhi", "(trackDphi+zPhi[0])");
	TreeSig->SetAlias("TrackPhi", "(RawTrackPhi+2*3.14159*(RawTrackPhi<-3.14159)-2*3.14159*(RawTrackPhi>3.14159))");
	TreeBkg->SetAlias("RawTrackPhi", "(trackDphi+zPhi[0])");
	TreeBkg->SetAlias("TrackPhi", "(RawTrackPhi+2*3.14159*(RawTrackPhi<-3.14159)-2*3.14159*(RawTrackPhi>3.14159))");
	TreePP0->SetAlias("RawTrackPhi", "(trackDphi+zPhi[0])");
	TreePP0->SetAlias("TrackPhi", "(RawTrackPhi+2*3.14159*(RawTrackPhi<-3.14159)-2*3.14159*(RawTrackPhi>3.14159))");
	TreeSgG->SetAlias("RawTrackPhi", "(trackDphi+genZPhi[0])");
	TreeSgG->SetAlias("TrackPhi", "(RawTrackPhi+2*3.14159*(RawTrackPhi<-3.14159)-2*3.14159*(RawTrackPhi>3.14159))");
	TreeBgG->SetAlias("RawTrackPhi", "(trackDphi+genZPhi[0])");
	TreeBgG->SetAlias("TrackPhi", "(RawTrackPhi+2*3.14159*(RawTrackPhi<-3.14159)-2*3.14159*(RawTrackPhi>3.14159))");
	
	TreeSig->SetAlias("ZPhi", "(zPhi+2*3.14159*(zPhi<-3.14159)-2*3.14159*(zPhi>3.14159))");
	TreeBkg->SetAlias("ZPhi", "(zPhi+2*3.14159*(zPhi<-3.14159)-2*3.14159*(zPhi>3.14159))");
	TreePP0->SetAlias("ZPhi", "(zPhi+2*3.14159*(zPhi<-3.14159)-2*3.14159*(zPhi>3.14159))");
	TreeSgG->SetAlias("ZPhi", "(genZPhi+2*3.14159*(genZPhi<-3.14159)-2*3.14159*(genZPhi>3.14159))");
	TreeBgG->SetAlias("ZPhi", "(genZPhi+2*3.14159*(genZPhi<-3.14159)-2*3.14159*(genZPhi>3.14159))");

*/

	TreeSig->Draw("PhiRangeCorrelation(trackDphi)>>+d10","(0.5*NCollWeight*ZWeight*VZWeight*trackWeight*trackResidualWeight)"*(evtCut&&trkCut),"goff");
	TreeBkg->Draw("PhiRangeCorrelation(trackDphi)>>+d20","(0.5*NCollWeight*ZWeight*VZWeight*trackWeight*trackResidualWeight)"*(evtCut&&trkCut),"goff");
	TreeSgG->Draw("PhiRangeCorrelation(trackDphi)>>+d4","(0.5*NCollWeight*ZWeight*VZWeight*trackWeight*trackResidualWeight)"*(evtCutGen&&trkCut),"goff");
	TreeBgG->Draw("PhiRangeCorrelation(trackDphi)>>+d5","(0.5*NCollWeight*ZWeight*VZWeight*trackWeight*trackResidualWeight)"*(evtCutGen&&trkCut),"goff");

//TreePP0->Draw("PhiRangeCorrelation(trackDphi)>>+d3","(0.5*NCollWeight*ZWeight*VZWeight*trackWeight*trackResidualWeight)"*(evtCutPP&&trkCut&&"NPU==0"),"goff");
	TreeSgG->Draw("PhiRangeCorrelation(trackDphi)>>+d3","(0.5*NCollWeight*ZWeight*VZWeight*trackWeight*trackResidualWeight)"*(evtCutGen&&trkCut&&"subevent==0"),"goff");

	TreeSig->Draw("PhiRangeCorrelation(-trackDphi)>>+d10","(0.5*NCollWeight*ZWeight*VZWeight*trackWeight*trackResidualWeight)"*(evtCut&&trkCut),"goff");
	TreeBkg->Draw("PhiRangeCorrelation(-trackDphi)>>+d20","(0.5*NCollWeight*ZWeight*VZWeight*trackWeight*trackResidualWeight)"*(evtCut&&trkCut),"goff");
	TreeSgG->Draw("PhiRangeCorrelation(-trackDphi)>>+d4","(0.5*NCollWeight*ZWeight*VZWeight*trackWeight*trackResidualWeight)"*(evtCutGen&&trkCut),"goff");
	TreeBgG->Draw("PhiRangeCorrelation(-trackDphi)>>+d5","(0.5*NCollWeight*ZWeight*VZWeight*trackWeight*trackResidualWeight)"*(evtCutGen&&trkCut),"goff");

//TreePP0->Draw("PhiRangeCorrelation(-trackDphi)>>+d3","(0.5*NCollWeight*ZWeight*VZWeight*trackWeight*trackResidualWeight)"*(evtCutPP&&trkCut&&"NPU==0"),"goff");
	TreeSgG->Draw("PhiRangeCorrelation(-trackDphi)>>+d3","(0.5*NCollWeight*ZWeight*VZWeight*trackWeight*trackResidualWeight)"*(evtCutGen&&trkCut&&"subevent==0"),"goff");


	//TreeSig->Draw("trackPt>>d10","(NCollWeight*trackWeight*trackResidualWeight)"*(evtCut&&trkCut),"goff");
	////TreeSig->Draw("trackPt>>d2","(NCollWeight*trackWeight*trackResidualWeight)"*(evtCut&&trkCut&&"hiBin<100"),"goff");
	//TreeBkg->Draw("trackPt>>d20","(NCollWeight*trackWeight*trackResidualWeight)"*(evtCut&&trkCut),"goff");
	//TreeSgG->Draw("trackPt>>d4","(NCollWeight*trackWeight*trackResidualWeight)"*(evtCutGen&&trkCut),"goff");
	//TreePP0->Draw("trackPt>>d3","(NCollWeight*trackWeight*trackResidualWeight)"*(evtCutPP&&trkCut&&"NPU==0"),"goff");
	//TreeBgG->Draw("trackPt>>d5","(NCollWeight*trackWeight*trackResidualWeight)"*(evtCutGen&&trkCut),"goff");


	//TreeSig->SetScanField(0);
	//TreeSig->Scan("trackPt:TrackEta:TrackPhi:NCollWeight:trackWeight:trackResidualWeight:ZWeight:zPt[0]:zEta[0]:zMass[0]",evtCut&&trkCut);


	//TreeSig->Draw("zPt>>d1","NCollWeight"*(evtCut&&"hiBin>100"),"goff");
	//TreeSig->Draw("zPt>>d2","NCollWeight"*(evtCut&&"hiBin<100"),"goff");
	////TreeBkg->Draw("zPt>>d2","NCollWeight"*evtCut,"goff");
	//TreePP0->Draw("zPt>>d3","NCollWeight"*(evtCut&&"subevent==0"),"goff");
	
	TH1D HNSig("HNSig","Normalization", 1, 0, 1);
	TH1D HNBkg("HNBkg","Normalization", 1, 0, 1);
	TH1D HNPP0("HNPP0","Normalization", 1, 0, 1);
	TH1D HNSgG("HNSgG","Normalization", 1, 0, 1);
	TH1D HNBgG("HNBgG","Normalization", 1, 0, 1);

	//TreeSig->Draw("0>>HNSig", "NCollWeight"*evtCut,"goff");
	////TreeSig->Draw("0>>HNBkg", "NCollWeight"*(evtCut&&"hiBin<100"),"goff");
	//TreeBkg->Draw("0>>HNBkg", "NCollWeight"*evtCut,"goff");
	//TreeSgG->Draw("0>>HNSgG", "NCollWeight"*evtCutGen,"goff");
	//TreePP0->Draw("0>>HNPP0", "NCollWeight"*(evtCutPP&&"NPU==0"),"goff");
	//TreeBgG->Draw("0>>HNBgG", "NCollWeight"*evtCutGen,"goff");

	TreeSig->Draw("0>>HNSig", "(NCollWeight*ZWeight*VZWeight)"*(evtCut&&trkCutSum),"goff");
	TreeBkg->Draw("0>>HNBkg", "(NCollWeight*ZWeight*VZWeight)"*(evtCut&&trkCutSum),"goff");
	TreeSgG->Draw("0>>HNSgG", "(NCollWeight*ZWeight*VZWeight)"*(evtCutGen&&trkCutSum),"goff");
	TreeBgG->Draw("0>>HNBgG", "(NCollWeight*ZWeight*VZWeight)"*(evtCutGen&&trkCutSum),"goff");

	//TreePP0->Draw("0>>HNPP0", "(NCollWeight*ZWeight*VZWeight)"*(evtCutPP&&"NPU==0"),"goff");
	TreeSgG->Draw("0>>HNPP0", "(NCollWeight*ZWeight*VZWeight)"*(evtCutGen&&trkCutSum&&"Sum$(subevent==0)>0"),"goff");

	double t1N = HNSig.GetBinContent(1);
	double t2N = HNBkg.GetBinContent(1);
	double t3N = HNPP0.GetBinContent(1);
	double t4N = HNSgG.GetBinContent(1);
	double t5N = HNBgG.GetBinContent(1);

	double t1NError = HNSig.GetBinError(1);
	double t2NError = HNBkg.GetBinError(1);
	double t3NError = HNPP0.GetBinError(1);
	double t4NError = HNSgG.GetBinError(1);
	double t5NError = HNBgG.GetBinError(1);

	std::cout<<"t1N = "<<t1N<<", t2N = "<<t2N<<", t3N = "<<t3N<<", t4N = "<<t4N<<", t5N = "<<t5N<<std::endl;

	std::cout<<"=============================="<<std::endl;

	std::cout<<"Number of Weighted Reco Z events in sig = "<<t1N<<"+-"<<t1NError<<std::endl;
	std::cout<<"Number of Weighted Reco Z events in bkg = "<<t2N<<"+-"<<t2NError<<std::endl;
	std::cout<<"Sig-Bkg = "<<t1N-t2N<<std::endl;
	std::cout<<"Number of Weighted Z events in pp0 = "<<t3N<<"+-"<<t3NError<<std::endl;
	std::cout<<"Number of Weighted Gen Z events in sig = "<<t4N<<"+-"<<t4NError<<std::endl;
	std::cout<<"Number of Weighted Gen Z events in bkg = "<<t5N<<"+-"<<t5NError<<std::endl;
	std::cout<<"Sig-Bkg = "<<t4N-t5N<<std::endl;

	std::cout<<"=============================="<<std::endl;

	std::cout<<"Number of Unweighted Reco Z events in sig = "<<HNSig.GetEntries()<<std::endl;
	std::cout<<"Number of Unweighted Reco Z events in bkg = "<<HNBkg.GetEntries()<<std::endl;
	std::cout<<"Sig-Bkg = "<<HNSig.GetEntries()-HNBkg.GetEntries()<<std::endl;
	std::cout<<"Number of Unweighted Z events in pp0 = "<<HNPP0.GetEntries()<<std::endl;
	std::cout<<"Number of Unweighted Gen Z events in sig = "<<HNSgG.GetEntries()<<std::endl;
	std::cout<<"Number of Unweighted Gen Z events in bkg = "<<HNBgG.GetEntries()<<std::endl;
	std::cout<<"Sig-Bkg = "<<HNSgG.GetEntries()-HNBgG.GetEntries()<<std::endl;

	std::cout<<"=============================="<<std::endl;

	std::cout<<"Number of Weighted Reco tracks in sig = "<<d10->Integral()<<std::endl;
	std::cout<<"Number of Weighted Reco tracks in bkg = "<<d20->Integral()<<std::endl;
	std::cout<<"Sig-Bkg = "<<d10->Integral()-d20->Integral()<<std::endl;
	std::cout<<"Number of Weighted tracks in pp0 = "<<d3->Integral()<<std::endl;
	std::cout<<"Number of Weighted Gen tracks in sig = "<<d4->Integral()<<std::endl;
	std::cout<<"Number of Weighted Gen tracks in bkg = "<<d5->Integral()<<std::endl;
	std::cout<<"Sig-Bkg = "<<d4->Integral()-d5->Integral()<<std::endl;

	std::cout<<"=============================="<<std::endl;

	std::cout<<"Number of Unweighted Reco tracks in sig = "<<d10->GetEntries()<<std::endl;
	std::cout<<"Number of Unweighted Reco tracks in bkg = "<<d20->GetEntries()<<std::endl;
	std::cout<<"Sig-Bkg = "<<d10->GetEntries()-d20->GetEntries()<<std::endl;
	std::cout<<"Number of Unweighted tracks in pp0 = "<<d3->GetEntries()<<std::endl;
	std::cout<<"Number of Unweighted Gen tracks in sig = "<<d4->GetEntries()<<std::endl;
	std::cout<<"Number of Unweighted Gen tracks in bkg = "<<d5->GetEntries()<<std::endl;
	std::cout<<"Sig-Bkg = "<<d4->GetEntries()-d5->GetEntries()<<std::endl;

	std::cout<<"=============================="<<std::endl;

	
	d10->Scale(1./t1N, "width");//sig reco
	d20->Scale(1./t2N, "width");//bkg reco
	d3 ->Scale(1./t3N, "width");//pp reco or sig gen sube=0
	d4 ->Scale(1./t4N, "width");//sig gen
	d5 ->Scale(1./t5N, "width");//bkg gen

/*
	auto updateErrors = [](TH1D* hist, double N, double NError) {
	  double relativeError = NError / N;
	  for (int i = 1; i <= hist->GetNbinsX(); i++) {
	      double content = hist->GetBinContent(i);
	      double error = hist->GetBinError(i);
	      double newRelativeError = sqrt(pow(error/content, 2) + pow(relativeError, 2));
	      hist->SetBinError(i, content * newRelativeError);
	  }
	};
	
	updateErrors(d10, t1N, t1NError);
	updateErrors(d20, t2N, t2NError);
	updateErrors(d3, t3N, t3NError);
	updateErrors(d4, t4N, t4NError);
	updateErrors(d5, t5N, t5NError);
*/
	TH1D *d1 = (TH1D*) d10->Clone("d1");//sig-bkg reco
	TH1D *d2 = (TH1D*) d4->Clone("d2");//sig-bkg gen

	d1->Add(d20,-1);
	d2->Add(d5,-1);

	d10->SetLineColor(kBlack);//sig reco
  d20->SetLineColor(kBlue);//bkg reco
  d1->SetLineColor(kRed);//sig-bkg reco
  d3->SetLineColor(kBlack);//pp reco or sig gen sube=0

  d10->SetMarkerStyle(kFullCircle);
  d20->SetMarkerStyle(kFullCircle);
  d1->SetMarkerStyle(kFullCircle);

  d10->SetMarkerColor(kBlack);
  d20->SetMarkerColor(kBlue);
  d1->SetMarkerColor(kRed);

  d10->SetMarkerSize(0.5);
  d20->SetMarkerSize(0.5);
  d1->SetMarkerSize(0.5);

  d4->SetLineColor(kGray+2);//grey //sig gen
	d5->SetLineColor(TColor::GetColor("#377eb8"));//blue //bkg gen
	d2->SetLineColor(kViolet);//sig-bkg gen

	d4->SetFillStyle(3345);
	d5->SetFillStyle(3354);
	d2->SetFillStyle(1001);//3395

	d4->SetFillColor(kGray+2);
	d5->SetFillColor(TColor::GetColor("#377eb8"));
	//d2->SetFillColor(kViolet);
	d2->SetFillColorAlpha(kViolet, 0.35);


	Double_t err1, err2, err3, err4, err5, err10, err20;

	t1N = d1->IntegralAndError(1,d1->GetNbinsX(),err1,"width");
	t2N = d2->IntegralAndError(1,d2->GetNbinsX(),err2,"width");
	t3N = d3->IntegralAndError(1,d3->GetNbinsX(),err3,"width");
	t4N = d4->IntegralAndError(1,d4->GetNbinsX(),err4,"width");
	t5N = d5->IntegralAndError(1,d5->GetNbinsX(),err5,"width");
	double t10N = d10->IntegralAndError(1,d10->GetNbinsX(),err10,"width");
	double t20N = d20->IntegralAndError(1,d20->GetNbinsX(),err20,"width");

	
/*
	if(lowpt==50){
		d1->Rebin(4); d2->Rebin(4); d3->Rebin(4);
	}else if(lowpt==10){
		d1->Rebin(2); d2->Rebin(2); d3->Rebin(2);
	}*/


	// == Start drawing == //

  

   double max1 = d10->GetMaximum();
   double max2 = d20->GetMaximum();
   double max3 = d3->GetMaximum();

   double min1 = d1->GetMinimum();
   double min2 = d2->GetMinimum();
   double min3 = d3->GetMinimum();

   if(min1>min2) min1=min2;

   if(min1>min3) min1=min3;

   Pad->cd();

   if(max1<max2) d20->Draw("ep");
   else d10->Draw("ep");
   d10->Draw("ep same");
   d20->Draw("ep same");

   d1->Draw("ep same");
   d3->Draw("hist same");

   d4->Draw("hist same");
   d5->Draw("hist same");

   if(max1<max2) max1=max2;

   d10->SetYTitle("dN/d#Delta#phi");
   d20->SetYTitle("dN/d#Delta#phi");

   d10->SetXTitle("#Delta#phi_{Z,track}");
   d20->SetXTitle("#Delta#phi_{Z,track}");

	TLegend leg1(0.54,0.62,0.98,0.92);
	leg1.AddEntry(d10 ,"Raw MC RECO","lep");
	leg1.AddEntry(d20 ,"Bkg MC RECO","lep");
	leg1.AddEntry(d1 ,"Raw-Bkg MC RECO","lep");
	leg1.AddEntry(d4 ,"Raw MC GEN","lf");
	leg1.AddEntry(d5 ,"Bkg MC GEN","lf");
	leg1.AddEntry(d2 ,"Raw-Bkg MC GEN","lf");
	leg1.AddEntry(d3 ,"Sig GEN (sube=0)","l");
	leg1.SetFillColorAlpha(kWhite,0);
	leg1.SetLineColor(kBlack);
	leg1.SetLineWidth(1);
	leg1.Draw();

	//TLatex *pt1 = new TLatex(0.1,0.97,Form("%d < Track p_{T} < %d",10,20));
	//TLatex *pt1 = new TLatex(0.18,0.82,Form("#splitline{%.1f < Z p_{T} < %.1f}{%.0f < Centrality < %.0f%}{%.1f < Track p_{T} < %.1f}",ptL,ptH,centL,centH,TptL,TptH));
  //pt1->SetTextFont(62);
  //pt1->SetTextSize(0.04);
  //pt1->SetNDC(kTRUE);
  //pt1->Draw();

  TLatex *pt0 = new TLatex(0.25,0.88,"Nominal MC Skim");
   pt0->SetTextFont(42);
   pt0->SetTextSize(0.03);
   pt0->SetNDC(kTRUE);
   pt0->Draw();

  TLatex *pt = new TLatex(0.25,0.82,Form("%.0f %%< Centrality < %.0f %%",centL,centH));
  pt->SetTextFont(42);
  pt->SetTextSize(0.03);
  pt->SetNDC(kTRUE);
  pt->Draw();

  TLatex *pt1 = new TLatex(0.25,0.76,Form("%.1f < Z p_{T} < %.1f GeV",ptL,ptH));
  pt1->SetTextFont(42);
  pt1->SetTextSize(0.03);
  pt1->SetNDC(kTRUE);
  pt1->Draw();

  TLatex *pt3 = new TLatex(0.25,0.70,Form("%.1f < Track p_{T} < %.1f GeV",TptL,TptH));
  pt3->SetTextFont(42);
  pt3->SetTextSize(0.03);
  pt3->SetNDC(kTRUE);
  pt3->Draw();

  // t1N  *= d1->GetBinWidth(1);
	// t2N  *= d2->GetBinWidth(1);
	// t3N  *= d3->GetBinWidth(1);
	// err1 *= d1->GetBinWidth(1);
	// err2 *= d2->GetBinWidth(1);
	// err3 *= d3->GetBinWidth(1);


	 TLatex *ptInt1 = new TLatex(0.25,0.58,Form("#Sigma raw RECO = %.2f #pm %.2f,  #Sigma bkg RECO = %.2f #pm %.2f", t10N , err10, t20N, err20));
   ptInt1->SetTextFont(42);
   ptInt1->SetTextSize(0.03);
   ptInt1->SetNDC(kTRUE);
   ptInt1->Draw();

   TLatex *ptInt11 = new TLatex(0.25,0.52,Form("#Sigma raw GEN = %.2f #pm %.2f,  #Sigma bkg GEN = %.2f #pm %.2f", t4N ,err4, t5N, err5));
   ptInt11->SetTextFont(42);
   ptInt11->SetTextSize(0.03);
   ptInt11->SetNDC(kTRUE);
   ptInt11->Draw();

   TLatex *ptInt2 = new TLatex(0.25,0.46,Form("#Sigma (raw-bkg) RECO = %.2f #pm %.2f, #Sigma (raw-bkg) GEN = %.2f #pm %.2f",t1N, err1, t2N, err2));
   ptInt2->SetTextFont(42);
   ptInt2->SetTextSize(0.03);
   ptInt2->SetNDC(kTRUE);
   ptInt2->Draw();

   TLatex *ptInt3 = new TLatex(0.25,0.40,Form("#Sigma sig GEN (sube=0) = %.2f #pm %.2f",t3N, err3));
   ptInt3->SetTextFont(42);
   ptInt3->SetTextSize(0.03);
   ptInt3->SetNDC(kTRUE);
   ptInt3->Draw();



  //TLatex *pt2 = new TLatex(0.1,0.97,Form("INT_{reco} = %.3f #pm %.3f, INT_{gen} = %.3f #pm %.3f, INT_{pp} = %.3f #pm %.3f",t1N,err1,t2N,err2,t3N,err3));
  //pt2->SetTextFont(42);
  //pt2->SetTextSize(0.03);
  //pt2->SetNDC(kTRUE);
  //pt2->Draw();

  d10->SetMaximum(2.5*max1);
  d20->SetMaximum(2.5*max1);

  if(min1<0){
  		d10->SetMinimum(2*min1);
			d20->SetMinimum(2*min1);
   }else{
      d10->SetMinimum(0);
			d20->SetMinimum(0);
   }

	

	RPad->cd();


	TH1D *PbPb_to_pp = (TH1D*) d1->Clone("PbPb_to_pp");
  PbPb_to_pp->Add(d3,-1);

  TH1D *Gen_to_pp = (TH1D*) d2->Clone("Gen_to_pp");
  Gen_to_pp->Add(d3,-1);

  TH1D *horiz_line = (TH1D*) d3->Clone("horiz_line");
  horiz_line->Add(d3,-1);

  horiz_line->SetLineColor(kBlack);
  PbPb_to_pp->SetLineColor(kRed);
  Gen_to_pp->SetLineColor(kViolet);
  Gen_to_pp->SetMarkerStyle(kFullSquare);
  Gen_to_pp->SetMarkerColor(kViolet);
  Gen_to_pp->SetMarkerSize(0.5);

  //Gen_to_pp->SetFillStyle(1001);
  //Gen_to_pp->SetFillColorAlpha(kViolet, 0.35);


  PbPb_to_pp->SetXTitle("");
  PbPb_to_pp->SetYTitle("Difference");

  PbPb_to_pp->SetMaximum(2);
  PbPb_to_pp->SetMinimum(-2);


  PbPb_to_pp->Draw("ep");
  Gen_to_pp->Draw("hist ep same");
  horiz_line->Draw("hist same");

  if(isRres)
  	gSystem->Exec(Form("mkdir -p /eos/user/p/pchou/figs/%s/%s",typeofdataRres,HistName.c_str()));
  else if(noZ)
  	gSystem->Exec(Form("mkdir -p /eos/user/p/pchou/figs/%s/%s",typeofdataNoZ,HistName.c_str()));
  else
  	gSystem->Exec(Form("mkdir -p /eos/user/p/pchou/figs/%s/%s",typeofdata,HistName.c_str()));

  if(isRres)
  	c->SaveAs(Form("/eos/user/p/pchou/figs/%s/%s/Ztrack_%s_com_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f.png",typeofdataRres,HistName.c_str(),typeofdataRres1,ptL,ptH,centL,centH,TptL,TptH,100*leptonveto));
  else if(noZ)
  	c->SaveAs(Form("/eos/user/p/pchou/figs/%s/%s/Ztrack_%s_com_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f.png",typeofdataNoZ,HistName.c_str(),typeofdataNoZ1,ptL,ptH,centL,centH,TptL,TptH,100*leptonveto)); 
  else
  	c->SaveAs(Form("/eos/user/p/pchou/figs/%s/%s/Ztrack_%s_com_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f.png",typeofdata,HistName.c_str(),typeofdata1,ptL,ptH,centL,centH,TptL,TptH,100*leptonveto)); 
   
  d4->SetFillStyle(0);
  d5->SetFillStyle(0);
  d10->SetMaximum(12);
  d20->SetMaximum(12);

  if(isRres)
  	c->SaveAs(Form("/eos/user/p/pchou/figs/%s/%s/Ztrack_%s_com_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f_low.png",typeofdataRres,HistName.c_str(),typeofdataRres1,ptL,ptH,centL,centH,TptL,TptH,100*leptonveto));
  else if(noZ)
  	c->SaveAs(Form("/eos/user/p/pchou/figs/%s/%s/Ztrack_%s_com_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f_low.png",typeofdataNoZ,HistName.c_str(),typeofdataNoZ1,ptL,ptH,centL,centH,TptL,TptH,100*leptonveto)); 
  else
  	c->SaveAs(Form("/eos/user/p/pchou/figs/%s/%s/Ztrack_%s_com_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f_low.png",typeofdata,HistName.c_str(),typeofdata1,ptL,ptH,centL,centH,TptL,TptH,100*leptonveto)); 

}

int main(int argc, char *argv[]){

	style();

	CommandLine CL(argc, argv);

	double ptL                   = CL.GetDouble("ptL", 40);
	double ptH                   = CL.GetDouble("ptH", 200);
	double centL                 = CL.GetDouble("centL", 0);
	double centH                 = CL.GetDouble("centH", 90);
	double TptL                  = CL.GetDouble("TptL", 1);
	double TptH                  = CL.GetDouble("TptH", 1000);

	bool isRres					 = CL.GetBool("isRres",false);
	bool noZ					 = CL.GetBool("noZ",false);

	double leptonveto			 = CL.GetDouble("leptonveto", 0.01);

	string filebase = "/eos/cms/store/group/phys_heavyions/pchou/SkimZHadron2024/";
	string cernbox = "/eos/home-p/pchou/SkimZHadron2024/";

	if(noZ)
		TreeSig->Add((filebase + "OutputMC_v4b_ee_NoZ_v3/Result*.root").c_str());
	else
		TreeSig->Add((cernbox + "OutputMC_v4d_ee_v3/Result*.root").c_str());

	if(isRres)
		TreeBkg->Add((cernbox + "OutputMCBkg_v4d_ee_Rres_v3/Result*.root").c_str());
	else if(noZ)
		TreeBkg->Add((cernbox + "OutputMCBkg_v4b_ee_NoZ_tol120_v3/Result*.root").c_str());
	else
		TreeBkg->Add((cernbox + "OutputMCBkg_v4b_ee_tol120/Result*.root").c_str());

	//TreePP0->Add((filebase + "OutputPPMC_v3c_ee/*.root").c_str());
	TreeSgG->Add((cernbox + "OutputMCGen_v4d_ee/Result*.root").c_str());
	TreeBgG->Add((cernbox + "OutputMCbkgGen_v4d_ee/Result*.root").c_str());

	//TreeSig->Add("/eos/cms/store/group/phys_heavyions/pchou/SkimMC_v14.root");
	////TreeBkg->Add("/eos/cms/store/group/phys_heavyions/pchou/SkimMCbkg_v14.root");
	//TreePP0->Add("/eos/cms/store/group/phys_heavyions/pchou/SkimPPMC_v14.root");
	//TreeSgG->Add("/eos/cms/store/group/phys_heavyions/pchou/SkimMCGen_v14.root");
	////TreeBgG->Add("/eos/cms/store/group/phys_heavyions/pchou/SkimMCbkgGen_v14.root");


	DrawSimple_single(ptL,ptH,centL,centH,TptL,TptH,isRres,noZ,leptonveto);


	//DrawSimple_single(40,200,0,10,1,2);

	//DrawSimple_single(20,2000,0,10,10,20);
	//DrawSimple_single(20,2000,0,10,20,50);
	//DrawSimple_single(20,2000,0,10,50,100);
	//DrawSimple_single(20,2000,0,90,10,20);
	//DrawSimple_single(20,2000,0,90,20,50);
	//DrawSimple_single(20,2000,0,90,50,100);
}
