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
#include <TLeaf.h>
#include <TLatex.h>
#include <TLine.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <glob.h>
using namespace std;

#include "CommandLine.h"
#include "ProgressBar.h"

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

const char *typeofdata = "ZHadron2024/DrawNTrk/ov1_v2fv3gen_sub0/20240726/";
const char *typeofdata1 = "ov1_v2fv3gen_sub0";


string filebase = "/eos/cms/store/group/phys_heavyions/pchou/SkimZHadron2024/";

string TreeSigNames = filebase + "OutputMC_v2f_ee/Result*.root";
string TreeBkgNames = filebase + "OutputMCBkg_v2f_ee/Result*.root";
string TreeSgGNames = filebase + "OutputMCGen_v3_ee/Result*.root";
string TreeBgGNames = filebase + "OutputMCbkgGen_v3_ee/Result*.root";

std::vector<std::string> my_glob(const char *pattern) {
    glob_t g;
    glob(pattern, GLOB_TILDE, nullptr, &g); // one should ensure glob returns 0!
    std::vector<std::string> filelist;
    filelist.reserve(g.gl_pathc);
    for (size_t i = 0; i < g.gl_pathc; ++i) {
        filelist.emplace_back(g.gl_pathv[i]);
    }
    globfree(&g);
    return filelist;
}

void FillNTrks(TH1D* d10, TTree *TreeSig, float ptL,float ptH,float centL,float centH,float TptL, float TptH, bool isSube0){

	int hiBin, bestZidx;
	float NCollWeight, ZWeight, VZWeight;

	vector<double> *zMass = nullptr;
	vector<double> *zPt = nullptr;

    vector<double> *trackPt = nullptr;
    vector<double> *trackMuTagged = nullptr;
    vector<double> *trackWeight = nullptr;
    vector<double> *trackResidualWeight = nullptr;
    vector<int> *subevent = nullptr;

    TreeSig->SetBranchAddress("zMass", &zMass);
    TreeSig->SetBranchAddress("zPt", &zPt);
    TreeSig->SetBranchAddress("hiBin", &hiBin);
    TreeSig->SetBranchAddress("bestZidx", &bestZidx);
    TreeSig->SetBranchAddress("NCollWeight", &NCollWeight);
    TreeSig->SetBranchAddress("ZWeight", &ZWeight);
    TreeSig->SetBranchAddress("VZWeight", &VZWeight);
    TreeSig->SetBranchAddress("trackPt", &trackPt);
    TreeSig->SetBranchAddress("trackMuTagged", &trackMuTagged);
    TreeSig->SetBranchAddress("trackWeight", &trackWeight);
    TreeSig->SetBranchAddress("trackResidualWeight", &trackResidualWeight);
    TreeSig->SetBranchAddress("subevent", &subevent);

    int EntryCount = TreeSig->GetEntries();
    ProgressBar Bar(cout, EntryCount);
    Bar.SetStyle(-1);

	for (int i = 0; i < EntryCount; i++) {
        TreeSig->GetEntry(i);

        if(EntryCount < 500 || (i % (EntryCount / 300)) == 0)
      	{
         Bar.Update(i);
         Bar.Print();
      	}

        // Apply event selection

        bool evtCut = false;

        if(bestZidx >= 0 && zMass->at(bestZidx) > 60 && zPt->at(bestZidx) > ptL && zPt->at(bestZidx) < ptH && hiBin > 2*centL && hiBin < 2*centH)
        	evtCut = true;

        if (!evtCut) 
        	continue;

        // Count tracks for this Z
        int nTracks = 0;

        // Loop over tracks
        for (int j = 0; j < trackPt->size(); j++) {
            
            if (trackPt->at(j) > TptL && trackPt->at(j) < TptH && trackMuTagged->at(j) == 0) {
            	if((isSube0 && subevent->at(j) == 0) || !isSube0)
                	nTracks+=trackWeight->at(j)*trackResidualWeight->at(j);
            }
        }

        // Fill the histogram
        float weight = NCollWeight * ZWeight * VZWeight;
        if(nTracks>100)
        	d10->Fill(nTracks, weight);
    }

    Bar.Update(EntryCount);
   	Bar.Print();
   	Bar.PrintLine();

}



void DrawNTrk_single(float ptL=20,float ptH=2000,float centL=0,float centH=90,float TptL=0,float TptH=10000){


	string HistName = "NTrk";

	TCanvas *c = new TCanvas("c", "", 800, 600);

	TH1D* d10 = new TH1D("d10","", 100, 0, 2000);
	TH1D* d20 = new TH1D("d20","", 100, 0, 2000);
	TH1D* d3 =  new TH1D( "d3","", 100, 0, 2000);
	TH1D* d4 =  new TH1D( "d4","", 100, 0, 2000);
	TH1D* d5 =  new TH1D( "d5","", 100, 0, 2000);

	TCut evtCut = Form("zMass[0]>60&&zPt[0]>%f&&zPt[0]<%f&&hiBin>%f&&hiBin<%f",ptL,ptH,2*centL,2*centH);
	TCut evtCutGen = Form("genZMass[0]>60&&genZPt[0]>%f&&genZPt[0]<%f&&hiBin>%f&&hiBin<%f",ptL,ptH,2*centL,2*centH);
	TCut trkCut = Form("trackMuTagged==0&&trackPt>%f&&trackPt<%f",TptL,TptH);

	TCut evtCutPP = Form("zMass[0]>60&&zPt[0]>%f&&zPt[0]<%f",ptL,ptH);
    
	cout<<"Filling TreeSig"<<endl;
    for (const auto &InputFileName : my_glob(TreeSigNames.c_str())){
    	TFile* input = TFile::Open(InputFileName.c_str(), "READ");

    	if (!input || input->IsZombie()) {
            cerr << "Error opening file: " << InputFileName << endl;
            if (input) delete input;
            continue;
        }

    	TTree* Tree = (TTree*)input->Get("Tree");

    	if (!Tree) {
            cerr << "Error getting tree from file: " << InputFileName << endl;
            delete input;
            continue;
        }

    	FillNTrks(d10, Tree, ptL,ptH,centL,centH,TptL,TptH, false);
    	input->Close();
    }

    cout<<"Filling TreeBkg"<<endl;

    for (const auto &InputFileName : my_glob(TreeBkgNames.c_str())){
    	TFile* input = TFile::Open(InputFileName.c_str(), "READ");
    	
    	if (!input || input->IsZombie()) {
            cerr << "Error opening file: " << InputFileName << endl;
            if (input) delete input;
            continue;
        }


    	TTree* Tree = (TTree*)input->Get("Tree");

    	if (!Tree) {
            cerr << "Error getting tree from file: " << InputFileName << endl;
            delete input;
            continue;
        }

    	FillNTrks(d20, Tree, ptL,ptH,centL,centH,TptL,TptH, false);
    	input->Close();
    }

    cout<<"Filling TreeSgG"<<endl;

    for (const auto &InputFileName : my_glob(TreeSgGNames.c_str())){
    	TFile* input = TFile::Open(InputFileName.c_str(), "READ");
    	
    	if (!input || input->IsZombie()) {
            cerr << "Error opening file: " << InputFileName << endl;
            if (input) delete input;
            continue;
        }

    	TTree* Tree = (TTree*)input->Get("Tree");

    	if (!Tree) {
            cerr << "Error getting tree from file: " << InputFileName << endl;
            delete input;
            continue;
        }

    	FillNTrks(d3,  Tree, ptL,ptH,centL,centH,TptL,TptH, true);
    	FillNTrks(d4,  Tree, ptL,ptH,centL,centH,TptL,TptH, false);
    	input->Close();
    }

    cout<<"Filling TreeBgG"<<endl;

    for (const auto &InputFileName : my_glob(TreeBgGNames.c_str())){
    	TFile* input = TFile::Open(InputFileName.c_str(), "READ");
    	
    	if (!input || input->IsZombie()) {
            cerr << "Error opening file: " << InputFileName << endl;
            if (input) delete input;
            continue;
        }

    	TTree* Tree = (TTree*)input->Get("Tree");

    	if (!Tree) {
            cerr << "Error getting tree from file: " << InputFileName << endl;
            delete input;
            continue;
        }

    	FillNTrks(d5,  Tree, ptL,ptH,centL,centH,TptL,TptH, false);
    	input->Close();
    }


	d10->SetLineColor(kBlack);//sig reco
  	d20->SetLineColor(kBlue);//bkg reco
  	d3->SetLineColor(kBlack);//pp reco or sig gen sube=0

  	d10->SetMarkerStyle(kFullCircle);
  	d20->SetMarkerStyle(kFullCircle);

  	d10->SetMarkerColor(kBlack);
  	d20->SetMarkerColor(kBlue);

  	d10->SetMarkerSize(0.5);
  	d20->SetMarkerSize(0.5);

  	d4->SetLineColor(kGray+2);//grey //sig gen
	d5->SetLineColor(TColor::GetColor("#377eb8"));//blue //bkg gen

	d4->SetFillStyle(3345);
	d5->SetFillStyle(3354);

	d4->SetFillColor(kGray+2);
	d5->SetFillColor(TColor::GetColor("#377eb8"));

	// == Start drawing == //

   	gSystem->Exec(Form("mkdir -p /eos/user/p/pchou/figs/%s/%s",typeofdata,HistName.c_str()));

   	d10->Draw("ep");
   	d20->Draw("ep same");
   	d3->Draw("hist same");

   	d4->Draw("hist same");
   	d5->Draw("hist same");
    d10->GetXaxis()->SetTitle("Number of Tracks per Z");
    d10->GetYaxis()->SetTitle("Number of Z Events");

   
	TLegend leg1(0.58,0.65,0.98,0.92);
	leg1.AddEntry(d10 ,"Raw MC RECO","lep");
	leg1.AddEntry(d20 ,"Bkg MC RECO","lep");
	leg1.AddEntry(d4 ,"Raw MC GEN","lf");
	leg1.AddEntry(d5 ,"Bkg MC GEN","lf");
	leg1.AddEntry(d3 ,"Sig GEN (sube=0)","l");
	leg1.SetFillColorAlpha(kWhite,0);
	leg1.SetLineColor(kBlack);
	leg1.SetLineWidth(1);
	leg1.Draw();

	double mean = d10->GetMean();
    double rms = d10->GetRMS();
    TLatex *latex = new TLatex(0.15, 0.85, Form("Reco Raw: Mean = %.2f, RMS = %.2f",  d10->GetMean(),  d10->GetRMS()));
    latex->SetNDC();
    latex->SetTextFont(42);
   	latex->SetTextSize(0.03);
    latex->Draw();
    latex = new TLatex(0.15, 0.80, Form("Reco Bkg: Mean = %.2f, RMS = %.2f",  d20->GetMean(),  d20->GetRMS()));
    latex->SetNDC();
    latex->SetTextFont(42);
   	latex->SetTextSize(0.03);
    latex->Draw();

  	c->SaveAs(Form("/eos/user/p/pchou/figs/%s/%s/Ztrack_%s_com_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f_low.png",typeofdata,HistName.c_str(),typeofdata1,ptL,ptH,centL,centH,TptL,TptH)); 

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

	DrawNTrk_single(ptL,ptH,centL,centH,TptL,TptH);

}
