#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TMath.h>
#include <TROOT.h>
#include <cmath>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TCut.h>
#include <TSystem.h>
#include <TLatex.h>
#include <vector>
#include <iostream>

using namespace std;

#include "CommandLine.h"

int main(int argc, char *argv[]);


TChain *TreeSig = new TChain("Tree"); 
TChain *TreeBkg = new TChain("Tree"); 
TChain *TreeSgG = new TChain("Tree"); 
TChain *TreeBgG = new TChain("Tree"); 

const char *typeofdata = "ZHadron2024/DrawNTrk/ov1_v4b_sub0_Rres_tol120/20240831/";
const char *typeofdata1 = "ov1_v4b_sub0_Rres_tol120";
    
string filebase = "/eos/cms/store/group/phys_heavyions/pchou/SkimZHadron2024/";
string cernbox = "/eos/home-p/pchou/SkimZHadron2024/";

int main(int argc, char *argv[]){
//void DrawNTrk_simple_single(float ptL=40,float ptH=200, float centL=0, float centH=10, float TptL=1, float TptH=2){

    CommandLine CL(argc, argv);

    double ptL                   = CL.GetDouble("ptL", 40);
    double ptH                   = CL.GetDouble("ptH", 200);
    double centL                 = CL.GetDouble("centL", 0);
    double centH                 = CL.GetDouble("centH", 10);
    double TptL                  = CL.GetDouble("TptL", 1);
    double TptH                  = CL.GetDouble("TptH", 2);
    double MinNtrk               = CL.GetDouble("MinNtrk", 0);
    double MaxNtrk               = CL.GetDouble("MaxNtrk", 4000);
    double Nbins                 = CL.GetDouble("Nbins", 100);


    gStyle->SetOptStat(0);

    
    TreeSig->Add((filebase + "OutputMC_v4b_ee_v3/Result*.root").c_str());
    TreeBkg->Add((cernbox  + "OutputMCBkg_v4b_ee_Rres_tol120_v3/Result*.root").c_str());
    TreeSgG->Add((filebase + "OutputMCGen_v3c_ee/Result*.root").c_str());
    TreeBgG->Add((filebase + "OutputMCbkgGen_v4_ee_tol120/Result*.root").c_str());
    
    //float ptL=40, ptH=200, centL=0, centH=10, TptL=1, TptH=2;
    
    TCut evtCut = Form("zMass[0]>60&&zPt[0]>%f&&zPt[0]<%f&&hiBin>%f&&hiBin<%f",ptL,ptH,2*centL,2*centH);
    TCut evtCutGen = Form("genZMass[0]>60&&genZPt[0]>%f&&genZPt[0]<%f&&hiBin>%f&&hiBin<%f",ptL,ptH,2*centL,2*centH);
    string trkCut = Form("trackMuTagged==0&&trackPt>%f&&trackPt<%f",TptL,TptH);
    string trkCutGen = Form("GenTrackMuTagged==0&&GenTrackPt>%f&&GenTrackPt<%f",TptL,TptH);
    
    TH1D* d1 =  new TH1D( "d1","",Nbins,MinNtrk,MaxNtrk);
    TH1D* d2 =  new TH1D( "d2","",Nbins,MinNtrk,MaxNtrk);
    TH1D* d3 =  new TH1D( "d3","",Nbins,MinNtrk,MaxNtrk);
    TH1D* d4 =  new TH1D( "d4","",Nbins,MinNtrk,MaxNtrk);
    TH1D* d5 =  new TH1D( "d5","",Nbins,MinNtrk,MaxNtrk);
    TH1D* d6 =  new TH1D( "d6","",Nbins,MinNtrk,MaxNtrk);
    TH1D* d7 =  new TH1D( "d7","",Nbins,MinNtrk,MaxNtrk);
    TH1D* d8 =  new TH1D( "d8","",Nbins,MinNtrk,MaxNtrk);
    TH1D* d9 =  new TH1D( "d9","",Nbins,MinNtrk,MaxNtrk);
    TH1D* d10= new TH1D( "d10","",Nbins,MinNtrk,MaxNtrk);
    
    
    TreeSig->Draw(Form("Sum$(trackWeight*trackResidualWeight*(%s))>>d1",trkCut.c_str()),"(NCollWeight*ZWeight*VZWeight)"*(evtCut),"goff");
    TreeSig->Draw(Form("Sum$((%s))>>d2",trkCut.c_str()),(evtCut),"goff");
    
    TreeBkg->Draw(Form("Sum$(trackWeight*trackResidualWeight*(%s))>>d3",trkCut.c_str()),"(NCollWeight*ZWeight*VZWeight)"*(evtCut),"goff");
    TreeBkg->Draw(Form("Sum$((%s))>>d4",trkCut.c_str()),(evtCut),"goff");
    
    TreeSgG->Draw(Form("Sum$(trackWeight*trackResidualWeight*(%s))>>d5",trkCut.c_str()),"(NCollWeight*ZWeight*VZWeight)"*(evtCutGen),"goff");
    TreeSgG->Draw(Form("Sum$((%s))>>d6",trkCut.c_str()),(evtCutGen),"goff");
    
    TreeBgG->Draw(Form("Sum$(trackWeight*trackResidualWeight*(%s))>>d7",trkCut.c_str()),"(NCollWeight*ZWeight*VZWeight)"*(evtCutGen),"goff");
    TreeBgG->Draw(Form("Sum$((%s))>>d8",trkCut.c_str()),(evtCutGen),"goff");

    TreeBkg->Draw(Form("Sum$((%s))>>d9",trkCutGen.c_str()),(evtCutGen),"goff");
    TreeBgG->Draw(Form("Sum$((%s))>>d10",trkCutGen.c_str()),(evtCutGen),"goff");

    
    d1->Sumw2();
    d2->Sumw2();
    d3->Sumw2();
    d4->Sumw2();
    d5->Sumw2();
    d6->Sumw2();
    d7->Sumw2();
    d8->Sumw2();
    d9->Sumw2();
    d10->Sumw2();

    
    
    d1->SetLineColor(kRed);
    d2->SetLineColor(kRed);
    d3->SetLineColor(kBlue);
    d4->SetLineColor(kBlue);
    d5->SetLineColor(kGreen);
    d6->SetLineColor(kGreen);
    d7->SetLineColor(kBlack);
    d8->SetLineColor(kBlack);
    d9->SetLineColor(kViolet);
    d10->SetLineColor(kViolet);
    
    d1->SetMarkerColor(kRed);
    d2->SetMarkerColor(kRed);
    d3->SetMarkerColor(kBlue);
    d4->SetMarkerColor(kBlue);
    d5->SetMarkerColor(kGreen);
    d6->SetMarkerColor(kGreen);
    d7->SetMarkerColor(kBlack);
    d8->SetMarkerColor(kBlack);
    d9->SetMarkerColor(kViolet);
    d10->SetMarkerColor(kViolet);
    

    
    d1->SetMarkerStyle(20);
    d2->SetMarkerStyle(20);
    d3->SetMarkerStyle(20);
    d4->SetMarkerStyle(20);
    d5->SetMarkerStyle(20);
    d6->SetMarkerStyle(20);
    d7->SetMarkerStyle(20);
    d8->SetMarkerStyle(20);
    d9->SetMarkerStyle(20);
    d10->SetMarkerStyle(20);
    

    
    d1->SetMarkerSize(0.5);
    d2->SetMarkerSize(0.5);
    d3->SetMarkerSize(0.5);
    d4->SetMarkerSize(0.5);
    d5->SetMarkerSize(0.5);
    d6->SetMarkerSize(0.5);
    d7->SetMarkerSize(0.5);
    d8->SetMarkerSize(0.5);
    d9->SetMarkerSize(0.5);
    d10->SetMarkerSize(0.5);
    

    
    d1->SetLineWidth(2);
    d2->SetLineWidth(2);
    d3->SetLineWidth(2);
    d4->SetLineWidth(2);
    d5->SetLineWidth(2);
    d6->SetLineWidth(2);
    d7->SetLineWidth(2);
    d8->SetLineWidth(2);
    d9->SetLineWidth(2);
    d10->SetLineWidth(2);
    


    d1->SetFillStyle(3345);
    d2->SetFillStyle(3354);
    d3->SetFillStyle(3345);
    d4->SetFillStyle(3354);
    d5->SetFillStyle(3345);
    d6->SetFillStyle(3354);
    d7->SetFillStyle(3345);
    d8->SetFillStyle(3354);
    d9->SetFillStyle(3345);
    d10->SetFillStyle(3354);
    


    d1->SetFillColor(kRed);
    d3->SetFillColor(kBlue);
    d2->SetFillColor(kRed);
    d4->SetFillColor(kBlue);

    d5->SetFillColor(kGreen);
    d6->SetFillColor(kGreen);
    d7->SetFillColor(kBlack);
    d8->SetFillColor(kBlack);
    d9->SetFillColor(kViolet);
    d10->SetFillColor(kViolet);

    
    d1->GetXaxis()->SetTitle("N_{trk}/Event");
    d2->GetXaxis()->SetTitle("N_{trk}/Event");
    d3->GetXaxis()->SetTitle("N_{trk}/Event");
    d4->GetXaxis()->SetTitle("N_{trk}/Event");
    d5->GetXaxis()->SetTitle("N_{trk}/Event");
    d6->GetXaxis()->SetTitle("N_{trk}/Event");
    d7->GetXaxis()->SetTitle("N_{trk}/Event");
    d8->GetXaxis()->SetTitle("N_{trk}/Event");
    d9->GetXaxis()->SetTitle("N_{trk}/Event");
    d10->GetXaxis()->SetTitle("N_{trk}/Event");

    
    d1->GetYaxis()->SetTitle("");
    d2->GetYaxis()->SetTitle("");
    d3->GetYaxis()->SetTitle("");
    d4->GetYaxis()->SetTitle("");
    d5->GetYaxis()->SetTitle("");
    d6->GetYaxis()->SetTitle("");
    d7->GetYaxis()->SetTitle("");
    d8->GetYaxis()->SetTitle("");
    d9->GetYaxis()->SetTitle("");
    d10->GetYaxis()->SetTitle("");



    double x0 = 0.5, y0 = 0.46, deltay = 0.04;

    TLatex *latex1 = new TLatex(x0, y0, Form("#mu_{1} = %.2f #pm %.2f, #sigma_{1} = %.2f #pm %.2f",   d1->GetMean(), d1->GetMeanError(),  d1->GetRMS(), d1->GetRMSError()));
    latex1->SetNDC();
    latex1->SetTextFont(42);
    latex1->SetTextSize(0.025);
    
    TLatex *latex2 = new TLatex(x0, y0-deltay, Form("#mu_{2} = %.2f #pm %.2f, #sigma_{2} = %.2f  #pm %.2f",  d2->GetMean(), d2->GetMeanError(),  d2->GetRMS(), d2->GetRMSError()));
    latex2->SetNDC();
    latex2->SetTextFont(42);
    latex2->SetTextSize(0.025);
    
    TLatex *latex3 = new TLatex(x0, y0-2*deltay, Form("#mu_{3} = %.2f #pm %.2f, #sigma_{3} = %.2f  #pm %.2f",  d3->GetMean(), d3->GetMeanError(),  d3->GetRMS(), d3->GetRMSError()));
    latex3->SetNDC();
    latex3->SetTextFont(42);
    latex3->SetTextSize(0.025);

    TLatex *latex4 = new TLatex(x0, y0-3*deltay, Form("#mu_{4} = %.2f #pm %.2f, #sigma_{4} = %.2f  #pm %.2f",  d4->GetMean(), d4->GetMeanError(),  d4->GetRMS(), d4->GetRMSError()));
    latex4->SetNDC();
    latex4->SetTextFont(42);
    latex4->SetTextSize(0.025);

    TLatex *latex5 = new TLatex(x0, y0-4*deltay, Form("#mu_{5} = %.2f #pm %.2f, #sigma_{5} = %.2f  #pm %.2f",  d5->GetMean(), d5->GetMeanError(),  d5->GetRMS(), d5->GetRMSError()));
    latex5->SetNDC();
    latex5->SetTextFont(42);
    latex5->SetTextSize(0.025);

    TLatex *latex6 = new TLatex(x0, y0-5*deltay, Form("#mu_{6} = %.2f #pm %.2f, #sigma_{6} = %.2f  #pm %.2f",  d6->GetMean(), d6->GetMeanError(),  d6->GetRMS(), d6->GetRMSError()));
    latex6->SetNDC();
    latex6->SetTextFont(42);
    latex6->SetTextSize(0.025);

    TLatex *latex7 = new TLatex(x0, y0-6*deltay, Form("#mu_{7} = %.2f #pm %.2f, #sigma_{7} = %.2f  #pm %.2f",  d7->GetMean(), d7->GetMeanError(),  d7->GetRMS(), d7->GetRMSError()));
    latex7->SetNDC();
    latex7->SetTextFont(42);
    latex7->SetTextSize(0.025);

    TLatex *latex8 = new TLatex(x0, y0-7*deltay, Form("#mu_{8} = %.2f #pm %.2f, #sigma_{8} = %.2f  #pm %.2f",  d8->GetMean(), d8->GetMeanError(),  d8->GetRMS(), d8->GetRMSError()));
    latex8->SetNDC();
    latex8->SetTextFont(42);
    latex8->SetTextSize(0.025);

    TLatex *latex9 = new TLatex(x0, y0-8*deltay, Form("#mu_{9} = %.2f #pm %.2f, #sigma_{9} = %.2f  #pm %.2f",  d9->GetMean(), d9->GetMeanError(),  d9->GetRMS(), d9->GetRMSError()));
    latex9->SetNDC();
    latex9->SetTextFont(42);
    latex9->SetTextSize(0.025);

    TLatex *latex10 = new TLatex(x0, y0-9*deltay, Form("#mu_{10} = %.2f #pm %.2f, #sigma_{10} = %.2f  #pm %.2f",  d10->GetMean(), d10->GetMeanError(),  d10->GetRMS(), d10->GetRMSError()));
    latex10->SetNDC();
    latex10->SetTextFont(42);
    latex10->SetTextSize(0.025);
    
    d1->Scale(1./d1->Integral());
    d2->Scale(1./d2->Integral());
    d3->Scale(1./d3->Integral());
    d4->Scale(1./d4->Integral());
    d5->Scale(1./d5->Integral());
    d6->Scale(1./d6->Integral());
    d7->Scale(1./d7->Integral());
    d8->Scale(1./d8->Integral());
    d9->Scale(1./d9->Integral());
    d10->Scale(1./d10->Integral());



    double max1 = max(d1->GetMaximum(),d2->GetMaximum());
    double max2 = max(d3->GetMaximum(),d4->GetMaximum());
    double max3 = max(d5->GetMaximum(),d6->GetMaximum());
    double max4 = max(d7->GetMaximum(),d8->GetMaximum());
    double max5 = max(d9->GetMaximum(),d10->GetMaximum());


    double max0 = max(max1,max(max2,max(max3,max(max4,max5))));
    d1->SetMaximum(max0*1.2);

    
    
    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    c1->cd();
    d1->Draw("hist");
    d2->Draw("hist same");
    d3->Draw("hist same");
    d4->Draw("hist same");
    d5->Draw("hist same");
    d6->Draw("hist same");
    d7->Draw("hist same");
    d8->Draw("hist same");
    d9->Draw("hist same");
    d10->Draw("hist same");


    
    TLegend *leg = new TLegend(0.5,0.5,0.9,0.9);
    leg->AddEntry(d1,"1: Raw (weighted)","lf");
    leg->AddEntry(d2,"2: Raw (unweighted)","lf");
    leg->AddEntry(d3,"3: Background (weighted)","lf");
    leg->AddEntry(d4,"4: Background (unweighted)","lf");
    leg->AddEntry(d5,"5: Raw Gen (weighted)","lf");
    leg->AddEntry(d6,"6: Raw Gen (unweighted)","lf");
    leg->AddEntry(d7,"7: Background Gen (weighted)","lf");
    leg->AddEntry(d8,"8: Background Gen (unweighted)","lf");
    leg->AddEntry(d9,"9: Background Reco Gen (unweighted)","lf");
    leg->AddEntry(d10,"10: Background Gen Gen (unweighted)","lf");
    leg->SetFillColorAlpha(kWhite,0);
    leg->SetLineColor(kBlack);
    leg->SetLineWidth(1);
    leg->Draw();

    latex1->Draw();
    latex2->Draw();
    latex3->Draw();
    latex4->Draw();
    latex5->Draw();
    latex6->Draw();
    latex7->Draw();
    latex8->Draw();
    latex9->Draw();
    latex10->Draw();



    TLatex *latex = new TLatex(0.1, 0.92, Form("p_{T}^{Z} = %.0f-%.0f GeV/c, Centrality = %.0f-%.0f %%, p_{T}^{trk} = %.0f-%.0f GeV/c",ptL,ptH,centL,centH,TptL,TptH));
    latex->SetNDC();
    latex->SetTextFont(42);
    latex->SetTextSize(0.025);
    latex->Draw();

    gSystem->mkdir(Form("/eos/user/p/pchou/figs/%s/%s",typeofdata,typeofdata1),kTRUE);
    //gSystem->Exec(Form("mkdir -p /eos/user/p/pchou/figs/%s/%s",typeofdata,typeofdata1));
    
    c1->SaveAs(Form("/eos/user/p/pchou/figs/%s/%s/Ntrk_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f.png",typeofdata,typeofdata1,ptL,ptH,centL,centH,TptL,TptH)); 

}