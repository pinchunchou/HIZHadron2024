#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <TROOT.h>
#include <cmath>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TLine.h>
#include <iostream>

using namespace std;

double TptL_min = 0.5;

int main(int argc, char *argv[]);
void ZBasicBkgSub_single(int binnum,float ptL,float ptH,float centL,float centH,float TptL,float TptH,
   string HistName, string XTitleName, string YTitleName, int rebinnum);
void ZBasicBkgSub_loop(int binnum,float ptL,float ptH,float centL,float centH,float TptL,float TptH);

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

TFile *file_sigMC;
TFile *file_bkgMC;
TFile *file_ppMC;
TFile *file_ppbkgMC;


const char *typeofdata = "ZHadron2024/BasicBkgSub/ov1_v1c_Reco_noZw_comparesub0/20240706/";
const char *typeofdata1 = "ov1_v1c_Reco_noZw_comparesub0";


bool selfmix = false;
bool isgen   = false;
bool drawlog = false;
bool drawrat = false;
bool drawlow = true;

bool comparesub0 = true;

void ZBasicBkgSub_single(int binnum=40,float ptL=20,float ptH=2000,float centL=0,float centH=90,float TptL=0,float TptH=10000, 
   string HistName="HPhi", string XTitleName = "#Delta#phi_{Z,track}", string YTitleName = "dN/d#Delta#phi", int rebinnum=1)
{
   std::cout<<"ptL = "<<ptL<<", ptH = "<<ptH<<", centL = "<<centL<<", centH = "<<centH<<", TptL = "<<TptL<<", TptH = "<<TptH<<std::endl;
   
   TCanvas *c = new TCanvas("c", "", 500, 600);
   TPad *Pad  = new TPad("Pad" , "", 0, 0.3, 1,   1);
   TPad *RPad = new TPad("RPad", "", 0, 0. , 1, 0.3);

   Pad->Draw();
   RPad->Draw();

   std::cout<<"Getting histograms..."<<std::endl;

   std::string FolderName = Form("Plot_ZPT_%.0f_%.0f_Cent_%.0f_%.0f_TrackPT_%.2f_%.2f",ptL,ptH,centL,centH,TptL,TptH);
   std::replace(FolderName.begin(), FolderName.end(), '.', 'p');

   std::string HistNameGen = HistName;
   if(isgen==false)
      HistNameGen.replace(HistNameGen.find('H'), 1, "HGen");
      //std::replace(HistNameGen.begin(), HistNameGen.end(), 'H', 'HGen');

   cout<<"HistNameGen = "<<HistNameGen<<endl;
   TH1D* hMC_phi = (TH1D*) file_sigMC->Get(Form("%s/%s", FolderName.c_str(), HistName.c_str()));
   cout<<"a"<<endl;

   std::string HistNamePP = HistName;
   if(comparesub0)
      HistNamePP.replace(HistNamePP.find('H'), 1, "HGen");

   TH1D* hpp_phi = (TH1D*) file_ppMC->Get(Form("%s/%s", FolderName.c_str(), HistNamePP.c_str()));//HistNameGen
   cout<<"b"<<endl;

   TH1D* hpp_bkg_phi;

   if(selfmix)
      hpp_bkg_phi = (TH1D*) file_ppbkgMC->Get(Form("%s/%s", FolderName.c_str(), HistName.c_str()));

   cout<<"c"<<endl;

   hMC_phi->SetName("hMC_phi");
   hpp_phi->SetName("hpp_phi");

   if(selfmix)
      hpp_bkg_phi->SetName("hpp_bkg_phi");

   TH1D* hMC_bkg_phi = (TH1D*) file_bkgMC->Get(Form("%s/%s", FolderName.c_str(), HistName.c_str()));
   hMC_bkg_phi->SetName("hMC_bkg_phi");

   std::cout<<"Getting Entries..."<<std::endl;

   //TH1D *nM_tN     = (TH1D *) file_sigMC->Get(Form("%s/HTotalEventCount",FolderName.c_str()));
   //TH1D *nMb_tN    = (TH1D *) file_bkgMC->Get(Form("%s/HTotalEventCount",FolderName.c_str()));
   //TH1D *npM_tN    = (TH1D *) file_ppMC ->Get(Form("%s/HTotalEventCount",FolderName.c_str()));
   //TH1D *npb_tN    = (TH1D *) file_ppbkgMC ->Get(Form("%s/HTotalEventCount",FolderName.c_str()));
   
   string NameEntryCount;
   if(isgen)
      NameEntryCount = "GenEntryCount";
   else
      NameEntryCount = "EntryCount";


   string NameEntryCountPP = NameEntryCount;
   if(comparesub0)
      NameEntryCountPP = "GenEntryCount";

   TNamed *nM_tN  = (TNamed *) file_sigMC->Get(Form("%s/%s",FolderName.c_str(),NameEntryCount.c_str()));
   TNamed *nMb_tN = (TNamed *) file_bkgMC->Get(Form("%s/%s",FolderName.c_str(),NameEntryCount.c_str()));
   TNamed *npM_tN = (TNamed *) file_ppMC->Get(Form("%s/%s",FolderName.c_str(),NameEntryCountPP.c_str()));//GenEntryCount
   TNamed *npb_tN;

   if(selfmix)
      npb_tN= (TNamed *) file_ppbkgMC->Get(Form("%s/%s",FolderName.c_str(),NameEntryCount.c_str()));

   std::cout<<"hi"<<std::endl;
   std::string sM_tN  = (std::string) nM_tN->GetTitle();
   std::string sMb_tN = (std::string) nMb_tN->GetTitle();
   std::string spM_tN = (std::string) npM_tN->GetTitle();
   std::string spb_tN = selfmix ? (std::string) npb_tN->GetTitle() : "0";

   //float tM_tN     =     nM_tN->GetBinContent(1);
   //float tMb_tN    =    nMb_tN->GetBinContent(1);
   //float tpM_tN    =    npM_tN->GetBinContent(1);  
   //float tpb_tN    =    npb_tN->GetBinContent(1);  

   std::cout<<"hi1"<<std::endl;
   float tM_tN  = std::stof(sM_tN);
   float tMb_tN = std::stof(sMb_tN);
   float tpM_tN = std::stof(spM_tN);
   float tpb_tN = std::stof(spb_tN);

   std::cout<<"tM_tN = "<<tM_tN<<std::endl;
   std::cout<<"tMb_tN = "<<tMb_tN<<std::endl;
   std::cout<<"tpM_tN = "<<tpM_tN<<std::endl;
   std::cout<<"tpb_tN = "<<tpb_tN<<std::endl;

   double bin_width = hMC_phi->GetBinWidth(1);

   hMC_phi->Sumw2();
   hMC_bkg_phi->Sumw2();
   hpp_phi->Sumw2();

   hMC_phi->Rebin(rebinnum);
   hMC_bkg_phi->Rebin(rebinnum);
   hpp_phi->Rebin(rebinnum);

   if(selfmix)
      hpp_bkg_phi->Rebin(rebinnum);

   bin_width*=rebinnum;

   hMC_phi->Scale(1./tM_tN/bin_width);
   hMC_bkg_phi->Scale(1./tMb_tN/bin_width);
   hpp_phi->Scale(1./tpM_tN/bin_width);

   if(selfmix)
      hpp_bkg_phi->Scale(1./tpb_tN/bin_width);


   TH1D *hMC_sb_phi = (TH1D*) hMC_phi->Clone("hMC_sb_phi");

   hMC_sb_phi->Add(hMC_bkg_phi,-1);

   if(selfmix)
      hpp_phi->Add(hpp_bkg_phi,-1); //SigBkg

   int countM0 = hMC_phi->GetEntries();
   std::cout<<"MC 0 = "<<countM0<<std::endl;

   int countMb = hMC_bkg_phi->GetEntries();
   std::cout<<"MC Bkg = "<<countMb<<std::endl;

   std::cout<<"Drawing..."<<std::endl;

   hMC_phi->SetMarkerStyle(24);
   hMC_bkg_phi->SetMarkerStyle(24);
   hMC_sb_phi->SetMarkerStyle(24);

   hMC_phi->SetMarkerColor(kRed);
   hMC_phi->SetLineColor(kRed);
   hMC_bkg_phi->SetMarkerColor(kRed);
   hMC_bkg_phi->SetLineColor(kRed);
   hMC_sb_phi->SetMarkerColor(kRed);
   hMC_sb_phi->SetLineColor(kRed);

   if(TptL==0) TptL=TptL_min;

   string typeofsample;

   if(selfmix)
      typeofsample = "Z - Z #times Z (MC)";
   else if(isgen)
      typeofsample = "Nominal MC Gen";
   else
      typeofsample = "Nominal MC";

   TLatex *pt0 = new TLatex(0.25,0.88,typeofsample.c_str());
   //TLatex *pt0 = new TLatex(0.15,0.82,"Nominal MC GEN (Pythia+Hydjet)");
   pt0->SetTextFont(42);
   pt0->SetTextSize(0.03);
   pt0->SetNDC(kTRUE);

   TLatex *pt = new TLatex(0.25,0.82,Form("%.0f %%< Centrality < %.0f %%",centL,centH));
   pt->SetTextFont(42);
   pt->SetTextSize(0.03);
   pt->SetNDC(kTRUE);

   TLatex *pt2 = new TLatex(0.25,0.76,Form("%.1f < Z p_{T} < %.1f GeV",ptL,ptH));
   pt2->SetTextFont(42);
   pt2->SetTextSize(0.03);
   pt2->SetNDC(kTRUE);

   TLatex *pt3 = new TLatex(0.25,0.70,Form("%.1f < Track p_{T} < %.1f GeV",TptL,TptH));
   pt3->SetTextFont(42);
   pt3->SetTextSize(0.03);
   pt3->SetNDC(kTRUE);

   // == Start drawing == //

   gSystem->Exec(Form("mkdir -p /eos/user/p/pchou/figs/%s/%s",typeofdata,HistName.c_str()));

   double max1 = hMC_phi->GetMaximum();
   double max2 = hMC_bkg_phi->GetMaximum();
   double max3 = hpp_phi->GetMaximum();

   hMC_phi->SetMarkerColor(kBlack);
   hMC_bkg_phi->SetMarkerColor(kBlue);
   hMC_sb_phi->SetMarkerColor(kRed);

   std::cout<< "bin_width = " << bin_width;
   std::cout<< ", rebinnum = " << rebinnum << std::endl;
   std::cout<< "hMC_phi->Integral() = "<< hMC_phi->Integral();
   std::cout<< ", hMC_bkg_phi->Integral() = "<< hMC_bkg_phi->Integral();
   std::cout<< ", hMC_sb_phi->Integral() = "<< hMC_sb_phi->Integral();
   std::cout<< ", hpp_phi->Integral() = "<< hpp_phi->Integral()<<std::endl;

   double sig_int=0, bkg_int=0, sb_int=0, pp_int=0;
   double sig_err=0, bkg_err=0, sb_err=0, pp_err=0;

   sig_int = hMC_phi->IntegralAndError(1,hMC_phi->GetNbinsX(),sig_err,"width");
   bkg_int = hMC_bkg_phi->IntegralAndError(1,hMC_bkg_phi->GetNbinsX(),bkg_err,"width");
   sb_int = hMC_sb_phi->IntegralAndError(1,hMC_sb_phi->GetNbinsX(),sb_err,"width");
   pp_int = hpp_phi->IntegralAndError(1,hpp_phi->GetNbinsX(),pp_err,"width");


   //TLatex *ptInt1 = new TLatex(0.15,0.56,Form("#Sigma raw = %.2f #pm %.2f,  #Sigma bkg = %.2f #pm %.2f",sig_int,sig_err,bkg_int,bkg_err));
   TLatex *ptInt1 = new TLatex(0.25,0.64,Form("#Sigma raw = %.2f,  #Sigma bkg = %.2f",sig_int,bkg_int));
   ptInt1->SetTextFont(42);
   ptInt1->SetTextSize(0.03);
   ptInt1->SetNDC(kTRUE);

   string pptext;

   if(selfmix)
      pptext = "#Sigma pp (raw-bkg)";
   else if(comparesub0)
      pptext = "sig pythia gen";
   else
      pptext = "#Sigma pp";//sig pythia gen

   //TLatex *ptInt2 = new TLatex(0.15,0.50,Form("#Sigma (raw-bkg) = %.2f #pm %.2f,  %s = %.2f #pm %.2f",sb_int,sb_err,pptext.c_str(),pp_int,pp_err));
   TLatex *ptInt2 = new TLatex(0.25,0.58,Form("#Sigma (raw-bkg) = %.2f,  %s = %.2f",sb_int,pptext.c_str(),pp_int));
   ptInt2->SetTextFont(42);
   ptInt2->SetTextSize(0.03);
   ptInt2->SetNDC(kTRUE);

   hMC_phi->SetLineColor(kBlack);
   hMC_bkg_phi->SetLineColor(kBlue);
   hMC_sb_phi->SetLineColor(kRed);
   hpp_phi->SetLineColor(kBlack);

   hMC_phi->SetMarkerStyle(kFullCircle);
   hMC_bkg_phi->SetMarkerStyle(kFullCircle);
   hMC_sb_phi->SetMarkerStyle(kFullCircle);

   hMC_phi->SetMarkerSize(0.5);
   hMC_bkg_phi->SetMarkerSize(0.5);
   hMC_sb_phi->SetMarkerSize(0.5);

   //hpp_phi->SetLineWidth(2);

   double min1 = hMC_sb_phi->GetMinimum();
   double min2 = hpp_phi->GetMinimum();
   
   Pad->cd();

   if(min1>min2) min1=min2;

   if(max1<max2) hMC_bkg_phi->Draw("ep");
   else hMC_phi->Draw("ep");
   hMC_phi->Draw("ep same");
   hMC_bkg_phi->Draw("ep same");

   hMC_sb_phi->Draw("ep same");
   hpp_phi->Draw("hist same");

   if(max1<max2) max1=max2;

   hMC_phi->SetYTitle(YTitleName.c_str());
   hMC_bkg_phi->SetYTitle(YTitleName.c_str());

   hMC_phi->SetXTitle(XTitleName.c_str());
   hMC_bkg_phi->SetXTitle(XTitleName.c_str());

   string genlegtxt = "";
   if(isgen)
      genlegtxt = " Gen";

   TLegend leg1(0.58,0.72,0.98,0.90);
   leg1.AddEntry(hMC_phi ,Form("raw%s",genlegtxt.c_str()),"lep");
   leg1.AddEntry(hMC_bkg_phi ,Form("bkg%s",genlegtxt.c_str()),"lep");
   //leg1.AddEntry(hMC_bkg_phi ,"bkg (subevt#neq0)","lep");
   //leg1.AddEntry(hMC_bkg_phi ,"raw Gen (subevt#neq0)","lep");
   leg1.AddEntry(hMC_sb_phi ,Form("raw-bkg%s",genlegtxt.c_str()),"lep");
   if(selfmix){
      leg1.AddEntry(hpp_phi ,Form("pp raw-bkg%s",genlegtxt.c_str()),"l");
   }else{
      if(comparesub0)
         leg1.AddEntry(hpp_phi ,"sig GEN (sube=0)","lep");
      else
         leg1.AddEntry(hpp_phi ,Form("pp%s (NPU=0)",genlegtxt.c_str()),"l");
      //leg1.AddEntry(hpp_phi ,"sig GEN","l");
   }
   leg1.SetFillColorAlpha(kWhite,0);
   leg1.SetLineColor(kBlack);
   leg1.SetLineWidth(1);
   leg1.Draw();

   pt0->Draw();
   pt->Draw();
   pt2->Draw();
   pt3->Draw();

   ptInt1->Draw();
   ptInt2->Draw();

   //std::cout<<"max1 = "<<max1<<std::endl;
   hMC_phi->SetMaximum(2*max1);
   hMC_bkg_phi->SetMaximum(2*max1);

   if(min1<0){
      hMC_phi->SetMinimum(2*min1);
      hMC_bkg_phi->SetMinimum(2*min1);
      hMC_sb_phi->SetMinimum(2*min1);
      hpp_phi->SetMinimum(2*min1);
   }else{
      hMC_phi->SetMinimum(0);
      hMC_bkg_phi->SetMinimum(0);
      hMC_sb_phi->SetMinimum(0);
      hpp_phi->SetMinimum(0);
   }

   RPad->cd();

   TH1D *PbPb_to_pp = (TH1D*) hMC_sb_phi->Clone("PbPb_to_pp");
   //PbPb_to_pp->Divide(hpp_phi);
   PbPb_to_pp->Add(hpp_phi,-1);

   TH1D *horiz_line = (TH1D*) hMC_sb_phi->Clone("horiz_line");
   //horiz_line->Divide(hMC_sb_phi);
   horiz_line->Add(hMC_sb_phi,-1);

   horiz_line->SetLineColor(kBlack);
   PbPb_to_pp->SetLineColor(kRed);

   PbPb_to_pp->SetXTitle("");

   if(selfmix)
      PbPb_to_pp->SetYTitle("Difference");
   else
      PbPb_to_pp->SetYTitle("Difference");
      //PbPb_to_pp->SetYTitle("(PbPb r-b) - pp");
   

   PbPb_to_pp->Draw("ep");
   horiz_line->Draw("hist same");

   c->SaveAs(Form("/eos/user/p/pchou/figs/%s/%s/Ztrack_%s_com_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f.png",typeofdata,HistName.c_str(),typeofdata1,ptL,ptH,centL,centH,TptL,TptH)); 
   
   if(drawrat){
      RPad->Clear();
      RPad->SetTopMargin(0);
      RPad->cd();

      TH1D *PbPb_rat_pp = (TH1D*) hMC_sb_phi->Clone("PbPb_rat_pp");
      PbPb_rat_pp->Divide(hpp_phi);
   
      TH1D *horir_line = (TH1D*) hMC_sb_phi->Clone("horir_line");
      horir_line->Divide(hMC_sb_phi);
   
      horir_line->SetLineColor(kBlack);
      PbPb_rat_pp->SetLineColor(kRed);
   
      PbPb_rat_pp->SetXTitle("");
   
      if(selfmix)
         PbPb_rat_pp->SetYTitle("Ratio");
      else
         PbPb_rat_pp->SetYTitle("Ratio");
         //PbPb_rat_pp->SetYTitle("(PbPb r-b) / pp");
      
   
      PbPb_rat_pp->Draw("ep");
      horir_line->Draw("hist same");
   
      c->SaveAs(Form("/eos/user/p/pchou/figs/%s/%s/Ztrack_%s_com_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f_rat.png",typeofdata,HistName.c_str(),typeofdata1,ptL,ptH,centL,centH,TptL,TptH)); 

   }
   if(drawlog){
      Pad->SetLogy();
   
      if(min1>0){
         hMC_phi->SetMinimum(0.5*min1);
         hMC_bkg_phi->SetMinimum(0.5*min1);
         hMC_sb_phi->SetMinimum(0.5*min1);
         hpp_phi->SetMinimum(0.5*min1);
      }else{
         hMC_phi->SetMinimum(0.01);
         hMC_bkg_phi->SetMinimum(0.01);
         hMC_sb_phi->SetMinimum(0.01);
         hpp_phi->SetMinimum(0.01);
      }
   
      hMC_phi->SetMaximum(1000*max1);
      hMC_bkg_phi->SetMaximum(1000*max1);
      hMC_sb_phi->SetMaximum(1000*max1);
      hpp_phi->SetMaximum(1000*max1);
   
      c->SaveAs(Form("/eos/user/p/pchou/figs/%s/%s/Ztrack_%s_comlog_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f_log.png",typeofdata,HistName.c_str(),typeofdata1,ptL,ptH,centL,centH,TptL,TptH)); 
      Pad->SetLogy(0);
   }
   
   
   if(drawlow){
      //hMC_phi->SetMaximum(1.2*max3);
      //hMC_bkg_phi->SetMaximum(1.2*max3);

      hMC_phi->SetMaximum(8);
      hMC_bkg_phi->SetMaximum(8);

      c->SaveAs(Form("/eos/user/p/pchou/figs/%s/%s/Ztrack_%s_com_%.0f_%.0f_%.0f_%.0f_%.0f_%.0f_low.png",typeofdata,HistName.c_str(),typeofdata1,ptL,ptH,centL,centH,TptL,TptH)); 
   }

   c->Clear(); 

}

void ZBasicBkgSub_loop(int binnum=40,float ptL=20,float ptH=2000,float centL=0,float centH=90,float TptL=0,float TptH=10000)
{

   double bineta = 0.032, binphi = M_PI/50;

   string gentxt="";
   if(isgen)
      gentxt="Gen";

   //string HistName[] = {Form("H%sZEta",gentxt.c_str()), Form("H%sZPhi",gentxt.c_str()), Form("H%sTrackEta",gentxt.c_str()), Form("H%sTrackPhi",gentxt.c_str()), Form("H%sEta",gentxt.c_str()), Form("H%sPhi",gentxt.c_str())};
   //string XTitleName[] = {"#eta_{Z}", "#phi_{Z}", "#eta_{track}", "#phi_{track}", "#Delta#eta_{Z,track}", "#Delta#phi_{Z,track}"};
   //string YTitleName[] = {"dN/d#eta", "dN/d#phi", "dN/d#eta", "dN/d#phi", "dN/d#Delta#eta", "dN/d#Delta#phi"};
   
   string HistName[] = {Form("H%sPhi",gentxt.c_str())};
   string XTitleName[] = {"#Delta#phi_{Z,track}"};
   string YTitleName[] = {"dN/d#Delta#phi"};
   
   int rebin_num[] = {2, 2, 2, 2, 1, 1};

   int i_draw = sizeof(HistName)/sizeof(HistName[0]);

   for(int i=0; i<i_draw; i++)
      ZBasicBkgSub_single(binnum, ptL, ptH, centL, centH, TptL, TptH, 
         HistName[i], XTitleName[i], YTitleName[i], rebin_num[i]);

}

int main(int argc, char *argv[]){

   style();
   string filebase = "/eos/home-p/pchou/BasicPlots/ZHadron2024/";

   if(isgen)
      file_sigMC = TFile::Open((filebase + "GraphMCSigGen_v1_ee.root").c_str(), "read");
   else
      file_sigMC = TFile::Open((filebase + "GraphMCSig_v1_ee_noZw.root").c_str(), "read");
   
   if(selfmix)
      file_bkgMC = TFile::Open("~/eos_base/BasicPlots/GraphMCSigBkg_v17d_PFmuon_350_10HF_ov20.root","read");
   else{
      if(isgen)
         file_bkgMC = TFile::Open((filebase + "GraphMCBkgGen_v1_ee.root").c_str(), "read");
      else
         file_bkgMC = TFile::Open((filebase + "GraphMCBkg_v1_ee_noZw.root").c_str(), "read");
   }

   if(isgen)
      file_ppMC  = TFile::Open((filebase + "GraphPPMCGen0NPU_v1_ee.root").c_str(), "read");
   else if(comparesub0)
      file_ppMC  = TFile::Open((filebase + "GraphMCSigGen0Sub_v1_ee.root").c_str(), "read");
   else
      file_ppMC  = TFile::Open((filebase + "GraphPPMC0NPU_v1_ee_noZw.root").c_str(), "read");

   if(selfmix)
      file_ppbkgMC  = TFile::Open("~/eos_base/BasicPlots/GraphPPMCSigBkg_v17d_PFmuon_143_10HF.root","read");
   else
      file_ppbkgMC  = TFile::Open("~/eos_base/BasicPlots/GraphPPMC_v18c_noW_v3.root","read");

   //file_ppbkgMC  = TFile::Open("~/eos_base/BasicPlots/GraphPPMCSigBkg_v17d_Cent10_143.root","read");


   ZBasicBkgSub_loop(40, 40, 200,  0, 100,  1,  2);
   //ZBasicBkgSub_loop(40, 40, 200,  0, 100,  2,  4);
   //ZBasicBkgSub_loop(40, 40, 200,  0, 100,  4, 10);
   ZBasicBkgSub_loop(40, 40, 200,  0,  10,  1,  2);
   //ZBasicBkgSub_loop(40, 40, 200,  0,  10,  2,  4);
   ZBasicBkgSub_loop(40, 40, 200,  0,  10,  4, 10);

   //ZBasicBkgSub_loop(40, 40, 200, 10,  30,  1,  2);
   //ZBasicBkgSub_loop(40, 40, 200, 30,  50,  1,  2);
   ZBasicBkgSub_loop(40, 40, 200, 50,  90,  1,  2);

   file_sigMC->Close();
   file_bkgMC->Close();
   file_ppMC->Close();

   if(selfmix)
      file_ppbkgMC->Close();

   //gSystem->Exec(Form("convert +append /eos/user/p/pchou/figs/track/%s/BasicBkgSub/*Eta/*40_200_0_10_1_2.png /eos/user/p/pchou/figs/track/%s/BasicBkgSub/Eta_merge.png",typeofdata,typeofdata));
   //gSystem->Exec(Form("convert +append /eos/user/p/pchou/figs/track/%s/BasicBkgSub/*Phi/*40_200_0_10_1_2.png /eos/user/p/pchou/figs/track/%s/BasicBkgSub/Phi_merge.png",typeofdata,typeofdata));
   //gSystem->Exec(Form("convert -append /eos/user/p/pchou/figs/track/%s/BasicBkgSub/*_merge.png /eos/user/p/pchou/figs/track/%s/BasicBkgSub/merge.png",typeofdata,typeofdata));

   return 0;

}
