#include <iostream>
using namespace std;

#include "TFile.h"
#include "TTree.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TF1.h"

#include "SetStyle.h"
#include "PlotHelper4.h"

int main();
void MakePlot(double XMin, double XMax, double YMin, double YMax, string Title,
   TProfile &P1, TProfile &P2, bool LogX, string Output);

int main()
{
   SetThesisStyle();

   TFile File("ZEfficiencyAll.root");
   TTree *Tree = (TTree *)File.Get("Tree");

   double PTBins[41] = {0};
   for(int i = 0; i <= 40; i++)
      PTBins[i] = exp(log(1) + (log(200) - log(1)) / 40 * i);

   Tree->SetAlias("TnP",          "(Mu1TnP1*Mu1TnP2*Mu1TnP3*Mu2TnP1*Mu2TnP2*Mu2TnP3)");
   Tree->SetAlias("AATnP_Sys1P",  "(Mu1TnP1P1*Mu1TnP2*Mu1TnP3*Mu2TnP1P1*Mu2TnP2*Mu2TnP3)");
   Tree->SetAlias("AATnP_Sys1N",  "(Mu1TnP1P2*Mu1TnP2*Mu1TnP3*Mu2TnP1P2*Mu2TnP2*Mu2TnP3)");
   Tree->SetAlias("AATnP_Sys2P",  "(Mu1TnP1*Mu1TnP2P1*Mu1TnP3*Mu2TnP1*Mu2TnP2P1*Mu2TnP3)");
   Tree->SetAlias("AATnP_Sys2N",  "(Mu1TnP1*Mu1TnP2P2*Mu1TnP3*Mu2TnP1*Mu2TnP2P2*Mu2TnP3)");
   Tree->SetAlias("AATnP_Sys3P",  "(Mu1TnP1*Mu1TnP2*Mu1TnP3P1*Mu2TnP1*Mu2TnP2*Mu2TnP3P1)");
   Tree->SetAlias("AATnP_Sys3N",  "(Mu1TnP1*Mu1TnP2*Mu1TnP3P2*Mu2TnP1*Mu2TnP2*Mu2TnP3P2)");
   Tree->SetAlias("AATnP_Stat1P", "(Mu1TnP1N1*Mu1TnP2*Mu1TnP3*Mu2TnP1N1*Mu2TnP2*Mu2TnP3)");
   Tree->SetAlias("AATnP_Stat1N", "(Mu1TnP1N2*Mu1TnP2*Mu1TnP3*Mu2TnP1N2*Mu2TnP2*Mu2TnP3)");
   Tree->SetAlias("AATnP_Stat2P", "(Mu1TnP1*Mu1TnP2N1*Mu1TnP3*Mu2TnP1*Mu2TnP2N1*Mu2TnP3)");
   Tree->SetAlias("AATnP_Stat2N", "(Mu1TnP1*Mu1TnP2N2*Mu1TnP3*Mu2TnP1*Mu2TnP2N2*Mu2TnP3)");
   Tree->SetAlias("AATnP_Stat3P", "(Mu1TnP1*Mu1TnP2*Mu1TnP3N1*Mu2TnP1*Mu2TnP2*Mu2TnP3N1)");
   Tree->SetAlias("AATnP_Stat3N", "(Mu1TnP1*Mu1TnP2*Mu1TnP3N2*Mu2TnP1*Mu2TnP2*Mu2TnP3N2)");

   // vector<string> Variations
   // {
   //    "AATnP_Sys1P", "AATnP_Sys1N",
   //    "AATnP_Sys2P", "AATnP_Sys2N",
   //    "AATnP_Sys3P", "AATnP_Sys3N",
   //    "AATnP_Stat1P", "AATnP_Stat1N",
   //    "AATnP_Stat2P", "AATnP_Stat2N",
   //    "AATnP_Stat3P", "AATnP_Stat3N"
   // };
   vector<string> Variations
   {
      "AATnP_Sys1P", "AATnP_Sys1N",
   };

   PdfFileHelper PdfFile("Variations.pdf");
   PdfFile.AddTextPage("Tag and probe variations");

   vector<int> HiBinBoundary{0, 20, 40, 80, 200};

   for(string Variation : Variations)
   {
      PdfFile.AddTextPage(Variation);

      TProfile HY    ("HY",     "Raw;;", 40, -2.4, 2.4);
      TProfile HPT   ("HPT",    "Raw;;", 40, PTBins);
      TProfile HHiBin("HHiBin", "Raw;;", 40, 0, 200);

      Tree->Draw(Form("%s/TnP:Y>>HY", Variation.c_str()),         "PT > 1", "prof");
      Tree->Draw(Form("%s/TnP:PT>>HPT", Variation.c_str()),       "PT > 1", "prof");
      Tree->Draw(Form("%s/TnP:HiBin>>HHiBin", Variation.c_str()), "PT > 1", "prof");
      
      PdfFile.AddPlot(HY, "");
      PdfFile.AddPlot(HPT, "", false, false, false, true);
      PdfFile.AddPlot(HHiBin, "");
      
      // Now we do things in the four centrality bins
      for(int iC = 0; iC < 4; iC++)
      {
         TProfile HCY ("HCY",  Form("Raw C%d;;", iC), 40, -2.4, 2.4);
         TProfile HCPT("HCPT", Form("Raw C%d;;", iC), 40, PTBins);

         string HiBinCut = Form("HiBin >= %d && HiBin < %d", HiBinBoundary[iC], HiBinBoundary[iC+1]);
         Tree->Draw(Form("%s/TnP:Y>>HCY", Variation.c_str()),   Form("PT > 1 && %s", HiBinCut.c_str()), "prof");
         Tree->Draw(Form("%s/TnP:PT>>HCPT", Variation.c_str()), Form("PT > 1 && %s", HiBinCut.c_str()), "prof");
         PdfFile.AddPlot(HCY);
         PdfFile.AddPlot(HCPT, "", false, false, false, true);
        
         TF1 FY("FY", "pol10", -3, 3);
         HCY.Fit(&FY);
         PdfFile.AddPlot(HCY);

         string YWeight = "0";
         for(int i = 10; i >= 0; i--)
            YWeight = "(" + YWeight + ")*Y+" + FY.GetParameter(i);
         YWeight = "(" + YWeight + ")";
         
         TProfile HCYP1 ("HCYP1",  Form("Pass 1 C%d;;", iC), 40, -2.4, 2.4);
         TProfile HCPTP1("HCPTP1", Form("Pass 1 C%d;;", iC), 40, PTBins);

         Tree->Draw(Form("%s/TnP/%s:Y>>HCYP1", Variation.c_str(), YWeight.c_str()),
            Form("PT > 1 && %s", HiBinCut.c_str()), "prof");
         Tree->Draw(Form("%s/TnP/%s:PT>>HCPTP1", Variation.c_str(), YWeight.c_str()),
            Form("PT > 1 && %s", HiBinCut.c_str()), "prof");
         PdfFile.AddPlot(HCYP1);
         PdfFile.AddPlot(HCPTP1, "", false, false, false, true);
      }
   }

   PdfFile.AddTimeStampPage();
   PdfFile.Close();

   File.Close();

   return 0;
}

