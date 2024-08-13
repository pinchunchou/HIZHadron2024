#include <iostream>

#include "TFile.h"
#include "TH1D.h"

class TrackResidualCorrector
{
public:
   TrackResidualCorrector(std::string filename = "ResidualCorrector_20230508/totalCorrection.root")
   {
      f = new TFile(filename.c_str());
      hPtCorrTotal  = (TH1D *)f->Get("hPtCorrTotal");
      hEtaCorrTotal = (TH1D *)f->Get("hEtaCorrTotal");
      hPhiCorrTotal = (TH1D *)f->Get("hPhiCorrTotal");
   }

   ~TrackResidualCorrector()
   {
      f->Close();
      delete f;
   }

   double GetCorrectionFactor(double pt, double eta, double phi, int type = 0)
   {
      if(phi < 0)
         phi += 2 * M_PI;

      double PTMax = hPtCorrTotal->GetXaxis()->GetBinUpEdge(hPtCorrTotal->GetNbinsX());

      int bin_pt = hPtCorrTotal->GetXaxis()->FindBin(pt);
      int bin_eta = hEtaCorrTotal->GetXaxis()->FindBin(eta);
      int bin_phi = hPhiCorrTotal->GetXaxis()->FindBin(phi);
     
      // last bin in PT is overflow
      if(pt >= PTMax)
         bin_pt = hPtCorrTotal->GetNbinsX();
      

      double ptCorr, etaCorr, phiCorr;

      // Handle different types for PT correction
      switch(type) {
        case 1:
            ptCorr = hPtCorrTotal->GetBinContent(bin_pt) + hPtCorrTotal->GetBinError(bin_pt);
            etaCorr = hEtaCorrTotal->GetBinContent(bin_eta) + hEtaCorrTotal->GetBinError(bin_eta);
            phiCorr = hPhiCorrTotal->GetBinContent(bin_phi) + hPhiCorrTotal->GetBinError(bin_phi);
            break;
        case -1:
            ptCorr = hPtCorrTotal->GetBinContent(bin_pt) - hPtCorrTotal->GetBinError(bin_pt);
            etaCorr = hEtaCorrTotal->GetBinContent(bin_eta) - hEtaCorrTotal->GetBinError(bin_eta);
            phiCorr = hPhiCorrTotal->GetBinContent(bin_phi) - hPhiCorrTotal->GetBinError(bin_phi);
            break;
        default:
            ptCorr = hPtCorrTotal->GetBinContent(bin_pt);
            etaCorr = hEtaCorrTotal->GetBinContent(bin_eta);
            phiCorr = hPhiCorrTotal->GetBinContent(bin_phi);
      }
  
      double corr = ptCorr * etaCorr * phiCorr;

      if(isnan(corr))
      {
         std::cerr << "Error!  nan efficiency! " << bin_pt << " " << bin_eta << " " << bin_phi << std::endl;
         corr = 1;
      }

      return corr;
   }

private:
   TFile* f;
   TH1D *hPtCorrTotal;
   TH1D *hEtaCorrTotal;
   TH1D *hPhiCorrTotal;
};

class TrackResidualCentralityCorrector
{
public:
   TrackResidualCentralityCorrector(std::string F1, std::string F2, std::string F3, std::string F4)
   {
      TRC1 = new TrackResidualCorrector(F1);
      TRC2 = new TrackResidualCorrector(F2);
      TRC3 = new TrackResidualCorrector(F3);
      TRC4 = new TrackResidualCorrector(F4);
   }
   
   TrackResidualCentralityCorrector(std::vector<std::string> F)
   {
      // cout << F.size() << " " << F[0] << endl;
      if(F.size() == 4)
      {
         TRC1 = new TrackResidualCorrector(F[0]);
         TRC2 = new TrackResidualCorrector(F[1]);
         TRC3 = new TrackResidualCorrector(F[2]);
         TRC4 = new TrackResidualCorrector(F[3]);
      }
      else if(F.size() == 1)
      {
         TRC1 = new TrackResidualCorrector(F[0]);
         TRC2 = new TrackResidualCorrector(F[0]);
         TRC3 = new TrackResidualCorrector(F[0]);
         TRC4 = new TrackResidualCorrector(F[0]);
      }
      else
      {
         TRC1 = nullptr;
         TRC2 = nullptr;
         TRC3 = nullptr;
         TRC4 = nullptr;
      }
   }

   ~TrackResidualCentralityCorrector()
   {
      if(TRC1 != nullptr)   delete TRC1;
      if(TRC2 != nullptr)   delete TRC2;
      if(TRC3 != nullptr)   delete TRC3;
      if(TRC4 != nullptr)   delete TRC4;
   }

   double GetCorrectionFactor(double pt, double eta, double phi, int hiBin, int type = 0)
   {
      if(hiBin < 20)         return TRC1->GetCorrectionFactor(pt, eta, phi, type);
      else if(hiBin < 60)    return TRC2->GetCorrectionFactor(pt, eta, phi, type);
      else if(hiBin < 100)   return TRC3->GetCorrectionFactor(pt, eta, phi, type);
      else                   return TRC4->GetCorrectionFactor(pt, eta, phi, type);
      return 0;
   }

private:
   TrackResidualCorrector *TRC1;
   TrackResidualCorrector *TRC2;
   TrackResidualCorrector *TRC3;
   TrackResidualCorrector *TRC4;
};


