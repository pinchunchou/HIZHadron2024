#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <glob.h>

using namespace std;

#include "TH1D.h"
#include "TH3D.h"
#include "TFile.h"
#include "TTree.h"

#include "CommandLine.h"
#include "CommonFunctions.h"

#include "PlotHelper4.h"

#define MAX 10000

std::vector<std::string> my_glob(const char *pattern);

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
void CalculateEfficiencyAndError(TH1D* hNum, TH1D* hDen, TH1D* hResult) {
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
}

int main(int argc, char *argv[])
{
	CommandLine CL(argc, argv);

	string InputFileNames 		  = CL.Get("Input", "PbPb/*.root");
	string OutputFileName         = CL.Get("Output", "TrackEfficiencyPlots.pdf");
   	string RootOutputFileName     = CL.Get("RootOutput", "TrkEffPbPb.root");
   	vector<string> TrackResidualPaths   = CL.GetStringVector("TrackResidualPath",vector<string>{"", "", "", ""});
	bool DoIteration              = CL.GetBool("DoIteration", false);
	int TrackResIdx 			  = CL.GetInt("TrackResIdx", 0);

	string TrackResidualPath 	  = TrackResidualPaths[TrackResIdx];

	std::cout<<"InputFileNames = "<<InputFileNames<<std::endl;
	std::cout<<"OutputFileName = "<<OutputFileName<<std::endl;
	std::cout<<"RootOutputFileName = "<<RootOutputFileName<<std::endl;
	std::cout<<"TrackResidualPath = "<<TrackResidualPath<<std::endl;

	PdfFileHelper PdfFile(OutputFileName);
   	PdfFile.AddTextPage("Track efficiency derivation");

	TH3D *HRecoTrack = nullptr;
	TH3D *HRecoTrackCorrected = nullptr;
   	TH3D *HRecoTrackResCorrected = nullptr;
   	TH3D *HGenTrack = nullptr;

   	TFile* ResFile = TFile::Open(TrackResidualPath.c_str(), "READ");
    TH1D *hPtCorrTotal_old = (TH1D*)ResFile->Get("hPtCorrTotal");
    TH1D *hEtaCorrTotal_old = (TH1D*)ResFile->Get("hEtaCorrTotal");
    TH1D *hPhiCorrTotal_old = (TH1D*)ResFile->Get("hPhiCorrTotal");

    hPtCorrTotal_old->SetName("hPtCorrTotal_old");
    hEtaCorrTotal_old->SetName("hEtaCorrTotal_old");
    hPhiCorrTotal_old->SetName("hPhiCorrTotal_old");

   std::vector<double> PTs(49);
   std::generate_n(PTs.begin(), 49, [n = 0.2]() mutable { return n += 0.2; });

   std::vector<double> Etas(51);
   std::generate_n(Etas.begin(), 51, [n = -2.4-0.096]() mutable { return n += 0.096; });

   std::vector<double> Phis(51);
   std::generate_n(Phis.begin(), 51, [n = -0.04*M_PI]() mutable { return n += 0.04*M_PI; });

   sort(Etas.begin(), Etas.end());
   sort(Phis.begin(), Phis.end());
   sort(PTs.begin(), PTs.end());

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

   	
    TH1D *hPtCorrTotal = new TH1D("hPtCorrTotal", "p_{T} Correction Factor", NPT, PTBins);
   TH1D *hEtaCorrTotal = new TH1D("hEtaCorrTotal", "#eta Correction Factor", NEta, EtaBins);
   TH1D *hPhiCorrTotal = new TH1D("hPhiCorrTotal", "#phi Correction Factor", NPhi, PhiBins);

   	



    for (const auto &InputFileName : my_glob(InputFileNames.c_str())){
    	std::cout << "Processing " << InputFileName << std::endl;


        TFile* input = TFile::Open(InputFileName.c_str(), "READ");        
        if (!input) continue;

        TH3D* hReco = (TH3D*)input->Get("HRecoTrack");
        TH3D* hRecoCorr = (TH3D*)input->Get("HRecoTrackCorrected");
        TH3D* hRecoResCorr = (TH3D*)input->Get("HRecoTrackResCorrected");
        TH3D* hGen = (TH3D*)input->Get("HGenTrack");


        if (!hReco || !hRecoCorr || !hRecoResCorr || !hGen) {
            input->Close();
            continue;
        }

        if (!HRecoTrack) {
            HRecoTrack = (TH3D*)hReco->Clone();
            HRecoTrack->SetDirectory(0);
        }else{
            HRecoTrack->Add(hReco);
        }

        if (!HRecoTrackCorrected) {
            HRecoTrackCorrected = (TH3D*)hRecoCorr->Clone();
            HRecoTrackCorrected->SetDirectory(0);
        }else{
            HRecoTrackCorrected->Add(hRecoCorr);
        }

        if (!HRecoTrackResCorrected) {
			HRecoTrackResCorrected = (TH3D*)hRecoResCorr->Clone();
			HRecoTrackResCorrected->SetDirectory(0);
		}else{
			HRecoTrackResCorrected->Add(hRecoResCorr);
		}

		if (!HGenTrack) {
			HGenTrack = (TH3D*)hGen->Clone();
			HGenTrack->SetDirectory(0);
		}else{
			HGenTrack->Add(hGen);
		}

        input->Close();
    }

    if (HRecoTrack && HRecoTrackCorrected && HRecoTrackResCorrected && HGenTrack) {


    	// Calculate efficiencies and errors
    	TH1D* HEfficiency_eta = new TH1D("HEfficiency_eta", "#eta Efficiency", NEta, EtaBins);
    	TH1D* HEfficiency_phi = new TH1D("HEfficiency_phi", "#phi Efficiency", NPhi, PhiBins);
    	TH1D* HEfficiency_pt = new TH1D("HEfficiency_pt", "p_{T} Efficiency", NPT, PTBins);
	
    	CalculateEfficiencyAndError(ProjectX(HRecoTrack), ProjectX(HGenTrack), HEfficiency_eta);
    	CalculateEfficiencyAndError(ProjectY(HRecoTrack), ProjectY(HGenTrack), HEfficiency_phi);
    	CalculateEfficiencyAndError(ProjectZ(HRecoTrack), ProjectZ(HGenTrack), HEfficiency_pt);

    	TH1D* HEfficiencyCorrected_eta = new TH1D("HEfficiencyCorrected_eta", "#eta Corrected Efficiency", NEta, EtaBins);
    	TH1D* HEfficiencyCorrected_phi = new TH1D("HEfficiencyCorrected_phi", "#phi Corrected Efficiency", NPhi, PhiBins);
    	TH1D* HEfficiencyCorrected_pt = new TH1D("HEfficiencyCorrected_pt", "p_{T} Corrected Efficiency", NPT, PTBins);

		CalculateEfficiencyAndError(ProjectX(HRecoTrackCorrected), ProjectX(HGenTrack), HEfficiencyCorrected_eta);
		CalculateEfficiencyAndError(ProjectY(HRecoTrackCorrected), ProjectY(HGenTrack), HEfficiencyCorrected_phi);
		CalculateEfficiencyAndError(ProjectZ(HRecoTrackCorrected), ProjectZ(HGenTrack), HEfficiencyCorrected_pt);

		TH1D* HEfficiencyResCorrected_eta = new TH1D("HEfficiencyResCorrected_eta", "#eta Residual Corrected Efficiency", NEta, EtaBins);
		TH1D* HEfficiencyResCorrected_phi = new TH1D("HEfficiencyResCorrected_phi", "#phi Residual Corrected Efficiency", NPhi, PhiBins);
		TH1D* HEfficiencyResCorrected_pt = new TH1D("HEfficiencyResCorrected_pt", "p_{T} Residual Corrected Efficiency", NPT, PTBins);

		CalculateEfficiencyAndError(ProjectX(HRecoTrackResCorrected), ProjectX(HGenTrack), HEfficiencyResCorrected_eta);
		CalculateEfficiencyAndError(ProjectY(HRecoTrackResCorrected), ProjectY(HGenTrack), HEfficiencyResCorrected_phi);
		CalculateEfficiencyAndError(ProjectZ(HRecoTrackResCorrected), ProjectZ(HGenTrack), HEfficiencyResCorrected_pt);


   		TH3D *HEfficiencyResCorrected = (TH3D*)HRecoTrackResCorrected->Clone("HEfficiencyResCorrected");
   		HEfficiencyResCorrected->Divide(HGenTrack);
   		 // Calculate correction factors  

   		double overallEfficiency = HRecoTrackResCorrected->Integral() / HGenTrack->Integral();
		std::cout << "Overall efficiency: " << overallEfficiency << std::endl;
		std::cout << "Overall efficiency with overflow: " << HRecoTrackResCorrected->Integral(0,-1,0,-1,0,-1) / HGenTrack->Integral(0,-1,0,-1,0,-1) << std::endl;


		for (int i = 1; i <= HEfficiencyResCorrected_eta->GetNbinsX(); ++i) {
			bool DoEtaIteration = DoIteration && HEfficiencyResCorrected_eta->GetNbinsX() == hEtaCorrTotal_old->GetNbinsX();
    	    double eff = DoEtaIteration ? HEfficiencyResCorrected_eta->GetBinContent(i)/hEtaCorrTotal_old->GetBinContent(i) : HEfficiencyCorrected_eta->GetBinContent(i);
    	    double err = (ProjectX(HGenTrack)->GetBinContent(i) > 0 && eff > 0 && eff < 1) ? sqrt(eff * (1-eff) / (ProjectX(HGenTrack)->GetBinContent(i))) : 0;
    	    hEtaCorrTotal->SetBinContent(i, eff > 0 ? 1.0 / eff : 1.0);
    	    hEtaCorrTotal->SetBinError(i, eff > 0 ? err / (eff * eff) : 0.0);
    	}

   		
   		for(int i = 1; i <= HEfficiencyResCorrected_phi->GetNbinsX(); i++){
   			bool DoPhiIteration = DoIteration && HEfficiencyResCorrected_phi->GetNbinsX() == hPhiCorrTotal_old->GetNbinsX();
   		    double eff = DoPhiIteration ? HEfficiencyResCorrected_phi->GetBinContent(i)/hPhiCorrTotal_old->GetBinContent(i) : HEfficiencyCorrected_phi->GetBinContent(i);
   		    double err = (ProjectY(HGenTrack)->GetBinContent(i) > 0 && eff > 0 && eff < 1) ? sqrt(eff * (1-eff) / (ProjectY(HGenTrack)->GetBinContent(i))) : 0;
   		    hPhiCorrTotal->SetBinContent(i, eff > 0 ? 1.0 / eff : 1.0);
   		    hPhiCorrTotal->SetBinError(i, eff > 0 ? err / (eff * eff) : 0.0);
   		}


   		for(int i = 1; i <= HEfficiencyResCorrected_pt->GetNbinsX(); i++){
   			bool DoPtIteration = DoIteration && HEfficiencyResCorrected_pt->GetNbinsX() == hPtCorrTotal_old->GetNbinsX();
   		    double eff = DoPtIteration ? HEfficiencyResCorrected_pt->GetBinContent(i)/hPtCorrTotal_old->GetBinContent(i) : HEfficiencyCorrected_pt->GetBinContent(i);
   		    double err = (ProjectZ(HGenTrack)->GetBinContent(i) > 0 && eff > 0 && eff < 1) ? sqrt(eff * (1-eff) / (ProjectZ(HGenTrack)->GetBinContent(i))) : 0;
   		    hPtCorrTotal->SetBinContent(i, eff > 0 ? 1.0 / eff : 1.0);
   		    hPtCorrTotal->SetBinError(i, eff > 0 ? err / (eff * eff) : 0.0);
   		}

   		//double ptIntegral_i = HEfficiencyResCorrected_pt->Integral();
   		//double etaIntegral_i = HEfficiencyResCorrected_eta->Integral();
   		//double phiIntegral_i = HEfficiencyResCorrected_phi->Integral();
		
   		double ptIntegral = hPtCorrTotal->Integral();
   		double etaIntegral = hEtaCorrTotal->Integral();
   		double phiIntegral = hPhiCorrTotal->Integral();

   		std::cout<<"ptIntegral = "<<ptIntegral<<std::endl;
   		std::cout<<"etaIntegral = "<<etaIntegral<<std::endl;
   		std::cout<<"phiIntegral = "<<phiIntegral<<std::endl;

   		std::cout<<"ptIntegral overflow = "<<hPtCorrTotal->Integral(1,-1)<<std::endl;
   		std::cout<<"etaIntegral overflow = "<<hEtaCorrTotal->Integral(1,-1)<<std::endl;
   		std::cout<<"phiIntegral overflow = "<<hPhiCorrTotal->Integral(1,-1)<<std::endl;


   		std::cout<<"hEtaCorrTotal->GetNbinsX() = "<<hEtaCorrTotal->GetNbinsX()<<std::endl;
   		std::cout<<"hPhiCorrTotal->GetNbinsX() = "<<hPhiCorrTotal->GetNbinsX()<<std::endl;
   		std::cout<<"hPtCorrTotal->GetNbinsX() = "<<hPtCorrTotal->GetNbinsX()<<std::endl;


   		double recoEffIntegral = HRecoTrackCorrected->Integral();
   		double recoResIntegral = HRecoTrackResCorrected->Integral();
   		double genIntegral = HGenTrack->Integral();

   		double CorrTotal = DoIteration ? genIntegral/recoResIntegral : genIntegral/recoEffIntegral ;

   		std::cout<<"CorrTotal = "<<CorrTotal<<std::endl;

   		//CorrTotal = 1;
/*
now
w1 /1
w2 /w2
w3 /w3

w1x = w1 * pow(w1*w2*w3,-2/9)
w2x = w2 * pow(w1*w2*w3,-2/9)
w3x = w3 * pow(w1*w2*w3,-2/9)

w1*w2*w3*pow(w1*w2*w3,-2/3) = pow(w1*w2*w3,1/3)
replace w2 to pow(w1*w2*w3, 2/9)


w1 * pow(w2*w3,x)
w2 * pow(w1*w3,x)
w3 * pow(w1*w2,x)


w1*w2*w3*pow(w1*w2*w3,2x) = pow(w1*w2*w3,1/3)

1+2x = 1/3
2x = -2/3
x = -1/3

w2


*/
   		CorrTotal = (hPtCorrTotal->GetNbinsX() / ptIntegral) * (hEtaCorrTotal->GetNbinsX() / etaIntegral) * (hPhiCorrTotal->GetNbinsX() / phiIntegral);
   		std::cout<<"CorrTotal = "<<1/CorrTotal<<std::endl;

   		//hPtCorrTotal ->Scale(pow(CorrTotal,2/9.));
   		//hEtaCorrTotal->Scale(pow(CorrTotal,2/9.));
   		//hPhiCorrTotal->Scale(pow(CorrTotal,2/9.));

   		//hPtCorrTotal ->Scale(pow((hEtaCorrTotal->GetNbinsX() / etaIntegral) * (hPhiCorrTotal->GetNbinsX() / phiIntegral),1/3.));
   		//hEtaCorrTotal->Scale(pow((hPtCorrTotal->GetNbinsX() / ptIntegral) * (hPhiCorrTotal->GetNbinsX() / phiIntegral),1/3.));
   		//hPhiCorrTotal->Scale(pow((hPtCorrTotal->GetNbinsX() / ptIntegral) * (hEtaCorrTotal->GetNbinsX() / etaIntegral),1/3.));

   		hEtaCorrTotal->Scale(hEtaCorrTotal->GetNbinsX() / etaIntegral);
   		hPhiCorrTotal->Scale(hPhiCorrTotal->GetNbinsX() / phiIntegral);
   		//hPtCorrTotal->Scale(hPtCorrTotal->GetNbinsX() / ptIntegral);


   		//HEfficiencyResCorrected_eta->Scale(HEfficiencyResCorrected_eta->GetNbinsX()/HEfficiencyResCorrected_eta->Integral());
   		//HEfficiencyResCorrected_phi->Scale(HEfficiencyResCorrected_phi->GetNbinsX()/HEfficiencyResCorrected_phi->Integral());
   		//HEfficiencyResCorrected_pt->Scale(HEfficiencyResCorrected_pt->GetNbinsX()/HEfficiencyResCorrected_pt->Integral());


   		PdfFile.AddTextPage("Track efficiency plots: #eta");

   		PdfFile.AddPlot(HGenTrack->ProjectionX(), "hist text00", true);
   		PdfFile.AddPlot(HRecoTrack->ProjectionX(), "hist text00", true);
   		PdfFile.AddPlot(HEfficiency_eta, "hist e text00", false);
   		PdfFile.AddPlot(HEfficiencyCorrected_eta, "hist e text00", false);
   		PdfFile.AddPlot(HEfficiencyResCorrected_eta, "hist e text00", false);
		
   		PdfFile.AddTextPage("Track efficiency plots: #phi");
		
   		PdfFile.AddPlot(HGenTrack->ProjectionY(), "hist text00", true);
   		PdfFile.AddPlot(HRecoTrack->ProjectionY(), "hist text00", true);
   		PdfFile.AddPlot(HEfficiency_phi, "hist e text00", false);
   		PdfFile.AddPlot(HEfficiencyCorrected_phi, "hist e text00", false);
   		PdfFile.AddPlot(HEfficiencyResCorrected_phi, "hist e text00", false);
		
   		PdfFile.AddTextPage("Track efficiency plots: p_{T}");
		
   		PdfFile.AddPlot(HGenTrack->ProjectionZ(), "hist text00", true);
   		PdfFile.AddPlot(HRecoTrack->ProjectionZ(), "hist text00", true);
   		PdfFile.AddPlot(HEfficiency_pt, "hist e text00", false);
   		PdfFile.AddPlot(HEfficiencyCorrected_pt, "hist e text00", false);
   		PdfFile.AddPlot(HEfficiencyResCorrected_pt, "hist e text00", false);
		
   		PdfFile.AddTimeStampPage();
   		PdfFile.Close();

        TFile* output = TFile::Open(RootOutputFileName.c_str(), "RECREATE");
        
        output->cd();

        //hPtCorrTotal->SetName("hPtCorrTotal");
    	//hEtaCorrTotal->SetName("hEtaCorrTotal");
    	//hPhiCorrTotal->SetName("hPhiCorrTotal");
//
    	//hPtCorrTotal->SetTitle("hPtCorrTotal");
    	//hEtaCorrTotal->SetTitle("hEtaCorrTotal");
    	//hPhiCorrTotal->SetTitle("hPhiCorrTotal");
        
        hPtCorrTotal->Write();
   		hEtaCorrTotal->Write();
   		hPhiCorrTotal->Write();
        
        output->Close();
        delete HRecoTrack;
        delete HRecoTrackCorrected;
        delete HRecoTrackResCorrected;
        delete HGenTrack;
    }
}

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
