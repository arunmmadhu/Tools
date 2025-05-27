#include "RooGlobalFunc.h"
#include "RooPolynomial.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooAbsPdf.h"
#include "RooFitResult.h"
//#include "RooMCStudy.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "TH1F.h"
#include "TH2F.h"
#include <cmath>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <iostream>
#include <algorithm>

using namespace RooFit;

void bbcc_norm_main () 
{
    gROOT->SetBatch(true);
    
    TFile *TreeFile_h = new TFile("LOCAL_COMBINED_ztau3mutauh_skimmer_default.root","READ");//tauh,mu,e doesn't matter in skimmer
    
    TH1D * nStop_values_Type1_bb  = (TH1D*)TreeFile_h->Get("ztau3mutauh_skimmer_default_Selection_bb_nStop_DistributionMC8");
    TH1D * nStop_values_Type3_bb  = (TH1D*)TreeFile_h->Get("ztau3mutauh_skimmer_default_Selection_bb_nStop_DistributionMC9");
    TH1D * nStop_values_Type5_bb  = (TH1D*)TreeFile_h->Get("ztau3mutauh_skimmer_default_Selection_bb_nStop_DistributionMC10");
    TH1D * nStop_values_Type9_bb  = (TH1D*)TreeFile_h->Get("ztau3mutauh_skimmer_default_Selection_bb_nStop_DistributionMC11");
    
    TH1D * nStop_values_Type1_cc  = (TH1D*)TreeFile_h->Get("ztau3mutauh_skimmer_default_Selection_cc_nStop_DistributionMC8");
    TH1D * nStop_values_Type3_cc  = (TH1D*)TreeFile_h->Get("ztau3mutauh_skimmer_default_Selection_cc_nStop_DistributionMC9");
    TH1D * nStop_values_Type5_cc  = (TH1D*)TreeFile_h->Get("ztau3mutauh_skimmer_default_Selection_cc_nStop_DistributionMC10");
    TH1D * nStop_values_Type9_cc  = (TH1D*)TreeFile_h->Get("ztau3mutauh_skimmer_default_Selection_cc_nStop_DistributionMC11");
    
    //Getting the mean
    double mean_type1_bb = nStop_values_Type1_bb->GetMean();
    double mean_type3_bb = nStop_values_Type3_bb->GetMean();
    double mean_type5_bb = nStop_values_Type5_bb->GetMean();
    double mean_type9_bb = nStop_values_Type9_bb->GetMean();
    
    std::cout << "The mean for type 1 bbbar is: "<< mean_type1_bb << std::endl;
    std::cout << "The mean for type 3 bbbar is: "<< mean_type3_bb << std::endl;
    std::cout << "The mean for type 5 bbbar is: "<< mean_type5_bb << std::endl;
    std::cout << "The mean for type 9 bbbar is: "<< mean_type9_bb << std::endl;
    
    //Getting the mean
    double mean_type1_cc = nStop_values_Type1_cc->GetMean();
    double mean_type3_cc = nStop_values_Type3_cc->GetMean();
    double mean_type5_cc = nStop_values_Type5_cc->GetMean();
    double mean_type9_cc = nStop_values_Type9_cc->GetMean();
    
    std::cout << "The mean for type 1 ccbar is: "<< mean_type1_cc << std::endl;
    std::cout << "The mean for type 3 ccbar is: "<< mean_type3_cc << std::endl;
    std::cout << "The mean for type 5 ccbar is: "<< mean_type5_cc << std::endl;
    std::cout << "The mean for type 9 ccbar is: "<< mean_type9_cc << std::endl;
    
    std::cout << "Bin content for type 1 bbbar is: "<< nStop_values_Type1_bb->Integral() << " and for ccbar is: " << nStop_values_Type1_cc->Integral() << " with k1: " << (43861.) / (89978.)*(nStop_values_Type1_bb->Integral()/(nStop_values_Type1_bb->Integral()+nStop_values_Type1_cc->Integral()))*(9.517/4.596) << " and k2: " << (43861.0) / (89978.0)*(nStop_values_Type1_cc->Integral()/(nStop_values_Type1_bb->Integral()+nStop_values_Type1_cc->Integral()))*(9.517/4.968) << std::endl;
    
    std::cout << "Bin content for type 3 bbbar is: "<< nStop_values_Type3_bb->Integral() << " and for ccbar is: " << nStop_values_Type3_cc->Integral() << " with k1: " << (35040.) / (699758.)*(nStop_values_Type3_bb->Integral()/(nStop_values_Type3_bb->Integral()+nStop_values_Type3_cc->Integral()))*(9.517/4.596) << " and k2: " << (35040.) / (699758.)*(nStop_values_Type3_cc->Integral()/(nStop_values_Type3_bb->Integral()+nStop_values_Type3_cc->Integral()))*(9.517/4.968) << std::endl;
    
    std::cout << "Bin content for type 5 bbbar is: "<< nStop_values_Type5_bb->Integral() << " and for ccbar is: " << nStop_values_Type5_cc->Integral() << " with k1: " << (21935.) / (839656.)*(nStop_values_Type5_bb->Integral()/(nStop_values_Type5_bb->Integral()+nStop_values_Type5_cc->Integral()))*(9.517/4.596) << " and k2: " << (21935.) / (839656.)*(nStop_values_Type5_cc->Integral()/(nStop_values_Type5_bb->Integral()+nStop_values_Type5_cc->Integral()))*(9.517/4.968) << std::endl;
    
    std::cout << "Bin content for type 9 bbbar is: "<< nStop_values_Type9_bb->Integral() << " and for ccbar is: " << nStop_values_Type9_cc->Integral() << " with k1: " << (49028.) / (1.47945e+06)*(nStop_values_Type9_bb->Integral()/(nStop_values_Type9_bb->Integral()+nStop_values_Type9_cc->Integral()))*(9.517/4.596) << " and k2: " << (49028.) / (1.47945e+06)*(nStop_values_Type9_cc->Integral()/(nStop_values_Type9_bb->Integral()+nStop_values_Type9_cc->Integral()))*(9.517/4.968) << std::endl;
    


    // Vector of histogram names and pointers
    std::vector<std::pair<std::string, TH1D*> > hists;
    hists.push_back(std::make_pair("nStop_values_Type1_bb", nStop_values_Type1_bb));
    hists.push_back(std::make_pair("nStop_values_Type3_bb", nStop_values_Type3_bb));
    hists.push_back(std::make_pair("nStop_values_Type5_bb", nStop_values_Type5_bb));
    hists.push_back(std::make_pair("nStop_values_Type9_bb", nStop_values_Type9_bb));
    hists.push_back(std::make_pair("nStop_values_Type1_cc", nStop_values_Type1_cc));
    hists.push_back(std::make_pair("nStop_values_Type3_cc", nStop_values_Type3_cc));
    hists.push_back(std::make_pair("nStop_values_Type5_cc", nStop_values_Type5_cc));
    hists.push_back(std::make_pair("nStop_values_Type9_cc", nStop_values_Type9_cc));

    for (size_t i = 0; i < hists.size(); ++i) {
        std::string name = hists[i].first;
        TH1D* hist = hists[i].second;

        if (!hist) {
            std::cerr << "Histogram not found: " << name << std::endl;
            continue;
        }

        double xmin = hist->GetXaxis()->GetXmin();
        double xmax = hist->GetXaxis()->GetXmax();

        // Single-component geometric fit
        TF1* singleGeo = new TF1(("singleGeo_" + name).c_str(),
            "[0]*[1]*pow(1-[1], x-1)", xmin, xmax);
        singleGeo->SetParameters(1000, 0.01);
        singleGeo->SetParNames("A", "p");
        singleGeo->SetParLimits(0, 0, 1e9);
        singleGeo->SetParLimits(1, 1e-7, 0.5);
        hist->Fit(singleGeo, "QRN");

        // Two-component geometric fit
        TF1* twoGeo = new TF1(("twoGeo_" + name).c_str(),
            "[0]*[1]*pow(1-[1], x-1) + [2]*[3]*pow(1-[3], x-1)", xmin, xmax);
        twoGeo->SetParameters(1000, 0.01, 500, 0.001);
        twoGeo->SetParNames("A1", "p1", "A2", "p2");
        twoGeo->SetParLimits(0, 0, 1e9);
        twoGeo->SetParLimits(1, 1e-7, 0.5);
        twoGeo->SetParLimits(2, 0, 1e9);
        twoGeo->SetParLimits(3, 1e-7, 0.5);
        hist->Fit(twoGeo, "QRN");

        // Output results
        std::cout << "\n====================" << std::endl;
        std::cout << "Fits for: " << name << std::endl;

        std::cout << "Single Geometric Fit:\n";
        std::cout << "  A  = " << singleGeo->GetParameter(0)
                  << ", p  = " << singleGeo->GetParameter(1) << std::endl;

        std::cout << "Two-Component Geometric Fit:\n";
        std::cout << "  A1 = " << twoGeo->GetParameter(0)
                  << ", p1 = " << twoGeo->GetParameter(1) << std::endl;
        std::cout << "  A2 = " << twoGeo->GetParameter(2)
                  << ", p2 = " << twoGeo->GetParameter(3) << std::endl;
                  
                  
                  
        // Create canvas and draw
        TCanvas* c = new TCanvas(("canvas_" + name).c_str(), name.c_str(), 800, 600);
        hist->SetStats(0);
        hist->SetLineColor(kBlack);
        hist->SetMarkerStyle(20);
        hist->Draw("E");

        singleGeo->SetLineColor(kBlue);
        singleGeo->Draw("same");

        twoGeo->SetLineColor(kRed);
        twoGeo->Draw("same");

        TLegend* leg = new TLegend(0.6, 0.7, 0.89, 0.89);
        leg->AddEntry(hist, "Data", "lep");
        leg->AddEntry(singleGeo, "Single geometric fit", "l");
        leg->AddEntry(twoGeo, "Two-component fit", "l");
        leg->Draw();

        // Save the plot
        std::string outfile = "plots/" + name + ".png";
        c->SaveAs(outfile.c_str());

        delete c; // cleanup
    
    
    
    }

        /*
        // Define the two-component geometric function
        TF1* fitFunc = new TF1("fitFunc",
            // f(x) = A1*p1*(1-p1)^(x-1) + A2*p2*(1-p2)^(x-1)
            "[0]*[1]*pow(1-[1], x-1) + [2]*[3]*pow(1-[3], x-1)",
            nStop_values_Type1_bb->GetXaxis()->GetXmin(),
            nStop_values_Type1_bb->GetXaxis()->GetXmax()
        );
        
        // Parameters: [0]=A1, [1]=p1, [2]=A2, [3]=p2
        fitFunc->SetParameters(1000, 0.01, 500, 0.001);  // initial guesses
        fitFunc->SetParNames("A1", "p1", "A2", "p2");
        
        // Optional: constrain probabilities to [1e-7, 0.5] and amplitudes positive
        fitFunc->SetParLimits(1, 1e-7, 0.5);  // p1
        fitFunc->SetParLimits(3, 1e-7, 0.5);  // p2
        fitFunc->SetParLimits(0, 0, 1e9);     // A1
        fitFunc->SetParLimits(2, 0, 1e9);     // A2
        
        // Fit
        nStop_values_Type1_bb->Fit(fitFunc, "R");
        
        fitFunc->Print();              // Print formula and parameter values
        fitFunc->GetParameter(1);      // Extract best-fit p
        nStop_values_Type1_bb->Draw();
        fitFunc->Draw("same");
        */
    
}