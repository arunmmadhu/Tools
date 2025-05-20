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
    TFile *TreeFile_h = new TFile("LOCAL_COMBINED_ztau3mutauh_skimmer_default.root","READ");//tauh,mu,e doesn't matter in skimmer
    
    TH1D * nStop_values_Type1  = (TH1D*)TreeFile_h->Get("ztau3mutauh_skimmer_default_Selection_nStop_DistributionMC8");
    TH1D * nStop_values_Type3  = (TH1D*)TreeFile_h->Get("ztau3mutauh_skimmer_default_Selection_nStop_DistributionMC9");
    TH1D * nStop_values_Type5  = (TH1D*)TreeFile_h->Get("ztau3mutauh_skimmer_default_Selection_nStop_DistributionMC10");
    TH1D * nStop_values_Type9  = (TH1D*)TreeFile_h->Get("ztau3mutauh_skimmer_default_Selection_nStop_DistributionMC11");
    
    //Getting the mean
    double mean_type1 = nStop_values_Type1->GetMean();
    double mean_type3 = nStop_values_Type3->GetMean();
    double mean_type5 = nStop_values_Type5->GetMean();
    double mean_type9 = nStop_values_Type9->GetMean();
    
    std::cout << "The mean for type 1 is: "<< mean_type1 << std::endl;
    std::cout << "The mean for type 3 is: "<< mean_type3 << std::endl;
    std::cout << "The mean for type 5 is: "<< mean_type5 << std::endl;
    std::cout << "The mean for type 9 is: "<< mean_type9 << std::endl;
    
    
    
    
}