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

void calc_multiplicity () 
{
    TFile *TreeFile_h = new TFile("LOCAL_COMBINED_ztau3mutauh_default_LumiScaled.root","READ");
    TFile *TreeFile_mu = new TFile("LOCAL_COMBINED_ztau3mutaumu_default_LumiScaled.root","READ");
    TFile *TreeFile_e = new TFile("LOCAL_COMBINED_ztau3mutaue_default_LumiScaled.root","READ");
    
    TH1D * taue_no_of_triplet  = (TH1D*)TreeFile_e->Get("ztau3mutaue_default_PostSelection_No_of_candidates_TripletMC1");
    TH1D * taumu_no_of_triplet  = (TH1D*)TreeFile_mu->Get("ztau3mutaumu_default_PostSelection_No_of_candidates_TripletMC2");
    TH1D * tauh_no_of_triplet_ha  = (TH1D*)TreeFile_h->Get("ztau3mutauh_default_PostSelection_A_No_of_candidates_TripletMC3");
    TH1D * tauh_no_of_triplet_hb  = (TH1D*)TreeFile_h->Get("ztau3mutauh_default_PostSelection_B_No_of_candidates_TripletMC3");
    
    TH1D * taue_no_of_triplet_Data  = (TH1D*)TreeFile_e->Get("ztau3mutaue_default_PostSelection_No_of_candidates_TripletData");
    TH1D * taumu_no_of_triplet_Data  = (TH1D*)TreeFile_mu->Get("ztau3mutaumu_default_PostSelection_No_of_candidates_TripletData");
    TH1D * tauh_no_of_triplet_ha_Data  = (TH1D*)TreeFile_h->Get("ztau3mutauh_default_PostSelection_A_No_of_candidates_TripletData");
    TH1D * tauh_no_of_triplet_hb_Data  = (TH1D*)TreeFile_h->Get("ztau3mutauh_default_PostSelection_B_No_of_candidates_TripletData");
    
    TH1D * taue_no_of_e  = (TH1D*)TreeFile_e->Get("ztau3mutaue_default_PostSelection_No_of_candidates_OS_sideMC1");
    TH1D * taumu_no_of_mu  = (TH1D*)TreeFile_mu->Get("ztau3mutaumu_default_PostSelection_No_of_candidates_OS_sideMC2");
    TH1D * tauh_no_of_ha  = (TH1D*)TreeFile_h->Get("ztau3mutauh_default_PostSelection_A_No_of_candidates_OS_sideMC3");
    TH1D * tauh_no_of_hb  = (TH1D*)TreeFile_h->Get("ztau3mutauh_default_PostSelection_B_No_of_candidates_OS_sideMC3");
    
    TH1D * taue_no_of_e_Data  = (TH1D*)TreeFile_e->Get("ztau3mutaue_default_PostSelection_No_of_candidates_OS_sideData");
    TH1D * taumu_no_of_mu_Data  = (TH1D*)TreeFile_mu->Get("ztau3mutaumu_default_PostSelection_No_of_candidates_OS_sideData");
    TH1D * tauh_no_of_ha_Data  = (TH1D*)TreeFile_h->Get("ztau3mutauh_default_PostSelection_A_No_of_candidates_OS_sideData");
    TH1D * tauh_no_of_hb_Data  = (TH1D*)TreeFile_h->Get("ztau3mutauh_default_PostSelection_B_No_of_candidates_OS_sideData");
    
    
    //For tau e
    std::cout << "For tau_e: "<<std::endl;
    
    std::cout << "Percent of taue MC events with 1 triplet: " << taue_no_of_triplet->GetBinContent(2)/( taue_no_of_triplet->Integral(1, 11) ) << " and with OS stuff: " << taue_no_of_e->GetBinContent(2)/( taue_no_of_e->Integral(1, 11) ) << " total MC events: "<< taue_no_of_e->Integral(1, 11) << std::endl;
    
    std::cout << "Percent of taue Data events with 1 triplet: " << taue_no_of_triplet_Data->GetBinContent(2)/( taue_no_of_triplet_Data->Integral(1, 11) ) << " and with OS stuff: " << taue_no_of_e_Data->GetBinContent(2)/( taue_no_of_e_Data->Integral(1, 11) ) << " total Data events: "<< taue_no_of_e_Data->Integral(1, 11) << std::endl;
    
    std::cout << "--- --- --- ---"<< std::endl; 
    std::cout << "--- --- --- ---"<< std::endl;
    std::cout << "--- --- --- ---"<< std::endl;
    
    //For tau mu
    std::cout << "For tau_mu: "<<std::endl;
    
    std::cout << "Percent of taumu MC events with 1 triplet: " << taumu_no_of_triplet->GetBinContent(2)/( taumu_no_of_triplet->Integral(1, 11) ) << " and with OS stuff: " << taumu_no_of_mu->GetBinContent(2)/( taumu_no_of_mu->Integral(1, 11) ) << " total MC events: "<< taumu_no_of_mu->Integral(1, 11) << std::endl;
    
    std::cout << "Percent of taumu Data events with 1 triplet: " << taumu_no_of_triplet_Data->GetBinContent(2)/( taumu_no_of_triplet_Data->Integral(1, 11) ) << " and with OS stuff: " << taumu_no_of_mu_Data->GetBinContent(2)/( taumu_no_of_mu_Data->Integral(1, 11) ) << " total Data events: "<< taumu_no_of_mu_Data->Integral(1, 11) << std::endl;
    
    std::cout << "--- --- --- ---"<< std::endl; 
    std::cout << "--- --- --- ---"<< std::endl;
    std::cout << "--- --- --- ---"<< std::endl;
    
    //For tau ha
    std::cout << "For tau_ha: "<<std::endl;
    
    std::cout << "Percent of tauha MC events with 1 triplet: " << tauh_no_of_triplet_ha->GetBinContent(2)/( tauh_no_of_triplet_ha->Integral(1, 11) ) << " and with OS stuff: " << tauh_no_of_ha->GetBinContent(2)/( tauh_no_of_ha->Integral(1, 11) ) << " total MC events: "<< tauh_no_of_ha->Integral(1, 11) << std::endl;
    
    std::cout << "Percent of tauha Data events with 1 triplet: " << tauh_no_of_triplet_ha_Data->GetBinContent(2)/( tauh_no_of_triplet_ha_Data->Integral(1, 11) ) << " and with OS stuff: " << tauh_no_of_ha_Data->GetBinContent(2)/( tauh_no_of_ha_Data->Integral(1, 11) ) << " total Data events: "<< tauh_no_of_ha_Data->Integral(1, 11) << std::endl;
    
    std::cout << "--- --- --- ---"<< std::endl; 
    std::cout << "--- --- --- ---"<< std::endl;
    std::cout << "--- --- --- ---"<< std::endl;
    
    //For tau ha
    std::cout << "For tau_hb: "<<std::endl;
    
    std::cout << "Percent of tauhb MC events with 1 triplet: " << tauh_no_of_triplet_hb->GetBinContent(2)/( tauh_no_of_triplet_hb->Integral(1, 11) ) << " and with OS stuff: " << tauh_no_of_hb->GetBinContent(2)/( tauh_no_of_hb->Integral(1, 11) ) << " total MC events: "<< tauh_no_of_hb->Integral(1, 11) << std::endl;
    
    std::cout << "Percent of tauhb Data events with 1 triplet: " << tauh_no_of_triplet_hb_Data->GetBinContent(2)/( tauh_no_of_triplet_hb_Data->Integral(1, 11) ) << " and with OS stuff: " << tauh_no_of_hb_Data->GetBinContent(2)/( tauh_no_of_hb_Data->Integral(1, 11) ) << " total Data events: "<< tauh_no_of_hb_Data->Integral(1, 11) << std::endl;
    
    std::cout << "--- --- --- ---"<< std::endl; 
    std::cout << "--- --- --- ---"<< std::endl;
    std::cout << "--- --- --- ---"<< std::endl;
    
    
    
    
}