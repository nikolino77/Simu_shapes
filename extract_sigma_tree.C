#include <time.h>
#include <string>
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <cmath>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <algorithm>
#include "TMath.h"
#include "TCanvas.h"
#include "TBranch.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TFile.h"

#ifdef __MAKECINT__
#pragma link C++ class <vector<double>+;
#endif

using namespace std;

int main()
{
    
    TFile* hfile = new TFile("/media/Elements/temp/simu_0_10000.root","OPEN");
    TFile* hfile_out = new TFile("/media/Elements/temp/simu_0_10000_out.root","RECREATE");
    
    TF1* shao = (TF1*)hfile->Get("shao");
    //TF1* shao_smear = (TF1*)hfile->Get("shao_smear");
    //TF1* cer = (TF1*)hfile->Get("cer");
    //TF1* cer_smear = (TF1*)hfile->Get("cer_smear");
    TF1* sum = (TF1*)hfile->Get("sum");
    TF1* conv_last = (TF1*)hfile->Get("conv_last");
    TF1* irf = (TF1*)hfile->Get("irf");
    
    TTree *Singles 	= (TTree*)hfile->Get("tree");
    cout << "Tree " << Singles -> GetName() << endl;
    
    vector<double> *shao_stamps = new vector<double>();
    Singles->SetBranchAddress("shao_stamps",&shao_stamps);
    
    vector<double >* shao_smear_stamps = new vector<double>();
    Singles->SetBranchAddress("shao_smear_stamps",&shao_smear_stamps);
    
    vector<double >* conv_stamps = new vector<double>();
    Singles->SetBranchAddress("conv_stamps",&conv_stamps);
    
    vector<double >* sum_stamps = new vector<double>();
    Singles->SetBranchAddress("sum_stamps",&sum_stamps);
    
    vector<double >* sum_smear_stamps = new vector<double>();
    Singles->SetBranchAddress("sum_smear_stamps",&sum_smear_stamps);
    
    Int_t nentries = Singles->GetEntries();
    cout << "Num entries = " << nentries << endl;
    
    double min = 0.;
    double max = 1000e-9;
    int nbins = 1000000;
    int order = 0;
    
    TTree *tree_out = new TTree("tree_out", "simu_crao_out");
    
    vector<TH1D>* n_photon_shao		= new vector<TH1D>();
    tree_out -> Branch("n_photon_shao",&n_photon_shao);
    
    
    // vector<TH1D* >* n_photon_shao		= new vector<TH1D* >();
    //   vector<TH1D* >* n_photon_shao_smear		= new vector<TH1D* >();
    //  vector<TH1D* >* n_photon_conv		= new vector<TH1D* >();
    //  vector<TH1D* >* n_photon_sum		= new vector<TH1D* >();
    //  vector<TH1D* >* n_photon_sum_smear		= new vector<TH1D* >();
    
    for(int j = 0; j < order; j++)
    {
        stringstream ss;
        ss << j + 1;
        string num = ss.str();
        string histo_name_shao 		= "shao_photon_" + num;
        string histo_name_shao_smear 	= "shao_smear_photon_" + num;
        string histo_name_conv 		= "conv_photon_" + num;
        string histo_name_sum 		= "sum_photon_" + num;
        string histo_name_sum_smear 	= "sum_smear_photon_" + num;
        
        n_photon_shao		-> push_back(new TH1D(histo_name_shao.c_str(),histo_name_shao.c_str(),nbins,min,max));
        // n_photon_shao_smear	-> push_back(new TH1D(histo_name_shao_smear.c_str(),histo_name_shao_smear.c_str(),nbins,min,max));
        // n_photon_conv		-> push_back(new TH1D(histo_name_conv.c_str(),histo_name_conv.c_str(),nbins,min,max));
        // n_photon_sum		-> push_back(new TH1D(histo_name_sum.c_str(),histo_name_sum.c_str(),nbins,min,max));
        // n_photon_sum_smear	-> push_back(new TH1D(histo_name_sum_smear.c_str(),histo_name_sum_smear.c_str(),nbins,min,max));
    }
    
    
    for(int i = 0; i < 10; i++)
    {
        if(i%500 == 0)
        {
            std::cout << "Evento " << i << std::endl;
        }
        
        Singles->GetEntry(i);
        
        for(int j = 0; j < order; j++)
        {
            n_photon_shao		-> at(j) -> Fill(shao_stamps -> at (j));
            //n_photon_shao_smear	-> at(j) -> Fill(shao_smear_stamps -> at (j));
            //n_photon_conv		-> at(j) -> Fill(conv_stamps -> at (j));
            //n_photon_sum		-> at(j) -> Fill(sum_stamps -> at (j));
            //n_photon_sum_smear	-> at(j) -> Fill(sum_smear_stamps -> at (j));
        }
    }
    tree_out -> Fill();
    
    hfile -> Close();
    
    
    
    
    hfile_out -> cd();

    tree_out -> Write();
    
//    shao 	-> Write();
    //shao_smear  -> Write();
    //cer	  -> Write();
    //cer_smear	  -> Write();
  //  sum		-> Write();
  //  conv_last	-> Write();
   // irf		-> Write();
    
    //n_photon_shao -> at(j) -> GetRMS());
    //n_photon_shao_smear	-> at(j) -> GetRMS());
    //n_photon_conv-> at(j) -> GetRMS());
    //n_photon_sum-> at(j) -> GetRMS());
    //n_photon_sum_smear-> at(j) -> GetRMS());   
   // hfile_out -> Write();
    hfile_out -> Close();
    
    return 0;
    
}