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
    int order = 100;
    
    TTree *tree_sigma = new TTree("tree_sigma", "extract_sigma");
    
    TH1D* n_photon_shao_ex = new TH1D("n_photon_shao_ex","n_photon_shao_ex",nbins, min, max);
    tree_sigma -> Branch("n_photon_shao_ex",&n_photon_shao_ex);
    
    TH1D* n_photon_shao_smear_ex = new TH1D("n_photon_shao_ex","n_photon_shao_smear_ex",nbins, min, max);
    tree_sigma -> Branch("n_photon_shao_smear_ex",&n_photon_shao_smear_ex);
    
    TH1D* n_photon_conv_ex = new TH1D("n_photon_conv_ex","n_photon_conv_ex",nbins, min, max);
    tree_sigma -> Branch("n_photon_conv_ex",&n_photon_conv_ex);
    
    TH1D* n_photon_sum_ex = new TH1D("n_photon_sum_ex","n_photon_sum_ex",nbins, min, max);
    tree_sigma -> Branch("n_photon_sum_ex",&n_photon_sum_ex);
    
    TH1D* n_photon_sum_smear_ex = new TH1D("n_photon_sum_smear_ex","n_photon_sum_smear_ex",nbins, min, max);
    tree_sigma -> Branch("n_photon_sum_smear_ex",&n_photon_sum_smear_ex);
    
    for(int j = 0; j < order; j++)
    {
            std::cout << "Evento " << j << std::endl;
      for(int i = 0; i < nentries; i++)
      {
        Singles->GetEntry(i);
          n_photon_shao_ex	->	Fill(shao_stamps -> at (j));
        n_photon_shao_smear_ex->	Fill(shao_smear_stamps -> at (j));
        n_photon_conv_ex	->	Fill(conv_stamps -> at (j));
        n_photon_sum_ex		-> Fill(sum_stamps -> at (j));
        n_photon_sum_smear_ex	->  Fill(sum_smear_stamps -> at (j));
      }
    
      tree_sigma -> Fill();
    }

    hfile -> Close();
    
    hfile_out -> cd();
    tree_sigma -> Write();
    hfile_out -> Close();
     
    return 0;

}
