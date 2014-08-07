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

    TFile* hfile = new TFile("simu_crao.root","OPEN");
    TFile* hfile_out = new TFile("simu_crao_out.root","RECREATE");
    
    TTree *Singles 	= (TTree*)hfile->Get("tree");
    cout << "Tree " << Singles -> GetName() << endl;
    
    vector<double> *shao_stamps = new vector<double>();
    Singles->SetBranchAddress("shao_stamps",&shao_stamps);
     
    vector<double >* shao_smear_stamps = new vector<double>();
    Singles->SetBranchAddress("shao_smear_stamps",&shao_smear_stamps);

    vector<double >* conv_stamps = new vector<double>();
    Singles->SetBranchAddress("conv_stamps",&conv_stamps);
   
    //vector<double >* sum_stamps = new vector<double>();
    //Singles->SetBranchAddress("s_stamps",&shao_stamps);
    
    vector<double >* sum_stamps = new vector<double>();
    Singles->SetBranchAddress("sum_stamps",&sum_stamps);

    Int_t nentries = Singles->GetEntries(); 
    cout << "Num entries = " << nentries << endl;	
    
    double min = 0.;
    double max = 200e-09;
    int nbins = 200000;
    int order = 100;
    
    vector<TH1D* >* n_photon_shao		= new vector<TH1D* >();
    vector<TH1D* >* n_photon_shao_smear	= new vector<TH1D* >();
    vector<TH1D* >* n_photon_conv		= new vector<TH1D* >();
    vector<TH1D* >* n_photon_sum		= new vector<TH1D* >();
		
    for(int j = 0; j < order; j++)
    {
      stringstream ss;
      ss << j + 1;
      string num = ss.str();
      string histo_name_shao 		= "shao_photon_" + num;
      string histo_name_shao_smear 	= "shao_smear_photon_" + num;
      string histo_name_conv 		= "conv_photon_" + num;
      string histo_name_sum 		= "sum_photon_" + num;

      n_photon_shao		-> push_back(new TH1D(histo_name_shao.c_str(),histo_name_shao.c_str(),nbins,min,max));    
      n_photon_shao_smear	-> push_back(new TH1D(histo_name_shao_smear.c_str(),histo_name_shao_smear.c_str(),nbins,min,max));    
      n_photon_conv		-> push_back(new TH1D(histo_name_conv.c_str(),histo_name_conv.c_str(),nbins,min,max));    
      n_photon_sum		-> push_back(new TH1D(histo_name_sum.c_str(),histo_name_sum.c_str(),nbins,min,max));    
    }

    
    for(int i = 0; i < nentries; i++)
    {
      if(i%500 == 0)
      {
        std::cout << "Evento " << i << std::endl;
      }
      
      Singles->GetEntry(i);
      
      for(int j = 0; j < order; j++)
      {
	n_photon_shao		-> at(j) -> Fill(shao_stamps -> at (j));
	n_photon_shao_smear	-> at(j) -> Fill(shao_smear_stamps -> at (j));
	n_photon_conv		-> at(j) -> Fill(conv_stamps -> at (j));
	n_photon_sum		-> at(j) -> Fill(sum_stamps -> at (j));
      }
    }
    
    hfile -> Close();
    
    
    TMultiGraph *mg = new TMultiGraph();
    
    TGraph* var_shao = new TGraph(order);
    var_shao -> SetLineColor(kRed);
    TGraph* var_shao_smear = new TGraph(order);
    var_shao_smear -> SetLineColor(kBlack);
    TGraph* var_conv = new TGraph(order);
    var_conv -> SetLineColor(kBlue);
    TGraph* var_sum = new TGraph(order);
    var_sum -> SetLineColor(kGreen);
    
    for(int j = 0; j < order; j++)
    {   
	var_shao -> SetPoint(j, j, n_photon_shao -> at(j) -> GetMean());
	var_shao_smear -> SetPoint(j, j, n_photon_shao_smear	-> at(j) -> GetMean());
	var_conv -> SetPoint(j, j, n_photon_conv-> at(j) -> GetMean());
	var_sum -> SetPoint(j, j, n_photon_sum-> at(j) -> GetMean());
    }
    
    mg->Add(var_shao,"lp");
    mg->Add(var_shao_smear,"lp");
    mg->Add(var_conv,"lp");
    mg->Add(var_sum,"lp");
    mg -> Draw("a");
    
    
    hfile_out -> cd();

    n_photon_shao		-> at(0) -> Write();
    n_photon_shao_smear		-> at(0) -> Write();
    n_photon_conv		-> at(0) -> Write();
    n_photon_sum		-> at(0) -> Write();
    mg -> Write();
    var_shao -> Write();
    var_shao_smear -> Write();
    var_conv -> Write();
    var_sum -> Write();
    
//    n_photon_shao	-> Write();
    hfile_out -> Close();
     
    return 0;

}