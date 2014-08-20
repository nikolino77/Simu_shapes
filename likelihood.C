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
  
  
    TFile* hfile = new TFile("/media/Elements/temp/simu_0_4000.root","OPEN");
    TFile* hfile_out = new TFile("/media/Elements/temp/simu_0_4000_out.root","RECREATE");
    
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
    
    double LY = 5000;
    double CY = 20;
    double a = LY / CY * (1. / (1 + LY / CY));
    double b = 1 / (1 + LY / CY);
    
    double theta = 1000e-12;
    
    double t_d = 30.3e-9; // rise time 2
    double t_r = 70e-12; //decay time 5

    double s = 66e-12; // sigma trans 1 / 3
    double tsk = 2e-9; // 3 / 0
    
    double l = 10e-12; // cerenkov sigma 2   
    double cer_min = theta; 
    double cer_mean = theta + 2* l;
    
    double norm_c 	= 1. / (l * sqrt(TMath::Pi() / 2) * (1 + TMath::Erf((cer_mean - cer_min) / sqrt(2) / l))); 
    double norm_irf	= 1. / (s * sqrt(TMath::Pi() / 2) * (1 + TMath::Erf(tsk / sqrt(2) / s)));
    double norm_s 	= 1. / (t_d-t_r) ; //normalization shao 4

    double norm_irf_easy = 1. / (s * sqrt(TMath::Pi() / 2) * (1 + TMath::Erf((tsk + cer_mean) / sqrt(2) / s)));
        
    //----------------- Cumulative total ------------------------//
    
    TF1 *CTot = new TF1("CTot", "[7] * [0] * sqrt(TMath::Pi() / 2) * [3] * (TMath::Erf((x-[2]-[6])/([3]*sqrt(2))) + TMath::Erf(([2]+[6])/([3]*sqrt(2)))) + [8] * [0]*[3]*sqrt(TMath::Pi()/2)*((TMath::Erf((x-[1]-[2])/[3]/sqrt(2))+TMath::Erf([2]/[3]/sqrt(2)))+([4]/([5]-[4])*TMath::Exp((-2*x*[4]+2*[4]*[1]+2*[4]*[2]+[3]*[3])/2/[4]/[4])*(TMath::Erf((x-[1]-[2]-[3]*[3]/[4])/([3]*sqrt(2)))+TMath::Erf(([2]+[3]*[3]/[4])/[3]/sqrt(2))))-([5]/([5]-[4])*TMath::Exp((-2*x*[5]+2*[5]*[1]+2*[5]*[2]+[3]*[3])/2/[5]/[5]) * (TMath::Erf((x-[1]-[2]-[3]*[3]/[5])/[3]/sqrt(2))+TMath::Erf(([2]+[3]*[3]/[5])/[3]/sqrt(2)))))", min, max);
    CTot -> SetNpx(10000000);    
    CTot -> SetParameters(norm_irf, theta, tsk, s, t_r, t_d, cer_mean, b, a);
    CTot -> SetLineColor(kGreen);
    
    vector<TH1D* >* n_photon_shao		= new vector<TH1D* >();
    vector<TH1D* >* n_photon_shao_smear	= new vector<TH1D* >();
    vector<TH1D* >* n_photon_conv		= new vector<TH1D* >();
    vector<TH1D* >* n_photon_sum		= new vector<TH1D* >();
    vector<TH1D* >* n_photon_sum_smear		= new vector<TH1D* >();
    
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
      n_photon_shao_smear	-> push_back(new TH1D(histo_name_shao_smear.c_str(),histo_name_shao_smear.c_str(),nbins,min,max));    
      n_photon_conv		-> push_back(new TH1D(histo_name_conv.c_str(),histo_name_conv.c_str(),nbins,min,max));    
      n_photon_sum		-> push_back(new TH1D(histo_name_sum.c_str(),histo_name_sum.c_str(),nbins,min,max));    
      n_photon_sum_smear	-> push_back(new TH1D(histo_name_sum_smear.c_str(),histo_name_sum_smear.c_str(),nbins,min,max));         
    }

    
    for(int i = 0; i < 100000; i++)
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
	n_photon_sum_smear	-> at(j) -> Fill(sum_smear_stamps -> at (j));	
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
    TGraph* var_sum_smear = new TGraph(order);
    var_sum_smear -> SetLineColor(kGreen);
    
    for(int j = 0; j < order; j++)
    {   
	var_shao 	-> SetPoint(j, j, n_photon_shao -> at(j) -> GetRMS());
	var_shao_smear 	-> SetPoint(j, j, n_photon_shao_smear	-> at(j) -> GetRMS());
	var_conv 	-> SetPoint(j, j, n_photon_conv-> at(j) -> GetRMS());
	var_sum 	-> SetPoint(j, j, n_photon_sum-> at(j) -> GetRMS());
   	var_sum_smear 	-> SetPoint(j, j, n_photon_sum_smear-> at(j) -> GetRMS());   
    }
    
    mg	-> Add(var_shao,"lp");
    mg	-> Add(var_shao_smear,"lp");
    mg	-> Add(var_conv,"lp");
    mg	-> Add(var_sum,"lp");
    mg	-> Add(var_sum_smear, "lp");
    mg 	-> Draw("a");
    
    hfile_out -> cd();

    n_photon_shao		-> at(0) -> Write();
    n_photon_shao_smear		-> at(0) -> Write();
    n_photon_conv		-> at(0) -> Write();
    n_photon_sum		-> at(0) -> Write();
    n_photon_sum_smear		-> at(0) -> Write();
    mg 			-> Write();
    var_shao 		-> Write();
    var_shao_smear 	-> Write();
    var_conv 		-> Write();
    var_sum 		-> Write();
    var_sum_smear 	-> Write();

    hfile_out -> Close();
     
    return 0;

}