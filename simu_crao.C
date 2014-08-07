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
#include "TRandom3.h"
#include "TFile.h"

#ifdef __MAKECINT__
#pragma link C++ class <vector<double>+;
#endif

using namespace std;

int main()
{

    TFile* hfile = new TFile("simu_crao.root","RECREATE");
   
    double min = 0.;
    double max = 1000e-009;
    
    double LY = 5000;
    double CY = 20;
    double a = LY / CY * (1. / (1 + LY / CY));
    double b = 1 / (1 + LY / CY);
    
    double theta = 200e-12;
    
    double t_d = 30.3e-9; // rise time 2
    double t_r = 70e-12; //decay time 5
    double norm_s = 1.0 / (t_d-t_r) ; //normalization shao 4
    double s = 66e-12; // sigma trans 1 / 3
    double tsk = 2e-9; // 3 / 0
    
    double norm_irf = 2. / s / TMath::Sqrt(2) / TMath::Sqrt(TMath::Pi()) / TMath::Erfc(- tsk / s / TMath::Sqrt(2));
    double l = 10e-12; // cerenkov sigma 2
    double cer_min = theta - 2 * l; 
    double norm_c = 2 / (l * TMath::Sqrt(2 * TMath::Pi()) * TMath::Erfc(- TMath::Sqrt(2))); 
    
    TF1 *shao = new TF1("shao", "(x > [3]) * [0] * (TMath::Exp(- (x - [3]) / [1]) - TMath::Exp(- (x - [3]) / [2]))", min, max);
    shao -> SetParameters(norm_s, t_d, t_r, theta);
    shao -> SetNpx(1000000);
    shao -> SetLineColor(kYellow);

    TF1 *shao_smear = new TF1("shao_smear", "TMath::Sqrt(TMath::Pi() / 2) * [0] * [1] * [2] * (-TMath::Exp(([0]*[0]-2*x*[3]+2*[3]*[4]+2*[3]*[6])/(2*[3]*[3]))*(TMath::Erf(([3]*([4]-x+[6]) + [0]*[0])/(TMath::Sqrt(2)*[0]*[3]))-TMath::Erf(([3]*[4] +[0]*[0])/(TMath::Sqrt(2)*[0]*[3])))+TMath::Exp(([0]*[0]-2*x*[5]+2*[5]*[4]+2*[5]*[6])/(2*[5]*[5]))*(TMath::Erf(([5]*([4]-x+[6]) + [0]*[0])/(TMath::Sqrt(2)*[0]*[5]))-TMath::Erf(([5]*[4] +[0]*[0])/(TMath::Sqrt(2)*[0]*[5]))))", min, max);
    shao_smear -> SetParameters(s, norm_irf, norm_s, t_d, tsk, t_r, theta);
    shao_smear -> SetNpx(1000000);
    shao_smear -> SetLineColor(kBlack);
   
    TF1 *irf = new TF1("irf", "(x > 0.0) * [0] * exp(- (x-[1]) * (x-[1]) / 2. / [2] / [2])", min, max);
    irf -> SetNpx(1000000);
    irf -> SetParameters(norm_irf, tsk, s);
    irf -> SetLineColor(kGreen);
   
    TF1 *sum = new TF1("sum", "(x > [3]) * [0] * (TMath::Exp(- (x - [3]) / [1]) - TMath::Exp(- (x - [3]) / [2])) +  (x > [4]) * [5] * exp(- (x-[3]) * (x-[3]) / 2. / [6] / [6])", min, max);
    sum -> SetParameters(a*norm_s, t_d, t_r, theta, cer_min, b*norm_c, l);
    sum -> SetNpx(1000000);
    sum -> SetLineColor(kYellow);

    TF1 *conv_last = new TF1("conv_last", "TMath::Sqrt(TMath::Pi() / 2) * [0] * [1] * [2] * (-TMath::Exp(([0]*[0]-2*x*[3]+2*[3]*[4]+2*[3]*[6])/(2*[3]*[3]))*(TMath::Erf(([3]*([4]-x+[6]) + [0]*[0])/(TMath::Sqrt(2)*[0]*[3]))-TMath::Erf(([3]*[4] +[0]*[0])/(TMath::Sqrt(2)*[0]*[3])))+TMath::Exp(([0]*[0]-2*x*[5]+2*[5]*[4]+2*[5]*[6])/(2*[5]*[5]))*(TMath::Erf(([5]*([4]-x+[6]) + [0]*[0])/(TMath::Sqrt(2)*[0]*[5]))-TMath::Erf(([5]*[4] +[0]*[0])/(TMath::Sqrt(2)*[0]*[5]))))+1.0 / (TMath::Sqrt([7]*[7]+[0]*[0])) * TMath::Sqrt(TMath::Pi() / 2) * [7] * [0] * [8] * [1] * (TMath::Exp(-([4]+[6]-x)*([4]+[6]-x)/2/([7]*[7]+[0]*[0]))*TMath::Erf((-[4]*[7]*[7]+[6]*[0]*[0]+[7]*[7]*(x-[9])-[0]*[0]*x+[0]*[0]*(x-[9]))/(TMath::Sqrt(2)*[0]*[7]*TMath::Sqrt([7]*[7]+[0]*[0])))-TMath::Exp(-([4]+[6]-x)*([4]+[6]-x)/2/([7]*[7]+[0]*[0]))*TMath::Erf((-[4]*[7]*[7]+[6]*[0]*[0]-[0]*[0]*x)/(TMath::Sqrt(2)*[0]*[7]*TMath::Sqrt([7]*[7]+[0]*[0]))))", min, max);
    conv_last -> SetParameters(s, norm_irf, norm_s*a, t_d, tsk, t_r, theta, l, norm_c*b, cer_min);
    conv_last -> SetNpx(1000000);
    conv_last -> SetLineColor(kBlack);
  
    TF1 *pois = new TF1("pois", "TMath::PoissonI(x, [0])", 0, 100000);
    pois -> SetNpx(1000000);
    pois -> SetParameter(0, LY+CY);

    TTree *tree = new TTree("tree", "simu_crao");
    
    vector<double >* shao_stamps = new vector<double>();
    tree -> Branch("shao_stamps",shao_stamps);

    vector<double >* shao_smear_stamps = new vector<double>();
    tree -> Branch("shao_smear_stamps",shao_smear_stamps);
    
    vector<double >* conv_stamps = new vector<double>();
    tree -> Branch("conv_stamps",conv_stamps);
    
    vector<double >* sum_stamps = new vector<double>();
    tree -> Branch("sum_stamps",sum_stamps);
    
    TRandom3* rand = new TRandom3();
    int n_trials = 1000000;
    
    for(int i = 0; i < n_trials; i++)
    {
      if(i%100 == 0)
      {
        std::cout << "Evento " << i << std::endl;
      }
      
      int mi_1 = pois -> GetRandom();
      for(int j = 0; j < mi_1; j++)
      {
	double shao_extr = shao -> GetRandom();
	double sum_extr = sum -> GetRandom();
	double irf_extr = irf -> GetRandom();
	double conv_extr = conv_last  -> GetRandom();
	//double conv_extr = 0;
	conv_stamps -> push_back(conv_extr);
	shao_smear_stamps -> push_back(shao_extr + irf_extr);
	shao_stamps -> push_back(shao_extr);
	sum_stamps -> push_back(sum_extr);
      }
      
      std::sort(shao_stamps -> begin(), shao_stamps -> end()); 
      std::sort(shao_smear_stamps -> begin(), shao_smear_stamps -> end());
      std::sort(conv_stamps -> begin(), conv_stamps -> end());
      std::sort(sum_stamps -> begin(), sum_stamps -> end());
      
      tree -> Fill();
      shao_stamps -> clear();
      conv_stamps -> clear();
      sum_stamps -> clear();
      shao_smear_stamps -> clear();
    }
    
    conv_last -> Write();
    sum -> Write();
    shao -> Write();
    irf -> Write();
    tree -> Write();
    
    hfile -> Close();
    return 0;

}