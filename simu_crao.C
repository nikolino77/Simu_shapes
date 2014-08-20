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
#include <string>
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
    string name("simu_0_10000");
    string path("/media/Elements/temp/");
    string root = path + name + ".root";
    string dat = path + name + ".txt";
	
    TFile* hfile = new TFile(root.c_str(),"RECREATE");
    ofstream out;
    out.open(dat.c_str());
    
    double min = 0.;
    double max = 1000e-9;
    int n_trials = 100000;
    
    double LY = 10000;
    double CY = 15;
    //double a = LY / CY * (1. / (1 + LY / CY));
    //double b = 1 / (1 + LY / CY);
    
    double a = 1.;
    double b = 0;
    
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
  
    out << "File name: "<< root	<< endl;
    out << "Min: " 	 << min 	<< endl;
    out << "Max: " 	 << max 	<< endl; 
    out << "N_trials: " << n_trials 	<< endl;
    out << "LY: " 	 << LY 		<< endl;
    out << "CY: " 	 << CY 		<< endl;
    out << "theta: " 	 << theta 	<< endl;
    out << "t_d: " 	 << t_d 	<< endl;
    out << "t_r: " 	 << t_r 	<< endl;
    out << "s: " 	 << s 		<< endl;
    out << "tsk: " 	 << tsk 	<< endl;
    out << "l: " 	 << l		<< endl;
    out << "cer_min: " 	 << cer_min	<< endl;
    out << "cer_mean: " << cer_mean	<< endl;
    out.close();
    
    TF1 *shao = new TF1("shao", "(x>[3])*[0]*(exp(-(x-[3])/[1])-exp(-(x-[3])/[2]))", min, max);
    shao -> SetParameters( norm_s, t_d, t_r, theta);
    shao -> SetNpx(1000000);
    shao -> SetLineColor(kYellow);
    cout << "Integral shao = " << shao -> Integral(0., 400e-009) << endl;

    TF1 *shao_smear = new TF1("shao_smear", "sqrt(TMath::Pi()/2)*[0]*[1]*[2]*(exp(([0]*[0]-2*x*[3]+2*[3]*[4]+2*[3]*[6])/(2*[3]*[3]))*(TMath::Erf(([3]*(x-[6]-[4])-[0]*[0])/(sqrt(2)*[0]*[3]))+TMath::Erf(([3]*[4] +[0]*[0])/(sqrt(2)*[0]*[3])))+exp(([0]*[0]-2*x*[5]+2*[5]*[4]+2*[5]*[6])/(2*[5]*[5]))*(TMath::Erf(([5]*(x-[6]-[4])-[0]*[0])/(sqrt(2)*[0]*[5]))+TMath::Erf(([5]*[4] +[0]*[0])/(sqrt(2)*[0]*[5]))))", min, max);
    shao_smear -> SetParameters(s, norm_irf, norm_s, t_d, tsk, t_r, theta);
    shao_smear -> SetNpx(1000000);
    shao_smear -> SetLineColor(kBlack);
    cout << "Integral shao smear = " << shao_smear -> Integral(0., 400e-009) << endl;
    
    TF1 *cer = new TF1("cer", "(x > [0])*[1]*exp(-(x-[2])*(x-[2])/2./[3]/[3])", min, max);
    cer -> SetNpx(1000000);
    cer -> SetParameters(cer_min, norm_c, cer_mean, l);
    cer -> SetLineColor(kRed);
    //cout << "Integral cer = " << cer -> Integral(1.,1.01e-009) << endl;
    
    TF1 *cer_smear = new TF1("cer_smear", "(1.0/(sqrt([0]*[0]+[1]*[1])))*sqrt(TMath::Pi()/2)*[0]*[1]*[2]*[6]*(exp(-([3]+[7]-x)*([3]+[7]-x)/2/([0]*[0]+[1]*[1]))*(TMath::Erf((-[3]*[0]*[0]+[7]*[1]*[1]+[0]*[0]*(x-[4])-[1]*[1]*x+[1]*[1]*(x-[4]))/(sqrt(2)*[0]*[1]*sqrt([0]*[0]+[1]*[1])))-TMath::Erf((-[3]*[0]*[0]+[7]*[1]*[1]-[1]*[1]*x)/(sqrt(2)*[0]*[1]*sqrt([0]*[0]+[1]*[1])))))", min, max);
    cer_smear -> SetNpx(1000000);
    cer_smear -> SetParameters(l, s, norm_c, tsk, cer_min, theta, norm_irf, cer_mean);
    cer_smear -> SetLineColor(kBlack);
    cout << "Integral cer_smear = " << cer_smear -> Integral(0.e-009, 30e-009) << endl;
    
    TF1 *sum = new TF1("sum", "(x > [3])*[0]*(exp(-(x-[3])/[1])-exp(-(x-[3])/[2]))+ (x > [4])*[5]*exp(-(x-[6])*(x-[6])/2./[7]/[7])", min, max);
    sum -> SetNpx(1000000);
    sum -> SetParameters(a*norm_s, t_d, t_r, theta, cer_min, b*norm_c, cer_mean , l);
    sum -> SetLineColor(kRed);
    cout << "Integral sum = " << sum -> Integral(0.,400e-009) << endl;
     
    TF1 *irf = new TF1("irf", "(x > 0.0) * [0] * exp(- (x-[1]) * (x-[1]) / 2. / [2] / [2])", min, max);
    irf -> SetNpx(1000000);
    irf -> SetParameters(norm_irf, tsk, s);
    irf -> SetLineColor(kGreen);
    cout << "Integral irf = " << irf -> Integral(0, 10e-009) << endl;
    
    TF1 *conv_last = new TF1("conv_last", "(1.0/(sqrt([7]*[7]+[0]*[0])))*sqrt(TMath::Pi()/2)*[7]*[0]*[8]*[1]*(exp(-([4]+[10]-x)*([4]+[10]-x)/2/([7]*[7]+[0]*[0]))*(TMath::Erf((-[4]*[7]*[7]+[10]*[0]*[0]+[7]*[7]*(x-[9])-[0]*[0]*x+[0]*[0]*(x-[9]))/(sqrt(2)*[7]*[0]*sqrt([7]*[7]+[0]*[0])))-TMath::Erf((-[4]*[7]*[7]+[10]*[0]*[0]-[0]*[0]*x)/(sqrt(2)*[7]*[0]*sqrt([7]*[7]+[0]*[0])))))+sqrt(TMath::Pi()/2)*[0]*[1]*[2]*(exp(([0]*[0]-2*x*[3]+2*[3]*[4]+2*[3]*[6])/(2*[3]*[3]))*(TMath::Erf(([3]*(x-[6]-[4])-[0]*[0])/(sqrt(2)*[0]*[3]))+TMath::Erf(([3]*[4] +[0]*[0])/(sqrt(2)*[0]*[3])))-exp(([0]*[0]-2*x*[5]+2*[5]*[4]+2*[5]*[6])/(2*[5]*[5]))*(TMath::Erf(([5]*(x-[6]-[4])-[0]*[0])/(sqrt(2)*[0]*[5]))+TMath::Erf(([5]*[4] +[0]*[0])/(sqrt(2)*[0]*[5]))))", min, max);
    conv_last -> SetParameters(s, norm_irf, norm_s*a, t_d, tsk, t_r, theta, l, norm_c*b, cer_min, cer_mean);  
    conv_last -> SetNpx(1000000);
    conv_last -> SetLineColor(kBlack);
    cout << "Integral conv_last = " << conv_last -> Integral(min, 400e-009) << endl;
  
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
    
    vector<double >* sum_smear_stamps = new vector<double>();
    tree -> Branch("sum_smear_stamps",sum_smear_stamps);
    
    TRandom3* rand = new TRandom3();
    
    for(int i = 0; i < n_trials; i++)
    {
      if(i%100 == 0)
      {
        std::cout << "Evento " << i << std::endl;
      }
      
      //int mi_1 = pois -> GetRandom();
      int mi_1 = LY+CY;
      for(int j = 0; j < mi_1; j++)
      {
	double shao_extr = shao 	-> GetRandom();
	double sum_extr  = sum 		-> GetRandom();
	double irf_extr  = irf 		-> GetRandom();
	double conv_extr = conv_last  	-> GetRandom();

	conv_stamps 		-> push_back(conv_extr);
	shao_smear_stamps 	-> push_back(shao_extr + irf_extr);
	shao_stamps 		-> push_back(shao_extr);
	sum_stamps 		-> push_back(sum_extr);
	sum_smear_stamps 	-> push_back(sum_extr + irf_extr);
      }
      
      std::sort(shao_stamps -> begin(), shao_stamps -> end()); 
      std::sort(shao_smear_stamps -> begin(), shao_smear_stamps -> end());
      std::sort(conv_stamps -> begin(), conv_stamps -> end());
      std::sort(sum_stamps -> begin(), sum_stamps -> end());
      std::sort(sum_smear_stamps -> begin(), sum_smear_stamps -> end());

      tree -> Fill();

      shao_stamps 	-> clear();
      conv_stamps 	-> clear();
      sum_stamps 	-> clear();
      shao_smear_stamps	-> clear();
      sum_smear_stamps 	-> clear();
    }
    
    conv_last	-> Write();
    sum		-> Write();
    shao	-> Write();
    irf		-> Write();
    tree	-> Write();
      
//     TH1F* hshao = new TH1F("hshao", "hshao", 200000, min,max);
//     TH1F* hshao_smear = new TH1F("hshao_smear", "hshao_smear", 200000,min,max);
//     TH1F* hconv = new TH1F("hconv", "hconv", 200000, min,max);
//     TH1F* hsum = new TH1F("hsum", "hsum", 200000, min,max);
//     TH1F* hsum_smear = new TH1F("hsum_smear", "hsum_smear", 200000,min,max);
//     
//     for(int i = 0; i < 10000000; i++)
//     {
//       if(i%100000 == 0)
//       {
//         std::cout << "Evento " << i << std::endl;
//       }
//       
//       double irf_extr = irf -> GetRandom();
//       double conv_extr = conv_last  -> GetRandom();
//       double shao_extr = shao -> GetRandom();
//       double sum_extr = sum  -> GetRandom();
//       double shao_smear_extr = shao_smear->GetRandom();
//       
// 	hshao ->Fill(shao_smear_extr);
// 	hsum -> Fill(sum_extr);
// 	hconv-> Fill(conv_extr);
// 	hsum_smear->Fill(sum_extr+irf_extr);
// 	hshao_smear->Fill(shao_extr+irf_extr);
// 	
//     }
//     
//     	hshao -> Write();
// 	hsum -> Write();
// 	hconv -> Write();
// 	hsum_smear -> Write();
// 	hshao_smear -> Write();
//  	shao_smear->Write();
//  	irf->Write();
//  	sum->Write();
//  	conv_last->Write();
	
    hfile -> Close();
    return 0;

}