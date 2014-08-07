#include "TF1.h"
#include "Math/GSLIntegrator.h"
#include "Math/WrappedTF1.h"
#include <iostream>



void testIntegralGaus() { 

  double min = 0.;
  double max = 200000;
    
  double LY = 4577;
  double CY = 20;
  double a = LY / CY * (1. / (1 + LY / CY));
  double b = 1 / (1 + LY / CY);
    
  double theta = 200;
    
  double t_d = 30000; // rise time 2
  double t_r = 70; //decay time 5
  double norm_s = a * 1.0 / (t_d-t_r) ; //normalization shao 4
  double s = 66; // sigma trans 1 / 3
  double tsk = 2000; // 3 / 0
    
  double norm_irf = 2. / s / TMath::Sqrt(2) / TMath::Sqrt(TMath::Pi()) / TMath::Erfc(- tsk / s / TMath::Sqrt(2));
  double l = 10; // cerenkov sigma 2
  double cer_min = theta - 2 * l; 
  double norm_c = b * 2 / (l * TMath::Sqrt(2 * TMath::Pi()) * TMath::Erfc(- TMath::Sqrt(2))); 
    
  TF1 *conv_last = new TF1("conv_last", "TMath::Sqrt(TMath::Pi() / 2) * [0] * [1] * [2] * (-TMath::Exp(([0]*[0]-2*x*[3]+2*[3]*[4]+2*[3]*[6])/(2*[3]*[3]))*(TMath::Erf(([3]*([4]-x+[6]) + [0]*[0])/(TMath::Sqrt(2)*[0]*[3]))-TMath::Erf(([3]*[4] +[0]*[0])/(TMath::Sqrt(2)*[0]*[3])))+TMath::Exp(([0]*[0]-2*x*[5]+2*[5]*[4]+2*[5]*[6])/(2*[5]*[5]))*(TMath::Erf(([5]*([4]-x+[6]) + [0]*[0])/(TMath::Sqrt(2)*[0]*[5]))-TMath::Erf(([5]*[4] +[0]*[0])/(TMath::Sqrt(2)*[0]*[5]))))+1.0 / (TMath::Sqrt([7]*[7]+[0]*[0])) * TMath::Sqrt(TMath::Pi() / 2) * [7] * [0] * [8] * [1] * (TMath::Exp(-([4]+[6]-x)*([4]+[6]-x)/2/([7]*[7]+[0]*[0]))*TMath::Erf((-[4]*[7]*[7]+[6]*[0]*[0]+[7]*[7]*(x-[9])-[0]*[0]*x+[0]*[0]*(x-[9]))/(TMath::Sqrt(2)*[0]*[7]*TMath::Sqrt([7]*[7]+[0]*[0])))-TMath::Exp(-([4]+[6]-x)*([4]+[6]-x)/2/([7]*[7]+[0]*[0]))*TMath::Erf((-[4]*[7]*[7]+[6]*[0]*[0]-[0]*[0]*x)/(TMath::Sqrt(2)*[0]*[7]*TMath::Sqrt([7]*[7]+[0]*[0]))))", min, max);
  conv_last -> SetParameters(s, norm_irf, norm_s, t_d, tsk, t_r, theta, l, norm_c, cer_min);
  conv_last -> SetNpx(1000000);
  conv_last -> SetLineColor(kBlack);
  
  cout << "Integral irf = " << conv_last -> Integral(min, max) << endl;
  cout << "norm_irf = " << norm_irf << endl;
  TCanvas *c6 = new TCanvas("c6");
  conv_last -> Draw();

}
