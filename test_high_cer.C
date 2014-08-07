{

  // A macro to check the ratio between cerenkov peaks and scintillation
   gSystem->Load("libMathMore");
  gROOT->SetStyle("Plain");
  gStyle->SetOptFit(1);
 
  double min = 0.;
  double max = 200e-9;
  
  double t_0 = 100e-12;
  double cer_sigma = 5e-12;
  
  double decay = 10.3e-9;
  double rise = 50e-12;
  double factor = 1. / (decay - rise);
  
  TF1 *shao = new TF1("shao1", "(x > [2]) * [0] * (exp(- (x-[2]) / [1]) - exp(- (x-[2]) / [3])) ", min, max);
  shao -> SetNpx(1000000);
  shao -> SetParameters(factor, decay, t_0, rise);
  cout << "Integral shao = " << shao -> Integral(min, max) << endl;
  cout << "factor = " << factor << endl;
  
  cout << "Maximum shao = " << shao -> GetMaximum() << endl;
  cout << "Maximum extracted = " << factor*(TMath::Exp(TMath::Log(rise/decay)/(1/rise-1/decay)/decay)-TMath::Exp(TMath::Log(rise/decay)/(1/rise-1/decay)/rise)) << endl;
  
  //double factor_2 = factor/2.;
  double t_1 = t_0 + cer_sigma;
  double cer_min = t_1 - 4 * cer_sigma; 
  double factor_2 = 1. / cer_sigma / TMath::Sqrt(2. * TMath::Pi()) * (1. + cer_sigma / 2. * TMath::Sqrt(2. * TMath::Pi()) * TMath::Erfc(1. / TMath::Sqrt(2) / cer_sigma * (t_1 - cer_min)));
  //double factor_2 = 1.;
  TF1 *cer = new TF1("cer", "(x > [0]) * [1] * exp(- (x-[2]) * (x-[2]) / 2. / [3] / [3])", min, max);
  cer -> SetNpx(1000000);
  cer -> SetParameters(cer_min, factor_2, t_1, cer_sigma);
  cer -> SetLineColor(kRed);
  cout << "Integral cer = " << cer -> Integral(min, max) << endl;
  cout << "factor_2 = " << factor_2 << endl;
  
  int nbin = 10000000;
  TH1D* prova = new TH1D("prova", "prova", nbin, min, max);
  for(int i =0; i<nbin; i++)
  {
    prova -> SetBinContent(i, cer->Eval(min + i*(max-min)/nbin));
  }
  TCanvas *c5 = new TCanvas("c5");
  prova -> Draw();
  cout << "TH1 test: " << prova -> Integral(0, nbin, "width") << endl;

  TRandom3* rand = new TRandom3();
  int extraction = 5000000;
  int nbins = 1000000;
  double sigma_irf = 66e-12;
  double sigma_trans = 100e-11;
  double norm_irf = 2. / sigma_irf / TMath::Sqrt(2) / TMath::Sqrt(TMath::Pi()) / TMath::Erfc(- sigma_trans / sigma_irf / TMath::Sqrt(2));
  TH1D* prova_smear_n = new TH1D("prova_smear_n", "prova_smear_n", nbins, min, max);
  TH1D* prova_smear = new TH1D("prova_smear", "prova_smear", nbins, min, max);
  TH1D* prova_irf = new TH1D("prova_irf", "prova_irf", nbins, min, max);
  TF1 *irf = new TF1("irf", "(x > 0.0) *[0] * exp(- (x-[1]) * (x-[1]) / 2. / [2] / [2])", min, max);
  irf -> SetNpx(1000000);
  irf -> SetParameters(norm_irf, sigma_trans, sigma_irf);
  for(int i =0; i<extraction; i++)
  {
    double extr = cer -> GetRandom();
    double extr_n = irf -> GetRandom();
    prova_smear_n -> Fill(extr);
    prova_smear -> Fill(extr + extr_n);
    prova_irf -> Fill(extr_n);
  }
  cout << "Integral irf = " << prova_irf -> Integral(0, nbins, "width") << endl;
  cout << "norm_irf = " << norm_irf << endl;
  TCanvas *c6 = new TCanvas("c6");
  prova_smear_n -> Draw();
  prova_smear -> Draw("same");
  prova_irf -> Draw("same");
  
  TCanvas *c7 = new TCanvas("c7");
  irf -> Draw();
  
  double LY = 10000;
  double CY = 20;
  double a = LY / CY * (1. / (1 + LY / CY));
  double b = 1 / (1 + LY / CY);
    
  TF1 *sum = new TF1("sum", "(x > [2]) * [0] * [8] * (exp(- (x-[2]) / [1]) - exp(- (x-[2]) / [3])) +  (x > [4]) * [5] * [9] * exp(- (x-[6]) * (x-[6]) / 2. / [7] / [7])", min, max);
  sum -> SetNpx(1000000);
  sum -> SetParameters(factor, decay, t_0, rise, cer_min, factor_2, t_1, cer_sigma, a, b);
  
  TCanvas *c1 = new TCanvas("c1");
  shao -> Draw();
  
  TCanvas *c2 = new TCanvas("c2");
  cer -> Draw();

  TCanvas *c3 = new TCanvas("c3");
  sum -> Draw();
  
  TCanvas *c4 = new TCanvas("c4");
  shao -> Draw();
  cer -> Draw("same");
  sum -> Draw("same");
        
}
