{
    double min = -10.;
    double max = 1000;
    
    double LY = 5000;
    double CY = 20;
    double a = LY / CY * (1. / (1 + LY / CY));
    double b = 1 / (1 + LY / CY);
    
    double theta = 2;
    
    double t_d = 30.3; // rise time 2
    double t_r = 0.07; //decay time 5

    double s = 0.066; // sigma trans 1 / 3
    double tsk = 1.; // 3 / 0
    
    double l = 0.005; // cerenkov sigma 2   
    double cer_min = theta; 
    double cer_mean = theta + 2*l;
    
    double norm_c 	= 1. / (l * sqrt(TMath::Pi() / 2) * (1 + TMath::Erf((cer_mean-cer_min) / sqrt(2) / l))); 
    double norm_irf	= 1. / (s * sqrt(TMath::Pi() / 2) * (1 + TMath::Erf(tsk / sqrt(2) / s)));
    double norm_s 	= 1. / (t_d-t_r) ; //normalization shao 4
  
    
    TF1 *cer = new TF1("cer", " (x > [0]) * [1] * exp(- (x-[2]) * (x-[2]) / 2. / [3] / [3])", min, max);
    cer -> SetNpx(10000000);
    cer -> SetParameters(cer_min, norm_c, cer_mean, l);
    cer -> SetLineColor(kRed);
    cout << "Integral cer = " << cer -> Integral(cer_min,2.3) << endl;
    cer -> Draw();
        TF1 *irf = new TF1("irf", "(x > 0.0) * [0] * exp(- (x-[1]) * (x-[1]) / 2. / [2] / [2])", min, max);
    irf -> SetNpx(1000000);
    irf -> SetParameters(norm_irf, tsk, s);
    irf -> SetLineColor(kGreen);
    cout << "Integral irf = " << irf -> Integral(0,10) << endl;
   irf->Draw("same");
    TF1 *cer_smear = new TF1("cer_smear", "(1.0 / (sqrt([0]*[0]+[1]*[1]))) * sqrt(TMath::Pi() / 2) * [0] * [1] * [2] * [6] * (exp(-([3]+[7]-x)*([3]+[7]-x)/2/([0]*[0]+[1]*[1]))*(TMath::Erf((-[3]*[0]*[0]+[7]*[1]*[1]+[0]*[0]*(x-[4])-[1]*[1]*x+[1]*[1]*(x-[4]))/(sqrt(2)*[0]*[1]*sqrt([0]*[0]+[1]*[1])))-TMath::Erf((-[3]*[0]*[0]+[7]*[1]*[1]-[1]*[1]*x)/(sqrt(2)*[0]*[1]*sqrt([0]*[0]+[1]*[1])))))", min, max);
    cer_smear -> SetNpx(1000000);
    cer_smear -> SetParameters(l, s, norm_c, tsk, cer_min, theta, norm_irf, cer_mean);
    cer_smear -> SetLineColor(kBlack);
    cout << "Integral cer_smear = " << cer_smear -> Integral(2.5, 4.) << endl;
    cer_smear->Draw("same");
}