{
    double min = 0.;
    double max = 200e-009;
    
    double LY = 20000;
    double CY = 20;
    double a = LY / CY * (1. / (1 + LY / CY));
    double b = 1 / (1 + LY / CY);
    
    double theta = 200e-12;
    
    double t_d = 30.3e-9; // rise time 2
    double t_r = 70e-12; //decay time 5
    double norm_s = a * 1.0 / (t_d-t_r) ; //normalization shao 4
    double s = 66e-12; // sigma trans 1 / 3
    double tsk = 2e-9; // 3 / 0
    
    double norm_irf = 2. / s / TMath::Sqrt(2) / TMath::Sqrt(TMath::Pi()) / TMath::Erfc(- tsk / s / TMath::Sqrt(2));
    double l = 10e-12; // cerenkov sigma 2
    double tsk_c = tsk - l; // delay cerenkov-scint 1
    double cer_min = theta - 2 * l; 
    //double norm_c = b / l / TMath::Sqrt(2. * TMath::Pi()) * (1. + l / 2. * TMath::Sqrt(2. * TMath::Pi()) * TMath::Erfc(1. / TMath::Sqrt(2) / l * (l + 4.*tsk_c))); // cerenkov normalization 4
    double norm_c = b * 2 / (l * TMath::Sqrt(2 * TMath::Pi()) * TMath::Erfc(- TMath::Sqrt(2))); 
    
    TF1 *shao = new TF1("conv", "(x > [3]) * [0] * (TMath::Exp(- (x - [3]) / [1]) - TMath::Exp(- (x - [3]) / [2]))", min, max);
    shao -> SetParameters(norm_s, t_d, t_r, theta);
    shao -> SetNpx(1000000);
    shao -> SetLineColor(kYellow);
    //cout << "Integral shao = " << shao -> Integral(min, max) << endl;

    TF1 *shao_smear = new TF1("shao_smear", "TMath::Sqrt(TMath::Pi() / 2) * [0] * [1] * [2] * (-TMath::Exp(([0]*[0]-2*x*[3]+2*[3]*[4]+2*[3]*[6])/(2*[3]*[3]))*(TMath::Erf(([3]*([4]-x+[6]) + [0]*[0])/(TMath::Sqrt(2)*[0]*[3]))-TMath::Erf(([3]*[4] +[0]*[0])/(TMath::Sqrt(2)*[0]*[3])))+TMath::Exp(([0]*[0]-2*x*[5]+2*[5]*[4]+2*[5]*[6])/(2*[5]*[5]))*(TMath::Erf(([5]*([4]-x+[6]) + [0]*[0])/(TMath::Sqrt(2)*[0]*[5]))-TMath::Erf(([5]*[4] +[0]*[0])/(TMath::Sqrt(2)*[0]*[5]))))", min, max);
    shao_smear -> SetParameters(s, norm_irf, norm_s, t_d, tsk, t_r, theta);
    shao_smear -> SetNpx(1000000);
    shao_smear -> SetLineColor(kBlack);
    //cout << "Integral shao smear = " << shao_smear -> Integral(min, max) << endl;
    
    TF1 *cer = new TF1("cer", " (x > [0]) * [1] * exp(- (x-[2]) * (x-[2]) / 2. / [3] / [3])", min, max);
    cer -> SetNpx(1000000);
    cer -> SetParameters(cer_min, norm_c, theta, l);
    cer -> SetLineColor(kRed);
    //cout << "Integral cer = " << cer -> Integral(cer_min, max) << endl;
    
    TF1 *cer_smear = new TF1("cer_smear", "1.0 / (TMath::Sqrt([0]*[0]+[1]*[1])) * TMath::Sqrt(TMath::Pi() / 2) * [0] * [1] * [2] * [6] * (TMath::Exp(-([3]+[5]-x)*([3]+[5]-x)/2/([0]*[0]+[1]*[1]))*TMath::Erf((-[3]*[0]*[0]+[5]*[1]*[1]+[0]*[0]*(x-[4])-[1]*[1]*x+[1]*[1]*(x-[4]))/(TMath::Sqrt(2)*[1]*[0]*TMath::Sqrt([0]*[0]+[1]*[1])))-TMath::Exp(-([3]+[5]-x)*([3]+[5]-x)/2/([0]*[0]+[1]*[1]))*TMath::Erf((-[3]*[0]*[0]+[5]*[1]*[1]-[1]*[1]*x)/(TMath::Sqrt(2)*[1]*[0]*TMath::Sqrt([0]*[0]+[1]*[1]))))", min, max);
    cer_smear -> SetNpx(1000000);
    cer_smear -> SetParameters(l, s, norm_c, tsk, cer_min, theta, norm_irf);
    cer_smear -> SetLineColor(kBlack);
    //cout << "Integral cer_smear = " << cer_smear -> Integral(tsk+theta-cer_min, max) << endl;
    
    TF1 *irf = new TF1("irf", "(x > 0.0) * [0] * exp(- (x-[1]) * (x-[1]) / 2. / [2] / [2])", min, max);
    irf -> SetNpx(1000000);
    irf -> SetParameters(norm_irf, tsk, s);
    irf -> SetLineColor(kGreen);
    //cout << "Integral irf = " << irf -> Integral(-10, max) << endl;
    
    TF1 *conv_last = new TF1("conv_last", "TMath::Sqrt(TMath::Pi() / 2) * [0] * [1] * [2] * (-TMath::Exp(([0]*[0]-2*x*[3]+2*[3]*[4]+2*[3]*[6])/(2*[3]*[3]))*(TMath::Erf(([3]*([4]-x+[6]) + [0]*[0])/(TMath::Sqrt(2)*[0]*[3]))-TMath::Erf(([3]*[4] +[0]*[0])/(TMath::Sqrt(2)*[0]*[3])))+TMath::Exp(([0]*[0]-2*x*[5]+2*[5]*[4]+2*[5]*[6])/(2*[5]*[5]))*(TMath::Erf(([5]*([4]-x+[6]) + [0]*[0])/(TMath::Sqrt(2)*[0]*[5]))-TMath::Erf(([5]*[4] +[0]*[0])/(TMath::Sqrt(2)*[0]*[5]))))+1.0 / (TMath::Sqrt([7]*[7]+[0]*[0])) * TMath::Sqrt(TMath::Pi() / 2) * [7] * [0] * [8] * [1] * (TMath::Exp(-([4]+[6]-x)*([4]+[6]-x)/2/([7]*[7]+[0]*[0]))*TMath::Erf((-[4]*[7]*[7]+[6]*[0]*[0]+[7]*[7]*(x-[9])-[0]*[0]*x+[0]*[0]*(x-[9]))/(TMath::Sqrt(2)*[0]*[7]*TMath::Sqrt([7]*[7]+[0]*[0])))-TMath::Exp(-([4]+[6]-x)*([4]+[6]-x)/2/([7]*[7]+[0]*[0]))*TMath::Erf((-[4]*[7]*[7]+[6]*[0]*[0]-[0]*[0]*x)/(TMath::Sqrt(2)*[0]*[7]*TMath::Sqrt([7]*[7]+[0]*[0]))))", min, max);
    conv_last -> SetParameters(s, norm_irf, norm_s, t_d, tsk, t_r, theta, l, norm_c, cer_min);
    conv_last -> SetNpx(1000000);
    conv_last -> SetLineColor(kBlack);
    //cout << "Integral conv_last = " << conv_last -> Integral(min, max) << endl;

    TCanvas* c1 = new TCanvas("c1");
    shao -> Draw();
    cer -> Draw("same");

    TCanvas* c2 = new TCanvas("c2");
    cer -> Draw();
    cer_smear -> Draw("same");
    
    TCanvas* c3 = new TCanvas("c3");
    shao -> Draw();
    shao_smear -> Draw("same");
    
    TCanvas *c4 = new TCanvas("c4");
    irf -> Draw();
    
    TCanvas *c5 = new TCanvas("c5");
    conv_last-> Draw();
    shao_smear -> Draw("same");
    
//     
//     TCanvas *c6 = new TCanvas("c6");
//     cer		-> Draw();
//     shao	-> Draw("same");
//     irf		-> Draw("same");
//     conv_last	-> Draw("same");
//     shao_smear -> Draw("same");
}