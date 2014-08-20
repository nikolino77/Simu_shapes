{
    double min = 0.;
    double max = 1000e-9;
    
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
    
    TF1 *shao = new TF1("shao", "(x > [3]) * [0] * (TMath::Exp(- (x - [3]) / [1]) - TMath::Exp(- (x - [3]) / [2]))", min, max);
    shao -> SetParameters(norm_s, t_d, t_r, theta);
    shao -> SetNpx(1000000);
    shao -> SetLineColor(kYellow);
    cout << "Integral shao = " << shao -> Integral(0., 400e-009) << endl;

    TF1 *shao_smear = new TF1("shao_smear", "TMath::Sqrt(TMath::Pi() / 2) * [0] * [1] * [2] * (-TMath::Exp(([0]*[0]-2*x*[3]+2*[3]*[4]+2*[3]*[6])/(2*[3]*[3]))*(TMath::Erf(([3]*([4]-x+[6]) + [0]*[0])/(TMath::Sqrt(2)*[0]*[3]))-TMath::Erf(([3]*[4] +[0]*[0])/(TMath::Sqrt(2)*[0]*[3])))+TMath::Exp(([0]*[0]-2*x*[5]+2*[5]*[4]+2*[5]*[6])/(2*[5]*[5]))*(TMath::Erf(([5]*([4]-x+[6]) + [0]*[0])/(TMath::Sqrt(2)*[0]*[5]))-TMath::Erf(([5]*[4] +[0]*[0])/(TMath::Sqrt(2)*[0]*[5]))))", min, max);
    shao_smear -> SetParameters(s, norm_irf, norm_s, t_d, tsk, t_r, theta);
    shao_smear -> SetNpx(1000000);
    shao_smear -> SetLineColor(kBlack);
    cout << "Integral shao smear = " << shao_smear -> Integral(0., 400e-009) << endl;
    
    TF1 *cer = new TF1("cer", " (x > [0]) * [1] * exp(- (x-[2]) * (x-[2]) / 2. / [3] / [3])", min, max);
    cer -> SetNpx(1000000);
    cer -> SetParameters(cer_min, norm_c, cer_mean, l);
    cer -> SetLineColor(kRed);
    //cout << "Integral cer = " << cer -> Integral(1.,1.01e-009) << endl;
    
    TF1 *cer_smear = new TF1("cer_smear", "(1.0 / (sqrt([0]*[0]+[1]*[1]))) * sqrt(TMath::Pi() / 2) * [0] * [1] * [2] * [6] * (exp(-([3]+[7]-x)*([3]+[7]-x)/2/([0]*[0]+[1]*[1]))*(TMath::Erf((-[3]*[0]*[0]+[7]*[1]*[1]+[0]*[0]*(x-[4])-[1]*[1]*x+[1]*[1]*(x-[4]))/(sqrt(2)*[0]*[1]*sqrt([0]*[0]+[1]*[1])))-TMath::Erf((-[3]*[0]*[0]+[7]*[1]*[1]-[1]*[1]*x)/(sqrt(2)*[0]*[1]*sqrt([0]*[0]+[1]*[1])))))", min, max);
    cer_smear -> SetNpx(1000000);
    cer_smear -> SetParameters(l, s, norm_c, tsk, cer_min, theta, norm_irf, cer_mean);
    cer_smear -> SetLineColor(kBlack);
    cout << "Integral cer_smear = " << cer_smear -> Integral(0.e-009, 30e-009) << endl;
    
    TF1 *sum = new TF1("sum", "(x > [3]) * [0] * (TMath::Exp(- (x - [3]) / [1]) - TMath::Exp(- (x - [3]) / [2])) +  (x > [4]) * [5] * exp(- (x-[6]) * (x-[6]) / 2. / [7] / [7])", min, max);
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
    cout << "Integral conv_last = " << conv_last -> Integral(min, 1000e-009) << endl;

    /*
    TH1F* test_conv = new TH1F("test_conv", "test_conv", 200000, 0, 200e-009);
    
    TH1F* test_sum = new TH1F("test_sum", "test_sum", 200000, 0, 200e-009);
    
    for(int i = 0; i< 10000000; i++)
    {
      if(i%100000==0)
      {
	cout << "Evento " << i << endl;
      }
      test_conv -> Fill(conv_last->GetRandom());
      test_sum -> Fill(sum->GetRandom() + irf ->GetRandom());
    }
   
     TCanvas *c1 = new TCanvas("c1");
     test_conv->Draw();
     test_sum ->Draw("same");*/
    
//      TCanvas* c1 = new TCanvas("c1");
//      shao -> Draw();
//      cer -> Draw("same");
// 
     TCanvas* c2 = new TCanvas("c2");
     conv_last -> Draw();
     sum -> Draw("same");

//     
//     TCanvas* c3 = new TCanvas("c3");
//     shao -> Draw();
//     shao_smear -> Draw("same");
//     
//     TCanvas *c4 = new TCanvas("c4");
//     irf -> Draw();
     
     TCanvas *c5 = new TCanvas("c5");
     //conv_last-> Draw();
     cer_smear -> Draw();
     shao_smear -> Draw("same");
     //          
//      TCanvas *c6 = new TCanvas("c6");
//      sum		-> Draw();
//      shao	-> Draw("same");
//      cer        -> Draw("same");
     //     irf		-> Draw("same");
//     conv_last	-> Draw("same");
//     shao_smear -> Draw("same");
}