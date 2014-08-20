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
    double norm_c = b / l / TMath::Sqrt(2. * TMath::Pi()) * (1. + l / 2. * TMath::Sqrt(2. * TMath::Pi()) * TMath::Erfc(1. / TMath::Sqrt(2) / l * (l + 4.*tsk_c))); // cerenkov normalization 4
    //double norm_c = 2 / (l * TMath::Sqrt(2 * TMath::Pi()) * TMath::Erfc(- TMath::Sqrt(2))); 
    double norm_temp = b * 2 / TMath::Sqrt(2*(l*l+s*s)*TMath::Pi()) / TMath::Erfc((cer_min - tsk - theta) / TMath::Sqrt(2*(l*l+s*s)));
    


    //----------- IRF group ------------------
  
    //double sigma_mod = TMath::Sqrt(l*l + s*s);
    
    TF1 *irf = new TF1("irf", "(x > 0.0) * [0] * exp(- (x-[1]) * (x-[1]) / 2. / [2] / [2])", min, max);
    irf -> SetNpx(1000000);
     irf -> SetParameters(norm_irf*b, tsk+theta, s);
    irf -> SetLineColor(kRed);
    cout << "Integral irf = " << irf -> Integral(0, 10e-009) << endl;
    
    
    //----------- CONV group ------------------
    
    TF1 *conv_last = new TF1("conv_last", "TMath::Sqrt(TMath::Pi() / 2) * [0] * [1] * [2] * (-TMath::Exp(([0]*[0]-2*x*[3]+2*[3]*[4]+2*[3]*[6])/(2*[3]*[3]))*(TMath::Erf(([3]*([4]-x+[6]) + [0]*[0])/(TMath::Sqrt(2)*[0]*[3]))-TMath::Erf(([3]*[4] +[0]*[0])/(TMath::Sqrt(2)*[0]*[3])))+TMath::Exp(([0]*[0]-2*x*[5]+2*[5]*[4]+2*[5]*[6])/(2*[5]*[5]))*(TMath::Erf(([5]*([4]-x+[6]) + [0]*[0])/(TMath::Sqrt(2)*[0]*[5]))-TMath::Erf(([5]*[4] +[0]*[0])/(TMath::Sqrt(2)*[0]*[5]))))+1.0 / (TMath::Sqrt([7]*[7]+[0]*[0])) * TMath::Sqrt(TMath::Pi() / 2) * [7] * [0] * [8] * [1] * (TMath::Exp(-([4]+[6]-x)*([4]+[6]-x)/2/([7]*[7]+[0]*[0]))*TMath::Erf((-[4]*[7]*[7]+[6]*[0]*[0]+[7]*[7]*(x-[9])-[0]*[0]*x+[0]*[0]*(x-[9]))/(TMath::Sqrt(2)*[0]*[7]*TMath::Sqrt([7]*[7]+[0]*[0])))-TMath::Exp(-([4]+[6]-x)*([4]+[6]-x)/2/([7]*[7]+[0]*[0]))*TMath::Erf((-[4]*[7]*[7]+[6]*[0]*[0]-[0]*[0]*x)/(TMath::Sqrt(2)*[0]*[7]*TMath::Sqrt([7]*[7]+[0]*[0]))))", min, max);
    conv_last -> SetParameters(s, norm_irf, norm_s, t_d, tsk, t_r, theta, l, norm_c, cer_min);
    conv_last -> SetNpx(1000000);
    conv_last -> SetLineColor(kBlack);
    //cout << "Integral conv_last = " << conv_last -> Integral(min, max) << endl;
    
    TF1 *Cconv_last = new TF1("Cconv_last", "(x > [6]) * [7] * TMath::Exp(- (x-[2]-[1]) * (x-[2]-[1]) / 2. / ([3]*[3]+[8]*[8])) + [0]*(0.5 * (TMath::Erf((x-[1]-[2])/[3]/TMath::Sqrt(2))+TMath::Erf([2]/[3]/TMath::Sqrt(2)))+([4]/([5]-[4]) * 0.5 * TMath::Exp((-2*x*[4]+2*[4]*[1]+2*[4]*[2]-[3]*[3])/2/[4]/[4]) * (TMath::Erf((x-[1]-[2]-[3]*[3]/[4])/[3]/TMath::Sqrt(2))+TMath::Erf(([2]+[3]*[3]/[4])/[3]/TMath::Sqrt(2))))-([5]/([5]-[4]) * 0.5 * TMath::Exp((-2*x*[5]+2*[5]*[1]+2*[5]*[2]-[3]*[3])/2/[5]/[5]) * (TMath::Erf((x-[1]-[2]-[3]*[3]/[5])/[3]/TMath::Sqrt(2))+TMath::Erf(([2]+[3]*[3]/[5])/[3]/TMath::Sqrt(2)))))", min, max);    
    Cconv_last -> SetNpx(10000000);
    Cconv_last -> SetParameters(norm_s, theta, tsk, s, t_r, t_d, cer_min, norm_temp, l);
    Cconv_last -> SetLineColor(kRed);

    //cout << "Integral conv_last = " << conv_last -> Integral(min, max) << endl;

    //----------- SHAO group ------------------
    
    
    TF1 *shao = new TF1("shao", "(x > [3]) * [0] * (TMath::Exp(- (x - [3]) / [1]) - TMath::Exp(- (x - [3]) / [2]))", min, max);
    shao -> SetParameters(norm_s, t_d, t_r, theta);
    shao -> SetNpx(1000000);
    shao -> SetLineColor(kYellow);
    cout << "Integral shao = " << shao -> Integral(min, max) << endl;
    
    TF1 *shao_smear = new TF1("shao_smear", "TMath::Sqrt(TMath::Pi() / 2) * [0] * [1] * [2] * (-TMath::Exp(([0]*[0]-2*x*[3]+2*[3]*[4]+2*[3]*[6])/(2*[3]*[3]))*(TMath::Erf(([3]*([4]-x+[6]) + [0]*[0])/(TMath::Sqrt(2)*[0]*[3]))-TMath::Erf(([3]*[4] +[0]*[0])/(TMath::Sqrt(2)*[0]*[3])))+TMath::Exp(([0]*[0]-2*x*[5]+2*[5]*[4]+2*[5]*[6])/(2*[5]*[5]))*(TMath::Erf(([5]*([4]-x+[6]) + [0]*[0])/(TMath::Sqrt(2)*[0]*[5]))-TMath::Erf(([5]*[4] +[0]*[0])/(TMath::Sqrt(2)*[0]*[5]))))", min, max);
    shao_smear -> SetParameters(s, norm_irf, norm_s, t_d, tsk, t_r, theta);
    shao_smear -> SetNpx(1000000);
    shao_smear -> SetLineColor(kGreen);
    cout << "Integral shao smear = " << shao_smear -> Integral(min, max) << endl;
    
    TF1 *Cshao = new TF1("Cshao", "(x > [3]) * [0] * ([1] - [2] - [1]*TMath::Exp(- (x - [3]) / [1]) + [2]*TMath::Exp(- (x - [3]) / [2]))", min, max);
    Cshao -> SetParameters(norm_s, t_d, t_r, theta);
    Cshao -> SetNpx(1000000);
    Cshao -> SetLineColor(kBlack);
    //cout << "Integral shao = " << shao -> Integral(min, max) << endl;
    
    TF1 *Cshao_smear = new TF1("Cshao_smear", "[0] * (0.5 * (TMath::Erf((x-[1]-[2])/[3]/TMath::Sqrt(2))+TMath::Erf([2]/[3]/TMath::Sqrt(2)))+([4]/([5]-[4]) * 0.5 * TMath::Exp((-2*x*[4]+2*[4]*[1]+2*[4]*[2]+[3]*[3])/2/[4]/[4]) * (TMath::Erf((x-[1]-[2]-[3]*[3]/[4])/[3]/TMath::Sqrt(2))+TMath::Erf(([2]+[3]*[3]/[4])/[3]/TMath::Sqrt(2))))-([5]/([5]-[4]) * 0.5 * TMath::Exp((-2*x*[5]+2*[5]*[1]+2*[5]*[2]+[3]*[3])/2/[5]/[5]) * (TMath::Erf((x-[1]-[2]-[3]*[3]/[5])/[3]/TMath::Sqrt(2))+TMath::Erf(([2]+[3]*[3]/[5])/[3]/TMath::Sqrt(2)))))", min, max);
    shao_smear -> SetNpx(10000000);
    Cshao_smear -> SetParameters(norm_c*a, theta, tsk, s, t_r, t_d);
    shao_smear -> SetLineColor(kRed);
    
//      for(int i = 0; i < 100000000; i++)
//      {
//        if(i%1000000==0)
//        {
//          cout << i<< endl;
//        }
//        trys->Fill(Cshao_smear -> GetRandom());
//      }

    
    //----------- CER group ------------------
    
    TF1 *cer = new TF1("cer", " (x > [0]) * [1] * exp(- (x-[2]) * (x-[2]) / 2. / [3] / [3])", min, max);
    cer -> SetNpx(1000000);
    cer -> SetParameters(cer_min, norm_c, theta, l);
    cer -> SetLineColor(kYellow);
    cout << "Integral cer = " << cer -> Integral(cer_min, 10e-009) << endl;
    
    TF1 *cer_smear = new TF1("cer_smear", "1.0 / (TMath::Sqrt([0]*[0]+[1]*[1])) * TMath::Sqrt(TMath::Pi() / 2) * [0] * [1] * [2] * [6] * (TMath::Exp(-([3]+[5]-x)*([3]+[5]-x)/2/([0]*[0]+[1]*[1]))*TMath::Erf((-[3]*[0]*[0]+[5]*[1]*[1]+[0]*[0]*(x-[4])-[1]*[1]*x+[1]*[1]*(x-[4]))/(TMath::Sqrt(2)*[1]*[0]*TMath::Sqrt([0]*[0]+[1]*[1])))-TMath::Exp(-([3]+[5]-x)*([3]+[5]-x)/2/([0]*[0]+[1]*[1]))*TMath::Erf((-[3]*[0]*[0]+[5]*[1]*[1]-[1]*[1]*x)/(TMath::Sqrt(2)*[1]*[0]*TMath::Sqrt([0]*[0]+[1]*[1]))))", min, max);
    cer_smear -> SetNpx(1000000);
    cer_smear -> SetParameters(l, s, norm_c*b, tsk, cer_min, theta, norm_irf);
    cer_smear -> SetLineColor(kGreen);
    cout << "Integral cer_smear = " << cer_smear -> Integral(1.5e-009, 6e-009) << endl;
      
    TF1 *Ccer = new TF1("Ccer", " (x > [0]) * [1] * [3] * TMath::Sqrt(TMath::Pi()/2) * (TMath::Erf(([2]-[0])/[3]/TMath::Sqrt(2))-TMath::Erf(([2]-x)/[3]/TMath::Sqrt(2)))", min, max);
    Ccer -> SetNpx(1000000);
    Ccer -> SetParameters(cer_min, norm_c, theta, l);
    Ccer -> SetLineColor(kBlack);
    //cout << "Integral cer = " << cer -> Integral(cer_min, max) << endl;
    
    double norm_irf_easy = b / TMath::Sqrt(TMath::Pi()/2) / s / (1+TMath::Erf((tsk+theta)/TMath::Sqrt(2)/s));
    TF1 *Ccer_smear_easy = new TF1("cer_smear_easy", "[0] * TMath::Sqrt(TMath::Pi()/2) * (-[1]) * (TMath::Erf(([2]+[3]-x)/([1]*TMath::Sqrt(2)))-TMath::Erf(([2]+[3]-[4])/([1]*TMath::Sqrt(2))))", min, max);
    cer_smear_easy -> SetNpx(1000000);
    cer_smear_easy -> SetParameters(norm_irf_easy, s, tsk, theta, 0.0);
    cer_smear_easy -> SetLineColor(kBlack);
    //cout << "Integral cer_smear = " << cer_smear -> Integral(tsk+theta-cer_min, max) << endl;

        
    double C = a / TMath::Erfc(- tsk / s / TMath::Sqrt(2));
    
    TF1 *temp= new TF1("temp", "[0] * ((TMath::Erf((x-[1]-[2])/[3]/TMath::Sqrt(2))+TMath::Erf([2]/[3]/TMath::Sqrt(2)))+([4]/([5]-[4]) * TMath::Exp((-2*x*[4]+2*[4]*[1]+2*[4]*[2]+[3]*[3])/2/[4]/[4]) * (TMath::Erf((x-[1]-[2]-[3]*[3]/[4])/[3]/TMath::Sqrt(2))+TMath::Erf(([2]+[3]*[3]/[4])/[3]/TMath::Sqrt(2))))-([5]/([5]-[4]) * TMath::Exp((-2*x*[5]+2*[5]*[1]+2*[5]*[2]+[3]*[3])/2/[5]/[5]) * (TMath::Erf((x-[1]-[2]-[3]*[3]/[5])/[3]/TMath::Sqrt(2))+TMath::Erf(([2]+[3]*[3]/[5])/[3]/TMath::Sqrt(2)))))+[6] * TMath::Sqrt(TMath::Pi()/2) * (-[3]) * (TMath::Erf(([2]+[1]-x)/([3]*TMath::Sqrt(2)))-TMath::Erf(([2]+[1]-[7])/([3]*TMath::Sqrt(2))))", min, max);
    temp -> SetNpx(1000000);
    temp -> SetParameters(0., theta, tsk, s, t_r, t_d, norm_irf_easy,0.0);
    temp -> SetLineColor(kBlack);
    
    TF1 *temp1= new TF1("temp1", "[0] * ((TMath::Erf((x-[1]-[2])/[3]/TMath::Sqrt(2))+TMath::Erf([2]/[3]/TMath::Sqrt(2)))+([4]/([5]-[4]) * TMath::Exp((-2*x*[4]+2*[4]*[1]+2*[4]*[2]+[3]*[3])/2/[4]/[4]) * (TMath::Erf((x-[1]-[2]-[3]*[3]/[4])/[3]/TMath::Sqrt(2))+TMath::Erf(([2]+[3]*[3]/[4])/[3]/TMath::Sqrt(2))))-([5]/([5]-[4]) * TMath::Exp((-2*x*[5]+2*[5]*[1]+2*[5]*[2]+[3]*[3])/2/[5]/[5]) * (TMath::Erf((x-[1]-[2]-[3]*[3]/[5])/[3]/TMath::Sqrt(2))+TMath::Erf(([2]+[3]*[3]/[5])/[3]/TMath::Sqrt(2)))))+[6] * TMath::Sqrt(TMath::Pi()/2) * (-[3]) * (TMath::Erf(([2]+[1]-x)/([3]*TMath::Sqrt(2)))-TMath::Erf(([2]+[1]-[7])/([3]*TMath::Sqrt(2))))", min, max);
    temp1 -> SetNpx(1000000);
    temp1 -> SetParameters(C, theta, tsk, s, t_r, t_d, 0.,0.0);
    temp1 -> SetLineColor(kYellow);
    
    TF1 *CTot= new TF1("CTot", "[0] * ((TMath::Erf((x-[1]-[2])/[3]/TMath::Sqrt(2))+TMath::Erf([2]/[3]/TMath::Sqrt(2)))+([4]/([5]-[4]) * TMath::Exp((-2*x*[4]+2*[4]*[1]+2*[4]*[2]+[3]*[3])/2/[4]/[4]) * (TMath::Erf((x-[1]-[2]-[3]*[3]/[4])/[3]/TMath::Sqrt(2))+TMath::Erf(([2]+[3]*[3]/[4])/[3]/TMath::Sqrt(2))))-([5]/([5]-[4]) * TMath::Exp((-2*x*[5]+2*[5]*[1]+2*[5]*[2]+[3]*[3])/2/[5]/[5]) * (TMath::Erf((x-[1]-[2]-[3]*[3]/[5])/[3]/TMath::Sqrt(2))+TMath::Erf(([2]+[3]*[3]/[5])/[3]/TMath::Sqrt(2)))))+[6] * TMath::Sqrt(TMath::Pi()/2) * (-[3]) * (TMath::Erf(([2]+[1]-x)/([3]*TMath::Sqrt(2)))-TMath::Erf(([2]+[1]-[7])/([3]*TMath::Sqrt(2))))", min, max);
    CTot -> SetNpx(1000000);
    CTot -> SetParameters(C, theta, tsk, s, t_r, t_d, norm_irf_easy,0.0);
    CTot -> SetLineColor(kGreen);
    //cout << "Integral cer_smear = " << cer_smear -> Integral(tsk+theta-cer_min, max) << endl;
   
    int order = 100;

    double arraypam[14];
    arraypam[0] = s;
    arraypam[1] = norm_irf;
    arraypam[2] = norm_s;
    arraypam[3] = t_d;
    arraypam[4] = tsk;
    arraypam[5] = t_r;
    arraypam[6] = theta;
    arraypam[7] = l;
    arraypam[8] = norm_c*b;
    arraypam[9] = 0.0;
    arraypam[10] = LY+CY;
    arraypam[11] = order;
    arraypam[12] = C;
    arraypam[13] = norm_irf_easy;
    
    TF1 *n_phot = new TF1("n_phot", "TMath::Binomial([10], [11]) * [11] *(([12] * ((TMath::Erf((x-[6]-[4])/[0]/sqrt(2))+TMath::Erf([4]/[0]/sqrt(2)))+([5]/([3]-[5]) * exp((-2*x*[5]+2*[5]*[6]+2*[5]*[4]+[0]*[0])/2/[5]/[5]) * (TMath::Erf((x-[6]-[4]-[0]*[0]/[5])/[0]/sqrt(2))+TMath::Erf(([4]+[0]*[0]/[5])/[0]/sqrt(2))))-([3]/([3]-[5]) * TMath::Exp((-2*x*[3]+2*[3]*[6]+2*[3]*[4]+[0]*[0])/2/[3]/[3]) * (TMath::Erf((x-[6]-[4]-[0]*[0]/[3])/[0]/TMath::Sqrt(2))+TMath::Erf(([4]+[0]*[0]/[3])/[0]/sqrt(2)))))+[13] * sqrt(TMath::Pi()/2) * (-[0]) * (TMath::Erf(([4]+[6]-x)/([0]*sqrt(2)))-TMath::Erf(([4]+[6]-[9])/([0]*sqrt(2)))))  ^([11]-1) )*((1 - ( [12] * ((TMath::Erf((x-[6]-[4])/[0]/sqrt(2))+TMath::Erf([4]/[0]/sqrt(2)))+([5]/([3]-[5]) * exp((-2*x*[5]+2*[5]*[6]+2*[5]*[4]+[0]*[0])/2/[5]/[5]) * (TMath::Erf((x-[6]-[4]-[0]*[0]/[5])/[0]/sqrt(2))+TMath::Erf(([4]+[0]*[0]/[5])/[0]/sqrt(2))))-([3]/([3]-[5]) * exp((-2*x*[3]+2*[3]*[6]+2*[3]*[4]+[0]*[0])/2/[3]/[3]) * (TMath::Erf((x-[6]-[4]-[0]*[0]/[3])/[0]/sqrt(2))+TMath::Erf(([4]+[0]*[0]/[3])/[0]/sqrt(2)))))+[13] * sqrt(TMath::Pi()/2) * (-[0]) * (TMath::Erf(([4]+[6]-x)/([0]*sqrt(2)))-TMath::Erf(([4]+[6]-[9])/([0]*sqrt(2))))))^([10]-[11])) *(sqrt(TMath::Pi() / 2) * [0] * [1] * [2] * ( - exp(([0]*[0]-2*x*[3]+2*[3]*[4]+2*[3]*[6])/(2*[3]*[3]))*(TMath::Erf(([3]*([4]-x+[6]) + [0]*[0])/(sqrt(2)*[0]*[3]))-TMath::Erf(([3]*[4] +[0]*[0])/(sqrt(2)*[0]*[3])))+ TMath::Exp(([0]*[0]-2*x*[5]+2*[5]*[4]+2*[5]*[6])/(2*[5]*[5]))*(TMath::Erf(([5]*([4]-x+[6]) + [0]*[0])/(sqrt(2)*[0]*[5]))-TMath::Erf(([5]*[4] +[0]*[0])/(sqrt(2)*[0]*[5]))))+ 1.0 / (sqrt([7]*[7]+[0]*[0])) * sqrt(TMath::Pi() / 2) * [7] * [0] * [8] * [1] * (exp(-([4]+[6]-x)*([4]+[6]-x)/2/([7]*[7]+[0]*[0]))*(TMath::Erf((-[4]*[7]*[7]+[6]*[0]*[0]+[7]*[7]*(x-[9])-[0]*[0]*x+[0]*[0]*(x-[9]))/(sqrt(2)*[0]*[7]*sqrt([7]*[7]+[0]*[0])))-TMath::Erf((-[4]*[7]*[7]+[6]*[0]*[0]-[0]*[0]*x)/(sqrt(2)*[0]*[7]*sqrt([7]*[7]+[0]*[0]))))))", min, max);
    n_phot -> SetNpx(1000000);
    n_phot -> SetParameters(arraypam);
    n_phot -> SetLineColor(kGreen);
 
    TCanvas* c1 = new TCanvas("c1");
    n_phot -> Draw();

//  int qa = 2;
// TF1* test = new TF1("test", "sqrt(2.) * [0]", 0, 1000);
// test->SetParameter(0, qa);
// test->Draw();


//irf -> Draw();
//    cer_smear -> Draw("same");
     
     //       temp1 -> Draw();
//        Cshao_smear -> Draw("same");
//        shao_smear -> Draw("same");
//  
 //    Cshao -> Draw("same");
//     
 //TCanvas* c2 = new TCanvas("c2");
 //    Cshao_smear -> Draw();
 

  // TCanvas* c3 = new TCanvas("c3");
  // CTot -> Draw(); 
  // temp1 -> Draw("same");
  // temp -> Draw("same");
   //Cshao_smear -> Draw("same");
 // 
   //Cshao_smear -> Draw();
    //cer -> Draw();
    //shao_smear -> Draw("same");
  
   //cer_smear -> Draw("same"); 
   //irf->Draw("same");
// TCanvas* c2 = new TCanvas("c2");
// irf -> Draw();
//      cer_smear -> Draw("same");
//         temp -> Draw("same");

      //   TCanvas* c4 = new TCanvas("c4");
	// CTot -> Draw();
  // Cconv_last -> Draw();

    /*    
    TCanvas* c3 = new TCanvas("c3");
    shao -> Draw();
    shao_smear -> Draw("same");
    
    TCanvas *c4 = new TCanvas("c4");
    irf -> Draw();
    */
//     TCanvas *c5 = new TCanvas("c5");
//     conv_last-> Draw();
//     shao_smear -> Draw("same");
    
//     
//     TCanvas *c6 = new TCanvas("c6");
//     cer		-> Draw();
//     shao	-> Draw("same");
//     irf		-> Draw("same");
//     conv_last	-> Draw("same");
//     shao_smear -> Draw("same");
}