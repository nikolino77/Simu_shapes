{
    double min = 0.;
    double max = 5e-9;
    
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

    
    //----------- IRF group ------------------
      
    TF1 *irf = new TF1("irf", "(x > 0.0) * [0] * exp(- (x-[1]) * (x-[1]) / 2. / [2] / [2])", min, max);
    irf -> SetNpx(1000000);
    irf -> SetParameters(norm_irf, tsk, s);
    irf -> SetLineColor(kGreen);
    cout << "Integral irf = " << irf -> Integral(0, 10e-009) << endl;
    
    
    //----------- CONV group ------------------
    
    TF1 *conv_last = new TF1("conv_last", "(1.0/(sqrt([7]*[7]+[0]*[0])))*sqrt(TMath::Pi()/2)*[7]*[0]*[8]*[1]*(exp(-([4]+[10]-x)*([4]+[10]-x)/2/([7]*[7]+[0]*[0]))*(TMath::Erf((-[4]*[7]*[7]+[10]*[0]*[0]+[7]*[7]*(x-[9])-[0]*[0]*x+[0]*[0]*(x-[9]))/(sqrt(2)*[7]*[0]*sqrt([7]*[7]+[0]*[0])))-TMath::Erf((-[4]*[7]*[7]+[10]*[0]*[0]-[0]*[0]*x)/(sqrt(2)*[7]*[0]*sqrt([7]*[7]+[0]*[0])))))+sqrt(TMath::Pi()/2)*[0]*[1]*[2]*(exp(([0]*[0]-2*x*[3]+2*[3]*[4]+2*[3]*[6])/(2*[3]*[3]))*(TMath::Erf(([3]*(x-[6]-[4])-[0]*[0])/(sqrt(2)*[0]*[3]))+TMath::Erf(([3]*[4] +[0]*[0])/(sqrt(2)*[0]*[3])))-exp(([0]*[0]-2*x*[5]+2*[5]*[4]+2*[5]*[6])/(2*[5]*[5]))*(TMath::Erf(([5]*(x-[6]-[4])-[0]*[0])/(sqrt(2)*[0]*[5]))+TMath::Erf(([5]*[4] +[0]*[0])/(sqrt(2)*[0]*[5]))))", min, max);
    conv_last -> SetParameters(s, norm_irf, norm_s*a, t_d, tsk, t_r, theta, l, norm_c*b, cer_min, cer_mean);  
    conv_last -> SetNpx(1000000);
    conv_last -> SetLineColor(kBlack);
    cout << "Integral conv_last = " << conv_last -> Integral(min, 400e-009) << endl;
    
    //----------- SHAO group ------------------
    
    TF1 *shao = new TF1("conv", "(x > [3]) * [0] * (exp(- (x - [3]) / [1]) - exp(- (x - [3]) / [2]))", min, max);
    shao -> SetParameters(norm_s, t_d, t_r, theta);
    shao -> SetNpx(1000000);
    shao -> SetLineColor(kYellow);
    cout << "Integral shao = " << shao -> Integral(0., 400e-009) << endl;

    TF1 *shao_smear = new TF1("shao_smear", "sqrt(TMath::Pi() / 2) * [0] * [1] * [2] * (-exp(([0]*[0]-2*x*[3]+2*[3]*[4]+2*[3]*[6])/(2*[3]*[3]))*(TMath::Erf(([3]*([4]-x+[6]) + [0]*[0])/(sqrt(2)*[0]*[3]))-TMath::Erf(([3]*[4] +[0]*[0])/(sqrt(2)*[0]*[3])))+exp(([0]*[0]-2*x*[5]+2*[5]*[4]+2*[5]*[6])/(2*[5]*[5]))*(TMath::Erf(([5]*([4]-x+[6]) + [0]*[0])/(sqrt(2)*[0]*[5]))-TMath::Erf(([5]*[4] +[0]*[0])/(sqrt(2)*[0]*[5]))))", min, max);
    shao_smear -> SetParameters(s, norm_irf, norm_s, t_d, tsk, t_r, theta);
    shao_smear -> SetNpx(1000000);
    shao_smear -> SetLineColor(kBlack);
    cout << "Integral shao smear = " << shao_smear -> Integral(0., 400e-009) << endl;
    
    TF1 *Cshao = new TF1("Cshao", "(x > [3]) * [0] * ([1] - [2] - [1]*exp(- (x - [3]) / [1]) + [2]*exp(- (x - [3]) / [2]))", min, max);
    Cshao -> SetParameters(norm_s, t_d, t_r, theta);
    Cshao -> SetNpx(1000000);
    Cshao -> SetLineColor(kBlack);
    
    TF1 *Cshao_smear = new TF1("Cshao_smear", "[0]*[3]*sqrt(TMath::Pi()/2)*((TMath::Erf((x-[1]-[2])/[3]/sqrt(2))+TMath::Erf([2]/[3]/sqrt(2)))+([4]/([5]-[4])*TMath::Exp((-2*x*[4]+2*[4]*[1]+2*[4]*[2]+[3]*[3])/2/[4]/[4])*(TMath::Erf((x-[1]-[2]-[3]*[3]/[4])/([3]*sqrt(2)))+TMath::Erf(([2]+[3]*[3]/[4])/[3]/sqrt(2))))-([5]/([5]-[4])*TMath::Exp((-2*x*[5]+2*[5]*[1]+2*[5]*[2]+[3]*[3])/2/[5]/[5]) * (TMath::Erf((x-[1]-[2]-[3]*[3]/[5])/[3]/sqrt(2))+TMath::Erf(([2]+[3]*[3]/[5])/[3]/sqrt(2)))))", min, max);
    Cshao_smear -> SetNpx(10000000);
    Cshao_smear -> SetParameters(norm_irf, theta, tsk, s, t_r, t_d);
    Cshao_smear -> SetLineColor(kRed);
    
    
    //----------- CER group ------------------
    
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
    
    TF1 *cer_smear_easy = new TF1("cer_smear_easy", "(x > 0.) * [0] * exp(- (x - [1] - [2]) * (x - [1] - [2]) / 2. / [3] / [3])", min, max);
    cer_smear_easy -> SetNpx(1000000);
    cer_smear_easy -> SetParameters(norm_irf_easy, tsk, cer_mean, s);
    cer_smear_easy -> SetLineColor(kGreen);
    cout << "Integral cer_smear_easy = " << cer_smear_easy -> Integral(0.e-009, 30e-009) << endl;
          
    TF1 *Ccer_smear_easy = new TF1("cer_smear_easy", "[0] * sqrt(TMath::Pi() / 2) * [1] * (TMath::Erf((x-[2]-[3])/([1]*sqrt(2))) + TMath::Erf(([2]+[3])/([1]*sqrt(2)))) ", min, max);
    Ccer_smear_easy -> SetNpx(1000000);
    Ccer_smear_easy -> SetParameters(norm_irf_easy, s, tsk, cer_mean);
    Ccer_smear_easy -> SetLineColor(kBlack);
        
    //----------------- Cumulative total ------------------------//
    
    TF1 *CTot = new TF1("CTot", "[7] * [0] * sqrt(TMath::Pi() / 2) * [3] * (TMath::Erf((x-[2]-[6])/([3]*sqrt(2))) + TMath::Erf(([2]+[6])/([3]*sqrt(2)))) + [8] * [0]*[3]*sqrt(TMath::Pi()/2)*((TMath::Erf((x-[1]-[2])/[3]/sqrt(2))+TMath::Erf([2]/[3]/sqrt(2)))+([4]/([5]-[4])*TMath::Exp((-2*x*[4]+2*[4]*[1]+2*[4]*[2]+[3]*[3])/2/[4]/[4])*(TMath::Erf((x-[1]-[2]-[3]*[3]/[4])/([3]*sqrt(2)))+TMath::Erf(([2]+[3]*[3]/[4])/[3]/sqrt(2))))-([5]/([5]-[4])*TMath::Exp((-2*x*[5]+2*[5]*[1]+2*[5]*[2]+[3]*[3])/2/[5]/[5]) * (TMath::Erf((x-[1]-[2]-[3]*[3]/[5])/[3]/sqrt(2))+TMath::Erf(([2]+[3]*[3]/[5])/[3]/sqrt(2)))))", min, max);
    CTot -> SetNpx(10000000);    
    CTot -> SetParameters(norm_irf, theta, tsk, s, t_r, t_d, cer_mean, b, a);
    CTot -> SetLineColor(kGreen);
    //cout << "Integral cer_smear = " << cer_smear -> Integral(tsk+theta-cer_min, max) << endl;
   
    
//     int order = 1;
// 
//     double arraypam[15];
//     arraypam[0] = norm_irf;
//     arraypam[1] = theta;
//     arraypam[2] = tsk;
//     arraypam[3] = s;
//     arraypam[4] = t_r;
//     arraypam[5] = t_d;
//     arraypam[6] = cer_mean;
//     arraypam[7] = b;
//     arraypam[8] = a;
//     arraypam[9] = LY+CY;
//     arraypam[10] = order;
//     arraypam[11] = l;
//     arraypam[12] = cer_min;
//     arraypam[13] = norm_c;
//     arraypam[14] = norm_s;
//            
//     TF1 *n_phot1 = new TF1("n_phot1", "TMath::Binomial([9], [10]) * [10] *(([7] * [0] * sqrt(TMath::Pi() / 2) * [3] * (TMath::Erf((x-[2]-[6])/([3]*sqrt(2))) + TMath::Erf(([2]+[6])/([3]*sqrt(2)))) + [8]*[0]*[3]*sqrt(TMath::Pi()/2)*((TMath::Erf((x-[1]-[2])/[3]/sqrt(2))+TMath::Erf([2]/[3]/sqrt(2)))+([4]/([5]-[4])*TMath::Exp((-2*x*[4]+2*[4]*[1]+2*[4]*[2]+[3]*[3])/2/[4]/[4])*(TMath::Erf((x-[1]-[2]-[3]*[3]/[4])/([3]*sqrt(2)))+TMath::Erf(([2]+[3]*[3]/[4])/[3]/sqrt(2))))-([5]/([5]-[4])*TMath::Exp((-2*x*[5]+2*[5]*[1]+2*[5]*[2]+[3]*[3])/2/[5]/[5])*(TMath::Erf((x-[1]-[2]-[3]*[3]/[5])/[3]/sqrt(2))+TMath::Erf(([2]+[3]*[3]/[5])/[3]/sqrt(2))))))^([10]-1))*( (1 - ([7] * [0] * sqrt(TMath::Pi() / 2) * [3] * (TMath::Erf((x-[2]-[6])/([3]*sqrt(2))) + TMath::Erf(([2]+[6])/([3]*sqrt(2)))) + [8]*[0]*[3]*sqrt(TMath::Pi()/2)*((TMath::Erf((x-[1]-[2])/[3]/sqrt(2))+TMath::Erf([2]/[3]/sqrt(2)))+([4]/([5]-[4])*TMath::Exp((-2*x*[4]+2*[4]*[1]+2*[4]*[2]+[3]*[3])/2/[4]/[4])*(TMath::Erf((x-[1]-[2]-[3]*[3]/[4])/([3]*sqrt(2)))+TMath::Erf(([2]+[3]
// *[3]/[4])/[3]/sqrt(2))))-([5]/([5]-[4])*TMath::Exp((-2*x*[5]+2*[5]*[1]+2*[5]*[2]+[3]*[3])/2/[5]/[5])*(TMath::Erf((x-[1]-[2]-[3]*[3]/[5])/[3]/sqrt(2))+TMath::Erf(([2]+[3]*[3]/[5])/[3]/sqrt(2)))))))^([9]-[10]))*( (1.0/(sqrt([11]*[11]+[3]*[3])))*sqrt(TMath::Pi()/2)*[11]*[3]*[13]*[7]*[0]*(exp(-([2]+[6]-x)*([2]+[6]-x)/2/([11]*[11]+[3]*[3]))*(TMath::Erf((-[2]*[11]*[11]+[6]*[3]*[3]+[11]*[11]*(x-[12])-[3]*[3]*x+[3]*[3]*(x-[12]))/(sqrt(2)*[11]*[3]*sqrt([11]*[11]+[3]*[3])))-TMath::Erf((-[2]*[11]*[11]+[6]*[3]*[3]-[3]*[3]*x)/(sqrt(2)*[11]*[3]*sqrt([11]*[11]+[3]*[3])))))+sqrt(TMath::Pi()/2)*[3]*[0]*[14]*[8]*(exp(([3]*[3]-2*x*[5]+2*[5]*[2]+2*[5]*[1])/(2*[5]*[5]))*(TMath::Erf(([5]*(x-[1]-[2])-[3]*[3])/(sqrt(2)*[3]*[5]))+TMath::Erf(([5]*[2] +[3]*[3])/(sqrt(2)*[3]*[5])))+exp(([3]*[3]-2*x*[4]+2*[4]*[2]+2*[4]*[1])/(2*[4]*[4]))*(TMath::Erf(([4]*(x-[1]-[2])-[3]*[3])/(sqrt(2)*[3]*[4]))+TMath::Erf(([4]*[2] +[3]*[3])/(sqrt(2)*[3]*[4])))))", min, max);
//     n_phot1 -> SetNpx(1000000);
//     n_phot1 -> SetParameters(arraypam);
//     n_phot1 -> SetLineColor(kGreen);
//     
//     TF1 *n_phot2 = new TF1("n_phot2", "TMath::Binomial([9], [10]) * [10] *(([7] * [0] * sqrt(TMath::Pi() / 2) * [3] * (TMath::Erf((x-[2]-[6])/([3]*sqrt(2))) + TMath::Erf(([2]+[6])/([3]*sqrt(2)))) + [8]*[0]*[3]*sqrt(TMath::Pi()/2)*((TMath::Erf((x-[1]-[2])/[3]/sqrt(2))+TMath::Erf([2]/[3]/sqrt(2)))+([4]/([5]-[4])*TMath::Exp((-2*x*[4]+2*[4]*[1]+2*[4]*[2]+[3]*[3])/2/[4]/[4])*(TMath::Erf((x-[1]-[2]-[3]*[3]/[4])/([3]*sqrt(2)))+TMath::Erf(([2]+[3]*[3]/[4])/[3]/sqrt(2))))-([5]/([5]-[4])*TMath::Exp((-2*x*[5]+2*[5]*[1]+2*[5]*[2]+[3]*[3])/2/[5]/[5])*(TMath::Erf((x-[1]-[2]-[3]*[3]/[5])/[3]/sqrt(2))+TMath::Erf(([2]+[3]*[3]/[5])/[3]/sqrt(2))))))^([10]-1))*( (1 - ([7] * [0] * sqrt(TMath::Pi() / 2) * [3] * (TMath::Erf((x-[2]-[6])/([3]*sqrt(2))) + TMath::Erf(([2]+[6])/([3]*sqrt(2)))) + [8]*[0]*[3]*sqrt(TMath::Pi()/2)*((TMath::Erf((x-[1]-[2])/[3]/sqrt(2))+TMath::Erf([2]/[3]/sqrt(2)))+([4]/([5]-[4])*TMath::Exp((-2*x*[4]+2*[4]*[1]+2*[4]*[2]+[3]*[3])/2/[4]/[4])*(TMath::Erf((x-[1]-[2]-[3]*[3]/[4])/([3]*sqrt(2)))+TMath::Erf(([2]+[3]
// *[3]/[4])/[3]/sqrt(2))))-([5]/([5]-[4])*TMath::Exp((-2*x*[5]+2*[5]*[1]+2*[5]*[2]+[3]*[3])/2/[5]/[5])*(TMath::Erf((x-[1]-[2]-[3]*[3]/[5])/[3]/sqrt(2))+TMath::Erf(([2]+[3]*[3]/[5])/[3]/sqrt(2)))))))^([9]-[10]))*( (1.0/(sqrt([11]*[11]+[3]*[3])))*sqrt(TMath::Pi()/2)*[11]*[3]*[13]*[7]*[0]*(exp(-([2]+[6]-x)*([2]+[6]-x)/2/([11]*[11]+[3]*[3]))*(TMath::Erf((-[2]*[11]*[11]+[6]*[3]*[3]+[11]*[11]*(x-[12])-[3]*[3]*x+[3]*[3]*(x-[12]))/(sqrt(2)*[11]*[3]*sqrt([11]*[11]+[3]*[3])))-TMath::Erf((-[2]*[11]*[11]+[6]*[3]*[3]-[3]*[3]*x)/(sqrt(2)*[11]*[3]*sqrt([11]*[11]+[3]*[3])))))+sqrt(TMath::Pi()/2)*[3]*[0]*[14]*[8]*(exp(([3]*[3]-2*x*[5]+2*[5]*[2]+2*[5]*[1])/(2*[5]*[5]))*(TMath::Erf(([5]*(x-[1]-[2])-[3]*[3])/(sqrt(2)*[3]*[5]))+TMath::Erf(([5]*[2] +[3]*[3])/(sqrt(2)*[3]*[5])))+exp(([3]*[3]-2*x*[4]+2*[4]*[2]+2*[4]*[1])/(2*[4]*[4]))*(TMath::Erf(([4]*(x-[1]-[2])-[3]*[3])/(sqrt(2)*[3]*[4]))+TMath::Erf(([4]*[2] +[3]*[3])/(sqrt(2)*[3]*[4])))))", min, max);
//     n_phot2 -> SetNpx(1000000);
//     arraypam[10] = 4.;
//     n_phot2 -> SetParameters(arraypam);
//     n_phot2 -> SetLineColor(kRed);
//     
//     TF1 *n_phot3 = new TF1("n_phot3", "TMath::Binomial([9], [10]) * [10] *(([7] * [0] * sqrt(TMath::Pi() / 2) * [3] * (TMath::Erf((x-[2]-[6])/([3]*sqrt(2))) + TMath::Erf(([2]+[6])/([3]*sqrt(2)))) + [8]*[0]*[3]*sqrt(TMath::Pi()/2)*((TMath::Erf((x-[1]-[2])/[3]/sqrt(2))+TMath::Erf([2]/[3]/sqrt(2)))+([4]/([5]-[4])*TMath::Exp((-2*x*[4]+2*[4]*[1]+2*[4]*[2]+[3]*[3])/2/[4]/[4])*(TMath::Erf((x-[1]-[2]-[3]*[3]/[4])/([3]*sqrt(2)))+TMath::Erf(([2]+[3]*[3]/[4])/[3]/sqrt(2))))-([5]/([5]-[4])*TMath::Exp((-2*x*[5]+2*[5]*[1]+2*[5]*[2]+[3]*[3])/2/[5]/[5])*(TMath::Erf((x-[1]-[2]-[3]*[3]/[5])/[3]/sqrt(2))+TMath::Erf(([2]+[3]*[3]/[5])/[3]/sqrt(2))))))^([10]-1))*( (1 - ([7] * [0] * sqrt(TMath::Pi() / 2) * [3] * (TMath::Erf((x-[2]-[6])/([3]*sqrt(2))) + TMath::Erf(([2]+[6])/([3]*sqrt(2)))) + [8]*[0]*[3]*sqrt(TMath::Pi()/2)*((TMath::Erf((x-[1]-[2])/[3]/sqrt(2))+TMath::Erf([2]/[3]/sqrt(2)))+([4]/([5]-[4])*TMath::Exp((-2*x*[4]+2*[4]*[1]+2*[4]*[2]+[3]*[3])/2/[4]/[4])*(TMath::Erf((x-[1]-[2]-[3]*[3]/[4])/([3]*sqrt(2)))+TMath::Erf(([2]+[3]
// *[3]/[4])/[3]/sqrt(2))))-([5]/([5]-[4])*TMath::Exp((-2*x*[5]+2*[5]*[1]+2*[5]*[2]+[3]*[3])/2/[5]/[5])*(TMath::Erf((x-[1]-[2]-[3]*[3]/[5])/[3]/sqrt(2))+TMath::Erf(([2]+[3]*[3]/[5])/[3]/sqrt(2)))))))^([9]-[10]))*( (1.0/(sqrt([11]*[11]+[3]*[3])))*sqrt(TMath::Pi()/2)*[11]*[3]*[13]*[7]*[0]*(exp(-([2]+[6]-x)*([2]+[6]-x)/2/([11]*[11]+[3]*[3]))*(TMath::Erf((-[2]*[11]*[11]+[6]*[3]*[3]+[11]*[11]*(x-[12])-[3]*[3]*x+[3]*[3]*(x-[12]))/(sqrt(2)*[11]*[3]*sqrt([11]*[11]+[3]*[3])))-TMath::Erf((-[2]*[11]*[11]+[6]*[3]*[3]-[3]*[3]*x)/(sqrt(2)*[11]*[3]*sqrt([11]*[11]+[3]*[3])))))+sqrt(TMath::Pi()/2)*[3]*[0]*[14]*[8]*(exp(([3]*[3]-2*x*[5]+2*[5]*[2]+2*[5]*[1])/(2*[5]*[5]))*(TMath::Erf(([5]*(x-[1]-[2])-[3]*[3])/(sqrt(2)*[3]*[5]))+TMath::Erf(([5]*[2] +[3]*[3])/(sqrt(2)*[3]*[5])))+exp(([3]*[3]-2*x*[4]+2*[4]*[2]+2*[4]*[1])/(2*[4]*[4]))*(TMath::Erf(([4]*(x-[1]-[2])-[3]*[3])/(sqrt(2)*[3]*[4]))+TMath::Erf(([4]*[2] +[3]*[3])/(sqrt(2)*[3]*[4])))))", min, max);
//     n_phot3 -> SetNpx(1000000);
//     arraypam[10] = 20.;
//     n_phot3 -> SetParameters(arraypam);
//     n_phot3 -> SetLineColor(kBlue);
//     
//     TF1 *n_phot4 = new TF1("n_phot4", "TMath::Binomial([9], [10]) * [10] *(([7] * [0] * sqrt(TMath::Pi() / 2) * [3] * (TMath::Erf((x-[2]-[6])/([3]*sqrt(2))) + TMath::Erf(([2]+[6])/([3]*sqrt(2)))) + [8]*[0]*[3]*sqrt(TMath::Pi()/2)*((TMath::Erf((x-[1]-[2])/[3]/sqrt(2))+TMath::Erf([2]/[3]/sqrt(2)))+([4]/([5]-[4])*TMath::Exp((-2*x*[4]+2*[4]*[1]+2*[4]*[2]+[3]*[3])/2/[4]/[4])*(TMath::Erf((x-[1]-[2]-[3]*[3]/[4])/([3]*sqrt(2)))+TMath::Erf(([2]+[3]*[3]/[4])/[3]/sqrt(2))))-([5]/([5]-[4])*TMath::Exp((-2*x*[5]+2*[5]*[1]+2*[5]*[2]+[3]*[3])/2/[5]/[5])*(TMath::Erf((x-[1]-[2]-[3]*[3]/[5])/[3]/sqrt(2))+TMath::Erf(([2]+[3]*[3]/[5])/[3]/sqrt(2))))))^([10]-1))*( (1 - ([7] * [0] * sqrt(TMath::Pi() / 2) * [3] * (TMath::Erf((x-[2]-[6])/([3]*sqrt(2))) + TMath::Erf(([2]+[6])/([3]*sqrt(2)))) + [8]*[0]*[3]*sqrt(TMath::Pi()/2)*((TMath::Erf((x-[1]-[2])/[3]/sqrt(2))+TMath::Erf([2]/[3]/sqrt(2)))+([4]/([5]-[4])*TMath::Exp((-2*x*[4]+2*[4]*[1]+2*[4]*[2]+[3]*[3])/2/[4]/[4])*(TMath::Erf((x-[1]-[2]-[3]*[3]/[4])/([3]*sqrt(2)))+TMath::Erf(([2]+[3]
// *[3]/[4])/[3]/sqrt(2))))-([5]/([5]-[4])*TMath::Exp((-2*x*[5]+2*[5]*[1]+2*[5]*[2]+[3]*[3])/2/[5]/[5])*(TMath::Erf((x-[1]-[2]-[3]*[3]/[5])/[3]/sqrt(2))+TMath::Erf(([2]+[3]*[3]/[5])/[3]/sqrt(2)))))))^([9]-[10]))*( (1.0/(sqrt([11]*[11]+[3]*[3])))*sqrt(TMath::Pi()/2)*[11]*[3]*[13]*[7]*[0]*(exp(-([2]+[6]-x)*([2]+[6]-x)/2/([11]*[11]+[3]*[3]))*(TMath::Erf((-[2]*[11]*[11]+[6]*[3]*[3]+[11]*[11]*(x-[12])-[3]*[3]*x+[3]*[3]*(x-[12]))/(sqrt(2)*[11]*[3]*sqrt([11]*[11]+[3]*[3])))-TMath::Erf((-[2]*[11]*[11]+[6]*[3]*[3]-[3]*[3]*x)/(sqrt(2)*[11]*[3]*sqrt([11]*[11]+[3]*[3])))))+sqrt(TMath::Pi()/2)*[3]*[0]*[14]*[8]*(exp(([3]*[3]-2*x*[5]+2*[5]*[2]+2*[5]*[1])/(2*[5]*[5]))*(TMath::Erf(([5]*(x-[1]-[2])-[3]*[3])/(sqrt(2)*[3]*[5]))+TMath::Erf(([5]*[2] +[3]*[3])/(sqrt(2)*[3]*[5])))+exp(([3]*[3]-2*x*[4]+2*[4]*[2]+2*[4]*[1])/(2*[4]*[4]))*(TMath::Erf(([4]*(x-[1]-[2])-[3]*[3])/(sqrt(2)*[3]*[4]))+TMath::Erf(([4]*[2] +[3]*[3])/(sqrt(2)*[3]*[4])))))", min, max);
//     n_phot4 -> SetNpx(1000000);
//     arraypam[10] = 1000.;
//     n_phot4 -> SetParameters(arraypam);
//     n_phot4 -> SetLineColor(kYellow);
//  
//     
//     TCanvas c1;
//     c1.Divide(2,2);
//     c1.cd(1);
//     n_phot1 -> Draw();
//     c1.cd(2);
//     n_phot2 -> Draw();
//     c1.cd(3);
//     n_phot3 -> Draw();
//     c1.cd(4);
//     n_phot4 -> Draw();
//     
//     TCanvas* c2 = new TCanvas("c2");
//     n_phot1 -> Draw();
//     n_phot2 -> Draw("same");
//     n_phot3 -> Draw("same");
//     n_phot4 -> Draw("same");

    //  int qa = 2;
// TF1* test = new TF1("test", "sqrt(2.) * [0]", 0, 1000);
// test->SetParameter(0, qa);
// test->Draw();

   //  TCanvas* c1 = new TCanvas("c1");
     //cer -> Draw();
//      cer_smear -> Draw();
//      irf->Draw("same");
//      cer_smear_easy->Draw("same");
     //       temp1 -> Draw();
//        Cshao_smear -> Draw("same");
//        shao_smear -> Draw("same");
//  
 //    Cshao -> Draw("same");
//     
//  TCanvas* c2 = new TCanvas("c2");
//  CTot-> Draw();    
//  Ccer_smear_easy->Draw("same");
//  Cshao_smear -> Draw("same");
     
 

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