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

int main(int argc,char** argv)
{
    double LY 	= atoi(argv[1]);
    double CY 	= atoi(argv[2]);
    double t_d 	= atof(argv[3]);
    double t_r	= atof(argv[4]);
    double s	= atof(argv[5]);
    double l	= atof(argv[6]);
    double res	= atof(argv[7]);
    
//    cout << LY << endl;
     cout << CY << endl;
//     cout << t_d << endl;
//     cout << t_r << endl;
//     cout << s << endl;
//     cout << l << endl;
//     cout << res << endl;

    string sLY 		= argv[1];
    string sCY 		= argv[2];
    string st_d 	= argv[3];
    string st_r		= argv[4];
    string ss		= argv[5];
    string sl		= argv[6];
    string sres		= argv[7];
    
    string path_sim("./sim");
    string path_dat("./dat");
    string filename 	= "_" + sLY + "_" + sCY + "_" + st_d + "_" + st_r + "_" + ss + "_" + sl + "_" + sres;
    string root		= path_sim + filename + ".root";
    string dat 		= path_dat + filename + ".txt";

    TFile* hfile = new TFile(root.c_str(),"RECREATE");
    ofstream out;
    out.open(dat.c_str());

    TF1 *gaus = new TF1("gaus", "TMath::Gaus(x, [0], [1])", 0, 100000);
    gaus -> SetNpx(1000000);
    gaus -> SetParameters(LY + CY, (LY + CY) * res / 2.355);

    gaus -> Write();
  
    hfile -> Close();
    return 0;
}