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

    TFile* hfile = new TFile("simu_crao.root","OPEN");
    TFile* hfile_out = new TFile("simu_crao_out_bis.root","RECREATE");
    
    TTree *Singles 	= (TTree*)hfile->Get("tree");
    cout << "Tree " << Singles -> GetName() << endl;
    
    vector<double> *shao_stamps = new vector<double>();
    Singles->SetBranchAddress("shao_stamps",&shao_stamps);
     
    vector<double >* shao_smear_stamps = new vector<double>();
    Singles->SetBranchAddress("shao_smear_stamps",&shao_smear_stamps);

    vector<double >* conv_stamps = new vector<double>();
    Singles->SetBranchAddress("conv_stamps",&conv_stamps);
   
    //vector<double >* sum_stamps = new vector<double>();
    //Singles->SetBranchAddress("s_stamps",&shao_stamps);
    
    vector<double >* sum_stamps = new vector<double>();
    Singles->SetBranchAddress("sum_stamps",&sum_stamps);
    
    TH1D* first = new TH1D("first", "first", 200000, 0, 200e-09);
    TH1D* second = new TH1D("second", "second", 200000, 0, 200e-09);
    TH1D* third = new TH1D("third", "third", 200000, 0, 200e-09);
    TH1D* fourth = new TH1D("fourth", "fourth", 200000, 0, 200e-09);
    //TH1D* fifth = new TH1D("fifth", "fifth", 200000, 0, 200e-09);
   
    Int_t nentries = Singles->GetEntries(); 
    cout << "Num entries = " << nentries << endl;	
    
    for(int i = 0; i < nentries; i++)
    {
      if(i%500 == 0)
      {
        std::cout << "Evento " << i << std::endl;
      }
      
      Singles->GetEntry(i);
      
      first -> Fill(shao_stamps->at(0));
      second->Fill(shao_smear_stamps->at(0));
      third->Fill(conv_stamps->at(0));
      fourth->Fill(sum_stamps->at(0));
      //fifth->Fill(shao_stamps->at(4));
   
    }
    
    hfile -> Close();
    
    hfile_out -> cd();
    first ->Write();
    second ->Write();
    third -> Write();
    fourth ->Write();
    //fifth->Write();
    hfile_out -> Close();
     
    return 0;

}