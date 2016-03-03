#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TStyle.h"

#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TF1.h"
#include "TString.h"
#include "TLatex.h"
#include "TMath.h"
#include "TROOT.h"

#include "SiWEcalSSSimHit.hh"

//===============================================================================
//Declaration here definition after main
std::string IntToString (int number1);
Double_t langaufun(Double_t *x, Double_t *par);
Double_t langaufungaus(Double_t *x, Double_t *par);

int main(){//main

  //======================================================================================
  //The root files with the individual plots for all test beam configs 
  const int numberoffiles = 11;
  //======================================================================================
  //3 layers times 4 Si Pads
  const int numberofpads = 12; 
  //======================================================================================
  //Particle type and energy
  TString particle = "mu-"; //mu+
  int particleenergy = 15; // 15 30 50

  //======================================================================================
  //Read all files and tree
  TString filename[numberoffiles];
  TString plotBase[numberoffiles];
  for (int k=0; k<numberoffiles; k++){
    TString fp = "SiWEcal_" + IntToString(k+1);	
    TString sp = "_" + particle;
    TString tp = "_" + IntToString(particleenergy);
    TString fop = "_GeV";
    filename[k] = fp + sp + tp + fop;
    plotBase[k] = "PLOTS/" ;
    //std::cout << filename[k] << std::endl;
  }

  //The tree in the files that we want to get
  TTree *lTree[numberoffiles]; 

  TFile* files[numberoffiles];
  for (int k=0; k<numberoffiles; k++){
    // files[k]= new TFile(filename[k]);
    //    files[k] = TFile::Open("/afs/cern.ch/work/a/apsallid/CMS/Geant4/SiWEcal/Parallel/testingfile/"+filename[k]+".root");
    files[k] = TFile::Open("/tmp/apsallid/Configs/"+filename[k]+"*.root");
    lTree[k] = (TTree*) files[k]->Get("SiWEcalSSTree");
    if (!lTree[k]){
      std::cout << " -- Error, tree cannot be opened. Exiting..." << std::endl;
      return 1;
    }
  }
  //Set the branches here
  std::vector<SiWEcalSSSimHit> * simhitvec = 0;
  for (int k=0; k<numberoffiles; k++){
    lTree[k]->SetBranchAddress("SiWEcalSSSimHitVec",&simhitvec);
  }
  // TCanvas *myc = new TCanvas("myc","myc",500,500);

  //======================================================================================
  // Histos
  TH2F *h_nHits = new TH2F("h_nHits","Number of hits vs pads",numberofpads,0,numberofpads,20,0,20);
  TH1F *h_hitEnergy_si = new TH1F("h_hitEnergy_si",";E (MeV);SimHits",250,0,1);
  std::vector<TH1F*> h_hitEnergy_si_pad;
  h_hitEnergy_si_pad.clear();
  for (Int_t k=0; k<numberofpads; k++){
    h_hitEnergy_si_pad.push_back(new TH1F(("h_hitEnergy_si_pad_" + IntToString(k+1)).c_str(),";E (MeV);SimHits",250,0,1));
  }
  //======================================================================================
  //The files that we will store the results of this analysis
  TString res[numberoffiles];
  for (int k=0; k<numberoffiles; k++){
    res[k] = filename[k] + "_final.root";
    //std::cout << res[k] << std::endl;
  }
  
  TFile* results[numberoffiles];
  for (int k=0; k<numberoffiles; k++){
    results[k]= new TFile(res[k],"recreate");
    std::cout << "Results file " << res[k] << std::endl;

    //======================================================================================
    // Loop on entries
    for (Int_t ievt=0; ievt<lTree[k]->GetEntries(); ievt++) {
      lTree[k]->GetEntry(ievt);
      if (ievt%(1000)==0) std::cout << "Entry " << ievt << std::endl;
 
      std::vector<unsigned> nHits;
      nHits.resize(numberofpads,0);
      unsigned nHitsall = 0;
      //loop on hits
      for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){
	SiWEcalSSSimHit lHit = (*simhitvec)[iH];
	double energy = lHit.energy();
	double layer = (4*lHit.layer()) + lHit.silayer(); 
	//std::cout << layer << std::endl;
	
	if (energy>0) {
	  h_hitEnergy_si->Fill(energy);
	  nHits[layer]++;
	  //Fill histos in each pad
	  h_hitEnergy_si_pad[layer]->Fill(energy);
	}
      }//loop on hits
      
      // loop to plot total hits in each layer and count all hits of even
      for (unsigned iL(0); iL<nHits.size(); ++iL){ 
	h_nHits->Fill(iL,nHits[iL]);
	++nHitsall;
      }
    


    } //  Loop on entries


    //======================================================================================
    //Write histos 
    h_nHits->Write();
    h_hitEnergy_si->Write();
    for (Int_t l=0; l<numberofpads; l++){
      h_hitEnergy_si_pad[l]->Write();
    }
    
    results[k]->Close();

    //======================================================================================
    //We should here clear the histograms because we want them empty for the next file. 
    //Reset histos
    h_nHits->Reset();
    h_hitEnergy_si->Reset();
    for (Int_t l=0; l<numberofpads; l++){
      h_hitEnergy_si_pad[l]->Reset();
    }

  } // Loop on files

  //======================================================================================
  //Monitor peak evolution as function of the layer/pad.
  //Energy distributions we want to fit 
  TString energy[numberofpads];
  for (Int_t k=0; k<numberofpads; k++){
    energy[k] = ("h_hitEnergy_si_pad_" + IntToString(k+1)).c_str();
  }

  TGraphErrors *grMPV[numberoffiles];
  TString mpvname[numberoffiles];
  
  // TGraphErrors *grSigma[numberoffiles];
  TString sigmaname[numberoffiles];

  //For the mpv and sigma plot 
  //Si layers numbers to position : 0:-- 1: -+ 2:+- 3:++
  const char *pad[numberofpads] = {"Layer 1 - Pad --", "Layer 1 - Pad -+", "Layer 1 - Pad +-","Layer 1 - Pad ++",
			     "Layer 2 - Pad --", "Layer 2 - Pad -+", "Layer 2 - Pad +-","Layer 2 - Pad ++",
			     "Layer 3 - Pad --", "Layer 3 - Pad -+", "Layer 3 - Pad +-","Layer 3 - Pad ++"};
  std::vector<double> padd;
  padd.push_back(1.);padd.push_back(2.);padd.push_back(3.);padd.push_back(4.);padd.push_back(5.);padd.push_back(6.);
  padd.push_back(7.);padd.push_back(8.);padd.push_back(9.);padd.push_back(10.);padd.push_back(11.);padd.push_back(12.);
  std::vector<Double_t>  mpv;
  std::vector<double>  thesigma;
  std::vector<double>  padr;
  std::vector<Double_t>  mpvr;
  std::vector<double>  thesigmaerr;

  //The TH2F to save the histos with mpv and sigma for the different energies
  TH2F * mpvhistos[numberoffiles];
  TH2F * sigmahistos[numberoffiles];
 
  //======================================================================================
  //Canvas
  TCanvas *can[numberoffiles][numberofpads];
  TString canname;

  // TCanvas *canmpv[numberoffiles];
  double mpshift = 0;//-0.22278298;
  double histo_mean = 0.;
  double RMS = 0.;
  //======================================================================================
  //Read the files and histograms here
  TFile* rl[numberoffiles];
  for (int k=0; k<numberoffiles; k++){
    rl[k] = TFile::Open(res[k]);
  }
  TH1F * histo[numberoffiles][numberofpads];
  for (int i=0; i<numberoffiles; i++){
    gROOT->Reset();
    mpvname[i] = "MPV_" + filename[i]; 
    sigmaname[i] = "Sigma" + filename[i];
    mpvhistos[i] = new TH2F(mpvname[i],mpvname[i],numberofpads,0,numberofpads,300,0.7,1.0);
    sigmahistos[i] = new TH2F(sigmaname[i],sigmaname[i],numberofpads,0,numberofpads,300,0.7,1.0);
    //mpvhistos[i]->SetStats(0);
    sigmahistos[i]->SetStats(0);
    for (int k=0; k<numberofpads; k++){
      histo[i][k]= (TH1F*) rl[i]->Get(energy[k]);

      if (!histo[i][k]){std::cerr << "Could not get histogram " << energy[k] << "in file "
  				    <<  results[i]->GetName() << std::endl;}
      std::cout << "=========================================================================" << std::endl;
      std::cout << results[i]->GetName() << std::endl;

      canname = filename[i] + histo[i][k]->GetName();
      can[i][k] = new TCanvas(canname, canname,800,600);
      can[i][k]->cd();

      histo[i][k]->Draw();

      TF1 *fit = 0;
      double mpc = 0;
      histo_mean = histo[i][k]->GetMean();
      RMS = histo[i][k]->GetRMS();
      double minfit = histo_mean - RMS;
      double maxfit = histo_mean + RMS;

      histo[i][k]->Fit("landau","R+","same",minfit,maxfit);
      fit = (TF1*) histo[i][k]->GetFunction("landau");

      mpc = fit->GetParameter(1) - mpshift * fit->GetParameter(2); 

      mpvhistos[i]->GetXaxis()->SetBinLabel( (k+1)  , pad[k]);
      sigmahistos[i]->GetXaxis()->SetBinLabel( (k+1)  , pad[k]);

      mpvhistos[i]->Fill(k , mpc * 1000. );
      // mpvhistos[i]->SetBinError( (k+1) , fit->GetParError(1) ); 

      sigmahistos[i]->Fill( k , fit->GetParameter(2) * 1000.);
      // sigmahistos[i]->SetBinError( (k+1) , fit->GetParError(2) );

      mpv.push_back(mpc * 1000.);
      mpvr.push_back(fit->GetParError(1) * 1000.);
      padr.push_back(0.);

      thesigma.push_back(fit->GetParameter(2) * 1000.);
      thesigmaerr.push_back(fit->GetParError(2) * 1000.);

      
      // if(i==0) {mpvname = "grMPV_50GeV"; sigmaname = "grSigma_50GeV";}
      // if(i==1) {mpvname = "grMPV_100GeV";sigmaname = "grSigma_100GeV";}
      // if(i==2) {mpvname = "grMPV_150GeV";sigmaname = "grSigma_150GeV";}
      // if(i==3) {mpvname = "grMPV_200GeV";sigmaname = "grSigma_200GeV";}
       
      std::cout << "Pad " << k << " MPV "<< mpc << std::endl;

      can[i][k]->Print(canname + ".root",".root");
      can[i][k]->Close();
  
      
    }// end of loop over histos  

     //For the multi plot
    grMPV[i] = new TGraphErrors(mpv.size(),&padd[0],&mpv[0],&padr[0],&mpvr[0]);
    // grSigma[i] = new TGraphErrors(thesigma.size(),&padd[0],&thesigma[0],&thesigmaerr[0],&padr[0]);
    
    mpv.clear();
    mpvr.clear();
    thesigma.clear();
    thesigmaerr.clear();
    padr.clear();

  } // end of loop over files

  TFile f("mpvandsigma.root","recreate");

    // TCanvas *c1 = new TCanvas("c1","mpv",200,10,800,600);
    // mpvhistos[i]->GetXaxis()->SetTitle("Pad");
    // mpvhistos[i]->GetYaxis()->SetTitle("Landau MPV (MIPs) (keV) ");
    // // mpvhistos[i]->GetXaxis()->SetTitleOffset(1.2);
    // // mpvhistos[i]->GetYaxis()->SetTitleOffset(1.3);
    // mpvhistos[i]->GetXaxis()->SetTitleSize(0.045);
    // mpvhistos[i]->GetYaxis()->SetTitleSize(0.045);
    // mpvhistos[i]->GetXaxis()->SetTitleFont(42);
    // mpvhistos[i]->GetYaxis()->SetTitleFont(42);
    // mpvhistos[i]->SetMarkerStyle(20);
    // mpvhistos[i]->SetMarkerSize(1.2);
    // mpvhistos[i]->SetLineColor(2);

    // mpvhistos[i]->Draw();
    // c1->Update();


  //======================================================================================
  //For the multi mpv plot
  TCanvas *c1 = new TCanvas("c1","mpv",200,10,1200,1200); 
  // TGraphErrors *grMPV[0] = new TGraphErrors(mpv.size(),&padd[0],&mpv[0],&mpvr[0],&padr[0]);
  gStyle->SetStripDecimals(false);  
  //Change the bin labels
  grMPV[0]->GetXaxis()->SetBinLabel( grMPV[0]->GetXaxis()->FindBin(1.0) , pad[0]); 
  grMPV[0]->GetXaxis()->SetBinLabel( grMPV[0]->GetXaxis()->FindBin(2.0) , pad[1]); 
  grMPV[0]->GetXaxis()->SetBinLabel( grMPV[0]->GetXaxis()->FindBin(3.0) , pad[2]); 
  grMPV[0]->GetXaxis()->SetBinLabel( grMPV[0]->GetXaxis()->FindBin(4.0) , pad[3]); 
  grMPV[0]->GetXaxis()->SetBinLabel( grMPV[0]->GetXaxis()->FindBin(5.0) , pad[4]); 
  grMPV[0]->GetXaxis()->SetBinLabel( grMPV[0]->GetXaxis()->FindBin(6.0) , pad[5]); 
  grMPV[0]->GetXaxis()->SetBinLabel( grMPV[0]->GetXaxis()->FindBin(7.0) , pad[6]); 
  grMPV[0]->GetXaxis()->SetBinLabel( grMPV[0]->GetXaxis()->FindBin(8.0) , pad[7]); 
  grMPV[0]->GetXaxis()->SetBinLabel( grMPV[0]->GetXaxis()->FindBin(9.0) , pad[8]); 
  grMPV[0]->GetXaxis()->SetBinLabel( grMPV[0]->GetXaxis()->FindBin(10.0) , pad[9]); 
  grMPV[0]->GetXaxis()->SetBinLabel( grMPV[0]->GetXaxis()->FindBin(11.0) , pad[10]); 
  grMPV[0]->GetXaxis()->SetBinLabel( grMPV[0]->GetXaxis()->FindBin(12.0) , pad[11]); 
  grMPV[0]->GetXaxis()->SetTitleOffset(2.0);
  grMPV[0]->GetYaxis()->SetTitleOffset(1.5);
  grMPV[0]->GetXaxis()->SetLabelSize(0.04);
  
  grMPV[0]->SetMarkerStyle(20); 
  grMPV[0]->SetMarkerSize(1.2); 
  grMPV[0]->SetMarkerColor(1);  
  grMPV[0]->SetLineColor(1);
  grMPV[0]->Draw("APL");
  c1->Update();
    
  for (int k=1; k<numberoffiles-2; k++){
    grMPV[k]->SetMarkerStyle(20+k);
    grMPV[k]->SetMarkerSize(1.2); 
    grMPV[k]->SetMarkerColor(1+k);
    grMPV[k]->SetLineColor(1+k);
    grMPV[k]->Draw("PSL");
    c1->Update();

  }

  TLegend *leg = new TLegend(0.8,0.6,0.89,0.89);  //coordinates are fractions of pad dimensions
  leg->AddEntry(grMPV[0],"SiWECAL_B0X0_I0X0_A0","LP");  
  leg->AddEntry(grMPV[1],"SiWECAL_B42X0W_I0X0_A0","LP");  
  leg->AddEntry(grMPV[2],"SiWECAL_B84X0W_I0X0_A0","LP");  
  leg->AddEntry(grMPV[3],"SiWECAL_B84X0W200Fe_I0X0_A0","LP");  
  leg->AddEntry(grMPV[4],"SiWECAL_B84X0W300Fe_I0X0_A0","LP");  
  leg->AddEntry(grMPV[5],"SiWECAL_B0X0_I42X0W_A0","LP");  
  leg->AddEntry(grMPV[6],"SiWECAL_B18X0W_I42X0W_A0","LP");  
  leg->AddEntry(grMPV[7],"SiWECAL_B24X0W_I42X0W_A0","LP");  
  leg->AddEntry(grMPV[8],"SiWECAL_B42X0W_I42X0W_A0","LP");  
  // leg->AddEntry(grMPV[9],"SiWECAL_B24X0W_I42X0W_A48","LP");  
  // leg->AddEntry(grMPV[10],"SiWECAL_B0X0_I0X0_A90","LP");  
  leg->SetFillColor(18);
  leg->SetHeader(" Detector Configurations ");
                                            
  leg->Draw("PS");
  c1->Update();
  
  grMPV[0]-> SetTitle( " Landau Most Probable Value for different detector configurations and #mu^{-} energy 15 GeV  " );
  grMPV[0]->GetHistogram()->SetXTitle(" Pad ");
  grMPV[0]->GetHistogram()->SetYTitle(" Landau MPV (MIPs) (keV) ");
  grMPV[0]->GetYaxis()->SetRangeUser(90.0,120.0);
  //grMPV[0]->GetXaxis()->SetRangeUser(0.,5000.);
  c1->Update();

  //c1->Print("mpv.png",".png");

  c1->Write();
  c1->Close();

  std::cout << "==================================================" << std::endl;
  std::cout << "The Mean MPV value for all configurations for the MeVToMIP value" << std::endl;
  for (int k=0; k<numberoffiles-2; k++){
    std::cout << "Config " << k << ": " << grMPV[k]->GetMean(2)/1000. << std::endl;
  }


  return 0; 

}//main

//===============================================================================
//Definitions here
std::string IntToString (int number1)
{
  std::ostringstream oss1;
  oss1<< number1;
  return oss1.str();
}
//===============================================================================
Double_t langaufun(Double_t *x, Double_t *par) {
  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation),
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.
  // Numeric constants
  Double_t invsq2pi = 0.3989422804014; // (2 pi)^(-1/2)
  Double_t mpshift = 0;//-0.22278298; // Landau maximum location
  // Control constants
  Double_t np = 100.0; // number of convolution steps
  Double_t sc = 5.0; // convolution extends to +-sc Gaussian sigmas
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  // MP shift correction
  mpc = par[1] - mpshift * par[0];
  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  step = (xupp-xlow) / np;
  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }
  return (par[2] * step * sum * invsq2pi / par[3]);
};

//===============================================================================
Double_t langaufungaus(Double_t *x, Double_t *par) {
return langaufun(x,par)+par[4]*TMath::Gaus(x[0],par[5],par[6],true);
};
