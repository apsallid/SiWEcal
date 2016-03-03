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

#include "SiWEcalSSSimHit.hh"

int main(int argc, char** argv){//main


  TString lSuffix = "SiWEcal_15GeVmuon_100kevents"; //SiWEcal_2
  TString plotBase = "PLOTS/SiWEcal_2/mu-";

  TFile *inputFile = TFile::Open("/afs/cern.ch/work/a/apsallid/CMS/Geant4/SiWEcal/"+lSuffix+".root");
  if (!inputFile) {
    std::cout << " -- Error, input file cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  TTree *lTree = (TTree*)inputFile->Get("SiWEcalSSTree");
  if (!lTree){
    std::cout << " -- Error, tree cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  
  TCanvas *myc = new TCanvas("myc","myc",500,500);

  const unsigned nLayers = 12;//12 layers for the 4 Si sensors in each of the 3 layers

  TH2F *p_nHits = new TH2F("nHits","Number of hits vs layer",nLayers,0,nLayers,20,0,20);
  TH1F *p_hitEnergy_si = new TH1F("hitEnergy_si",";E (MeV);SimHits",250,0,1);
  TH1F *p_hitEnergySel_si = new TH1F("hitEnergySel_si",";E (MeV);SimHits",250,0,1);

  std::vector<SiWEcalSSSimHit> * simhitvec = 0;
  
  lTree->SetBranchAddress("SiWEcalSSSimHitVec",&simhitvec);

  const unsigned nEvts = lTree->GetEntries();

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries

    lTree->GetEntry(ievt);
    
    if (ievt%(nLayers*1000)==0) std::cout << "Entry " << ievt << std::endl;

    std::vector<unsigned> nHits;
    nHits.resize(nLayers,0);
    unsigned nHitsall = 0;
    for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
      SiWEcalSSSimHit lHit = (*simhitvec)[iH];
      double energy = lHit.energy();
      double layer = lHit.layer();
      if (energy>0) {
      	p_hitEnergy_si->Fill(energy);
	nHits[layer]++;
      }
    }//loop on hits

    // loop to plot total hits in each layer and count all hits of even
    for (unsigned iL(0); iL<nHits.size(); ++iL){ 
      p_nHits->Fill(iL,nHits[iL]);
      ++nHitsall;
    }
    

    if (nHitsall>0 && nHitsall<2) {
      for (unsigned iH(0); iH<(*simhitvec).size(); ++iH){//loop on hits
	SiWEcalSSSimHit lHit = (*simhitvec)[iH];
	double energy = lHit.energy();
	if (energy>0){
	  p_hitEnergySel_si->Fill(energy);
	}
      }
    }

  }//loop on entries

  
  // myc->cd();
  // gPad->SetLogz(1);
  // gStyle->SetOptStat(1111110);
  // p_nHits->Draw("colz");

  // myc->Update();
  // myc->Print(plotBase+"/mipHits.png");
  // myc->Print(plotBase+"/mipHits.pdf");
  // myc->Print(plotBase+"/mipHits.C");


  myc->cd();
  gPad->SetLogy(1);
  gStyle->SetOptStat(1111110);
  gStyle->SetOptFit(1111);
  p_hitEnergy_si->Draw();
  p_hitEnergy_si->Fit("landau","R+","",0.069,0.18);

  myc->Update();
  myc->Print(plotBase+"/mipDepositAll_si.png");
  myc->Print(plotBase+"/mipDepositAll_si.pdf");
  myc->Print(plotBase+"/mipDepositAll_si.C");

  myc->cd();
  gPad->SetLogy(1);
  gStyle->SetOptStat(1111110);
  gStyle->SetOptFit(1111);
  p_hitEnergySel_si->Draw();
  p_hitEnergySel_si->Fit("landau","LR+","",0.02,1);

  myc->Update();
  myc->Print(plotBase+"/mipDepositSel_si.png");
  myc->Print(plotBase+"/mipDepositSel_si.pdf");
  myc->Print(plotBase+"/mipDepositSel_si.C");
    
  myc->cd();
  gPad->SetLogy(0);
  gStyle->SetOptStat(1111110);
  gStyle->SetOptFit(1111);
  p_nHits->Draw();

  myc->Update();
  myc->Print(plotBase+"/HitsperLayer_si.png");
  myc->Print(plotBase+"/HitsperLayer_si.pdf");
  myc->Print(plotBase+"/HitsperLayer_si.C");


  return 0; 

}//main
