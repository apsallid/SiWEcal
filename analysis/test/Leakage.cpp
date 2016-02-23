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
#include "TChain.h"
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
#include "SiWEcalSSSamplingSection.hh"

//===============================================================================
//Declaration here of any helpful function definition after main
std::string IntToString (int number1);

int main(){//main

  //Now we only work with config 3: SiWECAL_B84X0W_I0X0_A0
  const int numberofconfigs = 1; // 11
  std::vector<int> configs;
  configs.clear();
  configs.push_back(3);//SiWECAL_B84X0W_I0X0_A0
  //======================================================================================
  //The root files with the individual plots for all test beam configs 
  const int numberoffiles = 1;
  //======================================================================================
  //3 layers times 4 Si Pads
  // const int numberofpads = 12; 
  //======================================================================================
  //Particle type and energy
  TString particle = "mu-"; //mu+
  const int numberofenergies = 6; // 15 30 50 80 100 150
  //const int numberofenergies = 1; // 15 30 50 80 100 150
 // int particleenergy = 15;
  std::vector<int> particleenergies;
  particleenergies.clear();
  particleenergies.push_back(15);
  particleenergies.push_back(30);
  particleenergies.push_back(50);
  particleenergies.push_back(80);
  particleenergies.push_back(100);
  particleenergies.push_back(150);

  //======================================================================================
  //Read all files and tree
  TString filename[numberofconfigs][numberofenergies][numberoffiles];
  for (int k=0; k<numberofconfigs; k++){
    for (int l=0; l<numberofenergies; l++){

      // //Here we need to know the exact number of the files 
      // std::string linux_command1 =  "/eos/cms/store/group/phys_b2g/apsallid/SiWEcal/DetectorConfigurations/Config_" + IntToString(configs[k]);
      // std::string linux_command2 =  "/" + particle; 
      // std::string linux_command3 =  "/" + IntToString(particleenergies[l]);
      // std::string linux_command4 = "/results";
      // std::string path =  linux_command1 + linux_command2 + linux_command3 + linux_command4;
      // std::string linux_command5 = " >> pp.txt ";

      // std::string linux_command6 = "eos ls " + path; 
      // std::string linux_command = linux_command6 + linux_command5; 
      // system(linux_command.c_str());
      // std::ifstream pp_file("pp.txt");
      // std::string buffer;
      // getline(pp_file, buffer);
      // std::istringstream is(buffer);
      // is >> numoffilesindataset;
      // pp_file.close();
      // system("rm -rf pp.txt");
 






      for (int j=0; j<numberoffiles; j++){
      TString fp = "SiWEcal_" + IntToString(configs[k]);	
      TString sp = "_" + particle;
      TString tp = "_" + IntToString(particleenergies[l]);
      // TString fop = "GeV_" + IntToString(j);
      TString fop = "GeV";
      filename[k][l][j] = fp + sp + tp + fop;
      //std::cout << filename[k] << std::endl;
      }
    }
  }

  //The tree in the files that we want to get
  // TTree *lTree[numberofconfigs][numberofenergies]; 
  // ...and chains in the cases of multiple files
  TChain *lChain[numberofconfigs][numberofenergies];

  // TFile* files[numberofconfigs][numberofenergies][numberoffiles];
  for (int k=0; k<numberofconfigs; k++){
    for (int l=0; l<numberofenergies; l++){
      lChain[k][l] = new TChain("SiWEcalSSTree"); 
      for (int j=0; j<numberoffiles; j++){
	// files[k]= new TFile(filename[k]);
	//    files[k] = TFile::Open("/afs/cern.ch/work/a/apsallid/CMS/Geant4/SiWEcal/Parallel/testingfile/"+filename[k]+".root");
	// files[k][l][j] = TFile::Open("/afs/cern.ch/work/a/apsallid/CMS/Geant4/SiWEcal/"+filename[k][l][j]+".root");
	// files[k][l] = TFile::Open("/tmp/apsallid/Configs/"+filename[k][l]+".root");
	// lTree[k][l] = (TTree*) files[k][l]->Get("SiWEcalSSTree");
	lChain[k][l]->Add("/tmp/apsallid/"+filename[k][l][j]+".root");
      }
    }
  }
  //Set the branches here
  std::vector<SiWEcalSSSimHit> * simhitvec = 0;
  std::vector<SiWEcalSSSamplingSection> * ssvec = 0;
  for (int k=0; k<numberofconfigs; k++){
    for (int l=0; l<numberofenergies; l++){
      lChain[k][l]->SetBranchAddress("SiWEcalSSSimHitVec",&simhitvec);
      lChain[k][l]->SetBranchAddress("SiWEcalSSSamplingSectionVec",&ssvec);  
    }
  }
  // TCanvas *myc = new TCanvas("myc","myc",500,500);

  //======================================================================================
  // Histos
  //For the loss at the back of the detector
  std::vector<TH1F*> h_Leakagefrombehind;
  h_Leakagefrombehind.clear();
  //For the loss before the calorimeter
  std::vector<TH1F*> h_lossbeforecalo;
  h_lossbeforecalo.clear();
  //Loss in anything but the sensors and the back of the detector
  std::vector<TH1F*> h_lossinpassivelayers;
  h_lossinpassivelayers.clear();
  //Loss in silicon
  std::vector<TH1F*> h_lossinsilicon;
  h_lossinsilicon.clear();

  for (Int_t k=0; k<numberofenergies; k++){
    h_Leakagefrombehind.push_back(new TH1F(("h_Leakagefrombehind_" + IntToString(particleenergies[k])).c_str(),";E (GeV);SimHits",160,0.,160.));
    h_lossbeforecalo.push_back(new TH1F(("h_lossbeforecalo_" + IntToString(particleenergies[k])).c_str(),";E (GeV);SimHits",160,0.,160.));
    h_lossinpassivelayers.push_back(new TH1F(("h_lossinpassivelayers_" + IntToString(particleenergies[k])).c_str(),";E (GeV);SimHits",160,0.,160.));
    h_lossinsilicon.push_back(new TH1F(("h_lossinsilicon_" + IntToString(particleenergies[k])).c_str(),";E (MeV);SimHits",250,0.,250.));
  }

  //======================================================================================
  //The files that we will store the results of this analysis
  TString res[numberofconfigs];
  for (int k=0; k<numberofconfigs; k++){
    TString fp = "SiWEcal_" + IntToString(configs[k]);	
    res[k] = fp + "_final.root";
    //std::cout << res[k] << std::endl;
  }
  
  TFile* results[numberofconfigs];
  for (int k=0; k<numberofconfigs; k++){
    results[k]= new TFile(res[k],"recreate");
    std::cout << "Results file " << res[k] << std::endl;

    //======================================================================================
    //Loop on energies
    for (int l=0; l<numberofenergies; l++){
      //======================================================================================
      // Loop on entries
      for (Int_t ievt=0; ievt<lChain[k][l]->GetEntries(); ievt++) {
	// 	// if (ievt != 0) {continue;}
      	lChain[k][l]->GetEntry(ievt);
	if (ievt%(10000)==0) std::cout << "Entry " << ievt << std::endl;
	// std::cout << "Entry " << ievt << std::endl;
	double losseverywherebuttheendandsensors = 0.;
	double lossinsensors = 0.;
	//Loop on sampling sections
	for(unsigned iL(0); iL<(*ssvec).size(); iL++){
	  // SiWEcalSSSamplingSection lSamSec = (*ssvec)[iL];
	  // double energyabsorber = lSamSec.absorberE();

	  // std::cout << "Sampling Section " << iL << " energyabsorber " << energyabsorber  << std::endl; 
	  
 
	  if ( (*ssvec)[iL].volNb() == 4 ){
	    // std::cout << " Success " <<(*ssvec)[iL].volNb()  <<  " (*ssvec)[iL].absorberE()  " << (*ssvec)[iL].absorberE() << std::endl; 
	    h_Leakagefrombehind[l]->Fill( (*ssvec)[iL].absorberE()/1000. ); //in GeV
	  }
	  if ( (*ssvec)[iL].volNb() == 0 ){
	    h_lossbeforecalo[l]->Fill( (*ssvec)[iL].absorberE()/1000. ); //in GeV
	  }
	  losseverywherebuttheendandsensors = losseverywherebuttheendandsensors + (*ssvec)[iL].passiveE()/1000.;//in GeV
	  //4th layer is dead
	  if (iL!=3){
	    lossinsensors = lossinsensors + (*ssvec)[iL].measuredE(); //in MeV
	  }
	  if (iL==4){
	    h_lossinpassivelayers[l]->Fill( losseverywherebuttheendandsensors  );//in GeV
	    h_lossinsilicon[l]->Fill( lossinsensors );//in MeV
	  }

	} //loop on sampling sections
      
      } //  Loop on entries

      //======================================================================================
      //Write histos 
      h_Leakagefrombehind[l]->Write();
      h_lossbeforecalo[l]->Write();
      h_lossinpassivelayers[l]->Write();
      h_lossinsilicon[l]->Write();

      //======================================================================================
      //We should here clear the histograms because we want them empty for the next file. 
      //Reset histos
      h_Leakagefrombehind[l]->Reset();
      h_lossbeforecalo[l]->Reset();
      h_lossinpassivelayers[l]->Reset();
      h_lossinsilicon[l]->Reset();

    }//Loop on energies
   
    results[k]->Close();

  } // Loop on files

  //======================================================================================
  //Make one plot for all different leakage energies
  TCanvas* c1 = new TCanvas("c1", "  ");
  TCanvas* c2 = new TCanvas("c2", "  ");
  TCanvas* c3 = new TCanvas("c3", "  ");
  TCanvas* c4 = new TCanvas("c4", "  ");
  //For the legend
  // TLegend* leg[numberoffiles];
  TLegend* leg1 = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg1->SetHeader("Energy");
  leg1->SetFillColor(17);
  TLegend* leg2 = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg2->SetHeader("Energy");
  leg2->SetFillColor(17);
  TLegend* leg3 = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg3->SetHeader("Energy");
  leg3->SetFillColor(17);
  TLegend* leg4 = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg4->SetHeader("Energy");
  leg4->SetFillColor(17);

  TString procNa[numberofenergies];
  procNa[0] = "15 GeV";
  procNa[1] = "30 GeV";
  procNa[2] = "50 GeV";
  procNa[3] = "80 GeV";
  procNa[4] = "100 GeV";
  procNa[5] = "150 GeV";

  //======================================================================================
  //The files that we will store the results of this analysis for the combined plot
  TString res_com[numberofconfigs];
  TFile* results_com[numberofconfigs];
  for (int k=0; k<numberofconfigs; k++){
    TString fp = "SiWEcal_" + IntToString(configs[k]);	
    res_com[k] = fp + "_combinedplots.root";
    //std::cout << res[k] << std::endl;
  }
  TH1F* hist1,* hist2,* hist3,* hist4;

  TString titleofplot1,titleofplot2,titleofplot; 
  for (int k=0; k<numberofconfigs; k++){
    results[k]= TFile::Open(res[k],"read");
    std::cout << "Results file " << res[k] << std::endl;
    
    //======================================================================================
    //Loop on energies
    for (int l=0; l<numberofenergies; l++){
      
      //For the leakage
      //------------------------------------------------------------------------------------------------
      hist1 = (TH1F*)results[k]->Get(("h_Leakagefrombehind_" + IntToString(particleenergies[l])).c_str());
      titleofplot1 = "Leakage energy behind calorimeter for "; 
      titleofplot2 =  particle +  " beam particle gun"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist1->SetTitle(titleofplot); 
      hist1->GetXaxis()->SetTitle("Leakage energy (GeV)"); 
      hist1->GetYaxis()->SetTitle("Events/1 GeV");
      c1->SetLogy();
      if ( l == 0){hist1->SetLineColor(4);}
      if ( l == 1){hist1->SetLineColor(2);}
      if ( l == 2){hist1->SetLineColor(1);}
      if ( l == 3){hist1->SetLineColor(3);}
      if ( l == 4){hist1->SetLineColor(5);}
      if ( l == 5){hist1->SetLineColor(6);}
      //hist->SetLineColor(l+1);
      leg1->AddEntry(hist1, procNa[l], "L");
      //hist->GetYaxis()->SetRangeUser(0.,900.);//  21: 100. 22: 900.
      c1->cd();
      c1->Update(); 
      l == 0 ? hist1->Draw("HIST") : hist1->Draw("HISTsame");
      l == 0 ? leg1->Draw() : leg1->Draw("same");
      c1->Update(); 

      //For the loss in sensors
      //------------------------------------------------------------------------------------------------
      c2->cd();
      hist2 = (TH1F*)results[k]->Get(("h_lossinsilicon_" + IntToString(particleenergies[l])).c_str());
      titleofplot1 = "Energy loss in silicon sensors for "; 
      titleofplot2 =  particle +  " beam particle gun"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist2->SetTitle(titleofplot); 
      hist2->GetXaxis()->SetTitle("Energy loss in silicon sensors (MeV)"); 
      hist2->GetYaxis()->SetTitle("Events/1 MeV");
      c2->SetLogy();
      if ( l == 0){hist2->SetLineColor(4);}
      if ( l == 1){hist2->SetLineColor(2);}
      if ( l == 2){hist2->SetLineColor(1);}
      if ( l == 3){hist2->SetLineColor(3);}
      if ( l == 4){hist2->SetLineColor(5);}
      if ( l == 5){hist2->SetLineColor(6);}
      //hist->SetLineColor(l+1);
      leg2->AddEntry(hist2, procNa[l], "L");
      //hist->GetYaxis()->SetRangeUser(0.,900.);//  21: 100. 22: 900.
      c2->Update(); 
      l == 0 ? hist2->Draw("HIST") : hist2->Draw("HISTsame");
      l == 0 ? leg2->Draw() : leg2->Draw("same");
      c2->Update(); 

      //For the loss before calo
      //------------------------------------------------------------------------------------------------
      c3->cd();
      hist3 = (TH1F*)results[k]->Get(("h_lossbeforecalo_" + IntToString(particleenergies[l])).c_str());
      titleofplot1 = "Energy loss before calorimeter for "; 
      titleofplot2 =  particle +  " beam particle gun"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist3->SetTitle(titleofplot); 
      hist3->GetXaxis()->SetTitle("Energy loss before calorimeter (GeV)"); 
      hist3->GetYaxis()->SetTitle("Events/1 GeV");
      c3->SetLogy();
      if ( l == 0){hist3->SetLineColor(4);}
      if ( l == 1){hist3->SetLineColor(2);}
      if ( l == 2){hist3->SetLineColor(1);}
      if ( l == 3){hist3->SetLineColor(3);}
      if ( l == 4){hist3->SetLineColor(5);}
      if ( l == 5){hist3->SetLineColor(6);}
      //hist->SetLineColor(l+1);
      leg3->AddEntry(hist3, procNa[l], "L");
      //hist->GetYaxis()->SetRangeUser(0.,900.);//  21: 100. 22: 900.
      c3->Update(); 
      l == 0 ? hist3->Draw("HIST") : hist3->Draw("HISTsame");
      l == 0 ? leg3->Draw() : leg3->Draw("same");
      c3->Update(); 
      //For the loss in passive layers
      //------------------------------------------------------------------------------------------------
      c4->cd();
      hist4 = (TH1F*)results[k]->Get(("h_lossinpassivelayers_" + IntToString(particleenergies[l])).c_str());
      titleofplot1 = "Energy loss in passive layers for "; 
      titleofplot2 =  particle +  " beam particle gun"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist4->SetTitle(titleofplot); 
      hist4->GetXaxis()->SetTitle("Energy loss in passive layers (GeV)"); 
      hist4->GetYaxis()->SetTitle("Events/1 GeV");
      c4->SetLogy();
      if ( l == 0){hist4->SetLineColor(4);}
      if ( l == 1){hist4->SetLineColor(2);}
      if ( l == 2){hist4->SetLineColor(1);}
      if ( l == 3){hist4->SetLineColor(3);}
      if ( l == 4){hist4->SetLineColor(5);}
      if ( l == 5){hist4->SetLineColor(6);}
      //hist->SetLineColor(l+1);
      leg4->AddEntry(hist4, procNa[l], "L");
      //hist->GetYaxis()->SetRangeUser(0.,900.);//  21: 100. 22: 900.
      c4->Update(); 
      l == 0 ? hist4->Draw("HIST") : hist4->Draw("HISTsame");
      l == 0 ? leg4->Draw() : leg4->Draw("same");
      c4->Update(); 

      
    }//Loop on energies
    
    results_com[k]= new TFile(res_com[k],"recreate");
    
    //======================================================================================
    //Write canvas with combined plot 
    // c1->Print("Leakage.pdf",".pdf");
    c1->Write();
    c2->Write();
    c3->Write();
    c4->Write();
    //Reset histos
    // hist->Reset();

    results_com[k]->Close();
  
  }//Loop on configs




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
