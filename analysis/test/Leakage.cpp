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
  TString particle = "pi+"; //remaining e+, pi+
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
	// lChain[k][l]->Add("/tmp/apsallid/"+filename[k][l][j]+".root");
	lChain[k][l]->Add("root://eoscms//eos/cms/store/group/phys_b2g/apsallid/SiWEcal/DetectorConfigurations/Config_3/"+filename[k][l][j]+".root");
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
  //Total Loss
  std::vector<TH1F*> h_totalloss;
  h_totalloss.clear();
  //======================================================================================
  // For the fit
  std::vector<TGraphErrors*> ElossUpstreamvsTotalMeasuredMIPs;
  ElossUpstreamvsTotalMeasuredMIPs.clear();

  TGraphErrors *corr = new TGraphErrors();
  for (Int_t k=0; k<numberofenergies; k++){
    h_Leakagefrombehind.push_back(new TH1F(("h_Leakagefrombehind_" + IntToString(particleenergies[k])).c_str(),";E (GeV);SimHits",160,0.,160.));
    h_lossbeforecalo.push_back(new TH1F(("h_lossbeforecalo_" + IntToString(particleenergies[k])).c_str(),";E (MeV);SimHits",100,0.,100.));
    h_lossinpassivelayers.push_back(new TH1F(("h_lossinpassivelayers_" + IntToString(particleenergies[k])).c_str(),";E (GeV);SimHits",160,0.,160.));
    h_lossinsilicon.push_back(new TH1F(("h_lossinsilicon_" + IntToString(particleenergies[k])).c_str(),";E (MeV);SimHits",500,0.,500.));
    h_totalloss.push_back(new TH1F(("h_totalloss_" + IntToString(particleenergies[k])).c_str(),";E (GeV);SimHits",200,0.,200.));

    ElossUpstreamvsTotalMeasuredMIPs.push_back( (TGraphErrors *) corr->Clone( ("ElossUpstreamvsTotalMeasuredMIPs_" + IntToString(particleenergies[k])).c_str() )  );
    
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
	// if (ievt==100){break;}
	// std::cout << "Entry " << ievt << std::endl;
	double losseverywherebuttheendandsensors = 0.;
	double lossinsensors = 0.;
	double totalloss = 0.;
	double lossbeforecalo = 0.;
	//Loop on sampling sections
	for(unsigned iL(0); iL<(*ssvec).size(); iL++){
	  // SiWEcalSSSamplingSection lSamSec = (*ssvec)[iL];
	  // double energyabsorber = lSamSec.absorberE();

	  // std::cout << "Sampling Section " << iL << " energyabsorber " << energyabsorber  << std::endl; 
	  
	  //In the 4th dead layer everything is leakage.
	  //Since in our numbering scheme 0 is before calo, 4 is the dead layer
	  if ( iL == 4 ){
	    h_Leakagefrombehind[l]->Fill( (*ssvec)[iL].totalE()/1000. ); //in GeV
	  }
	  //In the last layer is the leakage also
	  if ( iL == 5 ){
	    // std::cout << " Success " <<(*ssvec)[iL].volNb()  <<  " (*ssvec)[iL].absorberE()  " << (*ssvec)[iL].absorberE() << std::endl; 
	    h_Leakagefrombehind[l]->Fill( (*ssvec)[iL].absorberE()/1000. ); //in GeV
	  }
	  //The loss before calo, that is the 0 sampling section.
	  if ( iL == 0 ){
	    lossbeforecalo = (*ssvec)[iL].absorberE();//in MeV
	    h_lossbeforecalo[l]->Fill( (*ssvec)[iL].absorberE() ); //in MeV
	  }
	  losseverywherebuttheendandsensors = losseverywherebuttheendandsensors + (*ssvec)[iL].passiveE()/1000.;//in GeV
	  //4th layer is dead
	  if (iL!=4){
	    lossinsensors = lossinsensors + (*ssvec)[iL].measuredE(); //in MeV
	  }
	  //Total loss
	  totalloss = totalloss + ((*ssvec)[iL].measuredE()/1000.) + ((*ssvec)[iL].absorberE()/1000.); 
	  if (iL==5){
	    h_lossinpassivelayers[l]->Fill( losseverywherebuttheendandsensors  );//in GeV
	    h_lossinsilicon[l]->Fill( lossinsensors );//in MeV
	    h_totalloss[l]->Fill( totalloss  );//in GeV
	    //For the correction in Ereco
	    Int_t np=ElossUpstreamvsTotalMeasuredMIPs[l]->GetN();
	    // std::cout << "np " << np << std::endl;
	    ElossUpstreamvsTotalMeasuredMIPs[l]->SetPoint(np, lossinsensors/(0.095)  , lossbeforecalo  ); //y: ElossUp in MeV, x: TotalMeasuredMIPs in MIPs
	  }


	} //loop on sampling sections
      
      } //  Loop on entries

      //======================================================================================
      //Write histos 
      h_Leakagefrombehind[l]->Write();
      h_lossbeforecalo[l]->Write();
      h_lossinpassivelayers[l]->Write();
      h_lossinsilicon[l]->Write();
      h_totalloss[l]->Write();
      ElossUpstreamvsTotalMeasuredMIPs[l]->Write();

      //======================================================================================
      //We should here clear the histograms because we want them empty for the next file. 
      //Reset histos
      h_Leakagefrombehind[l]->Reset();
      h_lossbeforecalo[l]->Reset();
      h_lossinpassivelayers[l]->Reset();
      h_lossinsilicon[l]->Reset();
      h_totalloss[l]->Reset();
      
      
    }//Loop on energies
   
    results[k]->Close();

  } // Loop on files

  //======================================================================================
  //Make one plot for all different leakage energies
  TCanvas* c1 = new TCanvas("c1", "  ");
  TCanvas* c2 = new TCanvas("c2", "  ");
  TCanvas* c3 = new TCanvas("c3", "  ");
  TCanvas* c4 = new TCanvas("c4", "  ");
  TCanvas* c5 = new TCanvas("c5", "  ");
  TCanvas *c6[numberofenergies];
  for (int l=0; l<numberofenergies; l++){
    c6[l] = new TCanvas( ("c6_"+IntToString(l)) .c_str(), "  ");
  }
  TCanvas* c7 = new TCanvas("c7", "  ");
  TCanvas* c8 = new TCanvas("c8", "  ");

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
  TLegend* leg5 = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg5->SetHeader("Energy");
  leg5->SetFillColor(17);

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
    res_com[k] = fp + "_combinedplots_closure.root";
    //std::cout << res[k] << std::endl;
  }
  TH1F* hist1,* hist2,* hist3,* hist4,* hist5;
  TGraphErrors *gr[numberofenergies];
  //For accessing the fit results
  TF1 *fit[numberofenergies];
  TGraphErrors *gr_dEUpstream_totalmeasuredMIPs_fit;
  TGraphErrors *gr_dEUpstream_totalmeasuredMIPs_fit_const;
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TString titleofplot1,titleofplot2,titleofplot3,titleofplot; 
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
      c1->SaveAs("Leakage_" + particle + ".png"); 
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
      c2->SaveAs("EnergyLossInSilicon_" + particle + ".png"); 
      c2->Update(); 

      //For the loss before calo
      //------------------------------------------------------------------------------------------------
      c3->cd();
      hist3 = (TH1F*)results[k]->Get(("h_lossbeforecalo_" + IntToString(particleenergies[l])).c_str());
      titleofplot1 = "Energy loss before calorimeter for "; 
      titleofplot2 =  particle +  " beam particle gun"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist3->SetTitle(titleofplot); 
      hist3->GetXaxis()->SetTitle("Energy loss before calorimeter (MeV)"); 
      hist3->GetYaxis()->SetTitle("Events/1 MeV");
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
      c3->SaveAs("EnergyLossBeforeCalo_" + particle + ".png"); 
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
      c4->SaveAs("EnergyLossInPassiveLayers_" + particle + ".png"); 
      c4->Update(); 
      //For the total loss
      //------------------------------------------------------------------------------------------------
      c5->cd();
      hist5 = (TH1F*)results[k]->Get(("h_totalloss_" + IntToString(particleenergies[l])).c_str());
      titleofplot1 = "Total Energy loss for "; 
      titleofplot2 =  particle +  " beam particle gun"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist5->SetTitle(titleofplot); 
      hist5->GetXaxis()->SetTitle("Total Energy loss (GeV)"); 
      hist5->GetYaxis()->SetTitle("Events/1 GeV");
      c5->SetLogy();
      if ( l == 0){hist5->SetLineColor(4);}
      if ( l == 1){hist5->SetLineColor(2);}
      if ( l == 2){hist5->SetLineColor(1);}
      if ( l == 3){hist5->SetLineColor(3);}
      if ( l == 4){hist5->SetLineColor(5);}
      if ( l == 5){hist5->SetLineColor(6);}
      //hist->SetLineColor(l+1);
      leg5->AddEntry(hist5, procNa[l], "L");
      //hist->GetYaxis()->SetRangeUser(0.,900.);//  21: 100. 22: 900.
      c5->Update(); 
      l == 0 ? hist5->Draw("HIST") : hist5->Draw("HISTsame");
      l == 0 ? leg5->Draw() : leg5->Draw("same");
      c5->SaveAs("TotalEnergyLoss_" + particle + ".png"); 
      c5->Update(); 

      //For the fit 
      //------------------------------------------------------------------------------------------------
      gr[l] = (TGraphErrors*)results[k]->Get( ("ElossUpstreamvsTotalMeasuredMIPs_" + IntToString(particleenergies[l])).c_str()  );
      gr[l]->GetXaxis()->SetRangeUser(0.,5000.);
      gr[l]->Fit("pol1");
	
      // prof[(int) iL] = hist[(int) iL]->ProfileX();
      // prof[(int) iL]->Fit("pol1");
      fit[l] = gr[l]->GetFunction("pol1");
      //Results of the fit
      std::cout << "First parameter: " << fit[l]->GetParameter(0) << " Second parameter: " << fit[l]->GetParameter(1) << std::endl;
	
      titleofplot1 = "Energy lost in upstream materials before calorimeter vs measured energy in all layer for "; 
      titleofplot2 =  particle +  " "; 
      titleofplot3 =  IntToString(particleenergies[l]) +  " GeV beam particle gun"; 
      titleofplot = titleofplot1 + titleofplot2 + titleofplot3;
      gr[l]->SetTitle(titleofplot); 
      gr[l]->GetHistogram()->SetXTitle("Measured energy in all layers (MIPs)"); 
      gr[l]->GetHistogram()->SetYTitle("dE_lossUp (MeV)");
      if ( l == 0){fit[l]->SetLineColor(4);gr[l]->SetMarkerColor(4);gr[l]->SetMarkerStyle(20);}
      if ( l == 1){fit[l]->SetLineColor(2);gr[l]->SetMarkerColor(2);gr[l]->SetMarkerStyle(21);}
      if ( l == 2){fit[l]->SetLineColor(1);gr[l]->SetMarkerColor(1);gr[l]->SetMarkerStyle(22);}
      if ( l == 3){fit[l]->SetLineColor(3);gr[l]->SetMarkerColor(3);gr[l]->SetMarkerStyle(23);}
      if ( l == 5){fit[l]->SetLineColor(6);gr[l]->SetMarkerColor(6);gr[l]->SetMarkerStyle(29);}
      if ( l == 4){fit[l]->SetLineColor(5);gr[l]->SetMarkerColor(5);gr[l]->SetMarkerStyle(28);}
      //hist[(int) iL]->SetLineColor(l+1);
      // leg_devsmipswithfit[(int) iL]->AddEntry(hist[(int) iL], procNa[l], "P");
      // leg_devsmipswithfit[(int) iL]->AddEntry(gr[(int) iL], procNa[l], "P");
      //hist[(int) iL]->GetYaxis()->SetRangeUser(0.,900.);//  21: 100. 22: 900.
      c6[l]->cd();
      c6[l]->Update(); 
      gr[l]->Draw("AP");
      // l == 0 ? gr[(int) iL]->Draw("AP") : gr[(int) iL]->Draw("APS");
      // leg_devsmipswithfit[(int) iL]->Draw("PS");
      // l == 0 ? leg_devsmipswithfit[(int) iL]->Draw() : leg_devsmipswithfit[(int) iL]->Draw("same");
      // l == 0 ? fit[(int) iL]->Draw() : fit[(int) iL]->Draw("same");
      gr[l]->GetHistogram()->GetXaxis()->SetRangeUser(0.,5000.);

      char buf[500];
      TLatex lat;
      double latx = gr[l]->GetHistogram()->GetXaxis()->GetXmin()+(gr[l]->GetHistogram()->GetXaxis()->GetXmax()-gr[l]->GetHistogram()->GetXaxis()->GetXmin())/20.;
      double laty = gr[l]->GetHistogram()->GetMaximum();
      sprintf(buf,"dE_upstream = p0 + p1 * TotMeasMIPs");
      lat.DrawLatex(latx,laty*0.8,buf);
      sprintf(buf,"p0 = %3.3f +/- %3.3f",fit[l]->GetParameter(0),fit[l]->GetParError(0));
      lat.DrawLatex(latx,laty*0.7,buf);
      sprintf(buf,"p1 = %3.3f +/- %3.3f",fit[l]->GetParameter(1),fit[l]->GetParError(1));
      lat.DrawLatex(latx,laty*0.6,buf);
      sprintf(buf,"#chi^{2}/N = %3.3f/%d = %3.3f",fit[l]->GetChisquare(),fit[l]->GetNDF(),fit[l]->GetChisquare()/fit[l]->GetNDF());
      lat.DrawLatex(latx,laty*0.5,buf);
 

      c6[l]->Update(); 
      TString cantopng1 = "dEUpstreamVStotalmeasuredMIPs_plusfit_"+ particle;
      TString cantopng2 = "_" + IntToString(particleenergies[l]);
      TString cantopng3 = "GeV";
      TString cantopng = cantopng1 + cantopng2 + cantopng3 + ".png";
	
      c6[l]->SaveAs(cantopng);

      
    }//Loop on energies
    
    results_com[k]= new TFile(res_com[k],"recreate");

    //======================================================================================
    //Write canvas with combined plot 
    // c1->Print("Leakage.pdf",".pdf");
    c1->Write();
    c2->Write();
    c3->Write();
    c4->Write();
    c5->Write();
    for (int l=0; l<numberofenergies; l++){
      c6[l]->Write();
    }


    //Reset histos
    // hist->Reset();

    //Here we will make the plot for the fit results
    c7->cd();
    gr_dEUpstream_totalmeasuredMIPs_fit = (TGraphErrors *) corr->Clone( "ElossUpstreamvsTotalMeasuredMIPs_fit");
    gr_dEUpstream_totalmeasuredMIPs_fit->SetMarkerStyle( 20 ); 
    // gr_dEUpstream_totalmeasuredMIPs_fit->SetMarkerSize(0.2); 
    gr_dEUpstream_totalmeasuredMIPs_fit->SetMarkerColor(1);  
    gr_dEUpstream_totalmeasuredMIPs_fit->SetLineColor(1);
    c7->Update(); 

    // Int_t np=gr_dEUpstream_totalmeasuredMIPs_fit->GetN();
    // gr_dEUpstream_totalmeasuredMIPs_fit->SetPoint(numberofenergies, particleenergies[l] , fit[l]->GetParameter(0) );
    // gr_dEUpstream_totalmeasuredMIPs_fit->SetPointError(numberofenergies, 0. , fit[l]->GetParError(0) );
    for (int l=0; l<numberofenergies; l++){
      gr_dEUpstream_totalmeasuredMIPs_fit->SetPoint( l , particleenergies[l] , fit[l]->GetParameter(1) );
      gr_dEUpstream_totalmeasuredMIPs_fit->SetPointError(l , 0. , fit[l]->GetParError(1) );
    }
    gr_dEUpstream_totalmeasuredMIPs_fit->Draw("APL");
    c7->Update(); 

    gr_dEUpstream_totalmeasuredMIPs_fit->SetTitle( "Energy lost Upstream before calorimeter vs. Total measured MIPs fit results for " + particle + " for all beam energies" );
    gr_dEUpstream_totalmeasuredMIPs_fit->GetHistogram()->SetXTitle("Beam Energy (GeV)"); 
    gr_dEUpstream_totalmeasuredMIPs_fit->GetHistogram()->SetYTitle("dE_lossUpstream/Total Measured MIPs slope");
    gr_dEUpstream_totalmeasuredMIPs_fit->GetXaxis()->SetRangeUser(0.,160.);  
    gr_dEUpstream_totalmeasuredMIPs_fit->GetYaxis()->SetRangeUser(-5.,5.);  

    c7->Update(); 
    TString cnpng = "dEUpstreamVStotalmeasuredMIPs_plusfit_slope"+ particle + ".png";
    c7->SaveAs(cnpng);

    c7->Write();

   //Here we will make the plot for the fit results constant term
    c8->cd();
    gr_dEUpstream_totalmeasuredMIPs_fit_const = (TGraphErrors *) corr->Clone( "ElossUpstreamvsTotalMeasuredMIPs_fit_const");
    gr_dEUpstream_totalmeasuredMIPs_fit_const->SetMarkerStyle( 20 ); 
    // gr_dEUpstream_totalmeasuredMIPs_fit_const->SetMarkerSize(0.2); 
    gr_dEUpstream_totalmeasuredMIPs_fit_const->SetMarkerColor(1);  
    gr_dEUpstream_totalmeasuredMIPs_fit_const->SetLineColor(1);
    c8->Update(); 

    // Int_t np=gr_dEUpstream_totalmeasuredMIPs_fit_const->GetN();
    // gr_dEUpstream_totalmeasuredMIPs_fit_const->SetPoint(numberofenergies, particleenergies[l] , fit[l]->GetParameter(0) );
    // gr_dEUpstream_totalmeasuredMIPs_fit_const->SetPointError(numberofenergies, 0. , fit[l]->GetParError(0) );
    for (int l=0; l<numberofenergies; l++){
      gr_dEUpstream_totalmeasuredMIPs_fit_const->SetPoint( l , particleenergies[l] , fit[l]->GetParameter(0) );
      gr_dEUpstream_totalmeasuredMIPs_fit_const->SetPointError(l , 0. , fit[l]->GetParError(0) );
    }
    gr_dEUpstream_totalmeasuredMIPs_fit_const->Draw("APL");
    c8->Update(); 

    gr_dEUpstream_totalmeasuredMIPs_fit_const->SetTitle( "Energy lost Upstream before calorimeter vs. Total measured MIPs fit results for " + particle + " for all beam energies" );
    gr_dEUpstream_totalmeasuredMIPs_fit_const->GetHistogram()->SetXTitle("Beam Energy (GeV)"); 
    gr_dEUpstream_totalmeasuredMIPs_fit_const->GetHistogram()->SetYTitle("dE_lossUpstream/Total Measured MIPs Constant Term");
    gr_dEUpstream_totalmeasuredMIPs_fit_const->GetXaxis()->SetRangeUser(0.,160.);  
    gr_dEUpstream_totalmeasuredMIPs_fit_const->GetYaxis()->SetRangeUser(-15.,20.);  

    c8->Update(); 
    TString cnpng1 = "dEUpstreamVStotalmeasuredMIPs_plusfit_const"+ particle + ".png";
    c8->SaveAs(cnpng1);
    c8->Write();

    results_com[k]->Close();

    std::ofstream myfile;
    TString titleoffile = "dE_lossUpvsMIPs_" + particle + ".txt";
    myfile.open(titleoffile);

    //Here print the w0 and w1 that are going to be used for Ereco_corrected
    for (int l=0; l<numberofenergies; l++){
      std::cout << "---------------------------------------------------------" << std::endl;
      std::cout << "Particle " <<  particle << std::endl;
      std::cout << "Energy of the beam " <<  particleenergies[l] << std::endl;
      std::cout << "LossUpstreamvsmeasuredMIPs slope: " <<  fit[l]->GetParameter(1) << " error slope " << fit[l]->GetParError(1) << std::endl;
      std::cout << "LossUpstreamvsmeasuredMIPs constant: " <<  fit[l]->GetParameter(0) << " error constant " << fit[l]->GetParError(0) << std::endl;
      
      //particle particleenergy constant slope errorconstant errorslope
      myfile << particle << " " << particleenergies[l] << " " << fit[l]->GetParameter(0) << " " << fit[l]->GetParameter(1) << " " << fit[l]->GetParError(0) <<  " " << fit[l]->GetParError(1) << "\n";

    }

    myfile.close();

  

    











 
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
