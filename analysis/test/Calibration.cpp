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
#include "TProfile.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
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
  TString particle = "pi+"; //mu+
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
	// lChain[k][l]->Add("/tmp/apsallid/"+filename[k][l][j]+".root");

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
  //dE vs MIPs
  std::map< int, std::vector<TH2F*> > h_dEvsMIPs;
  h_dEvsMIPs.clear();
  //SF for layers defined as: SF = (dE + Esilicon)/dE = (Epassive+Eactive)/Eactive
  std::map< int, std::vector<TH1F*> > h_SF;
  h_SF.clear();
  //======================================================================================
  // For the fit
  std::map< int, std::vector<TGraphErrors*> > dEvsMIPs_calib;
  dEvsMIPs_calib.clear();
 
  lChain[0][0]->GetEntry(0);
  TGraphErrors *calib = new TGraphErrors();
  for(unsigned iL(0); iL<((*ssvec).size()-2); iL++){
    for (Int_t k=0; k<numberofenergies; k++){
      std::string auxnam1 = "h_dEvsMIPs_" + IntToString((int) iL);
      std::string auxnam1_sf = "h_SF_" + IntToString((int) iL);
      std::string auxnam1_tgr = "dEvsMIPs_" + IntToString((int) iL);
      std::string auxnam2 = "_" + IntToString(particleenergies[k]);
      std::string auxnam = auxnam1 + auxnam2;
      std::string auxnam_sf = auxnam1_sf + auxnam2;
      std::string auxnam_tgr = auxnam1_tgr + auxnam2;
      h_dEvsMIPs[(int) iL].push_back(new TH2F( auxnam.c_str(),";E (MIPs);SimHits;E (GeV);SimHits",2000,0.,2000.,160,0.,160.));
      h_SF[(int) iL].push_back(new TH1F( auxnam_sf.c_str(),";Scale factors;Events",200,0.9,1.1));

      dEvsMIPs_calib[(int) iL].push_back( (TGraphErrors *) calib->Clone( auxnam_tgr.c_str() )  );
      
    }
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
	if (ievt%(1000)==0) std::cout << "Entry " << ievt << std::endl;
	if (ievt==100){break;}
	// std::cout << "Entry " << ievt << std::endl;
	// double lossbefcalo = 0.;
	// double scalefactor = 0.;
	//Loop on sampling sections
	for(unsigned iL(0); iL<(*ssvec).size(); iL++){
	  // SiWEcalSSSamplingSection lSamSec = (*ssvec)[iL];
	  // double energyabsorber = lSamSec.absorberE();

	  // std::cout << "Sampling Section " << iL << " energyabsorber " << energyabsorber  << std::endl; 
	  
 	  //dE vs MIPs 
	  //For the energy lost start with material just after the previous Silicon layer (e.g. layer1) and 
	  //stop just before the specific layer (layer 2). Be careful with layer 0 before calorimeter. 
	  //For the conversion we use the simulation value that we calculated: 0.095 MeV/MIP
	  //In the 4th dead layer everything is leakage.
	  //Since in our numbering scheme 0 is before calo, 4 is the dead layer
	  //In the last layer is remaining leakage. 
	  //Also, before calorimeter (wall, air gap) we should not add the energy lost.
	  if ( (iL!=0) && (iL!=4) && (iL!=5)  ){
	    h_dEvsMIPs[(int) (*ssvec)[iL].volNb()][l]->Fill( ( (*ssvec)[iL].measuredE()/(0.095) ) , (*ssvec)[iL].absorberE()/1000. ); //x is MIPs, y is dE in GeV. 
	    Int_t np=dEvsMIPs_calib[(int) iL][l]->GetN();
	    // std::cout << "np " << np << std::endl;
	    dEvsMIPs_calib[(int) iL][l]->SetPoint(np, ( (*ssvec)[iL].measuredE()/(0.095) )  , (*ssvec)[iL].absorberE()/1000.  );
	    // dEvsMIPs_calib[(int) iL][l]->SetPointError();
	    
	    //SF = (dE + Esilicon)/dE = (Epassive+Eactive)/Eactive
	    // scalefactor = ( ( (lossbefcalo * 1000.) + (*ssvec)[iL].measuredE() ) / (lossbefcalo * 1000.) ); //both in MeV here
	    // h_SF[(int) (*ssvec)[iL].volNb()][l]->Fill( scalefactor ); 
	  }	  


	} //loop on sampling sections
      
      } //  Loop on entries

      //======================================================================================
      //Write histos 
      lChain[0][0]->GetEntry(0);
      for(unsigned iL(0); iL<((*ssvec).size()-2); iL++){
	h_dEvsMIPs[(int) iL][l]->Write();
	h_SF[(int) iL][l]->Write();
	dEvsMIPs_calib[(int) iL][l]->Write();
      }


      //======================================================================================
      //We should here clear the histograms because we want them empty for the next file. 
      //Reset histos
      for(unsigned iL(0); iL<((*ssvec).size()-2); iL++){
	h_dEvsMIPs[(int) iL][l]->Reset();
	h_SF[(int) iL][l]->Reset();
	// dEvsMIPs_calib[(int) iL][l]->Reset();
      }
      
    }//Loop on energies
   
    results[k]->Close();

  } // Loop on files




  //======================================================================================
  //Make one plot for all different leakage energies
  TCanvas *c1[((*ssvec).size()-2)];
  TCanvas *c2[((*ssvec).size()-2)];
  TCanvas *c4[((*ssvec).size()-2)][numberofenergies];
  for(unsigned iL(0); iL<((*ssvec).size()-2); iL++){
    c1[(int) iL] = new TCanvas(("c1_"+ IntToString((int) iL)).c_str(), "  ");
    c2[(int) iL] = new TCanvas(("c2_"+ IntToString((int) iL)).c_str(), "  ");
    for (int l=0; l<numberofenergies; l++){
      std::string can_title = "c4_"+ IntToString((int) iL);
      c4[(int) iL][l] = new TCanvas( (can_title+IntToString(l)) .c_str(), "  ");
    }
  }
  TCanvas *c3 = new TCanvas("c3", "  ");
  //For the legend
  // TLegend* leg[numberoffiles];
  TLegend* leg[((*ssvec).size()-2)];
  for(unsigned iL(0); iL<((*ssvec).size()-2); iL++){
    leg[(int) iL] = new TLegend(0.5, 0.7, 0.8, 0.9);
    leg[(int) iL]->SetHeader("Energy");
    leg[(int) iL]->SetFillColor(17);
  }
  TLegend* leg_sf[((*ssvec).size()-2)];
  for(unsigned iL(0); iL<((*ssvec).size()-2); iL++){
    leg_sf[(int) iL] = new TLegend(0.5, 0.7, 0.8, 0.9);
    leg_sf[(int) iL]->SetHeader("Energy");
    leg_sf[(int) iL]->SetFillColor(17);
  }
  TLegend* leg_devsmipswithfit[((*ssvec).size()-2)];
  for(unsigned iL(0); iL<((*ssvec).size()-2); iL++){
    leg_devsmipswithfit[(int) iL] = new TLegend(0.5, 0.7, 0.8, 0.9);
    leg_devsmipswithfit[(int) iL]->SetHeader("Energy");
    leg_devsmipswithfit[(int) iL]->SetFillColor(17);
  }
  TLegend* leg_fit = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg_fit->SetHeader("Layers");
  leg_fit->SetFillColor(17);

  TString procNa[numberofenergies];
  procNa[0] = "15 GeV";
  procNa[1] = "30 GeV";
  procNa[2] = "50 GeV";
  procNa[3] = "80 GeV";
  procNa[4] = "100 GeV";
  procNa[5] = "150 GeV";

  TString procNa_layers[((*ssvec).size()-2)];
  procNa_layers[0] = "Before calo";
  procNa_layers[1] = "Layer 1";
  procNa_layers[2] = "Layer 2";
  procNa_layers[3] = "Layer 3";

  //======================================================================================
  //The files that we will store the results of this analysis for the combined plot
  TString res_com[numberofconfigs];
  TFile* results_com[numberofconfigs];
  for (int k=0; k<numberofconfigs; k++){
    TString fp = "SiWEcal_" + IntToString(configs[k]);	
    res_com[k] = fp + "_combinedplots.root";
    //std::cout << res[k] << std::endl;
  }
  TH2F *hist[((*ssvec).size()-2)];
  TH1F *hist_sf[((*ssvec).size()-2)];
  TGraphErrors *gr[((*ssvec).size()-2)];
  // //use a TProfile to convert the 2-d to 1-d fit problem
  // TProfile *prof[((*ssvec).size()-2)];
  //For accessing the fit results
  TF1 *fit[((*ssvec).size()-2)][numberofenergies];
  // TF1 *line0 = new TF1("line0","[0]*x+[1]",0.,20);
  // line0->SetParameter(0, 1);
  //line0->SetParameter(1, 0);
  TGraphErrors *gr_dEvsMIPs_fit[((*ssvec).size()-2)];
  // TMultiGraph* mg = new TMultiGraph();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TString titleofplot1,titleofplot2,titleofplot3,titleofplot4,titleofplot; 
  for (int k=0; k<numberofconfigs; k++){
    results[k]= TFile::Open(res[k],"read");
    std::cout << "Results file " << res[k] << std::endl;
    
    //======================================================================================
    //Loop on energies
    for (int l=0; l<numberofenergies; l++){
      
      //======================================================================================
      //Loop on layers
      //Starting from 1 since 0 is before calorimeter
      for(unsigned iL(1); iL<((*ssvec).size()-2); iL++){

	std::string auxnam1 = "h_dEvsMIPs_" + IntToString((int) iL);
	std::string auxnam1_tgr = "dEvsMIPs_" + IntToString((int) iL);
	std::string auxnam2 = "_" + IntToString(particleenergies[l]);
	std::string auxnam = auxnam1 + auxnam2;
	std::string auxnam_tgr = auxnam1_tgr + auxnam2;
	

	//For the dE vs MIPs
	//------------------------------------------------------------------------------------------------
	std::cout << auxnam << std::endl;
	hist[(int) iL] = (TH2F*)results[k]->Get(auxnam.c_str());
	titleofplot1 = "Energy lost in materials in front of sensor vs measured energy for layer "; 
	titleofplot2 = IntToString( ((int) iL) ) + " and ";
	titleofplot3 =  particle +  " beam particle gun"; 
	titleofplot = titleofplot1 + titleofplot2 + titleofplot3;
	hist[(int) iL]->SetTitle(titleofplot); 
	hist[(int) iL]->GetXaxis()->SetTitle("Measured energy (MIPs)"); 
	hist[(int) iL]->GetYaxis()->SetTitle("dE in passive material (GeV)");
	// if (l==5) hist[1]->GetXaxis()->SetRangeUser(0.,(hist[(int) iL]->FindLastBinAbove(0.000001,1)) + 100);
	// hist[(int) iL]->SetMarkerSize(1.5);
	// c1[(int) iL]->SetLogy();
	if ( l == 0){hist[(int) iL]->SetMarkerColor(4);hist[(int) iL]->SetMarkerStyle(20);}
	if ( l == 1){hist[(int) iL]->SetMarkerColor(2);hist[(int) iL]->SetMarkerStyle(21);}
	if ( l == 2){hist[(int) iL]->SetMarkerColor(1);hist[(int) iL]->SetMarkerStyle(22);}
	if ( l == 3){hist[(int) iL]->SetMarkerColor(3);hist[(int) iL]->SetMarkerStyle(23);}
	if ( l == 4){hist[(int) iL]->SetMarkerColor(5);hist[(int) iL]->SetMarkerStyle(28);}
	if ( l == 5){hist[(int) iL]->SetMarkerColor(6);hist[(int) iL]->SetMarkerStyle(29);}
	//hist[(int) iL]->SetLineColor(l+1);
	leg[(int) iL]->AddEntry(hist[(int) iL], procNa[l], "P");
	// leg[(int) iL]->AddEntry(gr[(int) iL], procNa[l], "P");
	//hist[(int) iL]->GetYaxis()->SetRangeUser(0.,900.);
	c1[(int) iL]->cd();
	c1[(int) iL]->Update(); 
	// if (l==5){mg->Draw("ap same");}
	l == 0 ? hist[(int) iL]->Draw("HIST") : hist[(int) iL]->Draw("HISTsame");
	l == 0 ? leg[(int) iL]->Draw() : leg[(int) iL]->Draw("same");
	// l == 0 ? fit[(int) iL]->Draw() : fit[(int) iL]->Draw("same");
	c1[(int) iL]->Update(); 
	c1[(int) iL]->SaveAs("dEvsMIPs_"+particle+"_layer"+ IntToString( ((int) iL) ) + ".png");
	
	//For the scale factors
	//------------------------------------------------------------------------------------------------
	std::string auxnam1_sf = "h_SF_" + IntToString((int) iL);
	std::string auxnam_sf = auxnam1_sf + auxnam2;
 	hist_sf[(int) iL] = (TH1F*)results[k]->Get(auxnam_sf.c_str());
	titleofplot1 = "Scale factors for layer "; 
	titleofplot2 = IntToString( ((int) iL) ) + " and ";
	titleofplot3 =  particle +  " beam particle gun"; 
	titleofplot = titleofplot1 + titleofplot2 + titleofplot3;
	hist_sf[(int) iL]->SetTitle(titleofplot); 
	hist_sf[(int) iL]->GetXaxis()->SetTitle("Scale factor"); 
	hist_sf[(int) iL]->GetYaxis()->SetTitle("Events/0.001");
	// c1->SetLogy();
	if ( l == 0){hist_sf[(int) iL]->SetLineColor(4);}
	if ( l == 1){hist_sf[(int) iL]->SetLineColor(2);}
	if ( l == 2){hist_sf[(int) iL]->SetLineColor(1);}
	if ( l == 3){hist_sf[(int) iL]->SetLineColor(3);}
	if ( l == 4){hist_sf[(int) iL]->SetLineColor(5);}
	if ( l == 5){hist_sf[(int) iL]->SetLineColor(6);}
	//hist->SetLineColor(l+1);
	leg_sf[(int) iL]->AddEntry(hist_sf[(int) iL], procNa[l], "L");
	//hist->GetYaxis()->SetRangeUser(0.,900.);//  21: 100. 22: 900.
	c2[(int) iL]->cd();
	c2[(int) iL]->Update(); 
	l == 0 ? hist_sf[(int) iL]->Draw("HIST") : hist_sf[(int) iL]->Draw("HISTsame");
	l == 0 ? leg_sf[(int) iL]->Draw() : leg_sf[(int) iL]->Draw("same");
	c2[(int) iL]->Update(); 

	//For the fit 
	//------------------------------------------------------------------------------------------------
	gr[(int) iL] = (TGraphErrors*)results[k]->Get(auxnam_tgr.c_str());
	gr[(int) iL]->GetXaxis()->SetRangeUser(0.,5000.);
	gr[(int) iL]->Fit("pol1");
	
	// prof[(int) iL] = hist[(int) iL]->ProfileX();
	// prof[(int) iL]->Fit("pol1");
	fit[(int) iL][l] = gr[(int) iL]->GetFunction("pol1");
	//Results of the fit
	std::cout << "First parameter: " << fit[(int) iL][l]->GetParameter(0) << " Second parameter: " << fit[(int) iL][l]->GetParameter(1) << std::endl;
	
	titleofplot1 = "Energy lost in materials in front of sensor vs measured energy for layer "; 
	titleofplot2 = IntToString( ((int) iL) ) + " and ";
	titleofplot3 =  particle +  " "; 
	titleofplot4 =  IntToString(particleenergies[l]) +  " GeV beam particle gun"; 
	titleofplot = titleofplot1 + titleofplot2 + titleofplot3 + titleofplot4;
	gr[(int) iL]->SetTitle(titleofplot); 
	gr[(int) iL]->GetHistogram()->SetXTitle("Measured energy (MIPs)"); 
	gr[(int) iL]->GetHistogram()->SetYTitle("dE in passive material (GeV)");
	if ( l == 0){fit[(int) iL][l]->SetLineColor(4);gr[(int) iL]->SetMarkerColor(4);gr[(int) iL]->SetMarkerStyle(20);}
	if ( l == 1){fit[(int) iL][l]->SetLineColor(2);gr[(int) iL]->SetMarkerColor(2);gr[(int) iL]->SetMarkerStyle(21);}
	if ( l == 2){fit[(int) iL][l]->SetLineColor(1);gr[(int) iL]->SetMarkerColor(1);gr[(int) iL]->SetMarkerStyle(22);}
	if ( l == 3){fit[(int) iL][l]->SetLineColor(3);gr[(int) iL]->SetMarkerColor(3);gr[(int) iL]->SetMarkerStyle(23);}
	if ( l == 4){fit[(int) iL][l]->SetLineColor(5);gr[(int) iL]->SetMarkerColor(5);gr[(int) iL]->SetMarkerStyle(28);}
	if ( l == 5){fit[(int) iL][l]->SetLineColor(6);gr[(int) iL]->SetMarkerColor(6);gr[(int) iL]->SetMarkerStyle(29);}
	//hist[(int) iL]->SetLineColor(l+1);
	// leg_devsmipswithfit[(int) iL]->AddEntry(hist[(int) iL], procNa[l], "P");
	// leg_devsmipswithfit[(int) iL]->AddEntry(gr[(int) iL], procNa[l], "P");
	//hist[(int) iL]->GetYaxis()->SetRangeUser(0.,900.);//  21: 100. 22: 900.
	c4[(int) iL][l]->cd();
	c4[(int) iL][l]->Update(); 
	gr[(int) iL]->Draw("AP");
	// l == 0 ? gr[(int) iL]->Draw("AP") : gr[(int) iL]->Draw("APS");
	// leg_devsmipswithfit[(int) iL]->Draw("PS");
	// l == 0 ? leg_devsmipswithfit[(int) iL]->Draw() : leg_devsmipswithfit[(int) iL]->Draw("same");
	// l == 0 ? fit[(int) iL]->Draw() : fit[(int) iL]->Draw("same");
	gr[(int) iL]->GetHistogram()->GetXaxis()->SetRangeUser(0.,5000.);

	char buf[500];
	TLatex lat;
	double latx = gr[(int) iL]->GetHistogram()->GetXaxis()->GetXmin()+(gr[(int) iL]->GetHistogram()->GetXaxis()->GetXmax()-gr[(int) iL]->GetHistogram()->GetXaxis()->GetXmin())/20.;
	double laty = gr[(int) iL]->GetHistogram()->GetMaximum();
	sprintf(buf,"dE = p0 + p1 * MIPs");
	lat.DrawLatex(latx,laty*0.8,buf);
	sprintf(buf,"p0 = %3.3f +/- %3.3f",fit[(int) iL][l]->GetParameter(0),fit[(int) iL][l]->GetParError(0));
	lat.DrawLatex(latx,laty*0.7,buf);
	sprintf(buf,"p1 = %3.3f +/- %3.3f",fit[(int) iL][l]->GetParameter(1),fit[(int) iL][l]->GetParError(1));
	lat.DrawLatex(latx,laty*0.6,buf);
	sprintf(buf,"#chi^{2}/N = %3.3f/%d = %3.3f",fit[(int) iL][l]->GetChisquare(),fit[(int) iL][l]->GetNDF(),fit[(int) iL][l]->GetChisquare()/fit[(int) iL][l]->GetNDF());
	lat.DrawLatex(latx,laty*0.5,buf);
 

 	c4[(int) iL][l]->Update(); 
	TString cantopng1 = "dEvsMIPs_plusfit_"+ particle;
	TString cantopng2 = "_" + IntToString(particleenergies[l]);
	TString cantopng3 = "GeV_layer" + IntToString((int) iL);
	TString cantopng = cantopng1 + cantopng2 + cantopng3 + ".png";
	
	c4[(int) iL][l]->SaveAs(cantopng);

      }//Loop on layers
      
    }//Loop on energies

    results_com[k]= new TFile(res_com[k],"recreate");
    
    //======================================================================================
    //Write canvas with combined plot 
    // c1->Print("Leakage.pdf",".pdf");
    for(unsigned iL(1); iL<((*ssvec).size()-2); iL++){
      c1[(int) iL]->Write();
      c2[(int) iL]->Write();
      for (int l=0; l<numberofenergies; l++){
	c4[(int) iL][l]->Write();
      }
      gr[(int) iL]->Write();
      
    }
    //Reset histos
    // hist->Reset();

    //Here we will make the plot for the fit results
    c3->cd();

    for(unsigned iL(1); iL<((*ssvec).size()-2); iL++){
      std::string auxnam1_tgr = "dEoverMIPs_" + IntToString((int) iL);
      // gr_dEvsMIPs_fit[(int) iL]= new TGraphErrors(auxnam1_tgr.c_str());
      gr_dEvsMIPs_fit[(int) iL]= (TGraphErrors *) calib->Clone( auxnam1_tgr.c_str() );

      
      gr_dEvsMIPs_fit[(int) iL]->SetMarkerStyle( 19+((int) iL)); 
      
      // gr_dEvsMIPs_fit[(int) iL]->SetMarkerSize(0.2); 
      gr_dEvsMIPs_fit[(int) iL]->SetMarkerColor(((int) iL));  
      gr_dEvsMIPs_fit[(int) iL]->SetLineColor(((int) iL));
 
      c3->Update(); 

      // Int_t np=gr_dEvsMIPs_fit[(int) iL]->GetN();
      for (int l=0; l<numberofenergies; l++){
 	gr_dEvsMIPs_fit[(int) iL]->SetPoint( l , particleenergies[l] , fit[(int) iL][l]->GetParameter(1) );
    	gr_dEvsMIPs_fit[(int) iL]->SetPointError(l , 0. , fit[(int) iL][l]->GetParError(1) );
     }
      leg_fit->AddEntry(gr_dEvsMIPs_fit[(int) iL], procNa_layers[(int) iL], "LP");
      ((int) iL) ==1 ? gr_dEvsMIPs_fit[(int) iL]->Draw("APL") : gr_dEvsMIPs_fit[(int) iL]->Draw("PSL");
    }
    leg_fit->Draw("PS");
    c3->Update(); 

    gr_dEvsMIPs_fit[1]->SetTitle( "Energy lost/MIPs fit results for " + particle + " for all beam energies" );
    gr_dEvsMIPs_fit[1]->GetHistogram()->SetXTitle("Beam Energy (GeV)"); 
    gr_dEvsMIPs_fit[1]->GetHistogram()->SetYTitle("dE_loss = f(Measured energy) slope (GeV/MIP)");
    gr_dEvsMIPs_fit[1]->GetXaxis()->SetRangeUser(0.,160.);  
    gr_dEvsMIPs_fit[1]->GetYaxis()->SetRangeUser(-1.,1.);  

    c3->Update(); 
    c3->SaveAs("dEvsMIPs_fitresults_" + particle + ".png"); 
    c3->Write();

    std::ofstream myfile;
    TString titleoffile = "dEvsMIPs_" + particle + ".txt";
    myfile.open(titleoffile);

    //Here print the p0 and p1 that are going to be used for Ereco
    for (int l=0; l<numberofenergies; l++){
      std::cout << "=========================================================" << std::endl;
      std::cout << "Particle " <<  particle << std::endl;
      std::cout << "Energy of the beam " <<  particleenergies[l] << std::endl;
      for(unsigned iL(1); iL<((*ssvec).size()-2); iL++){
	std::cout << "---------------------------------------------------------" << std::endl;
	std::cout << "Layer " << iL << std::endl;
	std::cout << "dEvsMIPs slope: " <<  fit[(int) iL][l]->GetParameter(1) << " error slope " << fit[(int) iL][l]->GetParError(1) << std::endl;
	std::cout << "dEvsMIPs constant: " <<  fit[(int) iL][l]->GetParameter(0) << " error constant " << fit[(int) iL][l]->GetParError(0) << std::endl;

	//particle particleenergy layer constant slope errorconstant errorslope
	myfile << particle << " " << particleenergies[l] << " " << iL << " " << fit[(int) iL][l]->GetParameter(0) << " " << fit[(int) iL][l]->GetParameter(1) << " " << fit[(int) iL][l]->GetParError(0) <<  " " << fit[(int) iL][l]->GetParError(1) << "\n";


      }
    }


















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
