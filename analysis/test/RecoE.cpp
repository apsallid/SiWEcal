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
  TString particle = "e-"; //mu+
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

  //======================================================================================
  // Histos
  //Reco energy
  std::vector<TH1F*> h_recoE_raw;
  h_recoE_raw.clear();
  std::vector<TH1F*> h_recoE;
  h_recoE.clear();
  std::vector<TH1F*> h_recoE_corr;
  h_recoE_corr.clear();

  for (Int_t k=0; k<numberofenergies; k++){
    h_recoE.push_back(new TH1F(("h_recoE_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",30.,0.,30.));
    h_recoE_corr.push_back(new TH1F(("h_recoE_corr_" + IntToString(particleenergies[k])).c_str(),";E (GeV);Events",300.,0.,300.));
    h_recoE_raw.push_back(new TH1F(("h_recoE_raw_" + IntToString(particleenergies[k])).c_str(),";E (MeV);Events",500.,0.,500.));
  }

  //Open file and save the scale factors
  TString titleoffile1 = "dE_lossUpvsMIPs_" + particle + ".txt";
  std::ifstream infile1(titleoffile1);
  std::string dummyforparticle; double dummyforenergy;
  double dE_lossUpvsMIPs_fitresults_slope[numberofenergies];
  double dE_lossUpvsMIPs_fitresults_const[numberofenergies];
  double dE_lossUpvsMIPs_fitresults_slope_error[numberofenergies];
  double dE_lossUpvsMIPs_fitresults_const_error[numberofenergies];

  std::string buffer;
  for (Int_t k=0; k<numberofenergies; k++){
    getline(infile1, buffer);
    std::istringstream is(buffer);
    //particle particleenergy constant slope errorconstant errorslope
    is >> dummyforparticle >> dummyforenergy >> dE_lossUpvsMIPs_fitresults_const[k] >> dE_lossUpvsMIPs_fitresults_slope[k] >> dE_lossUpvsMIPs_fitresults_const_error[k] >> dE_lossUpvsMIPs_fitresults_slope_error[k];
    //Test
    std::cout << dummyforparticle << " " << dummyforenergy << " " << dE_lossUpvsMIPs_fitresults_const[k] << " " << dE_lossUpvsMIPs_fitresults_slope[k] << " " << " " << dE_lossUpvsMIPs_fitresults_const_error[k] << " " << dE_lossUpvsMIPs_fitresults_slope_error[k] << std::endl;

  }
  infile1.close();

  TString titleoffile2 = "dEvsMIPs_" + particle + ".txt";
  std::ifstream infile2(titleoffile2);
  const int numoflayers = 4;
  double dEvsMIPs_fitresults_slope[numberofenergies][numoflayers];
  double dEvsMIPs_fitresults_const[numberofenergies][numoflayers];
  double dEvsMIPs_fitresults_slope_error[numberofenergies][numoflayers];
  double dEvsMIPs_fitresults_const_error[numberofenergies][numoflayers];

  while ( getline(infile2,buffer) ) {
    std::istringstream is(buffer);
    //particle particleenergy layer constant slope errorconstant errorslope
    double currfitconst, currfitslope, currfitconsterr, currfitslopeerr;
    int currlayer;
    is >> dummyforparticle >> dummyforenergy >> currlayer >> currfitconst >> currfitslope >> currfitconsterr >> currfitslopeerr;
    int k = 0;
    if (dummyforenergy == 15){k=0;}
    if (dummyforenergy == 30){k=1;}
    if (dummyforenergy == 50){k=2;}
    if (dummyforenergy == 80){k=3;}
    if (dummyforenergy == 100){k=4;}
    if (dummyforenergy == 150){k=5;}

    dEvsMIPs_fitresults_const[k][currlayer] = currfitconst;
    dEvsMIPs_fitresults_slope[k][currlayer] = currfitslope;
    dEvsMIPs_fitresults_const_error[k][currlayer] = currfitconsterr;
    dEvsMIPs_fitresults_slope_error[k][currlayer] = currfitslopeerr;

    std::cout << dummyforparticle << " " << dummyforenergy << " " << currlayer << " " << dEvsMIPs_fitresults_const[k][currlayer] << " " << dEvsMIPs_fitresults_slope[k][currlayer] << " " << " " << dEvsMIPs_fitresults_const_error[k][currlayer] << " " << dEvsMIPs_fitresults_slope_error[k][currlayer] << std::endl;


  }

  infile2.close();
   
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
	if (ievt==1000){break;}
	// std::cout << "Entry " << ievt << std::endl;
	double recE = 0.;
	double recE_inGeV = 0.;
	double recE_raw = 0.;
	double recE_sf = 0.;
	double recE_corr = 0.;
	double lossinsensors = 0.;
	//Loop on sampling sections
	for(unsigned iL(0); iL<(*ssvec).size(); iL++){
	  // std::cout << "Sampling Section " << iL << " energyabsorber " << energyabsorber  << std::endl; 
	  
	  //RecoE
	  //4th layer is dead
	  if ( (iL!=0) && (iL!=4) && (iL!=5)  ){
	    //dEvsMIPs_fitresults_slope[l][(int) iL]: GeV/MIP
	    //So, (MeV/1000.) * (1./GeV/MIP) 
	    recE_raw = recE_raw + (*ssvec)[iL].measuredE()/(0.095);//in MIPs
	    recE_inGeV = recE_inGeV + ((*ssvec)[iL].measuredE()/0.095) * (dEvsMIPs_fitresults_slope[l][(int) iL]);
	    // recE_sf = recE_sf + ;
	    recE = recE +  ((*ssvec)[iL].measuredE()/1000.) * (1./dEvsMIPs_fitresults_slope[l][(int) iL]) ; //in MIPs
	    lossinsensors = lossinsensors + (*ssvec)[iL].measuredE(); //in MeV
	  }
	  if (iL==5){
	    h_recoE[l]->Fill( recE );//in MIPs
	    h_recoE_raw[l]->Fill( recE_raw );//in MIPs
	    //The correction for the upstream material
	    //Erec_corr = Erec + w0 + w1 * EmeasuredIn3layers
	    //dE_lossUpvsMIPs_fitresults_const[l] : MeV
	    //dE_lossUpvsMIPs_fitresults_slope[l] : MeV/MIP
	    //So, (MeV/0.095) + ( (1./MeV/MIP)* MeV 
	    // recE_corr = recE + (fabs(dE_lossUpvsMIPs_fitresults_const[l])/0.095) + ( fabs(1./dE_lossUpvsMIPs_fitresults_slope[l])*(lossinsensors) ) ;
	    recE_corr = (recE_inGeV*1000.) + (fabs(dE_lossUpvsMIPs_fitresults_const[l])) + ( fabs(dE_lossUpvsMIPs_fitresults_slope[l])*(lossinsensors/0.095) ) ;
	    // std::cout << "recE_corr " << recE_corr << std::endl;
	    // h_recoE_corr[l]->Fill( recE_corr );//in MIPs
	    h_recoE_corr[l]->Fill( recE_corr/1000. );//in GeV
	  }
	
	} //loop on sampling sections


      } //  Loop on entries

      //======================================================================================
      //Write histos 
      h_recoE[l]->Write();
      h_recoE_corr[l]->Write();

      //======================================================================================
      //We should here clear the histograms because we want them empty for the next file. 
      //Reset histos
      h_recoE[l]->Reset();
      h_recoE_corr[l]->Reset();
     
      
    }//Loop on energies
  
    results[k]->Close();

  } // Loop on files




  //======================================================================================
  //Make the following plots
  //mean(recoE) = f (beam) 
  //mean(recoE_corr) = f (beam)
  
  TCanvas *c1 = new TCanvas("c1", "  ");
  TCanvas *c2 = new TCanvas("c2", "  ");
  TCanvas *c3 = new TCanvas("c3", "  ");
  TCanvas *c4 = new TCanvas("c4", "  ");

  //For the legend
  TLegend* leg;
  leg = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg->SetHeader("Energy");
  leg->SetFillColor(17);
  TLegend* leg2;
  leg2 = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg2->SetHeader("Energy");
  leg2->SetFillColor(17);
  TLegend* leg3;
  leg3 = new TLegend(0.5, 0.7, 0.8, 0.9);
  leg3->SetHeader("Energy");
  leg3->SetFillColor(17);

  TString procNa[numberofenergies];
  procNa[0] = "15 GeV";
  procNa[1] = "30 GeV";
  procNa[2] = "50 GeV";
  procNa[3] = "80 GeV";
  procNa[4] = "100 GeV";
  procNa[5] = "150 GeV";
 
  //======================================================================================
  TGraphErrors * recoEvsbeamenergy;
  TGraphErrors * recoE_corrvsbeamenergy;
  TGraphErrors *recoEvsbeam = new TGraphErrors();
  recoEvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoEvsbeamenergy" );
  recoE_corrvsbeamenergy = (TGraphErrors *) recoEvsbeam->Clone( "recoE_corrvsbeamenergy" );
  // TGraphErrors *gr_fit;
  // //For accessing the fit results
  // TF1 *fit[numberofenergies];

  //======================================================================================
  //The files that we will store the results of this analysis for the combined plot
  TString res_com[numberofconfigs];
  TFile* results_com[numberofconfigs];
  for (int k=0; k<numberofconfigs; k++){
    TString fp = "SiWEcal_" + IntToString(configs[k]);	
    res_com[k] = fp + "_combinedplots_recoE.root";
    //std::cout << res[k] << std::endl;
  }
  TH1F *hist,*hist2;
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TString titleofplot1,titleofplot2,titleofplot; 
  for (int k=0; k<numberofconfigs; k++){
    results[k]= TFile::Open(res[k],"read");
    std::cout << "Results file " << res[k] << std::endl;
    
    //======================================================================================
    //Loop on energies
    for (int l=0; l<numberofenergies; l++){
      
      hist = (TH1F*)results[k]->Get( ("h_recoE_" + IntToString(particleenergies[l])).c_str());
      titleofplot1 = "Reconstructed energy for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist->SetTitle(titleofplot); 
      hist->GetXaxis()->SetTitle("Reconstructed energy (MIPs)"); 
      hist->GetYaxis()->SetTitle("Events/MIP");
      // if (l==5) hist[1]->GetXaxis()->SetRangeUser(0.,(hist[(int) iL]->FindLastBinAbove(0.000001,1)) + 100);
      // hist->SetMarkerSize(1.5);
      // c1->SetLogy();
      if ( l == 0){hist->SetLineColor(4);}
      if ( l == 1){hist->SetLineColor(2);}
      if ( l == 2){hist->SetLineColor(1);}
      if ( l == 3){hist->SetLineColor(3);}
      if ( l == 4){hist->SetLineColor(5);}
      if ( l == 5){hist->SetLineColor(6);}
      //hist[(int) iL]->SetLineColor(l+1);
      leg->AddEntry(hist, procNa[l], "L");
      // leg[(int) iL]->AddEntry(gr[(int) iL], procNa[l], "P");
      //hist->GetYaxis()->SetRangeUser(0.,900.);
      c1->cd();
      c1->Update(); 
      // if (l==5){mg->Draw("ap same");}
      l == 0 ? hist->Draw("HIST") : hist->Draw("HISTsame");
      l == 0 ? leg->Draw() : leg->Draw("same");
      // l == 0 ? fit->Draw() : fit->Draw("same");
      c1->Update(); 

      hist2 = (TH1F*)results[k]->Get( ("h_recoE_corr_" + IntToString(particleenergies[l])).c_str());
      titleofplot1 = "Reconstructed energy corrected for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      hist2->SetTitle(titleofplot); 
      hist2->GetXaxis()->SetTitle("Reconstructed energy corrected (GeV)"); 
      hist2->GetYaxis()->SetTitle("Events/GeV");
      // if (l==5) hist2[1]->GetXaxis()->SetRangeUser(0.,(hist2[(int) iL]->FindLastBinAbove(0.000001,1)) + 100);
      // hist2->SetMarkerSize(1.5);
      // c2->SetLogy();
      if ( l == 0){hist2->SetLineColor(4);}
      if ( l == 1){hist2->SetLineColor(2);}
      if ( l == 2){hist2->SetLineColor(1);}
      if ( l == 3){hist2->SetLineColor(3);}
      if ( l == 4){hist2->SetLineColor(5);}
      if ( l == 5){hist2->SetLineColor(6);}
      //hist2[(int) iL]->SetLineColor(l+1);
      leg2->AddEntry(hist2, procNa[l], "L");
      // leg[(int) iL]->AddEntry(gr[(int) iL], procNa[l], "P");
      //hist2->GetYaxis()->SetRangeUser(0.,900.);
      c2->cd();
      c2->Update(); 
      // if (l==5){mg->Draw("ap same");}
      l == 0 ? hist2->Draw("HIST") : hist2->Draw("HISTsame");
      l == 0 ? leg2->Draw() : leg2->Draw("same");
      // l == 0 ? fit->Draw() : fit->Draw("same");
      c2->Update(); 

      c3->cd();
      c3->Update(); 
      recoEvsbeamenergy->SetPoint( l , particleenergies[l] , hist->GetMean() );
      recoEvsbeamenergy->SetPointError(l , 0. , hist->GetRMS() );
      l==0 ? recoEvsbeamenergy->Draw("APL") : recoEvsbeamenergy->Draw("PSL");
      recoEvsbeamenergy->GetXaxis()->SetRangeUser(0.,160.);
      titleofplot1 = "Reconstructed energy for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      recoEvsbeamenergy->SetTitle(titleofplot); 
      recoEvsbeamenergy->GetHistogram()->SetYTitle("Reconstructed energy (MIPs)"); 
      recoEvsbeamenergy->GetHistogram()->SetXTitle("Beam Energy (GeV)");
      c3-> Update(); 

      c4->cd();
      c4->Update(); 
      recoE_corrvsbeamenergy->SetPoint( l , particleenergies[l] , hist2->GetMean() );
      recoE_corrvsbeamenergy->SetPointError(l , 0. , hist2->GetRMS() );
      l==0 ? recoE_corrvsbeamenergy->Draw("APL") : recoE_corrvsbeamenergy->Draw("PSL");
      recoE_corrvsbeamenergy->GetXaxis()->SetRangeUser(0.,160.);
      titleofplot1 = "Reconstructed energy corrected for "; 
      titleofplot2 =  particle +  " and all beam energies"; 
      titleofplot = titleofplot1 + titleofplot2;
      recoE_corrvsbeamenergy->SetTitle(titleofplot); 
      recoE_corrvsbeamenergy->GetHistogram()->SetYTitle("Reconstructed energy corrected (GeV)"); 
      recoE_corrvsbeamenergy->GetHistogram()->SetXTitle("Beam Energy (GeV)");
      c4-> Update(); 

      // char buf[500];
      // TLatex lat;
      // double latx = recoEvsbeamenergy->GetHistogram()->GetXaxis()->GetXmin()+(recoEvsbeamenergy->GetHistogram()->GetXaxis()->GetXmax()-recoEvsbeamenergy->GetHistogram()->GetXaxis()->GetXmin())/20.;
      // double laty = recoEvsbeamenergy->GetHistogram()->GetMaximum();
      // sprintf(buf,"dE = p0 + p1 * MIPs");
      // lat.DrawLatex(latx,laty*0.8,buf);
      // sprintf(buf,"p0 = %3.3f +/- %3.3f",fit[(int) iL][l]->GetParameter(0),fit[(int) iL][l]->GetParError(0));
      // lat.DrawLatex(latx,laty*0.7,buf);
      // sprintf(buf,"p1 = %3.3f +/- %3.3f",fit[(int) iL][l]->GetParameter(1),fit[(int) iL][l]->GetParError(1));
      // lat.DrawLatex(latx,laty*0.6,buf);
      // sprintf(buf,"#chi^{2}/N = %3.3f/%d = %3.3f",fit[(int) iL][l]->GetChisquare(),fit[(int) iL][l]->GetNDF(),fit[(int) iL][l]->GetChisquare()/fit[(int) iL][l]->GetNDF());
      // lat.DrawLatex(latx,laty*0.5,buf);
 

 	// c4[(int) iL][l]->Update(); 
	// TString cantopng1 = "dEvsMIPs_plusfit_"+ particle;
	// TString cantopng2 = "_" + IntToString(particleenergies[l]);
	// TString cantopng3 = "GeV_layer" + IntToString((int) iL);
	// TString cantopng = cantopng1 + cantopng2 + cantopng3 + ".png";
	
	// c4[(int) iL][l]->SaveAs(cantopng);

      
      


    }//Loop on energies

    results_com[k]= new TFile(res_com[k],"recreate");
    
    //======================================================================================
    //Write canvas with combined plot 
    c1->Write();
    c2->Write();
    c3->Write();
    c3->SaveAs("RecoEvsBeamEne.png");
    c4->Write();
    c4->SaveAs("RecoE_corrvsBeamEne.png");
    recoEvsbeamenergy->Write();
    recoE_corrvsbeamenergy->Write();

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
