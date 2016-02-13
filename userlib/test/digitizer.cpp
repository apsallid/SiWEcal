#include<string>
#include<set>
#include<iostream>
#include<fstream>
#include<sstream>
#include <boost/algorithm/string.hpp>

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TNtuple.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom3.h"
#include "Math/Vector4D.h"

#include "SiWEcalSSEvent.hh"
#include "SiWEcalSSInfo.hh"
#include "SiWEcalSSSimHit.hh"
#include "SiWEcalSSRecoHit.hh"
#include "SiWEcalSSCalibration.hh"
#include "SiWEcalSSDigitisation.hh"
#include "SiWEcalSSDetector.hh"
#include "SiWEcalSSGeometryConversion.hh"

template <class T>
void extractParameterFromStr(std::string aStr,T & vec){ 
  if (aStr == "") return;
  std::vector<std::string> layVec;
  boost::split( layVec, aStr, boost::is_any_of(","));

  for (unsigned iE(0); iE<layVec.size(); ++iE){//loop on elements
    std::vector<std::string> lPair;
    boost::split( lPair, layVec[iE], boost::is_any_of(":"));
    if (lPair.size() != 2) {
      std::cout << " -- Wrong string for parameter given as input:" << layVec[iE] << " Try again, expecting exactly one symbol \":\" between two \",\" ..." << std::endl;
      exit(1);
    }
    std::vector<std::string> lLay;
    boost::split( lLay, lPair[0], boost::is_any_of("-"));
    if (lLay.size() > 2) {
      std::cout << " -- Wrong string for granularities given as input:" << lPair[0] << " Try again, expecting at most one symbol \"-\"." << std::endl;
      exit(1);
    }
    unsigned beginIdx =  atoi(lLay[0].c_str());
    unsigned endIdx = lLay.size() == 1 ? beginIdx :  atoi(lLay[1].c_str());
    for (unsigned iL(beginIdx); iL<endIdx+1; ++iL){
      if (iL < vec.size())
	std::istringstream(lPair[1])>>vec[iL];
      else {
	std::cout << " -- WARNING! Input parameter has more layers: " << endIdx << " than detector : " 
		  << vec.size()
		  << ". Ignoring additional layer #" << iL << "... PLEASE CHECK SETTINGS ARE CORRECT FOR EXISTING LAYERS!!"
		  << std::endl;
      }
    }
  }//loop on elements
}

void processHist(const unsigned iL,
		 const bool doSmall,
		 const TH2D* histE,
		 SiWEcalSSDigitisation & myDigitiser,
		 TH1F* & p_noise,
		 const TH2D* histZ,
		 const double & meanZpos,
		 // const bool isTBsetup,
		 const SiWEcalSSSubDetector & subdet,
		 const std::vector<unsigned> & pThreshInADC,
		 const bool pSaveDigis,
		 SiWEcalSSRecoHitVec & lDigiHits,
		 SiWEcalSSRecoHitVec & lRecoHits
		 ){

  DetectorEnum adet = subdet.type;
  bool isScint = subdet.isScint;
  bool isSi = subdet.isSi;

  double rLim = subdet.radiusLim;
  
  for (int iX(1); iX<histE->GetNbinsX()+1;++iX){
    for (int iY(1); iY<histE->GetNbinsY()+1;++iY){

      double x = histE->GetXaxis()->GetBinCenter(iX);
      double y = histE->GetYaxis()->GetBinCenter(iY);
      
      double radius = sqrt(pow(x,2)+pow(y,2));

      //if (iY==1) std::cout << " Layer " << iL << " " << x << " " << y << std::endl;

      if ( (doSmall && radius >= rLim) ||
	   (!doSmall && radius < rLim && !isScint) ) continue;

      double digiE = 0;
      double simE = histE->GetBinContent(iX,iY);
      // double time = 0;
      // if (simE>0) time = histTime->GetBinContent(iX,iY)/simE;
      
      //bool passTime = myDigitiser.passTimeCut(adet,time);
      //if (!passTime) continue;
      
      double posz = meanZpos;
      //for noise only hits
      // double posz = 0.;
      if (simE>0) posz = histZ->GetBinContent(iX,iY)/simE;
      else posz = meanZpos;

      // if (posz==0.){std::cout << " x=" << x << ", y=" << y << std::endl;}
      
      //if (fabs(x) > 500 || fabs(y)>500) std::cout << " x=" << x << ", y=" << y << std::endl;
      

      //correct for particle angle in conversion to MIP
      //This is simE in our case. Just like previour lines are missing.
      // double simEcor = isTBsetup ? xtalkE : myDigitiser.mipCor(xtalkE,x,y,posz);
      digiE = simE;

      // std::cout << "simE " << simE << std::endl;
      
      // if (isScint && simEcor>0) {
      // 	digiE = myDigitiser.digiE(simEcor);
      // }

      myDigitiser.addNoise(digiE,iL,p_noise);
      
      // double noiseFrac = 1.0;
      // if (simEcor>0) noiseFrac = (digiE-simEcor)/simEcor;
      
      //for silicon-based Calo
      unsigned adc = 0;
      if (isSi){
	adc = myDigitiser.adcConverter(digiE,adet);
	digiE = myDigitiser.adcToMIP(adc,adet);
      }

      // if(simE>0.)std::cout << "simE " << simE << " adc " <<  adc << " energy in digis " << digiE << std::endl;

      double noiseFrac = 1.0;
      if (simE>0) noiseFrac = (digiE-simE)/simE;
      
      // bool aboveThresh = digiE > (pThreshInADC[iL]*myDigitiser.adcToMIP(1,adet,false));
      bool aboveThresh = (isSi && adc > pThreshInADC[iL]); 
      // (isScint && digiE > pThreshInADC[iL]*myDigitiser.adcToMIP(1,adet,false));
      //histE->SetBinContent(iX,iY,digiE);
      // if ((!pSaveDigis && aboveThresh) ||
      // if ((!pSaveDigis && aboveThresh) ||
      // if ( ( (!pSaveDigis) && aboveThresh)  || pSaveDigis)
      if ( aboveThresh )
	{//save hits
	  //double calibE = myDigitiser.MIPtoGeV(subdet,digiE);
	  SiWEcalSSRecoHit lRecHit;
	  lRecHit.layer(iL);
	  lRecHit.energy(digiE);
	  lRecHit.adcCounts(adc);
	  lRecHit.x(x);
	  lRecHit.y(y);
	  lRecHit.z(posz);
	  lRecHit.noiseFraction(noiseFrac);
	  //unsigned x_cell = static_cast<unsigned>(fabs(x)/(cellSize*granularity[iL]));
	  //unsigned y_cell = static_cast<unsigned>(fabs(y)/(cellSize*granularity[iL]));
	  //lRecHit.encodeCellId(x>0,y>0,x_cell,y_cell,granularity[iL]);
	  
	  if (pSaveDigis) lDigiHits.push_back(lRecHit);
	  
	  lRecoHits.push_back(lRecHit);
	  
	}//save hits
    }//loop on y
  }//loop on x
     
}//processHist



int main(int argc, char** argv){//main  

  /////////////////////////////////////////////////////////////
  //parameters
  /////////////////////////////////////////////////////////////
  const unsigned nReqA = 9;
  const unsigned nPar = static_cast<unsigned>(argc);
  if (nPar < nReqA-1) {
    std::cout << " Usage: "
              << argv[0] << " <nEvts to process (0=all)>"<< std::endl
              << "<full path to input file>"<< std::endl
              << "<full path to output file>"<< std::endl
              << "<granularities \"layer_i-layer_j:factor,layer:factor,...\">"<< std::endl
              << "<noise (in Mips) \"layer_i-layer_j:factor,layer:factor,...\">"<< std::endl
              << "<threshold (in ADC counts) \"layer_i-layer_j:factor,layer:factor,...\">"<< std::endl
              << "<intercalib factor in %>" << std::endl
	      << "<Number of si layers for TB setups>" << std::endl
              << std::endl
              << "<optional: randomSeed (default=0)> "  << std::endl
              << "<optional: debug (default=0)>" << std::endl
              << "<optional: save sim hits (default=0)> " << std::endl
              << "<optional: save digi hits (default=0)> " << std::endl
              << std::endl;
    return 1;
  }

  const unsigned pNevts = atoi(argv[1]);
  std::string inFilePath = argv[2];
  std::string outFilePath = argv[3];
  std::string granulStr = argv[4];
  std::string noiseStr = argv[5];
  std::string threshStr = argv[6];
  const unsigned interCalib = atoi(argv[7]);
  const unsigned nSiLayers = atoi(argv[8]);

  std::cout << " ----------------------------------------" << std::endl
            << " -- Input parameters: " << std::endl
            << " -- Input file path: " << inFilePath << std::endl
            << " -- Output file path: " << outFilePath << std::endl
            << " -- Processing " ;
  if (pNevts>0) std::cout << pNevts;
  else std::cout << "all";
  std::cout << " events." << std::endl
            << " -- Granularities: " << granulStr << std::endl
            << " -- noise: " << noiseStr << std::endl
            << " -- thresholds: " << threshStr << std::endl
            << " -- intercalibration factor (in %): " << interCalib << std::endl;

	    
  //std::string pModel = "model2";
  unsigned debug = 0;
  unsigned pSeed = 0;
  bool pSaveDigis = 0;
  bool pSaveSims = 0;

  if (nPar > nReqA) std::istringstream(argv[nReqA])>>pSeed;
  if (nPar > nReqA+1) {
    debug = atoi(argv[nReqA+1]);
    std::cout << " -- DEBUG output is set to " << debug << std::endl;
  }
  if (nPar > nReqA+2) std::istringstream(argv[nReqA+2])>>pSaveSims;
  if (nPar > nReqA+3) std::istringstream(argv[nReqA+3])>>pSaveDigis;
  
  std::cout<< " -- Random seed will be set to : " << pSeed << std::endl;
  if (pSaveDigis) std::cout << " -- DigiHits are saved." << std::endl;
  if (pSaveSims) std::cout << " -- SimHits are saved." << std::endl;
  std::cout << " ----------------------------------------" << std::endl;
  
  //input
  std::string inputStr = inFilePath ;
  TFile *inputFile = TFile::Open(inputStr.c_str());
  
  if (!inputFile) {
    std::cout << " -- Error, input file " << inputStr << " cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  TTree *inputTree = (TTree*)inputFile->Get("SiWEcalSSTree");
  if (!inputTree){
    std::cout << " -- Error, tree SiWEcalSSTree  cannot be opened. Exiting..." << std::endl;
    return 1;
  }

  SiWEcalSSInfo * info=(SiWEcalSSInfo*)inputFile->Get("Info");
  const double cellSize = info->cellSize();
  const unsigned versionNumber = info->version();
  const unsigned model = info->model();
  
  bool isTBsetup = true;
  if (isTBsetup) std::cout << " -- Number of Si layers: " << nSiLayers << std::endl;
  else std::cout << " -- Number of Si layers ignored: hardcoded as a function of radius in SiWEcalSSGeometryConversion class." << std::endl;

  SiWEcalSSEvent * event=0;
  std::vector<SiWEcalSSSimHit> * hitvec = 0;

  inputTree->SetBranchAddress("SiWEcalSSEvent",&event);
  inputTree->SetBranchAddress("SiWEcalSSSimHitVec",&hitvec);
    
  //initialise detector
  SiWEcalSSDetector & myDetector = theDetector();

  std::cout << " -- Version number is : " << versionNumber 
	    << ", model = " << model
	    << ", cellSize = " << cellSize
	    << std::endl;

  bool bypassR = false;
  if (isTBsetup) bypassR = true;
  myDetector.buildDetector(versionNumber,bypassR);

  //initialise calibration class
  SiWEcalSSCalibration mycalib(inFilePath,bypassR,nSiLayers);

  const unsigned nLayers = myDetector.nLayers();

  SiWEcalSSGeometryConversion geomConv(inFilePath,model,cellSize,bypassR,nSiLayers);
  geomConv.setXYwidth(88.48); // Single Si Sensor width.
  const double xWidth = geomConv.getXYwidth();
  std::cout << " -- XYwidth of silicon pads " << xWidth << std::endl;

  SiWEcalSSDigitisation myDigitiser;

  std::vector<unsigned> granularity;
  granularity.resize(nLayers,1); //keeping the virtual sim cell size for the rec hit cell size also
  std::vector<double> pNoiseInMips;
  pNoiseInMips.resize(nLayers,0.05556); //Initialize with 0.05556 MIPS for noise, but can change on input
  std::vector<unsigned> pThreshInADC;
  pThreshInADC.resize(nLayers,17); // 5 x noise ~ 17 ADC channels. Offline only. Not cut here. 

  extractParameterFromStr<std::vector<unsigned> >(granulStr,granularity);
  extractParameterFromStr<std::vector<double> >(noiseStr,pNoiseInMips);
  extractParameterFromStr<std::vector<unsigned> >(threshStr,pThreshInADC);

  //unsigned nbCells = 0;

  std::cout << " -- Granularities and noise are setup like this:" << std::endl;
  for (unsigned iL(0); iL<nLayers; ++iL){
    std::cout << "Layer " ;
    if (iL<10) std::cout << " ";
    std::cout << iL << " : " << granularity[iL] << ", " << pNoiseInMips[iL] << " mips, " << pThreshInADC[iL] << " adc - ";
    if (iL%5==4) std::cout << std::endl;
    myDigitiser.setNoise(iL,pNoiseInMips[iL]);
    //nbCells += N_CELLS_XY_MAX/(granularity[iL]*granularity[iL]);
  }
  std::cout << std::endl;     
  //std::cout << " -- Total number of cells = " << nbCells << std::endl;

  geomConv.setGranularity(granularity);
  geomConv.initialiseHistos();

  TRandom3 *lRndm = new TRandom3();
  lRndm->SetSeed(pSeed);
  myDigitiser.setRandomSeed(pSeed);

  std::cout << " -- Random3 seed = " << lRndm->GetSeed() << std::endl
	    << " ----------------------------------------" << std::endl;

  myDigitiser.setIntercalibrationFactor(interCalib);

  //output
  std::ostringstream outputStr;
  outputStr << outFilePath << "/DigiPFcal" ;
  if (pSaveDigis)  outputStr << "_withDigiHits";
  if (pSaveSims)  outputStr << "_withSimHits";
  outputStr << ".root";
  
  TFile *outputFile = TFile::Open(outputStr.str().c_str(),"RECREATE");

  if (!outputFile) {
    std::cout << " -- Error, output file " << outputStr.str() << " cannot be opened. Exiting..." << std::endl;
    return 1;
  }
  else {
    std::cout << " -- File will be saved as " << outputStr.str() << std::endl;
  }

  SiWEcalSSInfo *lInfo = new SiWEcalSSInfo();
  lInfo->cellSize(cellSize);
  lInfo->version(versionNumber);
  lInfo->model(model);

  TTree *outputTree = new TTree("RecoTree","SiWEcal Standalone simulation reco tree");
  SiWEcalSSSimHitVec lSimHits;
  SiWEcalSSRecoHitVec lDigiHits;
  SiWEcalSSRecoHitVec lRecoHits;

  unsigned maxSimHits = 0;
  unsigned maxRecHits = 0;

  SiWEcalSSEvent lEvent;
  outputTree->Branch("SiWEcalSSEvent",&lEvent);
  if (pSaveSims) outputTree->Branch("SiWEcalSSSimHitVec","std::vector<SiWEcalSSSimHit>",&lSimHits);
  if (pSaveDigis) outputTree->Branch("SiWEcalSSDigiHitVec","std::vector<SiWEcalSSRecoHit>",&lDigiHits);
  outputTree->Branch("SiWEcalSSRecoHitVec","std::vector<SiWEcalSSRecoHit>",&lRecoHits);

  TH1F * p_noise = new TH1F("noiseCheck",";noise (MIPs)",1000,-5,5);

  /////////////////////////////////////////////////////////////
  //Loop on events
  /////////////////////////////////////////////////////////////

  const unsigned nEvts = (pNevts > inputTree->GetEntries() || pNevts==0) ? static_cast<unsigned>(inputTree->GetEntries()) : pNevts;

  std::cout << "- Processing = " << nEvts  << " events out of " << inputTree->GetEntries() << std::endl;

  for (unsigned ievt(0); ievt<nEvts; ++ievt){//loop on entries

    inputTree->GetEntry(ievt);
    lEvent.eventNumber(event->eventNumber());
    lEvent.vtx_x(event->vtx_x());
    lEvent.vtx_y(event->vtx_y());
    lEvent.vtx_z(event->vtx_z());
    //unsigned layer = volNb;
    
    mycalib.setVertex(lEvent.vtx_x(),lEvent.vtx_y(),lEvent.vtx_z());

    if (debug>0) {
      std::cout << " **DEBUG** Processing evt " << ievt << std::endl;
    }
    else if (ievt%50 == 0) std::cout << "... Processing event: " << ievt << std::endl;
    

    unsigned prevLayer = 10000;
    DetectorEnum type = DetectorEnum::SiWEcal;
    unsigned subdetLayer=0;
    for (unsigned iH(0); iH<(*hitvec).size(); ++iH){//loop on hits
      SiWEcalSSSimHit lHit = (*hitvec)[iH];

      //do not save hits with 0 energy...
      if (lHit.energy()>0 && pSaveSims) lSimHits.push_back(lHit);
      double layer = (4*lHit.layer()) + lHit.silayer(); 
      //unsigned layer = lHit.layer();
      if (layer != prevLayer){
	const SiWEcalSSSubDetector & subdet = myDetector.subDetectorByLayer(layer);
	type = subdet.type;
	subdetLayer = layer-subdet.layerIdMin;
	prevLayer = layer;
	if (debug > 1) std::cout << " - layer " << layer << " " << subdet.name << " " << subdetLayer << std::endl;
      }

      //This is the position of the hit with respect to virtual cells that divide each layer 0,1,2
      //So the whole layer contains silayers 0,1,2,3 
      //We should make transformation to the real pixel cells taking into account the 0.1 mm gap in the center
      double posx = lHit.get_x(cellSize);
      double posy = lHit.get_y(cellSize);
       std::cout << "posx " << posx << " posy " << posy <<std::endl;
      double radius = sqrt(pow(posx,2)+pow(posy,2));
      double posz = lHit.get_z();
      // double energy = lHit.energy()*mycalib.MeVToMip(layer,radius);
      double energy = lHit.energy()*mycalib.MeVToMip(layer,false);
      // std::cout << "Uncalib energy " << lHit.energy() << " MeVToMIP " << mycalib.MeVToMip(layer,false) << " Calib energy " <<  energy << std::endl;
      double realtime = mycalib.correctTime(lHit.time(),posx,posy,posz);
      bool passTime = myDigitiser.passTimeCut(type,realtime);
      if (!passTime) continue;
      // std::cout << "*******" << (4*lHit.layer()) + lHit.silayer() << " " << geomConv.getNumberOfSiLayers(type,radius) << std::endl;
      if (energy>0 && ( layer < geomConv.getNumberOfSiLayers(type,radius) ) ){
	if (debug > 1) std::cout << " hit " << iH 
				 << " lay " << layer  
				 << " x " << posx 
				 << " y " << posy
				 << " z " << posz
				 << " t " << lHit.time() << " " << realtime
				 << std::endl;
	geomConv.fill(type,subdetLayer,energy,realtime,posx,posy,posz);
      }

    }//loop on input simhits

    if (debug>1) {
      std::cout << " **DEBUG** simhits = " << (*hitvec).size() << " " << lSimHits.size() << std::endl;
    }

    //create hits, everywhere to have also pure noise
    //digitise
    //do not apply threshold
    //save
    unsigned nTotBins = 0;
    for (unsigned iL(0); iL<nLayers; ++iL){//loop on layers
      TH2D *histE = geomConv.get2DHist(iL,"E");
      // TH2D *histEs = geomConv.get2DHist(iL,"ESmall");
      // TH2D *histTime = geomConv.get2DHist(iL,"Time");
      TH2D *histZ = geomConv.get2DHist(iL,"Z");
      const SiWEcalSSSubDetector & subdet = myDetector.subDetectorByLayer(iL);
      bool isScint = subdet.isScint;
      nTotBins += histE->GetNbinsX()*histE->GetNbinsY();
      if (pSaveDigis) lDigiHits.reserve(nTotBins);

      //Here is the z position
      //meanZpos = (E_i * Z_i)/E_i
      double meanZpos = geomConv.getAverageZ(iL);
      
      //std::cout << iL << " " << meanZpos << std::endl;

      if (debug>1){
	std::cout << " -- Layer " << iL << " " << subdet.name << " z=" << meanZpos
		  << " totbins = " << nTotBins << " histE entries = " << histE->GetEntries() << std::endl;
      }

      //cell-to-cell cross-talk for scintillator
      if (isScint){
	//2.5% per 30-mm edge
	myDigitiser.setIPCrossTalk(0.025*histE->GetXaxis()->GetBinWidth(1)/30.);
      }
      else {
	myDigitiser.setIPCrossTalk(0);
      }

      processHist(iL,false,histE,myDigitiser,p_noise,histZ,meanZpos,subdet,pThreshInADC,pSaveDigis,lDigiHits,lRecoHits);
 
    }//loop on layers

    if (debug) {
      std::cout << " **DEBUG** sim-digi-reco hits = " << (*hitvec).size() << "-" << lDigiHits.size() << "-" << lRecoHits.size() << std::endl;
    }
    
    
    outputTree->Fill();
    //reserve necessary space and clear vectors.
    if (lSimHits.size() > maxSimHits) {
      maxSimHits = 2*lSimHits.size();
      std::cout << " -- INFO: event " << ievt << " maxSimHits updated to " << maxSimHits << std::endl;
    }
    if (lRecoHits.size() > maxRecHits) {
      maxRecHits = 2*lRecoHits.size();
      std::cout << " -- INFO: event " << ievt << " maxRecHits updated to " << maxRecHits << std::endl;
    }
    lSimHits.clear();
    lDigiHits.clear();
    lRecoHits.clear();
    geomConv.initialiseHistos();
    if (pSaveSims) lSimHits.reserve(maxSimHits);
    lRecoHits.reserve(maxRecHits);
    
  }//loop on entries

  outputFile->cd();
  outputFile->WriteObjectAny(lInfo,"SiWEcalSSInfo","Info");
  outputTree->Write();
  p_noise->Write();
  outputFile->Close();

  return 0;

}//main
