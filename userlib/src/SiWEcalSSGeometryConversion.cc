#include "SiWEcalSSGeometryConversion.hh"
#include <sstream>
#include <iostream>
#include <cmath>

SiWEcalSSGeometryConversion::SiWEcalSSGeometryConversion(std::string filePath, const unsigned & model, const double & cellsize, const bool bypassR, const unsigned nSiLayers){

  dopatch_=false;
  width_ = 200;//mm
  model_ = model;
  if (model==1) width_ = 500;
  else if (model == 2) width_ = 1700*2;
  else if (model == 3) width_ = 1000;
  else if (model == 4) width_ = 1300;
  
  cellSize_ = cellsize;
  bypassRadius_ = bypassR;
  nSiLayers_ = nSiLayers;
}
/*
unsigned SiWEcalSSGeometryConversion::getNumberOfSiLayers(const DetectorEnum type,
						      const double & eta) const
{
  if (model_ != 2) return 3;
  if (type == DetectorEnum::FHCAL) return 3;
  unsigned etaBin = 0;
  if (fabs(eta)>=1.4 && fabs(eta)<=1.75) etaBin = 1;
  else if (fabs(eta)>1.75 && fabs(eta)<=2.15) etaBin = 2;
  else if (fabs(eta) > 2.15) etaBin = 3;
  if (etaBin==0){
    if (type == DetectorEnum::FECAL) return 2;
    else if (type == DetectorEnum::MECAL) return 2;
    else if (type == DetectorEnum::BECAL) return 2;
  }
  else {
    if (etaBin==1) return 3;
    else if (etaBin==2) return 2;
    else return 1;
  }
  return 3;
  }*/

unsigned SiWEcalSSGeometryConversion::getNumberOfSiLayers(const DetectorEnum type,
						      const double & radius) const
{
  if (theDetector().subDetectorByEnum(type).isScint) return 3;
  if (model_ != 2) return nSiLayers_;

  double r1 = 1200;
  if (type == DetectorEnum::FHCAL) {
    r1 = 1000;
  }
  double r2 = theDetector().subDetectorByEnum(type).radiusLim;
  if (radius>r1) return 3;
  else if (radius>r2) return 2;
  else return 1;
  return 3;
}

SiWEcalSSGeometryConversion::~SiWEcalSSGeometryConversion(){
  std::map<DetectorEnum,std::vector<TH2D *> >::iterator liter =
    HistMapE_.begin();
  for (; liter !=HistMapE_.end();++liter){
    deleteHistos(liter->second);
  }
  HistMapE_.clear();
  if (!bypassRadius_) {
    liter =  HistMapESmall_.begin();
    for (; liter !=HistMapESmall_.end();++liter){
      deleteHistos(liter->second);
    }
    HistMapESmall_.clear();
  }
  liter = HistMapTime_.begin();
  for (; liter !=HistMapTime_.end();++liter){
    deleteHistos(liter->second);
  }
  HistMapTime_.clear();
  if (!bypassRadius_) {
    liter = HistMapTimeSmall_.begin();
    for (; liter !=HistMapTimeSmall_.end();++liter){
      deleteHistos(liter->second);
    }
    HistMapTimeSmall_.clear();
  }
  liter = HistMapZ_.begin();
  for (; liter !=HistMapZ_.end();++liter){
    deleteHistos(liter->second);
  }
  HistMapZ_.clear();
}


void SiWEcalSSGeometryConversion::deleteHistos(std::vector<TH2D *> & aVec){
  if (aVec.size()!=0){
    for (unsigned iL(0); iL<aVec.size();++iL){
      //if (aVec[iL]) aVec[iL]->Delete();
      delete aVec[iL]; aVec[iL]=0;
    }
    aVec.clear();
  }
}

void SiWEcalSSGeometryConversion::setGranularity(const std::vector<unsigned> & granul){
  granularity_.reserve(granul.size());
  for (unsigned iL(0); iL<granul.size();++iL){
    granularity_.push_back(granul[iL]);
  }
}

double SiWEcalSSGeometryConversion::cellSize(const unsigned aLayer, const double aR) const{
  if (theDetector().subDetectorByLayer(aLayer).isScint || model_ != 2) return cellSize_*granularity_[aLayer];
  double r1 = theDetector().subDetectorByLayer(aLayer).radiusLim;
  if (aR<r1) return cellSize_*3;
  else return cellSize_*4;
}

double SiWEcalSSGeometryConversion::cellSizeInCm(const unsigned aLayer, const double aR) const{
  return cellSize(aLayer, aR)/10.;
}

void SiWEcalSSGeometryConversion::initialiseHistos(const bool recreate,
					       std::string uniqStr,
					       const bool print){

   for (unsigned iS(0); iS<theDetector().nSections();++iS){
     resetVector(HistMapE_[theDetector().detType(iS)],"EmipHits"+uniqStr,theDetector().detName(iS),theDetector().subDetectorBySection(iS),theDetector().nLayers(iS),recreate,print);
     if (!bypassRadius_) resetVector(HistMapESmall_[theDetector().detType(iS)],"EmipHitsSmall"+uniqStr,theDetector().detName(iS),theDetector().subDetectorBySection(iS),theDetector().nLayers(iS),recreate,print);
     //std::cout << " check: " << HistMapE_[theDetector().detType(iS)].size() << std::endl;
     
     std::vector<double> avgvecE;
     avgvecE.resize(theDetector().nLayers(iS),0);
     avgMapE_[theDetector().detType(iS)]=avgvecE;
     
     resetVector(HistMapTime_[theDetector().detType(iS)],"TimeHits"+uniqStr,theDetector().detName(iS),theDetector().subDetectorBySection(iS),theDetector().nLayers(iS),recreate,print);
     if (!bypassRadius_) resetVector(HistMapTimeSmall_[theDetector().detType(iS)],"TimeHitsSmall"+uniqStr,theDetector().detName(iS),theDetector().subDetectorBySection(iS),theDetector().nLayers(iS),recreate,print);
     //std::cout << " check: " << HistMapTime_[theDetector().detType(iS)].size() << std::endl;

     resetVector(HistMapZ_[theDetector().detType(iS)],"zHits"+uniqStr,theDetector().detName(iS),theDetector().subDetectorBySection(iS),theDetector().nLayers(iS),recreate,print);
     //std::cout << " check: " << HistMapZ_[theDetector().detType(iS)].size() << std::endl;
     
     std::vector<double> avgvecZ;
     avgvecZ.resize(theDetector().nLayers(iS),0);
     avgMapZ_[theDetector().detType(iS)]=avgvecZ;
   }
 }

void SiWEcalSSGeometryConversion::fill(const DetectorEnum type,
				   const unsigned newlayer,
				   const double & weightedE,
				   const double & aTime,
				   const double & posx,
				   const double & posy,
				   const double & posz)
{
  double radius = sqrt(posx*posx+posy*posy);
  double r1 = theDetector().subDetectorByEnum(type).radiusLim;
  //patch for crack regions
  if (dopatch_ && model_==4){
    double xcrack1 = -160;
    double xcrack2 = 150;
    if (type==DetectorEnum::MECAL){
    xcrack1 += 310/3.;
    xcrack2 += 310/3.;
    }
    else if (type == DetectorEnum::BECAL){
      xcrack1 += 2*310/3.;
      xcrack2 += 2*310/3.;
    }

    if ( fabs(posx-xcrack1)<5 || fabs(posx-xcrack2)<5 ) {
      //std::cout << " -- skipping hit " << posx << ", det " << type << " layer " << newlayer << " cracks at " << xcrack1 << " " << xcrack2 << std::endl;
      return;
    }

  }
  HistMapE_[type][newlayer]->Fill(posx,posy,weightedE);
  //if (radius >= r1) HistMapE_[type][newlayer]->Fill(posx,posy,weightedE);
  //else 
  if (!bypassRadius_ && theDetector().subDetectorByEnum(type).isSi && radius<r1) {
    HistMapESmall_[type][newlayer]->Fill(posx,posy,weightedE);
    HistMapTimeSmall_[type][newlayer]->Fill(posx,posy,weightedE*aTime);
  }
  HistMapTime_[type][newlayer]->Fill(posx,posy,weightedE*aTime);
  HistMapZ_[type][newlayer]->Fill(posx,posy,weightedE*posz);
  avgMapZ_[type][newlayer] += weightedE*posz;
  avgMapE_[type][newlayer] += weightedE;
}

double SiWEcalSSGeometryConversion::getAverageZ(const unsigned layer){
  const SiWEcalSSSubDetector & subdet = theDetector().subDetectorByLayer(layer);
  unsigned newlayer = layer-subdet.layerIdMin;
  double avg = 0;
  if (avgMapE_[subdet.type][newlayer]>0)
    avg =avgMapZ_[subdet.type][newlayer]/avgMapE_[subdet.type][newlayer];
  return avg;
}

TH2D * SiWEcalSSGeometryConversion::get2DHist(const unsigned layer,std::string name){
  const SiWEcalSSSubDetector & subdet = theDetector().subDetectorByLayer(layer);
  unsigned newlayer = layer-subdet.layerIdMin;
  if (name == "E") return HistMapE_[subdet.type][newlayer];
  else if (name == "ESmall") {
    if (bypassRadius_) return 0;
    return HistMapESmall_[subdet.type][newlayer];
  }
  else if (name == "Time") return HistMapTime_[subdet.type][newlayer];
  else if (name == "TimeSmall"){
    if (bypassRadius_) return 0;
    return HistMapTimeSmall_[subdet.type][newlayer];
  }
  else if (name == "Z") return HistMapZ_[subdet.type][newlayer];
  else {
    std::cerr << " ERROR !! Unknown histogram name. Exiting..." << std::endl;
    exit(1);
  }
}

unsigned SiWEcalSSGeometryConversion::getGranularity(const unsigned aLayer, const SiWEcalSSSubDetector & adet){
  unsigned idx = adet.layerIdMin+aLayer;
  return granularity_[idx];
}

void SiWEcalSSGeometryConversion::resetVector(std::vector<TH2D *> & aVec,
					  std::string aVar,
					  std::string aString,
					  const SiWEcalSSSubDetector & aDet,
					  const unsigned nLayers,
					  bool recreate, 
					  bool print)
{
  //std::cout << " vector size: " << aVar << " " << aString << " = " << aVec.size() << std::endl;
  if (recreate){
    for (unsigned iL(0); iL<aVec.size();++iL){
      aVec[iL]->Delete();
    }
    aVec.clear();
  }
  if (aVec.size()!=0){
    for (unsigned iL(0); iL<aVec.size();++iL){
      aVec[iL]->Reset();
    }
  }
  else {
    if (nLayers > 0){
      aVec.resize(nLayers,0);
      if (print) std::cout << " -- Creating " << nLayers << " 2D histograms for " << aVar << " " << aString 
	//<< " with " << nBins << " bins between " << min << " and " << max 
		<< std::endl;
      for (unsigned iL(0); iL<nLayers;++iL){
	std::ostringstream lname;
	lname << aVar << "_" << aString << "_" << iL ;
	/*unsigned nBins = 0;
	  double xylim = 525;//750/sqrt(2) and multiple of 7.5
	if (aDet.type == DetectorEnum::FHCAL) {
	  xylim = 435;//600/sqrt(2) and multiple of 7.5
	}
	double xy = -1695;
	while(1){
	  if (xy>=(-1.*xylim) && xy < xylim){
	    xyedges[nBins] = xy;
	    //if (iL==0) std::cout << " xyedges[" << nBins << "]=" << xy << std::endl;
	    xy += 7.5;
	    nBins++;
	  }
	  else if (fabs(xy)<=1700 && fabs(yy)<=1700){
	    xyedges[nBins] = xy;
	    //if (iL==0) std::cout << " xyedges[" << nBins << "]=" << xy << std::endl;
	    xy += 10;
	    nBins++;
	  }
	  else break;
	}

	aVec[iL] = new TH2D(lname.str().c_str(),";x(mm);y(mm)",nBins-1,xyedges,nBins-1,xyedges);
	if (print && aVar == "EmipHits") {
	  std::cout << " ---- Layer " << iL << " bins, min, max = " << nBins-1 << " " << xyedges[0] << " " << xyedges[nBins-1] << std::endl;
	  //for (unsigned ib(0);ib<nBins;++ib){
	  //std::cout << xyedges[ib] << " ";
	  //}
	  //std::cout << std::endl;
	}
	*/
	// double newcellsize = cellSize_*getGranularity(iL,aDet);
	// //take smallest pair integer to be sure to fit in the geometry
	// //even if small dead area introduced at the edge
	// //take 0,0,0 as center of new cell.
	unsigned nBins = (aVar.find("Small")!=aVar.npos) ? 452 : 339;
	double min=-1695;
	double max=1695;
	double minx = min;
	double maxx = max;
	// if (getGranularity(iL,aDet) == 1) {
	//   nBins = static_cast<unsigned>(width_*1./cellSize_);
	//   min = -1.0*nBins*newcellsize/2.;
	//   max = nBins*newcellsize/2.;
	// }
	// else {
	if (aDet.isScint || bypassRadius_ || model_!=2){
	  double newcellsize = cellSize_*getGranularity(iL,aDet);
	  nBins = static_cast<unsigned>(width_*1./(newcellsize*2.))*2-2;
	  min = -1.0*nBins*newcellsize/2.-newcellsize/2.;
	  max = nBins*newcellsize/2.+newcellsize/2.;
	  nBins+=1;
	  minx = min;
	  maxx = max;

	  if (dopatch_ && model_==4){
	    if (aDet.type == DetectorEnum::MECAL) {
	      //set specific binning for displaced sections
	      minx += 310/3.;
	      maxx += 310/3.;
	    }
	    if (aDet.type == DetectorEnum::BECAL){
	      minx += 2*310/3.;
	      maxx += 2*310/3.;
	    }
	  }

	}
	aVec[iL] = new TH2D(lname.str().c_str(),";x(mm);y(mm)",nBins,minx,maxx,nBins,min,max);
	if (print && aVar.find("EmipHits")!=aVar.npos) std::cout << " ---- Layer " << iL << " bins, min, max = " << nBins << " " << min << " " << max << std::endl;
      }
    }
  }
  //std::cout << " vector size after: " << aVar << " " << aString << " = " << aVec.size() << std::endl;
}





