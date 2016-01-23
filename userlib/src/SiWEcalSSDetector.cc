#include "SiWEcalSSDetector.hh"
#include <iostream>

SiWEcalSSDetector & theDetector(){
  static SiWEcalSSDetector lDet;
  static bool firstDet=true;
  if (firstDet) std::cout << " -- Created detector static object." << std::endl;
  firstDet=false;
  return lDet;
}

void SiWEcalSSDetector::buildDetector(const unsigned versionNumber,
				  bool concept,
				  bool isCaliceHcal,
				  bool bypassR){
  
  bypassRadius_ = bypassR;
  reset();
  initialiseIndices(versionNumber);
  SiWEcalSSSubDetector SiWEcal;
  SiWEcal.type = DetectorEnum::SiWEcal;
  SiWEcal.name = "SiWEcal";
  SiWEcal.layerIdMin = indices_[0];
  SiWEcal.layerIdMax = indices_[1];
  SiWEcal.mipWeight = 1./0.0600;//Fill the MIP here from study 
  SiWEcal.absWeight = 1.;//ratio of abs dedx
  SiWEcal.gevWeight = 1.0; //MIPToGeV 
  SiWEcal.gevOffset = 0.0;
  SiWEcal.isSi = true;
  //if (versionNumber>=30) 
  if (!bypassRadius_) SiWEcal.radiusLim = 750;
  else SiWEcal.radiusLim = 0;
  if (SiWEcal.nLayers()>0) theDetector().addSubdetector(SiWEcal);
  
  finishInitialisation();
  
}


const SiWEcalSSSubDetector & SiWEcalSSDetector::subDetectorByLayer(const unsigned aLayer){
  unsigned section = getSection(aLayer);
  return subdets_[section];
}

unsigned SiWEcalSSDetector::getSection(const unsigned aLayer) const{
  if (aLayer>=nLayers_) {
    std::cerr << " -- Error ! Trying to access layer " << aLayer 
	      << " outside of range. nLayers = " << nLayers_
	      << std::endl;
    exit(1);
  }
  return section_[aLayer];
}

void SiWEcalSSDetector::addSubdetector(const SiWEcalSSSubDetector & adet){
  subdets_.push_back(adet);
  enumMap_[adet.type]=subdets_.size()-1;
  //indices_.push_back(adet.layerIdMin);
}
  
void SiWEcalSSDetector::finishInitialisation(){
  nSections_ = subdets_.size();
  //indices_.push_back(subdets_[nSections_-1].layerIdMax);
  const unsigned lastEle = indices_.size()-1;
  nLayers_ = indices_[lastEle];

  unsigned secIndex[lastEle];
  unsigned iS(0);
  for(unsigned i(0); i<lastEle ;i++){
     secIndex[i] = iS;
     if(indices_[i] < indices_[i+1])iS +=1;
  }
  //initialise layer-section conversion
  section_.resize(nLayers_,0);
  for (unsigned iL(0); iL<nLayers_;++iL){
    for (unsigned i(0); i<lastEle;++i){
      if (iL >= indices_[i] && iL < indices_[i+1]) section_[iL] = secIndex[i];
    }
  }
  printDetector(std::cout);
}

const SiWEcalSSSubDetector & SiWEcalSSDetector::subDetectorByEnum(DetectorEnum adet){
  if (enumMap_.find(adet) == enumMap_.end()){
    std::cerr << " -- Error ! Trying to access subdetector enum not present in this detector: "
	      << adet 
	      << std::endl;
    exit(1);
  } 
  return subdets_[enumMap_[adet]];
}

void SiWEcalSSDetector::reset() {
  subdets_.clear();
  enumMap_.clear();
  indices_.clear();
  section_.clear();
}

void SiWEcalSSDetector::printDetector(std::ostream & aOs) const{
  std::cout << " -------------------------- " << std::endl
	    << " -- Detector information -- " << std::endl
	    << " -------------------------- " << std::endl
	    << " - nSections = " << nSections_ << std::endl
	    << " - nLayers = " << nLayers_ << std::endl
	    << " - detNames = " ;
  for (unsigned i(0); i<nSections_;++i){
    std::cout << " " << detName(i);
  }
  std::cout << std::endl;
  std::cout << " -------------------------- " << std::endl;
}
