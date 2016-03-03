#ifndef SiWEcalSSDetector_h
#define SiWEcalSSDetector_h

#include <string>
#include <vector>
#include <map>
#include "TH2D.h"

enum DetectorEnum {
  SiWEcal
};

class SiWEcalSSSubDetector {

public:
  SiWEcalSSSubDetector():
    type(DetectorEnum::SiWEcal),
    name(""),
    layerIdMin(0),
    layerIdMax(0),
    mipWeight(1),
    absWeight(1),
    gevWeight(1),
    gevOffset(0),
    isSi(false),
    isScint(false),
    radiusLim(0)
  {};
  ~SiWEcalSSSubDetector(){};

  DetectorEnum type;
  std::string name;
  unsigned layerIdMin;
  unsigned layerIdMax;
  double mipWeight;
  double absWeight;
  double gevWeight;
  double gevOffset;
  bool isSi;
  bool isScint;
  double radiusLim;

  inline unsigned nLayers() const{
    return (layerIdMax-layerIdMin);
  };

private:

};

class SiWEcalSSDetector {

public:
  friend SiWEcalSSDetector & theDetector();

  inline void initialiseIndices(const unsigned versionNumber){
    
    indices_.clear();
    indices_.resize(2,0);
    //fill layer indices
    if (versionNumber!=0){
      indices_[0] = 0;
      indices_[1] = 16;
    }


  };

  void buildDetector(const unsigned versionNumber,
		     bool bypassR=false);

  const SiWEcalSSSubDetector & subDetectorByLayer(const unsigned aLayer);

  unsigned getSection(const unsigned aLayer) const;
  inline unsigned section(const DetectorEnum adet){
    if (enumMap_.find(adet) != enumMap_.end())
      return enumMap_[adet];
    return nSections();
  };

  void addSubdetector(const SiWEcalSSSubDetector & adet);
  
  void finishInitialisation();

  inline unsigned nLayers(const unsigned aSection) const{
    return subdets_[aSection].nLayers();
  };

  inline unsigned nLayers(DetectorEnum adet){
    return subdets_[enumMap_[adet]].nLayers();
  };

  const SiWEcalSSSubDetector & subDetectorByEnum(DetectorEnum adet);
  inline const SiWEcalSSSubDetector & subDetectorBySection(const unsigned aSection) const{
    return subdets_[aSection];
  };

  inline unsigned nLayers() const{
    return nLayers_;
  };

  inline unsigned nSections() const{
    return nSections_;
  };

  inline DetectorEnum detType(const unsigned aSection) const{
    return subdets_[aSection].type;
  };
  
  inline DetectorEnum detTypeLayer(const unsigned aLayer) const{
    return subdets_[getSection(aLayer)].type;
  };
  
  inline std::string detName(const unsigned aSection) const{
    return subdets_[aSection].name;
  };


  void reset();

  void printDetector(std::ostream & aOs) const ;

private:
  SiWEcalSSDetector(){
    bypassRadius_ = false;
  };

  ~SiWEcalSSDetector(){
    reset();
  };
  
  std::vector<SiWEcalSSSubDetector> subdets_;
  std::vector<unsigned> indices_;
  std::vector<unsigned> section_;
  std::map<DetectorEnum,unsigned> enumMap_;

  unsigned nLayers_;
  unsigned nSections_;
  bool bypassRadius_;
};

SiWEcalSSDetector & theDetector();




#endif
