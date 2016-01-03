#ifndef SiWEcalSSPUenergy_h
#define SiWEcalSSPUenergy_h

#include <string>
#include <fstream>
#include <vector>
#include "TF1.h"
#include "TMath.h"

class SiWEcalSSPUenergy{
  
public:
  SiWEcalSSPUenergy(){}; 
  SiWEcalSSPUenergy(const std::string filePath);
  ~SiWEcalSSPUenergy(); 
  double getDensity(const double & eta, const unsigned layer, const double & cellSize, const unsigned PU) const;
  
private:
  std::vector<double> p0_;
  std::vector<double> p1_;
};

#endif 
