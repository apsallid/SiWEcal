#ifndef _siwecalssinfo_hh_
#define _siwecalssinfo_hh_
#include <iomanip>
#include <vector>
#include "Rtypes.h"
#include <sstream>
#include <map>

class SiWEcalSSInfo{


public:
  SiWEcalSSInfo()
  {
    version_ = -1;
    model_ = -1;
    cellsize_ = 2.5;//mm
  };
  
  ~SiWEcalSSInfo(){};

  inline void version(const int & aVal){
    version_ = aVal;
  };
  
  inline int version() const{
    return version_;
  };

  inline void model(const int & aVal){
    model_ = aVal;
  };
  
  inline int model() const{
    return model_;
  };

  inline void cellSize(const double & aVal){
    cellsize_ = aVal;
  };
  
  inline double cellSize() const{
    return cellsize_;
  };

private:

  int version_;
  int model_;
  double cellsize_;


  ClassDef(SiWEcalSSInfo,1);



};


#endif
