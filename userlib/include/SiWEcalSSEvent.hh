#ifndef _siwecalssevent_hh_
#define _siwecalssevent_hh_
#include <iomanip>
#include <vector>
#include "Rtypes.h"
#include <sstream>
#include <map>

class SiWEcalSSEvent{


public:
  SiWEcalSSEvent():
    event_(0)
  {
    
  };

  ~SiWEcalSSEvent(){};

  inline unsigned eventNumber() const{
    return event_;
  };

  inline void eventNumber(const unsigned aNum){
    event_ = aNum;
  };

  inline double vtx_x() const{
    return xvtx_;
  };

  inline void vtx_x(const double x){
    xvtx_ = x;
  };

  inline double vtx_y() const{
    return yvtx_;
  };

  inline void vtx_y(const double y){
    yvtx_ = y;
  };

  inline double vtx_z() const{
    return zvtx_;
  };

  inline void vtx_z(const double z){
    zvtx_ = z;
  };


private:

  unsigned event_;
  double xvtx_;
  double yvtx_;
  double zvtx_;

  ClassDef(SiWEcalSSEvent,1);



};


#endif
