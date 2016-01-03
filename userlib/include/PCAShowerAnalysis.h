#ifndef PCAShowerAnalysis_h
#define PCAShowerAnalysis_h

#include "SiWEcalSSRecoHit.hh"
#include "SiWEcalSSCluster.hh"

#include "TPrincipal.h"
#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"
#include "Math/Point2D.h"
#include "Math/Point2Dfwd.h"

#include "TVector3.h"

class PCAShowerAnalysis
{

  public:

  PCAShowerAnalysis(bool segmented=true, bool logweighting=true, bool debug=false ) ;
  
  void showerParameters( const SiWEcalSSCluster & );

  ROOT::Math::XYZPoint showerBarycenter;
  ROOT::Math::XYZVector showerAxis;
  ROOT::Math::XYZVector showerEigenValues;
  ROOT::Math::XYZVector showerSigmas;

  ~PCAShowerAnalysis();
  
private:

  TPrincipal *principal_;
  
  double mip_;
  double entryz_;
  
  bool logweighting_;
  bool segmented_;
  
  bool alreadyfilled_;
  bool debug_;
  
};
#endif
