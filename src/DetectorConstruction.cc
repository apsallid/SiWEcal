#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include <boost/algorithm/string.hpp>

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4CSGSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4PhysicalConstants.hh"

using namespace std;

//
DetectorConstruction::DetectorConstruction(G4int ver, G4int mod) 
  : version_(ver), model_(mod), addPrePCB_(false)
{
  SetGapBetweenSensorPads(0.1*mm,0.1*mm);

  //radiation lengths: cf. http://pdg.lbl.gov/2012/AtomicNuclearProperties/
  //W 3.504 mm
  //Pb 5.612 mm
  //Cu 14.36 mm
  switch(version_)
    {
      //cf. http://arxiv.org/abs/0805.4833
    case v_CALICE:
      {
	G4cout << "[DetectorConstruction] starting v_CALICE (10x0.4+10x0.8+10x1.2)X_0 with Tungsten" << G4endl;
	std::vector<G4double> lThick;
	std::vector<std::string> lEle;
	lThick.push_back(0.4*3.504*mm);lEle.push_back("W");
	lThick.push_back(0.525*mm);lEle.push_back("Si");
	lThick.push_back(1.0*mm);lEle.push_back("PCB");
	lThick.push_back(2.5*mm);lEle.push_back("Air");

	for(unsigned i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	lThick[0] = 0.8*3.504*mm;
	for(unsigned i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	lThick[0] = 1.2*3.504*mm;
	for(unsigned i=0; i<10; i++) m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	break;
      }
    case v_SiWECAL_B0X0_I0X0_A0:
      {
	G4cout << "[DetectorConstruction] starting SiWECAL_B0X0W_I0X0_A0 0.0 X_0 with Tungsten" << G4endl;
	std::vector<G4double> lThick;
	std::vector<std::string> lEle;
	lThick.push_back(3.3*mm);lEle.push_back("Fe"); // This is for the front wall -- Check and correct 
	lThick.push_back(40.0*mm);lEle.push_back("Air"); // Gap between front wall and absorber
	lThick.push_back(0.0*mm);lEle.push_back("W");
	lThick.push_back(4.0*mm);lEle.push_back("Air"); // Gap between before absorber and plexiglass
	lThick.push_back(18.0*mm);lEle.push_back("Plexiglass"); // Plexiglass to press the connections of PCB to adaptor board
	
	//This is for the 16 SKIROC chips
	unsigned nChips = 16;
	for(unsigned i=0; i<nChips; i++) {
	  lThick.push_back(1.1*mm);lEle.push_back("PVC"); // Will correct this to be the material of SKIROC
	}
	lThick.push_back(1.6*mm);lEle.push_back("PCB");	
	lThick.push_back(0.2*mm);lEle.push_back("Air"); // this is for the glue dots
	lThick.push_back(0.325*mm);lEle.push_back("Si"); // For visualization increase the thickness of Si Sensors by a factor of 5
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	//plastic plate here 
	//kapton here
        lThick.push_back(25.0*mm);lEle.push_back("Al");
	lThick.push_back(23.0*mm);lEle.push_back("Air");

	unsigned nLay = 3;
	for(unsigned i=0; i<nLay; i++) {
	  if(i!=0){
	    //Front wall, gap, absorber, gap only in first layer
	    lThick[0] = 0.0*mm; 
	    lThick[1] = 0.0*mm; 
	    lThick[2] = 0.0*mm; 
	    lThick[3] = 0.0*mm; 
	  }
	  m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	}
	break;
      }
    case v_SiWECAL_B42X0W_I0X0_A0:
      {
	G4cout << "[DetectorConstruction] starting SiWECAL_B42X0W_I0X0_A0 4.2 X_0 with Tungsten" << G4endl;
	std::vector<G4double> lThick;
	std::vector<std::string> lEle;
	lThick.push_back(3.3*mm);lEle.push_back("Fe"); // This is for the front wall -- Check and correct 
	lThick.push_back(40.0*mm);lEle.push_back("Air"); // Gap between front wall and absorber
	lThick.push_back(4.2*3.504*mm);lEle.push_back("W");
	lThick.push_back(4.0*mm);lEle.push_back("Air"); // Gap between before absorber and plexiglass
	lThick.push_back(18.0*mm);lEle.push_back("Plexiglass"); // Plexiglass to press the connections of PCB to adaptor board
	
	//This is for the 16 SKIROC chips
	unsigned nChips = 16;
	for(unsigned i=0; i<nChips; i++) {
	  lThick.push_back(1.1*mm);lEle.push_back("PVC"); // Will correct this to be the material of SKIROC
	}
	lThick.push_back(1.6*mm);lEle.push_back("PCB");	
	lThick.push_back(0.2*mm);lEle.push_back("Air"); // this is for the glue dots
	lThick.push_back(0.325*mm);lEle.push_back("Si"); // For visualization increase the thickness of Si Sensors by a factor of 5
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	//plastic plate here 
	//kapton here
        lThick.push_back(25.0*mm);lEle.push_back("Al");
	lThick.push_back(23.0*mm);lEle.push_back("Air");

	unsigned nLay = 3;
	for(unsigned i=0; i<nLay; i++) {
	  if(i!=0){
	    //Front wall, gap, absorber, gap only in first layer
	    lThick[0] = 0.0*mm; 
	    lThick[1] = 0.0*mm; 
	    lThick[2] = 0.0*mm; 
	    lThick[3] = 0.0*mm; 
	  }
	  m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	}
	break;
      }
    case v_SiWECAL_B84X0W_I0X0_A0:
      {
	G4cout << "[DetectorConstruction] starting SiWECAL_B84X0W_I0X0_A0 8.4 X_0 with Tungsten" << G4endl;
	std::vector<G4double> lThick;
	std::vector<std::string> lEle;
	//Before calo here
	lThick.push_back(3.3*mm);lEle.push_back("Fe"); // This is for the front wall -- Check and correct 
	lThick.push_back(40.0*mm);lEle.push_back("Air"); // Gap between front wall and absorber
	m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	lThick.clear();
	lEle.clear();
	//Starting calo first layer here
	lThick.push_back(8.4*3.504*mm);lEle.push_back("W");
	lThick.push_back(4.0*mm);lEle.push_back("Air"); // Gap between before absorber and plexiglass
	lThick.push_back(18.0*mm);lEle.push_back("Plexiglass"); // Plexiglass to press the connections of PCB to adaptor board
	
	//This is for the 16 SKIROC chips
	unsigned nChips = 16;
	for(unsigned i=0; i<nChips; i++) {
	  lThick.push_back(1.1*mm);lEle.push_back("PVC"); // Will correct this to be the material of SKIROC
	}
	lThick.push_back(1.6*mm);lEle.push_back("PCB");	
	lThick.push_back(0.2*mm);lEle.push_back("Air"); // this is for the glue dots
	lThick.push_back(0.325*mm);lEle.push_back("Si"); // For visualization increase the thickness of Si Sensors by a factor of 5
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	lThick.clear();
	lEle.clear();
	
	//Next 3 layers here - 4th layer that wasn't read out
	//plastic plate here 
	//kapton here
        lThick.push_back(25.0*mm);lEle.push_back("Al");
	lThick.push_back(23.0*mm);lEle.push_back("Air");
	lThick.push_back(18.0*mm);lEle.push_back("Plexiglass"); // Plexiglass to press the connections of PCB to adaptor board
	
	//This is for the 16 SKIROC chips
	for(unsigned i=0; i<nChips; i++) {
	  lThick.push_back(1.1*mm);lEle.push_back("PVC"); // Will correct this to be the material of SKIROC
	}
	lThick.push_back(1.6*mm);lEle.push_back("PCB");	
	lThick.push_back(0.2*mm);lEle.push_back("Air"); // this is for the glue dots
	lThick.push_back(0.325*mm);lEle.push_back("Si"); // For visualization increase the thickness of Si Sensors by a factor of 5
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	lThick.push_back(0.325*mm);lEle.push_back("Si");

	//Putting also the 4th layer that wasn't read out just in case of back scattering
	//That is 4 - 1 that already put = 3 layers here
	unsigned nLay = 3;
	for(unsigned i=0; i<nLay; i++) {
	  m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	}
	//Will put an infinite volume here to collect the losses/leakage behind the end of the detector
	lThick.clear();
	lEle.clear();
	//plastic plate here 
	//kapton here
        lThick.push_back(25.0*mm);lEle.push_back("Al");
	lThick.push_back(23.0*mm);lEle.push_back("Air");
	lThick.push_back(1000.*m);lEle.push_back("SSteel");

	m_caloStruct.push_back( SamplingSection(lThick,lEle) );


	break;
      }
    case v_SiWECAL_B84X0W200Fe_I0X0_A0:
      {
	G4cout << "[DetectorConstruction] starting SiWECAL_B84X0W200Fe_I0X0_A0 8.4 X_0 with Tungsten and 200 mm of Fe" << G4endl;
	std::vector<G4double> lThick;
	std::vector<std::string> lEle;
	lThick.push_back(3.3*mm);lEle.push_back("Fe"); // This is for the front wall -- Check and correct 
	lThick.push_back(40.0*mm);lEle.push_back("Air"); // Gap between front wall and absorber
	lThick.push_back(8.4*3.504*mm);lEle.push_back("W");
	lThick.push_back(200.0*mm);lEle.push_back("Fe");
	lThick.push_back(4.0*mm);lEle.push_back("Air"); // Gap between before absorber and plexiglass
	lThick.push_back(18.0*mm);lEle.push_back("Plexiglass"); // Plexiglass to press the connections of PCB to adaptor board
	
	//This is for the 16 SKIROC chips
	unsigned nChips = 16;
	for(unsigned i=0; i<nChips; i++) {
	  lThick.push_back(1.1*mm);lEle.push_back("PVC"); // Will correct this to be the material of SKIROC
	}
	lThick.push_back(1.6*mm);lEle.push_back("PCB");	
	lThick.push_back(0.2*mm);lEle.push_back("Air"); // this is for the glue dots
	lThick.push_back(0.325*mm);lEle.push_back("Si"); // For visualization increase the thickness of Si Sensors by a factor of 5
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	//plastic plate here 
	//kapton here
        lThick.push_back(25.0*mm);lEle.push_back("Al");
	lThick.push_back(23.0*mm);lEle.push_back("Air");
	
	unsigned nLay = 3;
	for(unsigned i=0; i<nLay; i++) {
	  if(i!=0){
	    //Front wall, gap, absorber, gap only in first layer
	    lThick[0] = 0.0*mm; 
	    lThick[1] = 0.0*mm; 
	    lThick[2] = 0.0*mm; 
	    lThick[3] = 0.0*mm; 
	    lThick[4] = 0.0*mm; 
	  }
	  m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	}

	break;
      }
    case v_SiWECAL_B84X0W300Fe_I0X0_A0:
      {
	G4cout << "[DetectorConstruction] starting SiWECAL_B84X0W300Fe_I0X0_A0 8.4 X_0 with Tungsten and 300 mm of Fe" << G4endl;
	std::vector<G4double> lThick;
	std::vector<std::string> lEle;
	lThick.push_back(3.3*mm);lEle.push_back("Fe"); // This is for the front wall -- Check and correct 
	lThick.push_back(40.0*mm);lEle.push_back("Air"); // Gap between front wall and absorber
	lThick.push_back(8.4*3.504*mm);lEle.push_back("W");
	lThick.push_back(300.0*mm);lEle.push_back("Fe");
	lThick.push_back(4.0*mm);lEle.push_back("Air"); // Gap between before absorber and plexiglass
	lThick.push_back(18.0*mm);lEle.push_back("Plexiglass"); // Plexiglass to press the connections of PCB to adaptor board
	
	//This is for the 16 SKIROC chips
	unsigned nChips = 16;
	for(unsigned i=0; i<nChips; i++) {
	  lThick.push_back(1.1*mm);lEle.push_back("PVC"); // Will correct this to be the material of SKIROC
	}
	lThick.push_back(1.6*mm);lEle.push_back("PCB");	
	lThick.push_back(0.2*mm);lEle.push_back("Air"); // this is for the glue dots
	lThick.push_back(0.325*mm);lEle.push_back("Si"); // For visualization increase the thickness of Si Sensors by a factor of 5
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	//plastic plate here 
	//kapton here
        lThick.push_back(25.0*mm);lEle.push_back("Al");
	lThick.push_back(23.0*mm);lEle.push_back("Air");

	unsigned nLay = 3;
	for(unsigned i=0; i<nLay; i++) {
	  if(i!=0){
	    //Front wall, gap, absorber, gap only in first layer
	    lThick[0] = 0.0*mm; 
	    lThick[1] = 0.0*mm; 
	    lThick[2] = 0.0*mm; 
	    lThick[3] = 0.0*mm; 
	    lThick[4] = 0.0*mm; 
	  }
	  m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	}
	break;
      }
    case v_SiWECAL_B0X0_I42X0W_A0:
      {
	G4cout << "[DetectorConstruction] starting SiWECAL_B0X0W_I42X0_A0 4.2 X_0 with Tungsten between each plate" << G4endl;
	std::vector<G4double> lThick;
	std::vector<std::string> lEle;
	lThick.push_back(3.3*mm);lEle.push_back("Fe"); // This is for the front wall -- Check and correct 
	lThick.push_back(40.0*mm);lEle.push_back("Air"); // Gap between front wall and absorber
	lThick.push_back(0.0*mm);lEle.push_back("W");
	lThick.push_back(4.0*mm);lEle.push_back("Air"); // Gap between before absorber and plexiglass
	lThick.push_back(18.0*mm);lEle.push_back("Plexiglass"); // Plexiglass to press the connections of PCB to adaptor board
	
	//This is for the 16 SKIROC chips
	unsigned nChips = 16;
	for(unsigned i=0; i<nChips; i++) {
	  lThick.push_back(1.1*mm);lEle.push_back("PVC"); // Will correct this to be the material of SKIROC
	}
	lThick.push_back(1.6*mm);lEle.push_back("PCB");	
	lThick.push_back(0.2*mm);lEle.push_back("Air"); // this is for the glue dots
	lThick.push_back(0.325*mm);lEle.push_back("Si"); // For visualization increase the thickness of Si Sensors by a factor of 5
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	//plastic plate here 
	//kapton here
        lThick.push_back(25.0*mm);lEle.push_back("Al");
	lThick.push_back(15.0*mm);lEle.push_back("Air");

	unsigned nLay = 3;
	for(unsigned i=0; i<nLay; i++) {
	  if(i!=0){
	    //Front wall, gap, absorber, gap only in first layer
	    lThick[0] = 0.0*mm; 
	    lThick[1] = 0.0*mm; 
	    lThick[2] = 4.2*3.504*mm; // absorber between plates 
	  }
	  m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	}
	break;
      }
    case v_SiWECAL_B18X0W_I42X0W_A0:
      {
	G4cout << "[DetectorConstruction] starting SiWECAL_B18X0W_I42X0_A0 1.8 X_0 and 4.2 X_0 with Tungsten before and between each plate" << G4endl;
	std::vector<G4double> lThick;
	std::vector<std::string> lEle;
	lThick.push_back(3.3*mm);lEle.push_back("Fe"); // This is for the front wall -- Check and correct 
	lThick.push_back(40.0*mm);lEle.push_back("Air"); // Gap between front wall and absorber
	lThick.push_back(1.8*3.504*mm);lEle.push_back("W");
	lThick.push_back(4.0*mm);lEle.push_back("Air"); // Gap between before absorber and plexiglass
	lThick.push_back(18.0*mm);lEle.push_back("Plexiglass"); // Plexiglass to press the connections of PCB to adaptor board
	
	//This is for the 16 SKIROC chips
	unsigned nChips = 16;
	for(unsigned i=0; i<nChips; i++) {
	  lThick.push_back(1.1*mm);lEle.push_back("PVC"); // Will correct this to be the material of SKIROC
	}
	lThick.push_back(1.6*mm);lEle.push_back("PCB");	
	lThick.push_back(0.2*mm);lEle.push_back("Air"); // this is for the glue dots
	lThick.push_back(0.325*mm);lEle.push_back("Si"); // For visualization increase the thickness of Si Sensors by a factor of 5
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	//plastic plate here 
	//kapton here
        lThick.push_back(25.0*mm);lEle.push_back("Al");
	lThick.push_back(15.0*mm);lEle.push_back("Air");

	unsigned nLay = 3;
	for(unsigned i=0; i<nLay; i++) {
	  if(i!=0){
	    //Front wall, gap, absorber, gap only in first layer
	    lThick[0] = 0.0*mm; 
	    lThick[1] = 0.0*mm; 
	    lThick[2] = 4.2*3.504*mm; // absorber between plates 
	  }
	  m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	}
	break;
      }
    case v_SiWECAL_B24X0W_I42X0W_A0:
      {
	G4cout << "[DetectorConstruction] starting SiWECAL_B24X0W_I42X0_A0 2.4 X_0 and 4.2 X_0 with Tungsten before and between each plate" << G4endl;
	std::vector<G4double> lThick;
	std::vector<std::string> lEle;
	lThick.push_back(3.3*mm);lEle.push_back("Fe"); // This is for the front wall -- Check and correct 
	lThick.push_back(40.0*mm);lEle.push_back("Air"); // Gap between front wall and absorber
	lThick.push_back(2.4*3.504*mm);lEle.push_back("W");
	lThick.push_back(4.0*mm);lEle.push_back("Air"); // Gap between before absorber and plexiglass
	lThick.push_back(18.0*mm);lEle.push_back("Plexiglass"); // Plexiglass to press the connections of PCB to adaptor board
	
	//This is for the 16 SKIROC chips
	unsigned nChips = 16;
	for(unsigned i=0; i<nChips; i++) {
	  lThick.push_back(1.1*mm);lEle.push_back("PVC"); // Will correct this to be the material of SKIROC
	}
	lThick.push_back(1.6*mm);lEle.push_back("PCB");	
	lThick.push_back(0.2*mm);lEle.push_back("Air"); // this is for the glue dots
	lThick.push_back(0.325*mm);lEle.push_back("Si"); // For visualization increase the thickness of Si Sensors by a factor of 5
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	//plastic plate here 
	//kapton here
        lThick.push_back(25.0*mm);lEle.push_back("Al");
	lThick.push_back(15.0*mm);lEle.push_back("Air");

	unsigned nLay = 3;
	for(unsigned i=0; i<nLay; i++) {
	  if(i!=0){
	    //Front wall, gap, absorber, gap only in first layer
	    lThick[0] = 0.0*mm; 
	    lThick[1] = 0.0*mm; 
	    lThick[2] = 4.2*3.504*mm; // absorber between plates 
	  }
	  m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	}
	break;
      }
    case v_SiWECAL_B42X0W_I42X0W_A0:
      {
	G4cout << "[DetectorConstruction] starting SiWECAL_B42X0W_I42X0_A0 4.2 X_0 and 4.2 X_0 with Tungsten before and between each plate" << G4endl;
	std::vector<G4double> lThick;
	std::vector<std::string> lEle;
	lThick.push_back(3.3*mm);lEle.push_back("Fe"); // This is for the front wall -- Check and correct 
	lThick.push_back(40.0*mm);lEle.push_back("Air"); // Gap between front wall and absorber
	lThick.push_back(4.2*3.504*mm);lEle.push_back("W");
	lThick.push_back(4.0*mm);lEle.push_back("Air"); // Gap between before absorber and plexiglass
	lThick.push_back(18.0*mm);lEle.push_back("Plexiglass"); // Plexiglass to press the connections of PCB to adaptor board
	
	//This is for the 16 SKIROC chips
	unsigned nChips = 16;
	for(unsigned i=0; i<nChips; i++) {
	  lThick.push_back(1.1*mm);lEle.push_back("PVC"); // Will correct this to be the material of SKIROC
	}
	lThick.push_back(1.6*mm);lEle.push_back("PCB");	
	lThick.push_back(0.2*mm);lEle.push_back("Air"); // this is for the glue dots
	lThick.push_back(0.325*mm);lEle.push_back("Si"); // For visualization increase the thickness of Si Sensors by a factor of 5
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	//plastic plate here 
	//kapton here
        lThick.push_back(25.0*mm);lEle.push_back("Al");
	lThick.push_back(15.0*mm);lEle.push_back("Air");

	unsigned nLay = 3;
	for(unsigned i=0; i<nLay; i++) {
	  if(i!=0){
	    //Front wall, gap, absorber, gap only in first layer
	    lThick[0] = 0.0*mm; 
	    lThick[1] = 0.0*mm; 
	    lThick[2] = 4.2*3.504*mm; // absorber between plates 
	  }
	  m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	}
	break;
      }
    case v_SiWECAL_B24X0W_I42X0W_A48:
      {
	G4cout << "[DetectorConstruction] starting SiWECAL_B24X0W_I42X0_A0 2.4 X_0 and 4.2 X_0 with Tungsten before and between each plate and 48 degree angle" << G4endl;
	std::vector<G4double> lThick;
	std::vector<std::string> lEle;
	lThick.push_back(3.3*mm);lEle.push_back("Fe"); // This is for the front wall -- Check and correct 
	lThick.push_back(40.0*mm);lEle.push_back("Air"); // Gap between front wall and absorber
	lThick.push_back(2.4*3.504*mm);lEle.push_back("W");
	lThick.push_back(4.0*mm);lEle.push_back("Air"); // Gap between before absorber and plexiglass
	lThick.push_back(18.0*mm);lEle.push_back("Plexiglass"); // Plexiglass to press the connections of PCB to adaptor board
	
	//This is for the 16 SKIROC chips
	unsigned nChips = 16;
	for(unsigned i=0; i<nChips; i++) {
	  lThick.push_back(1.1*mm);lEle.push_back("PVC"); // Will correct this to be the material of SKIROC
	}
	lThick.push_back(1.6*mm);lEle.push_back("PCB");	
	lThick.push_back(0.2*mm);lEle.push_back("Air"); // this is for the glue dots
	lThick.push_back(0.325*mm);lEle.push_back("Si"); // For visualization increase the thickness of Si Sensors by a factor of 5
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	//plastic plate here 
	//kapton here
        lThick.push_back(25.0*mm);lEle.push_back("Al");
	lThick.push_back(15.0*mm);lEle.push_back("Air");

	unsigned nLay = 3;
	for(unsigned i=0; i<nLay; i++) {
	  if(i!=0){
	    //Front wall, gap, absorber, gap only in first layer
	    lThick[0] = 0.0*mm; 
	    lThick[1] = 0.0*mm; 
	    lThick[2] = 4.2*3.504*mm; // absorber between plates 
	  }
	  m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	}
	break;
      }
    case v_SiWECAL_B0X0_I0X0_A90:
      {
	G4cout << "[DetectorConstruction] starting SiWECAL_B0X0W_I0X0_A0 0.0 X_0 with Tungsten with 90 degree angle" << G4endl;
	std::vector<G4double> lThick;
	std::vector<std::string> lEle;
	lThick.push_back(3.3*mm);lEle.push_back("Fe"); // This is for the front wall -- Check and correct 
	lThick.push_back(40.0*mm);lEle.push_back("Air"); // Gap between front wall and absorber
	lThick.push_back(0.0*mm);lEle.push_back("W");
	lThick.push_back(4.0*mm);lEle.push_back("Air"); // Gap between before absorber and plexiglass
	lThick.push_back(18.0*mm);lEle.push_back("Plexiglass"); // Plexiglass to press the connections of PCB to adaptor board
	
	//This is for the 16 SKIROC chips
	unsigned nChips = 16;
	for(unsigned i=0; i<nChips; i++) {
	  lThick.push_back(1.1*mm);lEle.push_back("PVC"); // Will correct this to be the material of SKIROC
	}
	lThick.push_back(1.6*mm);lEle.push_back("PCB");	
	lThick.push_back(0.2*mm);lEle.push_back("Air"); // this is for the glue dots
	lThick.push_back(0.325*mm);lEle.push_back("Si"); // For visualization increase the thickness of Si Sensors by a factor of 5
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	lThick.push_back(0.325*mm);lEle.push_back("Si");
	//plastic plate here 
	//kapton here
        lThick.push_back(25.0*mm);lEle.push_back("Al");
	lThick.push_back(23.0*mm);lEle.push_back("Air");

	unsigned nLay = 3;
	for(unsigned i=0; i<nLay; i++) {
	  if(i!=0){
	    //Front wall, gap, absorber, gap only in first layer
	    lThick[0] = 0.0*mm; 
	    lThick[1] = 0.0*mm; 
	    lThick[2] = 0.0*mm; 
	    lThick[3] = 0.0*mm; 
	  }
	  m_caloStruct.push_back( SamplingSection(lThick,lEle) );
	}
	break;
      }

      
    }

  DefineMaterials();
  SetMagField(0);
  m_detectorMessenger = new DetectorMessenger(this);
  UpdateCalorSize();
}

DetectorConstruction::~DetectorConstruction() { delete m_detectorMessenger;}

void DetectorConstruction::DefineMaterials()
{ 
  G4NistManager* nistManager = G4NistManager::Instance();
  m_materials["Abs"] = nistManager->FindOrBuildMaterial("G4_W",false); 
  m_materials["Al"] = nistManager->FindOrBuildMaterial("G4_Al",false);
  m_dEdx["Al"] = 0.4358;
  m_materials["W"] = nistManager->FindOrBuildMaterial("G4_W",false); 
  m_dEdx["W"] = 2.210;
  m_materials["Pb"] = nistManager->FindOrBuildMaterial("G4_Pb",false); 
  m_dEdx["Pb"] = 1.274;
  m_materials["Cu"] = nistManager->FindOrBuildMaterial("G4_Cu",false); 
  m_dEdx["Cu"] = 1.257;
  m_materials["Si"] = nistManager->FindOrBuildMaterial("G4_Si",false);
  m_dEdx["Si"] = 0.3876;
  m_materials["Zn"] = nistManager->FindOrBuildMaterial("G4_Zn",false);
  m_dEdx["Zn"] = 1.007;
  m_materials["Air"]=nistManager->FindOrBuildMaterial("G4_AIR",false);
  m_dEdx["Air"] = 0;
  m_materials["Fe"] = nistManager->FindOrBuildMaterial("G4_Fe",false);
  m_dEdx["Fe"] = 1.143;
  m_materials["Mn"] = nistManager->FindOrBuildMaterial("G4_Mn",false);
  m_dEdx["Mn"] = 1.062 ;
  m_materials["C"] = nistManager->FindOrBuildMaterial("G4_C",false); 
  m_dEdx["C"] = 0.3952;
  m_materials["H"] = nistManager->FindOrBuildMaterial("G4_H",false); 
  m_dEdx["H"] =  0;
  m_materials["Cl"] = nistManager->FindOrBuildMaterial("G4_Cl",false); 
  m_dEdx["Cl"] = 0;
  m_materials["Cr"] = nistManager->FindOrBuildMaterial("G4_Cr",false); 
  m_dEdx["Cr"] = 1.046;
  m_materials["Ni"] = nistManager->FindOrBuildMaterial("G4_Ni",false); 
  m_dEdx["Ni"] = 1.307;
  m_materials["Plexiglass"] = nistManager->FindOrBuildMaterial("G4_PLEXIGLASS",false); 
  m_dEdx["Plexiglass"] = 0.6*m_dEdx["C"];

  m_materials["PCB"] = new G4Material("G10",1.700*g/cm3,4);
  m_materials["PCB"]->AddElement(nistManager->FindOrBuildElement(14), 1);
  m_materials["PCB"]->AddElement(nistManager->FindOrBuildElement(8) , 2);
  m_materials["PCB"]->AddElement(nistManager->FindOrBuildElement(6) , 3);
  m_materials["PCB"]->AddElement(nistManager->FindOrBuildElement(1) , 3);
  m_dEdx["PCB"] = 0;
  m_materials["Brass"]= new G4Material("Brass",8.5*g/cm3,2);
  m_materials["Brass"]->AddMaterial(m_materials["Cu"]  , 70*perCent);
  m_materials["Brass"]->AddMaterial(m_materials["Zn"]  , 30*perCent);
  m_dEdx["Brass"] = 0.7*m_dEdx["Cu"]+0.3*m_dEdx["Zn"];
  m_materials["Steel"]= new G4Material("Steel",7.87*g/cm3,3);
  m_materials["Steel"]->AddMaterial(m_materials["Fe"]  , 0.9843);
  m_materials["Steel"]->AddMaterial(m_materials["Mn"], 0.014);
  m_materials["Steel"]->AddMaterial(m_materials["C"], 0.0017);
  m_dEdx["Steel"] = 0.9843*m_dEdx["Fe"]+0.014*m_dEdx["Mn"]+0.0017*m_dEdx["C"];
  m_materials["SSteel"]= new G4Material("SSteel",8.02*g/cm3,4);
  m_materials["SSteel"]->AddMaterial(m_materials["Fe"]  , 0.70);
  m_materials["SSteel"]->AddMaterial(m_materials["Mn"], 0.01);
  m_materials["SSteel"]->AddMaterial(m_materials["Cr"], 0.19);
  m_materials["SSteel"]->AddMaterial(m_materials["Ni"], 0.10);
  m_dEdx["SSteel"] = 0.7*m_dEdx["Fe"]+0.01*m_dEdx["Mn"]+0.19*m_dEdx["Cr"]+0.1*m_dEdx["Ni"];
  m_materials["Scintillator"]= nistManager->FindOrBuildMaterial("G4_POLYSTYRENE",false); 
  m_dEdx["Scintillator"] = m_dEdx["C"];
  //m_materials["Scintillator"]= new G4Material("Scintillator",1.032*g/cm3,2);
  //m_materials["Scintillator"]->AddMaterial(m_materials["C"]  , 91.512109*perCent);
  //m_materials["Scintillator"]->AddMaterial(m_materials["H"]  , 8.4878906*perCent);
  G4cout << m_materials["Scintillator"] << G4endl;
  m_materials["Polystyrole"]= new G4Material("Polystyrole",1.065*g/cm3,2);
  m_materials["Polystyrole"]->AddMaterial(m_materials["H"]  , 50*perCent);
  m_materials["Polystyrole"]->AddMaterial(m_materials["C"]  , 50*perCent);
  m_dEdx["Polystyrole"] = 0.5*m_dEdx["C"];

  m_materials["PVC"]= new G4Material("PVC",1.350*g/cm3,3);
  m_materials["PVC"]->AddMaterial(m_materials["H"]  , 50*perCent);
  m_materials["PVC"]->AddMaterial(m_materials["C"]  , 33.33*perCent);
  m_materials["PVC"]->AddMaterial(m_materials["Cl"]  , 16.67*perCent);
  m_dEdx["PVC"] = 0.33*m_dEdx["C"];

  m_materials["CFMix"]= new G4Material("CFMix",0.120*g/cm3,3);
  m_materials["CFMix"]->AddMaterial(m_materials["Air"]  , 0.009);
  m_materials["CFMix"]->AddMaterial(m_materials["PVC"]  , 0.872);
  m_materials["CFMix"]->AddMaterial(m_materials["Polystyrole"]  , 0.119);
  m_dEdx["CFMix"] = 0;

  m_materials["Foam"]= new G4Material("Foam",0.0999*g/cm3,2);
  m_materials["Foam"]->AddMaterial(m_materials["C"]  , 0.856);
  m_materials["Foam"]->AddMaterial(m_materials["H"]  , 0.144);
  m_dEdx["Foam"] = 1.749*0.856*0.0999/10.;

  m_materials["WCu"]= new G4Material("WCu",14.979*g/cm3,2);
  m_materials["WCu"]->AddMaterial(m_materials["W"]  , 75*perCent);
  m_materials["WCu"]->AddMaterial(m_materials["Cu"]  , 25*perCent);
  m_dEdx["WCu"] = 0.75*m_dEdx["W"]+0.25*m_dEdx["Cu"];

  m_materials["NeutMod"]= new G4Material("NeutMod",0.950*g/cm3,2);
  m_materials["NeutMod"]->AddMaterial(m_materials["C"]  , 0.85628);
  m_materials["NeutMod"]->AddMaterial(m_materials["H"]  , 0.14372);
  m_dEdx["NeutMod"] = 1.749*0.86*0.950/10.;

}

//
void DetectorConstruction::UpdateCalorSize(){  

  m_CalorSizeZ=0;
  
  for(size_t i=0; i<m_caloStruct.size(); i++){
    m_CalorSizeZ=m_CalorSizeZ+m_caloStruct[i].Total_thick;
  }
 
  m_nSectors = 1;
  for(size_t i=0; i<m_caloStruct.size(); i++) m_caloStruct[i].setNumberOfSectors(m_nSectors);

  m_CalorSizeXY=500;
  
  m_WorldSizeZ=m_CalorSizeZ*3.1;  
  m_WorldSizeXY=m_CalorSizeXY*3.1;

  G4cout << "[DetectorConstruction][UpdateCalorSize] Z x XY = " 
	      << m_CalorSizeZ << " x " 
	      << m_CalorSizeXY << " mm " 
	      <<  G4endl;

}

//
G4VPhysicalVolume* DetectorConstruction::Construct()
{

  //clean old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  //world
  G4double expHall_z = 40*m;
  G4double expHall_x = 3*m;
  G4double expHall_y = 3*m;

  G4Box* experimentalHall_box = new G4Box("expHall_box",expHall_x,expHall_y,expHall_z);

  G4LogicalVolume* experimentalHall_log = new G4LogicalVolume(experimentalHall_box, m_materials["Air"],"expHall_log");
  G4VPhysicalVolume* experimentalHall_phys
    = new G4PVPlacement(0,                       // no rotation
			G4ThreeVector(0.,0.,0.), // translation position
			experimentalHall_log,    // its logical volume
			"expHall",               // its name
			0,                       // its mother volume
			false,                   // no boolean operations
			0);                      // its copy number

  //detector's World
  G4double pos_x = 0.;
  G4double pos_y = 0.;
  G4double pos_z = 0.;

  m_solidWorld = new G4Box("Wbox",m_WorldSizeXY/2,m_WorldSizeXY/2,m_WorldSizeZ/2);
  m_logicWorld = new G4LogicalVolume(m_solidWorld, m_materials["Air"], "Wlog");
  m_physWorld = new G4PVPlacement(0, G4ThreeVector(pos_x,pos_y,pos_z), m_logicWorld, "Wphys", experimentalHall_log, false, 0);

  for (unsigned iS(0); iS<m_nSectors; ++iS){
    buildStack(iS);
  }

  // Visualization attributes
  //
  m_logicWorld->SetVisAttributes(G4VisAttributes::Invisible);

  //return m_physWorld;
  return experimentalHall_phys;
}

void DetectorConstruction::buildStack(const unsigned sectorNum)
{


  //build the stack
  G4double zOffset(0.); //m_CalorSizeZ/2.
  G4double zOverburden(0.);

  char nameBuf[10];
  G4CSGSolid *solid;

  G4double totalLengthX0 = 0;
  G4double totalLengthL0 = 0;

  for(size_t i=0; i<m_caloStruct.size(); i++)
    {
      const unsigned nEle = m_caloStruct[i].n_elements;
      //index for counting Si sensitive layers 
      unsigned idx = 0;
      //index for counting SKIROC 
      unsigned chip_index = 0;     

      for (unsigned ie(0); ie<nEle;++ie){
	std::string eleName = m_caloStruct[i].ele_name[ie];
	
	sprintf(nameBuf,"%s%d",eleName.c_str(),int(i+1));
	//bool to split the Si sensor to 4 pads
	bool sivol = false;
	if (eleName=="Si") {
	  sprintf(nameBuf,"Si%d_%d",int(i+1),idx);
	  idx++;
	  sivol = true;
	}
	//for the chips geometry
	bool skirocvol = false;
	if (eleName=="PVC") {
	  sprintf(nameBuf,"SKIROC%d_%d",int(i+1),chip_index);
	  chip_index++;
	  skirocvol = true;
	}

	std::string baseName(nameBuf);
	G4double thick = m_caloStruct[i].ele_thick[ie];
	//One sensor has a size of 89.960 +/- 0.040 mm, and it is a matrix of
	//16x16 pixels with a pitch of 5.53 mm. So, 16 * 5.53 = 88.48 mm 
	//Since width is used for every volume we set with to 
	//2*88.48 mm + 2 * gap between sensors in center  = 177.16 mm 
	G4double width = 177.16*mm; //Hardcoded should fix these
	G4double height = 17.0*mm; //For SKIROC
	if (skirocvol){width = 17.0*mm;}

	if (eleName=="SSteel") {width = 1000.*m; height = 1000.*m;}

	if(thick>0){
	  solid = constructSolid(baseName,width,height,thick,sivol,skirocvol);
	  G4LogicalVolume *logi = new G4LogicalVolume(solid, m_materials[eleName], baseName+"log");
 	  m_caloStruct[i].ele_X0[ie] = m_materials[eleName]->GetRadlen();
	  m_caloStruct[i].ele_dEdx[ie] = m_dEdx[eleName];
	  m_caloStruct[i].ele_L0[ie] = m_materials[eleName]->GetNuclearInterLength();
	  G4cout << "************ " << eleName ;
	  G4cout << " layer " << i << " dEdx = " << m_caloStruct[i].ele_dEdx[ie] << " X0= " << m_caloStruct[i].ele_X0[ie] << " L0= " << m_caloStruct[i].ele_L0[ie] << " zpos=" << zOffset+zOverburden << " mm w= " << m_caloStruct[i].ele_thick[ie] << " mm";
	  //G4cout << " d=" << m_materials[eleName]->GetDensity();
	  //G4cout << G4endl;
	  //G4cout << *(m_materials[eleName]->GetMaterialTable()) << G4endl;

	  totalLengthX0 += m_caloStruct[i].ele_thick[ie]/m_caloStruct[i].ele_X0[ie];
	  G4cout << " TotX0 = " << totalLengthX0;// << G4endl;
	  totalLengthL0 += m_caloStruct[i].ele_thick[ie]/m_caloStruct[i].ele_L0[ie];
	  G4cout << " TotLambda = " << totalLengthL0 << G4endl;
	  
	  //G4double xpvpos = -m_CalorSizeXY/2.;
	  G4double xpvpos = 0.;
	  G4double ypvpos = 0.;
	  if (m_caloStruct[i].isSensitiveElement(ie)) m_logicSi.push_back(logi);
	  //std::cout << "Element of ele_vol is " << nEle*sectorNum  << std::endl;
	  if (eleName=="Si"){
	    switch (idx) {
	    case 1: xpvpos = -(width-0.2)/4. - gapx; ypvpos = -(width-0.2)/4. - gapy; break;
	    case 2: xpvpos = -(width-0.2)/4. - gapx; ypvpos = (width-0.2)/4. + gapy; break;
	    case 3: xpvpos = (width-0.2)/4. + gapx; ypvpos = -(width-0.2)/4. - gapy; break;
	    case 4: xpvpos = (width-0.2)/4. + gapx; ypvpos = (width-0.2)/4. + gapy; break;
	    }	    
	  }
	  if(skirocvol){
	    //The 4x4 structure of chips. I will assume for now that the 4 chips
	    //are symmetrical placed in every Si pad and will take the Si pad hardcoded for now. 
  	    G4double sipadwidth = 88.48*mm;//Should fix this !!!
	    G4double auxdisvarx =  gapx + (sipadwidth/4.);
	    G4double auxdisvary =  gapy + (sipadwidth/4.);
	    switch (chip_index) {
	    case 1: xpvpos = - auxdisvarx; ypvpos =  - auxdisvary - (sipadwidth/2.); break;
	    case 2: xpvpos = - auxdisvarx; ypvpos =  - auxdisvary; break;
	    case 3: xpvpos = - auxdisvarx - (sipadwidth/2.); ypvpos =  - auxdisvary - (sipadwidth/2.); break;
	    case 4: xpvpos = - auxdisvarx - (sipadwidth/2.); ypvpos =  - auxdisvary; break;
	      
	    case 5: xpvpos = - auxdisvarx; ypvpos =  auxdisvary + (sipadwidth/2.); break;
	    case 6: xpvpos = - auxdisvarx; ypvpos =  auxdisvary; break;
	    case 7: xpvpos = - auxdisvarx - (sipadwidth/2.); ypvpos =  auxdisvary + (sipadwidth/2.); break;
	    case 8: xpvpos = - auxdisvarx - (sipadwidth/2.); ypvpos =  auxdisvary; break;

	    case 9: xpvpos = auxdisvarx; ypvpos =  auxdisvary + (sipadwidth/2.); break;
	    case 10: xpvpos = auxdisvarx; ypvpos =  auxdisvary; break;
	    case 11: xpvpos = auxdisvarx + (sipadwidth/2.); ypvpos =  auxdisvary + (sipadwidth/2.); break;
	    case 12: xpvpos = auxdisvarx + (sipadwidth/2.); ypvpos =  auxdisvary; break;

	    case 13: xpvpos = auxdisvarx; ypvpos =  - auxdisvary - (sipadwidth/2.); break;
	    case 14: xpvpos = auxdisvarx; ypvpos =  - auxdisvary; break;
	    case 15: xpvpos = auxdisvarx + (sipadwidth/2.); ypvpos =  - auxdisvary - (sipadwidth/2.); break;
	    case 16: xpvpos = auxdisvarx + (sipadwidth/2.); ypvpos =  - auxdisvary; break;
	    }
	    
	  }


	  m_caloStruct[i].ele_vol[nEle*sectorNum+ie]= new G4PVPlacement(0, G4ThreeVector(xpvpos,ypvpos,zOffset+zOverburden+thick/2.), logi, baseName+"phys", m_logicWorld, false, 0);	
	  std::cout << " positionning layer at " << xpvpos << " and zpos " << zOffset+zOverburden << std::endl;
	  
	  G4VisAttributes *simpleBoxVisAtt= new G4VisAttributes(m_caloStruct[i].g4Colour(ie));
	  simpleBoxVisAtt->SetVisibility(true);
	  logi->SetVisAttributes(simpleBoxVisAtt);
	  //The z position of the Si pads and SKIROC must be the same
	  if ( ( (eleName=="Si") && idx==4) || ( (eleName!="Si") && (eleName!="PVC") ) || (eleName=="PVC" && chip_index==16) ){
	    zOverburden = zOverburden + thick;
	  }
	  //The z position of the SKIROC must be the same
	  // if ( (eleName=="PVC" && chip_index==16) || (eleName!="PVC") ){
	  //   zOverburden = zOverburden + thick;
	  // }
	  //for sensitive volumes
	  //add region to be able to set specific cuts for it
	  //just for Si
	  if (eleName=="Si"){
	    unsigned nlogicsi = m_logicSi.size();
	    G4Region* aRegion = new G4Region(baseName+"Reg");
	    m_logicSi[nlogicsi-1]->SetRegion(aRegion);
	    aRegion->AddRootLogicalVolume(m_logicSi[nlogicsi-1]);
	  }
	}

      }//loop on elements
    }//loop on layers

  
  


}//buildstack

void DetectorConstruction::SetMagField(G4double fieldValue)
{

  if(fieldValue<=0) return; 

  //apply a global uniform magnetic field along Z axis
  G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  if(m_magField) delete m_magField;                //delete the existing magn field
  m_magField = new G4UniformMagField(G4ThreeVector(0.,0.,fieldValue));
  fieldMgr->SetDetectorField(m_magField);
  fieldMgr->CreateChordFinder(m_magField);
  fieldMgr->SetDetectorField(m_magField);  
}

void DetectorConstruction::SetDetModel(G4int model)
{
  if (model <= 0) return;
  std::cout << " -- Setting detector model to " << model << std::endl;
  model_ = model;
}

//Gap between the Si sensors of single layer (2x2)
void DetectorConstruction::SetGapBetweenSensorPads(G4double gapx_, G4double gapy_)
{
  if ( (gapx_ <= 0) || (gapx_ <= 0) ){
    std::cout << " Error in the gap between sensor pads: Gap_x " << gapx_ << " Gap_y " << gapy_ << std::endl;
    std::cout << " These should be positive numbers " << std::endl;
    exit(1);
  } 
  gapx = gapx_;
  gapy = gapy_;
  std::cout << " Setting the gap between sensor pads: Gap_x: " << gapx << " mm Gap_y: " << gapy << " mm" << std::endl;
}

G4CSGSolid *DetectorConstruction::constructSolid (std::string baseName, const G4double & width, const G4double & height, G4double thick, bool sivol, bool skiroc){
  
  G4CSGSolid *solid;
  if(sivol){
    //This is for the 2 by 2 Si pads
    solid = new G4Box(baseName+"box", (width-0.2)/4., (width-0.2)/4., thick/2 );
  } else if (skiroc){
    //SKIROC
    solid = new G4Box(baseName+"box", width/2, height/2, thick/2 );
  } else {
    solid = new G4Box(baseName+"box", width/2, width/2, thick/2 );
  }

  return solid;
}
