#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "SamplingSection.hh"

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include <map>
#include <string>

class G4CSGSolid;
class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class DetectorMessenger;
class G4Colour;

/**
   @class DetectorConstruction
   @short builds a simple detector
 */
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

  enum DetectorVersion { 
    v_CALICE=0,
    v_SiWECAL_1=1, // Absorber before: 0.0 Absorber inside 0.0 	Angle 0 
    v_SiWECAL_2=2, // Absorber before: 4.2X0W  Absorber inside	0.0  	Angle 0
    v_SiWECAL_3=3, // Absorber before: 8.4X0W  Absorber inside 	0.0 	Angle 0
    v_SiWECAL_4=4, // Absorber before: 8.4X0W+200mmFe  Absorber inside 	0.0 	Angle 0
    v_SiWECAL_5=5, // Absorber before: 8.4X0W+300mmFe  Absorber inside 	0.0 	Angle 0
    v_SiWECAL_6=6, // Absorber before: 0.0   Absorber inside	4.2X0W between each plate 	Angle 0
    v_SiWECAL_7=7, // Absorber before: 1.8X0W  Absorber inside 	4.2X0W between each plate 	Angle 0
    v_SiWECAL_8=8, // Absorber before: 2.4X0W  Absorber inside 	4.2X0W between each plate 	Angle 0
    v_SiWECAL_9=9, // Absorber before: 4.2X0W  Absorber inside 	4.2X0W between each plate 	Angle 0
    v_SiWECAL_10=10, // Absorber before: 2.4X0W   Absorber inside	4.2X0W between each plate 	Angle 48
    v_SiWECAL_11=11 // Absorber before: 0.0   Absorber inside	0.0 	Angle 90 

  };

  enum DetectorModel {
    m_SiWECAL_0=0
  };

  /**
     @short CTOR
   */
  DetectorConstruction(G4int ver=DetectorConstruction::v_SiWECAL_1, 
		       G4int mod=DetectorConstruction::m_SiWECAL_0,
		       std::string absThickW="1.75,1.75,1.75,1.75,1.75,2.8,2.8,2.8,2.8,2.8,4.2,4.2,4.2,4.2,4.2");

  /**
     @short calorimeter structure (sampling sections)
   */
  std::vector<SamplingSection> m_caloStruct;
  std::vector<SamplingSection> *getStructure() { return &m_caloStruct; }

  int getModel() const { return model_; }
  int getVersion() const { return version_; }

  const std::vector<G4LogicalVolume*>  & getSiLogVol() {return m_logicSi; }
  const std::vector<G4LogicalVolume*>  & getAbsLogVol() {return m_logicAbs; }

  /**
     @short define the calorimeter materials
   */
  void DefineMaterials(); 
  std::map<std::string, G4Material *> m_materials;
  std::map<std::string, G4double > m_dEdx;
  std::map<std::string, G4Colour > m_colours;

  /**
     @short set magnetic field
   */
  void SetMagField(G4double fieldValue);
  G4UniformMagField* m_magField;      //pointer to the magnetic field

  /**
     @short set detector model
   */

  void SetDetModel(G4int model);

  void SetWThick(std::string thick);

  void SetGapBetweenSensorPads(G4double gapx_, G4double gapy_);

  /**
     @short DTOR
   */
  ~DetectorConstruction();
  

  /**
     @short getters
   */
  G4double GetCalorSizeXY() { return m_CalorSizeXY; }
  G4double GetCalorSizeZ()  { return m_CalorSizeZ; }
  G4double GetWorldSizeXY() { return m_WorldSizeXY; }
  G4double GetWorldSizeZ()  { return m_WorldSizeZ; }

  /**
     @short build the detector
   */

  G4VPhysicalVolume* Construct();

private:

  //detector version
  int version_;
  //integer to define detector model
  int model_;

  //add a pre PCB plate
  bool addPrePCB_;

  std::vector<G4double> absThickW_;

  //For the gap between Si sensor pad of a single layer (2x2)
  G4double gapx;
  G4double gapy;

  /**
     @short compute the calor dimensions
   */
  void UpdateCalorSize(); 

  /**
     @short build the calorimeter
   */
  G4VPhysicalVolume* ConstructCalorimeter();     

  void buildStack(const unsigned sectorNum);

  G4CSGSolid *constructSolid (std::string baseName, const G4double & width, const G4double & height, G4double thick, bool sivol, bool skirocvol);
  

  std::vector<G4Material* > m_SensitiveMaterial;
  
  G4double           m_CalorSizeXY, m_CalorSizeZ;
  G4double           m_WorldSizeXY, m_WorldSizeZ;
  G4double m_nSectors,m_sectorWidth,m_interSectorWidth;
            
  G4CSGSolid*        m_solidWorld;    //pointer to the solid World 
  G4LogicalVolume*   m_logicWorld;    //pointer to the logical World
  G4VPhysicalVolume* m_physWorld;     //pointer to the physical World  
  
  std::vector<G4LogicalVolume*>   m_logicSi;    //pointer to the logical Si volumes
  std::vector<G4LogicalVolume*>   m_logicAbs;    //pointer to the logical absorber volumes situated just before the si

  DetectorMessenger* m_detectorMessenger;  //pointer to the Messenger
};


#endif

