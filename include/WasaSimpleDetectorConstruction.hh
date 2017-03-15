// --------------------------------------------------------- 
// Definition of the WasaSimpleDetectorConstruction class
// Created by C.Rappold (c.rappold@gsi.de)
//----------------------------------------------------------


#ifndef WASADET_h
#define WASADET_h

//#include "WasaSimpleMaterialManager.hh"
//#include "WasaSimpleFieldSetup.hh"
//#include "WasaSimpleCommon.h"
#include "G4SolVDetectorConstruction.hh"

//#include "HypHIVDetectorConstruction.hh"
#include "globals.hh"
//#include "HypHIFrsMagneticField.hh"
#include "G4SolSimpleMagneticField.hh"
//#include "THypHi_Par.hh"
#include "G4FieldManager.hh"

#include "G4Color.hh"
#include "G4SolConfig.hh"
#include "G4VPhysicalVolume.hh"

//class HypHIFrsMagneticField;
//class THypHi_Par;
class G4VSolid;
class G4PVPlacement;
class G4VPhysicalVolume;
class G4Material;
class G4VSensitiveDetector;
class G4VisAttributes;
class G4GenericMessenger;

class G4SolSimpleMagneticField;
class G4FieldManager;



class WasaSimpleDetectorConstruction : public G4SolVDetectorConstruction
{
private:
  static G4ThreadLocal G4SolSimpleMagneticField* fMagneticField;
  static G4ThreadLocal G4FieldManager* fFieldMgr;

  const G4Colour  Blue        = {0.0, 0.0, 1.0} ;
  const G4Colour  Gray        = {0.5, 0.5, 0.5} ;
  const G4Colour  Red         = {1.0, 0.0, 0.0} ;
  const G4Colour  LightBlue   = {0.0, 0.0, 1.0,0.7} ;
  const G4Colour  Yellow      = {1.0, 1.0, 0.0} ; 
  const G4Colour  Purple      = {.5, 0.0, .5} ;
  const G4Colour  LightPurple = {106./255., 27./255., 189./255.} ;
  const G4Colour  Green       = {0.0, 1.0, 0.0};
  const G4Colour  Orange      = {1.0,0.647,0};
  const G4Colour  Pink        = {1.0,0.753,0.796};  

  const G4Colour PaletteOrange[5] = {{254/255.,237/255.,222/255.}, {253/255.,190/255.,133/255.}, {253/255.,141/255.,060/255.},
				     {230/255., 85/255.,13/255.}, {166/255.,54/255.,3/255.}};
  const G4Colour PaletteGreen[5]  = {{237/255.,248/255.,233/255.}, {186/255.,228/255.,179/255.}, {116/255.,196/255.,118/255.},
				     { 49/255.,163/255.,84/255.}, {0,109/255.,44/255.}};
  const G4Colour PaletteRed[5]    = {{254/255.,229/255.,217/255.}, {252/255.,174/255.,145/255.}, {251/255.,106/255.,74/255.}, 
				     {222/255.,45/255.,38/255.}, {165/255.,15/255.,21/255.}};
    
  
  const G4Colour ColorCDC[17] = {PaletteOrange[0],PaletteOrange[0],PaletteOrange[0],PaletteOrange[0],PaletteOrange[0],
				 PaletteOrange[2],PaletteOrange[2],PaletteOrange[2],PaletteOrange[2],PaletteOrange[2],PaletteOrange[2],
				 PaletteOrange[3],PaletteOrange[3],PaletteOrange[3],PaletteOrange[3],PaletteOrange[3],PaletteOrange[3]};
  

  void ConstructCD(G4double cd_rmax,G4double cd_z, G4double cdsPos_x, G4double cdsPos_y, G4double cdsPos_z);
  void ConstructPSB();
  void ConstructDownTracker(G4double cds_z, G4double RelativePos, G4double TargetPos_x, G4double TargetPos_y, G4double TargetPos_z);
  void ConstructEndFMF2();
  //G4LogicalVolume* FindVolume(const G4String& name);

  void DefineCommands();
  G4GenericMessenger* fMessenger;

  //const G4SolConfig& Par;
  
public:

  explicit WasaSimpleDetectorConstruction(const G4SolConfig& conf);
  ~WasaSimpleDetectorConstruction();

  void ConstructMaterials();
  G4VPhysicalVolume* Construct() override;
  void ConstructSDandField() override;
  //const std::vector<G4String>& GetNameDetectors() const { return NameDetectorsSD;}
  bool DoCD;
  bool DoPSB;
  bool DoEndFMF2;
  
  bool DoModHypHI;
  int DoOnlySense;

  bool DoForRoot;
  
  G4VSolid*        experimentalHall_box;
  G4LogicalVolume* experimentalHall_log;
  G4VPhysicalVolume*   experimentalHall_phys;

  //G4LogicalVolume* experimentalHall_logOutRoot;
  //G4VPhysicalVolume*   experimentalHall_physOutRoot;

  G4VSolid*        CD_tube;
  G4LogicalVolume* CD_log;
  G4PVPlacement*   CD_phys;

  G4VSolid*        Solenoid_tube;
  G4LogicalVolume* Solenoid_log;
  G4PVPlacement*   Solenoid_phys;

  
  G4LogicalVolume* CD_logOutRoot;
  G4PVPlacement*   CD_physOutRoot;

  G4LogicalVolume* CDS_endcap_log;

  G4VSolid*        HypHI_Target;
  G4LogicalVolume* HypHI_Target_log;
  G4PVPlacement*   HypHI_Target_phys;
  
  G4VSolid*        HypHI_InTracker;
  G4LogicalVolume* HypHI_InTracker_log;
  G4PVPlacement*   HypHI_InTracker_phys;

  G4VSolid*        HypHI_Endcap;
  G4LogicalVolume* HypHI_Endcap_log;
  G4PVPlacement*   HypHI_Endcap_phys;

  G4LogicalVolume* HypHI_TrackerFwd_log;
  G4LogicalVolume* HypHI_RPC_l_log;
  G4LogicalVolume* HypHI_RPC_h_log;
  
  std::vector<G4PVPlacement*> AllPlacements;
  
  //std::vector<G4String> NameDetectorsSD;
  
  double FieldInCD;         

  
  double TargetLength;       
  int TargetMaterialID;   
  std::string ChamberGasSelection;
  std::string CDCType;            

  double ReductionRadius;


  const double CD_Rmax = 300.;
  const double CD_Rmin = 30.;
  const double CD_Z = 560.;

  const double Solenoid_RadiusMin = 267.8 ; //mm
  const double Solenoid_RadiusMax = 288.8 ; //mm
  const double Solenoid_Length = 465. ; //mm

  const double MDC_RadiusMin = 41.;
  const double MDC_RadiusMax = 203.;
  const double MDC_LayerSpacing = 3.5;
  const int MDC_NbLayer = 17;
  const double MDC_Cell_Thinkness [17] = {4., 4., 4., 4., 4.,
					  6., 6., 6., 6., 6., 6.,
					  8., 8., 8., 8., 8., 8.}; //mm

  const int MDC_Cell_XXP[17] = {0, 1, 0, 1, 0,
				0, 1, 0, 1, 0, 1,
				0, 1, 0, 1, 0, 1};
  
  const double PSB_Length = 550.;
  const double PSB_Width = 38.;
  const double PSB_Thickness = 8.;
  const int PSB_Nb = 24;
  
  std::unique_ptr<G4Material> ArIsoButane;
  
  
};

#endif
