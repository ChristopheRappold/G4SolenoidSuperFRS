#ifndef KNUCLDET_h
#define KNUCLDET_h

//#include "KnuclMaterialManager.hh"
//#include "KnuclFieldSetup.hh"
//#include "KnuclCommon.h"
#include "G4VUserDetectorConstruction.hh"

//#include "HypHIVDetectorConstruction.hh"
#include "globals.hh"
//#include "HypHIFrsMagneticField.hh"
#include "G4SolSimpleMagneticField.hh"
//#include "THypHi_Par.hh"
#include "G4FieldManager.hh"

#include "G4Color.hh"

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



class KnuclDetectorConstruction : public G4VUserDetectorConstruction
{
private:
  static G4ThreadLocal G4SolSimpleMagneticField* fMagneticField;
  static G4ThreadLocal G4FieldManager* fFieldMgr;

  const G4Colour  Blue      = {0.0, 0.0, 1.0} ;
  const G4Colour  Gray      = {0.5, 0.5, 0.5} ;
  const G4Colour  Red       = {1.0, 0.0, 0.0} ;
  const G4Colour  LightBlue = {0.0, 0.0, 1.0,0.7} ;
  const G4Colour  Yellow    = {1.0, 1.0, 0.0} ; 
  const G4Colour  Purple    = {.5, 0.0, .5} ;
  const G4Colour  Green     = {0.0, 1.0, 0.0};
  const G4Colour  Orange    = {1.0,0.647,0};
  const G4Colour  Pink      = {1.0,0.753,0.796};  

  void ConstructKurama();
  void ConstructTOFn();
  void ConstructTOFp();
  void ConstructCDS(G4double cds_rmax,G4double cds_z, G4double cdsPos_x, G4double cdsPos_y, G4double cdsPos_z);
  void ConstructCDH(G4double cdsPos_x, G4double cdsPos_y, G4double cdsPos_z);
  void ConstructTargetChamber(G4double cdsPos_x, G4double cdsPos_y, G4double cdsPos_z, G4double tarCham_rmax, G4double tarCham_z, G4double tarChamPos_x, G4double tarChamPos_y, G4double tarChamPos_z);
  void ConstructKaonVeto(G4double cdsPos_x, G4double cdsPos_y, G4double tarCham_rmax);
  void ConstructAC(G4double cds_z, G4double CDS_AC_space, G4double AC_STC_space, G4double STC_BLC_space, G4double BLC_BLC_space);

  void DefineCommands();
  G4GenericMessenger* fMessenger;

public:
  //explicit KnuclDetectorConstruction();//const THypHi_Par& _par);
  KnuclDetectorConstruction();//KnuclAnaManager* ana);
  ~KnuclDetectorConstruction();

  void ConstructMaterials();
  G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();

  bool DoKurama;
  bool DoTOFn;
  bool DoTOFp;
  bool DoCDS;
  bool DoCDH;
  bool DoTargetChamber;
  bool DoAC;

  
  G4VSolid*        experimentalHall_box;
  G4LogicalVolume* experimentalHall_log;
  G4PVPlacement*   experimentalHall_phys;

  G4VSolid*        Kurama_box;
  G4LogicalVolume* Kurama_log;
  G4PVPlacement*   Kurama_phys;

  G4LogicalVolume* KuramaAperture_log;
  
  
  G4VSolid*        TOF_box;
  G4LogicalVolume* TOF_log;
  G4PVPlacement*   TOF_phys;
  
  G4VSolid*        CDS_tube;
  G4LogicalVolume* CDS_log;
  G4PVPlacement*   CDS_phys;

  G4LogicalVolume* CDS_endcap_log;

  std::vector<G4PVPlacement*> AllPlacements;
  
  double FieldInCDC;         
  double FieldInKurama;      
  double TofRefPos;          
  int Ncdh;               
  double ThicknessOfTgtCell; 
  double ThicknessOfTgtAl;   
  double ThicknessOfTgtCFRP; 
  double ThicknessOfCDCwCFRP;
  double TargetLength;       
  int TargetMaterialID;   
  std::string ChamberGasSelection;
  std::string CDCType;            
  double CDSLength;          
  int ChamberSmear;       
  double CDCResolution;      
  double ZVCLength;          
  int ZVertexChamber;       
  double ZVertexChamberCell;



  
  //### constant for CDS(solenoid) ###//
  const double CDS_RMAX = 590.0; // unit :: mm
  const double CDS_Z = 1170.0; // unit :: mm
  const double CDSYOKE_RMAX = 985.0; // unit :: mm
  const double CDSYOKE_RMIN = CDS_RMAX; // unit :: mm
  const double CDS_ENDCAP_RMIN = 150.0; // unit :: mm
  const double CDS_ENDCAP_Z = 155.0; // unit :: mm


  //### constant for target ###//
  const double TARGET_CHM_RMAX = 75.0; // unit :: mm
  const double TARGET_RMAX = 35.0; // unit :: mm
  
  //### constant for ZVC ###//
  const double ZVC_R = 130.0; // unit :: mm
  
  //### constant for CDC ###//
  const double CDC_RMIN = 150.0; // unit :: mm
  const double CDC_RMAX = 530.0; // unit :: mm
  const double CDC_RADIUS[15] = {190.5, 204.0, 217.5,
				 248.5, 262.0, 293.0, 306.5,
				 337.5, 351.0,
				 382.0, 395.5, 426.5, 440.0,
				 471.0, 484.5}; // unit :: mm
  const int N_CDC_CELL[15] = {72, 72, 72,
			      90, 90, 100, 100,
			      120, 120,
			      150, 150, 160, 160,
			      180, 180};
  const int CDC_CELL_XXP[15] = {0, 1, 0,
				0, 1, 0, 1,
				0, 1,
				0, 1, 0, 1,
				0, 1};
  const double CDC_OFFSET = 3.0; // unit :: cell
  const double CDC_CELL_OFFSET[15] = {0.0, 0.0, 0.0,
					CDC_OFFSET, CDC_OFFSET, -CDC_OFFSET,
					-CDC_OFFSET,0.0, 0.0,
					CDC_OFFSET, CDC_OFFSET, -CDC_OFFSET,
					-CDC_OFFSET,0.0, 0.0};
  
  
  //### constant for CDH ###//
  const double CDH_RMAX = 544.0; // unit :: mm
  const double CDH_RMIN = 574.0; // unit :: mm
  const double CDH_Z = 790.0; // unit :: mm



  
  
};

#endif
