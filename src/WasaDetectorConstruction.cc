//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: WasaDetectorConstruction.cc 77601 2013-11-26 17:08:44Z gcosmo $
// 
/// \file WasaDetectorConstruction.cc
/// \brief Implementation of the WasaDetectorConstruction class

#include "WasaDetectorConstruction.hh"

#include "globals.hh"

//#include "G4RunManager.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
//#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4NistManager.hh"
#include "G4PhysicalConstants.hh"

// VGM demo
#include "Geant4GM/volumes/Factory.h"
#include "RootGM/volumes/Factory.h"
#include "TGeoManager.h"
// end VGM demo

#include "G4SystemOfUnits.hh"

// #include "G4ChordFinder.hh"
// #include "G4Mag_UsualEqRhs.hh"

// #include "G4ClassicalRK4.hh"
// #include "G4HelixExplicitEuler.hh"
// #include "G4HelixImplicitEuler.hh"
// #include "G4HelixSimpleRunge.hh"
// #include "G4CashKarpRKF45.hh"
// #include "G4NystromRK4.hh"
// #include "G4SimpleHeum.hh"

#include "G4UserLimits.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"


#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include "G4SDManager.hh"

#include "G4SolSimpleMagneticField.hh"
#include "G4SolSensitiveD.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4AutoDelete.hh"

#include <iostream>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4SolSimpleMagneticField* WasaDetectorConstruction::fMagneticField = 0;
G4ThreadLocal G4FieldManager* WasaDetectorConstruction::fFieldMgr = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WasaDetectorConstruction::WasaDetectorConstruction(const G4SolConfig& _par) : G4SolVDetectorConstruction(_par),experimentalHall_log(nullptr),experimentalHall_phys(nullptr),MFLD_log(nullptr),MFLD_phys(nullptr), fCheckOverlaps(true)
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  // ---------------------------------------------------------------------------
  // VGM demo 
  //


WasaDetectorConstruction::~WasaDetectorConstruction()
{ 

}

G4VPhysicalVolume* WasaDetectorConstruction::Construct()
{


  G4VisAttributes VisDetectorSD(Blue);

  // 
  // Import geometry from Root
  //

  // Import geometry from the root file
  new TGeoManager("Wasa", "Geant4 basic example Wasa");

  const std::string nameGeometry = Par.Get<std::string>("Geometry_Namefile");
  gGeoManager->Import(nameGeometry.c_str());

  // Import geometry from Root to VGM
  RootGM::Factory rtFactory;
  rtFactory.SetDebug(0);
  rtFactory.Import(gGeoManager->GetTopNode());

  // Export VGM geometry to Geant4
  //
  Geant4GM::Factory g4Factory;
  g4Factory.SetDebug(0);
  rtFactory.Export(&g4Factory);
  
  G4VPhysicalVolume* world = g4Factory.World();

  experimentalHall_log = world->GetLogicalVolume();
  experimentalHall_phys = world;

  
  MFLD_log = FindVolume("MFLD");
  MFLD_phys = FindVolPhys("MFLD");

  std::cout<<"MFLD Phys:"<<MFLD_phys<<" "<<MFLD_phys->GetInstanceID()<<"\n";
  transMFLD = MFLD_phys->GetObjectTranslation();
  rotMFLD = MFLD_phys->GetObjectRotationValue();
  std::cout<<"Trans:"<<transMFLD<<"\n";
  std::cout<<"Rot:"<<rotMFLD<<"\n";

  const G4double Wasa_Zshift = Par.Get<double>("Wasa_ShiftZ");
  const G4double Systematic_shift =  Par.Get<double>("Systematic_Shift");
  const G4double TargetPosX = Par.Get<double>("Target_PosX");
  const G4double TargetPosY = Par.Get<double>("Target_PosY");
  const G4double TargetPosZ = Par.Get<double>("Target_PosZ");

  G4ThreeVector transMFLD_move(0.,0.,Wasa_Zshift+Systematic_shift);
  G4ThreeVector transMFLD_new = transMFLD+transMFLD_move;

  const int WasaSide = Par.Get<int>("Wasa_Side");
  if(WasaSide==1)
    rotMFLD.rotateY(180*degree);
  
  std::cout<<"Trans:"<<transMFLD_new<<"\n";
  std::cout<<"RotAfter:"<<rotMFLD<<"\n";

  MFLD_phys->SetTranslation(transMFLD_new);

  if(WasaSide==1)
    MFLD_phys->SetRotation(&rotMFLD);
  
  // Get volumes from logical volume store by name
  //G4LogicalVolume* calorLV = FindVolume("Calorimeter");
  
  G4LogicalVolume* WasaLV = FindVolume("WASA");
  G4LogicalVolume* InnerLV = FindVolume("INNER");
  if(!Par.IsAvailable("HypHI_InnerTrackerBox_Visible"))
    InnerLV->SetVisAttributes(G4VisAttributes::Invisible);

  INNER_log = FindVolume("INNER");
  INNER_phys = FindVolPhys("INNER");

  // G4LogicalVolume* DMagLV = FindVolume("DMag1");
  // G4LogicalVolume* VolGapLV = FindVolume("VolGap");
  // G4LogicalVolume* SecondMagFieldLV = FindVolume("SecondMagField");
  // G4LogicalVolume* SecondMagLV = FindVolume("SecondMag");
  // G4LogicalVolume* SetupLV = FindVolume("Setup");
  // G4LogicalVolume* IronQuadLV = FindVolume("IronQuad");

  // G4Region* aTargetRegion = new G4Region("TargetRegion");

  // G4LogicalVolume* TargetLV = FindVolume("Target");
  // G4NistManager* materialMgr = G4NistManager::Instance();
  // G4Material* Vacuum = materialMgr->FindOrBuildMaterial("G4_Galactic");
  // TargetLV->UpdateMaterial(Vacuum);
  
  // TargetLV->SetRegion(aTargetRegion);
  // aTargetRegion->AddRootLogicalVolume(TargetLV);

  //G4Region* aDetectorRegion = new G4Region("DetectorRegion");
  std::vector<G4String> NameCDC1_Invisible = { "ME01","ME02","ME03","ME04","ME05","ME06","ME07","ME08","ME09","ME10","ME11","ME12","ME13","ME14","ME15","ME16","ME17"};
  std::vector<G4String> NameCDC2_Invisible = { "MD01","MD02","MD03","MD04","MD05","MD06","MD07","MD08","MD09","MD10","MD11","MD12","MD13","MD14","MD15","MD16","MD17"};
  for(auto& CurrentName : NameCDC1_Invisible)
    {
      G4LogicalVolume* UTracker = FindVolume(CurrentName);
      UTracker->SetVisAttributes(G4VisAttributes::Invisible);
    }

  //int iColorT = 0;
  for(auto& CurrentName : NameCDC2_Invisible)
    {
      G4LogicalVolume* UTracker = FindVolume(CurrentName);
      //G4VisAttributes MG_Color(ColorCDC[iColorT]);
      //UTracker->SetVisAttributes(MG_Color);
      //++iColorT;
      UTracker->SetVisAttributes(G4VisAttributes::Invisible);
    }
  FindVolume("SOL")->SetVisAttributes(G4VisAttributes::Invisible);

					       
  std::vector<G4String> NameSD = { "MG01","MG02","MG03","MG04","MG05","MG06","MG07","MG08","MG09","MG10","MG11","MG12","MG13","MG14","MG15","MG16","MG17",
				   "PSCE","PSBE","PSFE"};

  NameDetectorsSD = NameSD;
  
  std::vector<G4String> NameSD_Color = { "MG01","MG02","MG03","MG04","MG05","MG06","MG07","MG08","MG09","MG10","MG11","MG12","MG13","MG14","MG15","MG16","MG17"};
 
  FindVolume("PSCE")->SetVisAttributes(G4VisAttributes(Blue));
  FindVolume("PSBE")->SetVisAttributes(G4VisAttributes(Blue));
  FindVolume("PSFE")->SetVisAttributes(G4VisAttributes(LightBlue));

  //     UTracker->SetRegion(aDetectorRegion);
  //     aDetectorRegion->AddRootLogicalVolume(UTracker);
      
  //UTracker->SetVisAttributes(VisDetectorSD);
    
  for(size_t iColor = 0; iColor<NameSD_Color.size();++iColor)
    {
      G4LogicalVolume* UTracker = FindVolume(NameSD_Color[iColor]);
      G4VisAttributes MG_Color(ColorCDC[iColor]);
      UTracker->SetVisAttributes(MG_Color);      
      //UTracker->SetVisAttributes(G4VisAttributes::Invisible);
    }

  

  // SetupLV->SetUserLimits( new G4UserLimits(DBL_MAX,Par.Get_CutLength_Track(),10*s,0.,0.) );
  
  // std::vector<double> cutsDet (4,Par.Get_CutValue_Plastic());
  // aDetectorRegion->SetProductionCuts(new G4ProductionCuts());
  // aDetectorRegion->GetProductionCuts()->SetProductionCuts(cutsDet);

  // std::vector<double> cutsTarget (4,Par.Get_CutValue_Target());
  // aTargetRegion->SetProductionCuts(new G4ProductionCuts());
  // aTargetRegion->GetProductionCuts()->SetProductionCuts(cutsTarget);

  
  // Visualization attributes
  //

  // G4LogicalVolume* worldLV = world->GetLogicalVolume();
  // worldLV->SetVisAttributes (G4VisAttributes::Invisible);

  // worldLV->SetUserLimits( new G4UserLimits(DBL_MAX,2*m,10*s,0.,0.) );  
  
  if(WasaLV)
    WasaLV->SetVisAttributes(G4VisAttributes::Invisible);
  if(MFLD_log)
    if(!Par.IsAvailable("HypHI_InnerTrackerBox_Visible"))
      MFLD_log->SetVisAttributes(G4VisAttributes::Invisible);
  // if(VolGapLV)
  //   VolGapLV->SetVisAttributes(G4VisAttributes::Invisible);
  // if(SecondMagFieldLV)
  //   SecondMagFieldLV->SetVisAttributes(G4VisAttributes::Invisible);
  // if(SecondMagLV)
  //   SecondMagLV->SetVisAttributes(G4VisAttributes::Invisible);
  // if(SetupLV)
  //   SetupLV->SetVisAttributes(G4VisAttributes::Invisible);
  // if(IronQuadLV)
  //   IronQuadLV->SetVisAttributes(G4VisAttributes::Invisible);

  

  // G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  // simpleBoxVisAtt->SetVisibility(true);
  // if (calorLV) calorLV->SetVisAttributes(simpleBoxVisAtt);
  double Sign = WasaSide == 0 ? 1. : -1.;

  G4NistManager* materialMgr = G4NistManager::Instance();
    
  //G4Material* Air = materialMgr->FindOrBuildMaterial("G4_AIR");
  G4Material* Air = materialMgr->FindOrBuildMaterial("G4_AIR");
  G4Material* Vacuum = materialMgr->FindOrBuildMaterial("G4_Galactic");
  G4Material* Si = materialMgr->FindOrBuildMaterial("G4_Si");//"Plastic");
  G4Material* Scinti = materialMgr->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");//G4_POLYETHYLENE");//"Plastic");
  G4Material* FiberCoreScinti =  materialMgr->FindOrBuildMaterial("G4_POLYSTYRENE");
  G4Material* Carbon = materialMgr->FindOrBuildMaterial("G4_C");//"Plastic");

  G4VisAttributes *Si_att = new G4VisAttributes(Pink);

  G4double TargetLength = Par.Get<double>("Target_Size");//1.0*cm;
  G4double BeamHoleSize = 0.*cm;
  if(Par.IsAvailable("HypHI_BeamHole"))
    BeamHoleSize = Par.Get<double>("HypHI_BeamHole");
  
  G4VSolid*        HypHI_Target = new G4Box("HypHI_Target", TargetLength, TargetLength, TargetLength);
  //G4LogicalVolume* HypHI_Target_log = new G4LogicalVolume(HypHI_Target, Carbon,"HypHI_Target_log", 0,0,0);
  G4LogicalVolume* HypHI_Target_log = new G4LogicalVolume(HypHI_Target, Vacuum,"HypHI_Target_log", 0,0,0);
  //G4PVPlacement*   HypHI_Target_phys = 
  G4ThreeVector TargetTrans = G4ThreeVector(TargetPosX, TargetPosY , TargetPosZ)-transMFLD_new;
  G4RotationMatrix* TargetRot = new G4RotationMatrix;
  // if(WasaSide==1)
  //   TargetRot->rotateY(-180*degree);
  G4Transform3D posTarget(*TargetRot, Sign*TargetTrans);
  AllPlacements.emplace_back(new G4PVPlacement(posTarget, HypHI_Target_log, "HypHI_Target_Phys",
					       MFLD_log, false,0));

// -------------------------- @@ -------------------------
/* 
// -------------------------- @@ -------------------------

  G4int nb_panel = 16;
  G4double HypHI_Si_minR = Par.Get<double>("HypHI_Si_minR");
  G4double HypHI_Si_maxR = Par.Get<double>("HypHI_Si_maxR");
  
  G4VSolid* HypHI_SiliciumSeg = new G4Tubs("HypHI_SiSeg",HypHI_Si_minR,HypHI_Si_maxR, 3*mm,
					   -CLHEP::twopi/static_cast<double>(2*nb_panel),2.*CLHEP::twopi/static_cast<double>(2*nb_panel)); 

  std::vector<double> posZ = {20.*cm, 25.*cm,30.*cm,40.*cm};
  if(Par.IsAvailable("HypHI_InnerTracker_Spec"))
    {
      double PosZInTracker = Par.Get<double>("HypHI_InnerTracker_PosZ");
      int NbInTracker = Par.Get<int>("HypHI_InnerTracker_Nb");
      double DistInTracker = Par.Get<double>("HypHI_InnerTracker_Spacing");
      posZ.resize(NbInTracker);
      for(size_t idL = 0; idL < posZ.size();++idL)
	posZ[idL] = TargetLength + TargetPosZ+PosZInTracker+static_cast<double>(idL)*DistInTracker;
    }

  for(size_t idLayer = 0;idLayer<posZ.size();++idLayer)
    {
      std::string name_Si("HypHI_InSi_log");
      name_Si+=std::to_string(idLayer);
      G4LogicalVolume* HypHI_InSi_log = new G4LogicalVolume(HypHI_SiliciumSeg,Si,name_Si, 0,0,0);
      NameDetectorsSD.push_back(HypHI_InSi_log->GetName());

      for(G4int IdSi = 0 ; IdSi<nb_panel;++IdSi)
	{
	  G4RotationMatrix* rotSi = new G4RotationMatrix;
	  double rotAngle = CLHEP::twopi/static_cast<double>(nb_panel)*static_cast<double>(IdSi);
	  rotSi->rotateZ(rotAngle);
	  std::string nameSi ("HypHI_Layer");
	  nameSi += std::to_string(idLayer);
	  nameSi += "_SiSeg_";
	  nameSi += std::to_string(IdSi);
	  G4ThreeVector SiTrans = G4ThreeVector(0., 0., posZ[idLayer])-transMFLD_new;
	  // if(WasaSide==1)
	  //   rotSi->rotateY(-180*degree);
	  G4Transform3D posSi(*rotSi, Sign*SiTrans);

	  AllPlacements.emplace_back(new G4PVPlacement(posSi,HypHI_InSi_log, nameSi, MFLD_log, false, IdSi));
	}

      Si_att->SetForceWireframe(false);
      HypHI_InSi_log->SetVisAttributes(Si_att);
    }


  const double TR1_posZ = Par.Get<double>("HypHI_TR1_posZ");
  G4VSolid* TR1_box = nullptr;
  if(Par.IsAvailable("HypHI_BeamHole"))
    {
      G4VSolid* TR1_box_init = new G4Box("TR1_box_init",15.*cm,15.*cm,1.*mm);
      G4VSolid* TR1_box_hole = new G4Box("TR1_box",BeamHoleSize*0.5,BeamHoleSize*0.5,1.*mm);
      TR1_box = new G4SubtractionSolid("TR1_box", TR1_box_init, TR1_box_hole);
    }
  else
    TR1_box = new G4Box("TR1_box",15.*cm,15.*cm,1.*mm);

  G4LogicalVolume* TR1_log = new G4LogicalVolume(TR1_box, Scinti, "TR1_log", 0, 0, 0);
  AllPlacements.emplace_back(new G4PVPlacement(0,Sign*(G4ThreeVector(0., 0., TR1_posZ+Systematic_shift)-transMFLD_new), TR1_log, "TR1_phys", MFLD_log,false,0));  
  NameDetectorsSD.push_back(TR1_log->GetName());

  const double TR2_posZ = Par.Get<double>("HypHI_TR2_posZ");
  
  G4VSolid* TR2_box = nullptr;
  if(Par.IsAvailable("HypHI_BeamHole"))
    {
      G4VSolid* TR2_box_init = new G4Box("TR2_box_init",15.*cm,15.*cm,1.*mm);
      G4VSolid* TR2_box_hole = new G4Box("TR2_box_hole",BeamHoleSize*0.5, BeamHoleSize*0.5,1.*mm);
      TR2_box = new G4SubtractionSolid("TR2_box", TR2_box_init, TR2_box_hole);
    }
  else
    TR2_box = new G4Box("TR2_box",15.*cm,15.*cm,1.*mm);

  G4LogicalVolume* TR2_log = new G4LogicalVolume(TR2_box, Scinti, "TR2_log", 0, 0, 0);
  AllPlacements.emplace_back(new G4PVPlacement(0,Sign*(G4ThreeVector(0., 0., TR2_posZ+Systematic_shift)-transMFLD_new), TR2_log, "TR2_phys", MFLD_log,false,0));  
  NameDetectorsSD.push_back(TR2_log->GetName());
// -------------------------- @@ -------------------------
*/
// -------------------------- @@ -------------------------



// ----------------------------- @@ -----------------------------
// 		     Fiber Detectors
// ----------------------------- @@ -----------------------------
	const double HypHI_FiberTracker1_posZ = Par.Get<double>("HypHI_FiberTracker1_posZ");
	const double HypHI_FiberTracker2_posZ = Par.Get<double>("HypHI_FiberTracker2_posZ");
	const double HypHI_FiberTracker3_posZ = Par.Get<double>("HypHI_FiberTracker3_posZ");
// -------------------------- --------------------------
	G4RotationMatrix* rotFib1 = new G4RotationMatrix;
	rotFib1->rotateX(90.*deg);
	G4RotationMatrix* rotFib2 = new G4RotationMatrix;
	rotFib2->rotateX(90.*deg); 
	rotFib2->rotateY(30.*deg);
	G4RotationMatrix* rotFib3 = new G4RotationMatrix;
	rotFib3->rotateX(90.*deg);
	rotFib3->rotateY(-30.*deg);
	G4double spacingX = 0.55 *mm;
	G4double startZ1 = 0.*mm;
	G4double startZ2 = 0.47631397*mm;
	std::vector<G4double> posZshift = {-4.*mm, 0.*mm, 4.*mm};
	G4ThreeVector posFib;


	G4VisAttributes* visAttributes_x = new G4VisAttributes(color_x);
	G4VisAttributes* visAttributes_u = new G4VisAttributes(color_u);
	G4VisAttributes* visAttributes_v = new G4VisAttributes(color_v);

// -------------------------- First Fiber Detector --------------------------
	G4VSolid* FiberD1_MothVol_solid = new G4Box("FiberDetector1",15.*cm,20.*cm,3.*cm);
	G4VSolid* FiberD1_layerX_MothVol_solid  = new G4Box("FiberD1_layerX_solid", 15.*cm, 20.*cm, 2.*mm);
	G4VSolid* FiberD1_layerU_MothVol_solid  = new G4Box("FiberD1_layerU_solid", 15.*cm, 20.*cm, 2.*mm);
	G4VSolid* FiberD1_layerV_MothVol_solid  = new G4Box("FiberD1_layerV_solid", 15.*cm, 20.*cm, 2.*mm);
	G4LogicalVolume* FiberD1_MothVol_log = new G4LogicalVolume(FiberD1_MothVol_solid, Air, "FiberD1_log", 0, 0, 0);
	G4LogicalVolume* FiberD1_MothVol_log_x  = new G4LogicalVolume(FiberD1_layerX_MothVol_solid, Air, "FiberD1_log_x", 0, 0, 0);
	G4LogicalVolume* FiberD1_MothVol_log_u  = new G4LogicalVolume(FiberD1_layerU_MothVol_solid, Air, "FiberD1_log_u", 0, 0, 0);
	G4LogicalVolume* FiberD1_MothVol_log_v  = new G4LogicalVolume(FiberD1_layerV_MothVol_solid, Air, "FiberD1_log_v", 0, 0, 0);

	AllPlacements.emplace_back(new G4PVPlacement(0,Sign*(G4ThreeVector(0. ,0. , HypHI_FiberTracker1_posZ + Systematic_shift)-transMFLD_new), FiberD1_MothVol_log, "FiberDetector1", MFLD_log,false,0)); 
	AllPlacements.emplace_back(new G4PVPlacement(0,Sign*(G4ThreeVector(0. ,0. , posZshift[0])), FiberD1_MothVol_log_x, "FiberD1_x", FiberD1_MothVol_log,false,0)); 
	AllPlacements.emplace_back(new G4PVPlacement(0,Sign*(G4ThreeVector(0. ,0. , posZshift[1])), FiberD1_MothVol_log_u, "FiberD1_u", FiberD1_MothVol_log,false,0)); 
	AllPlacements.emplace_back(new G4PVPlacement(0,Sign*(G4ThreeVector(0. ,0. , posZshift[2])), FiberD1_MothVol_log_v, "FiberD1_v", FiberD1_MothVol_log,false,0)); 	

	std::vector<G4double> FD1_startX1 = {-70.125*mm, -39.325*mm, -4.125*mm,  31.075*mm};
	std::vector<G4double> FD1_startX2 = {-69.85*mm,  -30.25*mm,   4.95*mm,   40.15*mm };
	std::vector<G4double> numFibersTopLayer1 = {56, 64, 64, 72};

	G4VSolid* FiberD1_Core_solid_x = new G4Tubs("FiberD1_Core_solid_x",0, 0.24*mm, 10.5*cm, 0.*deg, 360.*deg);
	G4VSolid* FiberD1_Core_solid_u = new G4Tubs("FiberD1_Core_solid_u",0, 0.24*mm, (10.5/sqrt(3)*2 )*cm, 0.*deg, 360.*deg);
	G4VSolid* FiberD1_Core_solid_v = new G4Tubs("FiberD1_Core_solid_v",0, 0.24*mm, (10.5/sqrt(3)*2 )*cm, 0.*deg, 360.*deg);
	G4LogicalVolume* FiberD1_Core_log_x = new G4LogicalVolume(FiberD1_Core_solid_x,FiberCoreScinti, "FiberD1_Core_log_x", 0, 0, 0);
	G4LogicalVolume* FiberD1_Core_log_u = new G4LogicalVolume(FiberD1_Core_solid_u,FiberCoreScinti, "FiberD1_Core_log_u", 0, 0, 0);	
	G4LogicalVolume* FiberD1_Core_log_v = new G4LogicalVolume(FiberD1_Core_solid_v,FiberCoreScinti, "FiberD1_Core_log_v", 0, 0, 0);	

	G4VSolid* FiberD1_Cladding_solid_x = new G4Tubs("FiberD1_Cladding_solid_x",0.24*mm, 0.25*mm, 10.5*cm, 0.*deg, 360.*deg);
	G4VSolid* FiberD1_Cladding_solid_u = new G4Tubs("FiberD1_Cladding_solid_u",0.24*mm, 0.25*mm, (10.5/sqrt(3)*2 )*cm, 0.*deg, 360.*deg);
	G4VSolid* FiberD1_Cladding_solid_v = new G4Tubs("FiberD1_Cladding_solid_v",0.24*mm, 0.25*mm, (10.5/sqrt(3)*2 )*cm, 0.*deg, 360.*deg);
	G4LogicalVolume* FiberD1_Cladding_log_x = new G4LogicalVolume(FiberD1_Cladding_solid_x, Scinti, "FiberD1_Cladding_log_x", 0, 0, 0);	
	G4LogicalVolume* FiberD1_Cladding_log_u = new G4LogicalVolume(FiberD1_Cladding_solid_u, Scinti, "FiberD1_Cladding_log_u", 0, 0, 0);
	G4LogicalVolume* FiberD1_Cladding_log_v = new G4LogicalVolume(FiberD1_Cladding_solid_v, Scinti, "FiberD1_Cladding_log_v", 0, 0, 0);

	for(G4int IdReadout = 0 ; IdReadout < 4; ++IdReadout){
		for(G4int IdFib = 0 ; IdFib < 128; ++IdFib){
			if(IdFib<numFibersTopLayer1[IdReadout]){
				posFib = G4ThreeVector(FD1_startX1[IdReadout] + IdFib*spacingX ,0. , -startZ1);
			}else{
				posFib = G4ThreeVector(FD1_startX2[IdReadout] + (IdFib-numFibersTopLayer1[IdReadout])*spacingX ,0. , -startZ2);
			}

			AllPlacements.emplace_back(new G4PVPlacement(rotFib1,Sign*posFib, FiberD1_Cladding_log_x, "FiberD1_Cladding_x", FiberD1_MothVol_log_x,false,IdReadout*128+IdFib)); 
			AllPlacements.emplace_back(new G4PVPlacement(rotFib1,Sign*posFib, FiberD1_Core_log_x, "FiberD1_Core_x", FiberD1_MothVol_log_x,false,IdReadout*128+IdFib));
			AllPlacements.emplace_back(new G4PVPlacement(rotFib2,Sign*posFib, FiberD1_Cladding_log_u, "FiberD1_Cladding_u", FiberD1_MothVol_log_u,false,IdReadout*128+IdFib));
			AllPlacements.emplace_back(new G4PVPlacement(rotFib2,Sign*posFib, FiberD1_Core_log_u, "FiberD1_Core_u", FiberD1_MothVol_log_u,false,IdReadout*128+IdFib));
			AllPlacements.emplace_back(new G4PVPlacement(rotFib3,Sign*posFib, FiberD1_Cladding_log_v, "FiberD1_Cladding_v", FiberD1_MothVol_log_v,false,IdReadout*128+IdFib));
			AllPlacements.emplace_back(new G4PVPlacement(rotFib3,Sign*posFib, FiberD1_Core_log_v, "FiberD1_Core_v", FiberD1_MothVol_log_v,false,IdReadout*128+IdFib));
		}
	}
	NameDetectorsSD.push_back(FiberD1_Core_log_x->GetName());
	NameDetectorsSD.push_back(FiberD1_Core_log_u->GetName());
	NameDetectorsSD.push_back(FiberD1_Core_log_v->GetName());

	FiberD1_Core_log_x->SetVisAttributes(visAttributes_x);
	FiberD1_Core_log_u->SetVisAttributes(visAttributes_u);
	FiberD1_Core_log_v->SetVisAttributes(visAttributes_v);
	FiberD1_Cladding_log_x->SetVisAttributes(visAttributes_x);
	FiberD1_Cladding_log_u->SetVisAttributes(visAttributes_u);
	FiberD1_Cladding_log_v->SetVisAttributes(visAttributes_v);
	FiberD1_MothVol_log->SetVisAttributes(G4VisAttributes::Invisible);
	FiberD1_MothVol_log_x->SetVisAttributes(G4VisAttributes::Invisible);
	FiberD1_MothVol_log_u->SetVisAttributes(G4VisAttributes::Invisible);
	FiberD1_MothVol_log_v->SetVisAttributes(G4VisAttributes::Invisible);

// -------------------------- Second Fiber Detector --------------------------
	G4VSolid* FiberD2_MothVol_solid = new G4Box("FiberDetector2",15.*cm,20.*cm,3.*cm);
	G4VSolid* FiberD2_layerX_MothVol_solid  = new G4Box("FiberD2_layerX_solid", 15.*cm, 20.*cm, 2.*mm);
	G4VSolid* FiberD2_layerU_MothVol_solid  = new G4Box("FiberD2_layerU_solid", 15.*cm, 20.*cm, 2.*mm);
	G4VSolid* FiberD2_layerV_MothVol_solid  = new G4Box("FiberD2_layerV_solid", 15.*cm, 20.*cm, 2.*mm);
	G4LogicalVolume* FiberD2_MothVol_log = new G4LogicalVolume(FiberD2_MothVol_solid, Air, "FiberD2_log", 0, 0, 0);
	G4LogicalVolume* FiberD2_MothVol_log_x  = new G4LogicalVolume(FiberD2_layerX_MothVol_solid, Air, "FiberD2_log_x", 0, 0, 0);
	G4LogicalVolume* FiberD2_MothVol_log_u  = new G4LogicalVolume(FiberD2_layerU_MothVol_solid, Air, "FiberD2_log_u", 0, 0, 0);
	G4LogicalVolume* FiberD2_MothVol_log_v  = new G4LogicalVolume(FiberD2_layerV_MothVol_solid, Air, "FiberD2_log_v", 0, 0, 0);

	AllPlacements.emplace_back(new G4PVPlacement(0,Sign*(G4ThreeVector(0. ,0. , HypHI_FiberTracker2_posZ + Systematic_shift)-transMFLD_new), FiberD2_MothVol_log, "FiberDetector2", MFLD_log,false,0)); 
	AllPlacements.emplace_back(new G4PVPlacement(0,Sign*(G4ThreeVector(0. ,0. , posZshift[0])), FiberD2_MothVol_log_x, "FiberD2_x", FiberD2_MothVol_log,false,0)); 
	AllPlacements.emplace_back(new G4PVPlacement(0,Sign*(G4ThreeVector(0. ,0. , posZshift[1])), FiberD2_MothVol_log_u, "FiberD2_u", FiberD2_MothVol_log,false,0)); 
	AllPlacements.emplace_back(new G4PVPlacement(0,Sign*(G4ThreeVector(0. ,0. , posZshift[2])), FiberD2_MothVol_log_v, "FiberD2_v", FiberD2_MothVol_log,false,0)); 	

	std::vector<G4double> FD2_startX1 = {-70.125*mm, -39.325*mm, -4.125*mm,  31.075*mm};
	std::vector<G4double> FD2_startX2 = {-69.85*mm,  -30.25*mm,   4.95*mm,   40.15*mm };
	std::vector<G4double> numFibersTopLayer2 = {56, 64, 64, 72};

	G4VSolid* FiberD2_Core_solid_x = new G4Tubs("FiberD2_Core_solid_x",0, 0.24*mm, 10.5*cm, 0.*deg, 360.*deg);
	G4VSolid* FiberD2_Core_solid_u = new G4Tubs("FiberD2_Core_solid_u",0, 0.24*mm, (10.5/sqrt(3)*2 )*cm, 0.*deg, 360.*deg);
	G4VSolid* FiberD2_Core_solid_v = new G4Tubs("FiberD2_Core_solid_v",0, 0.24*mm, (10.5/sqrt(3)*2 )*cm, 0.*deg, 360.*deg);
	G4LogicalVolume* FiberD2_Core_log_x = new G4LogicalVolume(FiberD2_Core_solid_x,FiberCoreScinti, "FiberD2_Core_log_x", 0, 0, 0);
	G4LogicalVolume* FiberD2_Core_log_u = new G4LogicalVolume(FiberD2_Core_solid_u,FiberCoreScinti, "FiberD2_Core_log_u", 0, 0, 0);	
	G4LogicalVolume* FiberD2_Core_log_v = new G4LogicalVolume(FiberD2_Core_solid_v,FiberCoreScinti, "FiberD2_Core_log_v", 0, 0, 0);	

	G4VSolid* FiberD2_Cladding_solid_x = new G4Tubs("FiberD2_Cladding_solid_x",0.24*mm, 0.25*mm, 10.5*cm, 0.*deg, 360.*deg);
	G4VSolid* FiberD2_Cladding_solid_u = new G4Tubs("FiberD2_Cladding_solid_u",0.24*mm, 0.25*mm, (10.5/sqrt(3)*2 )*cm, 0.*deg, 360.*deg);
	G4VSolid* FiberD2_Cladding_solid_v = new G4Tubs("FiberD2_Cladding_solid_v",0.24*mm, 0.25*mm, (10.5/sqrt(3)*2 )*cm, 0.*deg, 360.*deg);
	G4LogicalVolume* FiberD2_Cladding_log_x = new G4LogicalVolume(FiberD2_Cladding_solid_x, Scinti, "FiberD2_Cladding_log_x", 0, 0, 0);	
	G4LogicalVolume* FiberD2_Cladding_log_u = new G4LogicalVolume(FiberD2_Cladding_solid_u, Scinti, "FiberD2_Cladding_log_u", 0, 0, 0);
	G4LogicalVolume* FiberD2_Cladding_log_v = new G4LogicalVolume(FiberD2_Cladding_solid_v, Scinti, "FiberD2_Cladding_log_v", 0, 0, 0);

	for(G4int IdReadout = 0 ; IdReadout < 4; ++IdReadout){
		for(G4int IdFib = 0 ; IdFib < 128; ++IdFib){
			if(IdFib<numFibersTopLayer2[IdReadout]){
				posFib = G4ThreeVector(FD2_startX1[IdReadout] + IdFib*spacingX ,0. , -startZ1);
			}else{
				posFib = G4ThreeVector(FD2_startX2[IdReadout] + (IdFib-numFibersTopLayer2[IdReadout])*spacingX ,0. , -startZ2);
			}

			AllPlacements.emplace_back(new G4PVPlacement(rotFib1,Sign*posFib, FiberD2_Cladding_log_x, "FiberD2_Cladding_x", FiberD2_MothVol_log_x,false,IdReadout*128+IdFib)); 
			AllPlacements.emplace_back(new G4PVPlacement(rotFib1,Sign*posFib, FiberD2_Core_log_x, "FiberD2_Core_x", FiberD2_MothVol_log_x,false,IdReadout*128+IdFib));
			AllPlacements.emplace_back(new G4PVPlacement(rotFib2,Sign*posFib, FiberD2_Cladding_log_u, "FiberD2_Cladding_u", FiberD2_MothVol_log_u,false,IdReadout*128+IdFib));
			AllPlacements.emplace_back(new G4PVPlacement(rotFib2,Sign*posFib, FiberD2_Core_log_u, "FiberD2_Core_u", FiberD2_MothVol_log_u,false,IdReadout*128+IdFib));
			AllPlacements.emplace_back(new G4PVPlacement(rotFib3,Sign*posFib, FiberD2_Cladding_log_v, "FiberD2_Cladding_v", FiberD2_MothVol_log_v,false,IdReadout*128+IdFib));
			AllPlacements.emplace_back(new G4PVPlacement(rotFib3,Sign*posFib, FiberD2_Core_log_v, "FiberD2_Core_v", FiberD2_MothVol_log_v,false,IdReadout*128+IdFib));
		}
	}
	NameDetectorsSD.push_back(FiberD2_Core_log_x->GetName());
	NameDetectorsSD.push_back(FiberD2_Core_log_u->GetName());
	NameDetectorsSD.push_back(FiberD2_Core_log_v->GetName());

	FiberD2_Core_log_x->SetVisAttributes(visAttributes_x);
	FiberD2_Core_log_u->SetVisAttributes(visAttributes_u);
	FiberD2_Core_log_v->SetVisAttributes(visAttributes_v);
	FiberD2_Cladding_log_x->SetVisAttributes(visAttributes_x);
	FiberD2_Cladding_log_u->SetVisAttributes(visAttributes_u);
	FiberD2_Cladding_log_v->SetVisAttributes(visAttributes_v);
	FiberD2_MothVol_log->SetVisAttributes(G4VisAttributes::Invisible);
	FiberD2_MothVol_log_x->SetVisAttributes(G4VisAttributes::Invisible);
	FiberD2_MothVol_log_u->SetVisAttributes(G4VisAttributes::Invisible);
	FiberD2_MothVol_log_v->SetVisAttributes(G4VisAttributes::Invisible);


// -------------------------- Third Fiber Detector --------------------------
	G4VSolid* FiberD3_MothVol_solid = new G4Box("FiberDetector3",15.*cm,20.*cm,3.*cm);
	G4VSolid* FiberD3_layerX_MothVol_solid  = new G4Box("FiberD3_layerX_solid", 15.*cm, 20.*cm, 2.*mm);
	G4VSolid* FiberD3_layerU_MothVol_solid  = new G4Box("FiberD3_layerU_solid", 15.*cm, 20.*cm, 2.*mm);
	G4VSolid* FiberD3_layerV_MothVol_solid  = new G4Box("FiberD3_layerV_solid", 15.*cm, 20.*cm, 2.*mm);
	G4LogicalVolume* FiberD3_MothVol_log = new G4LogicalVolume(FiberD3_MothVol_solid, Air, "FiberD3_log", 0, 0, 0);
	G4LogicalVolume* FiberD3_MothVol_log_x  = new G4LogicalVolume(FiberD3_layerX_MothVol_solid, Air, "FiberD3_log_x", 0, 0, 0);
	G4LogicalVolume* FiberD3_MothVol_log_u  = new G4LogicalVolume(FiberD3_layerU_MothVol_solid, Air, "FiberD3_log_u", 0, 0, 0);
	G4LogicalVolume* FiberD3_MothVol_log_v  = new G4LogicalVolume(FiberD3_layerV_MothVol_solid, Air, "FiberD3_log_v", 0, 0, 0);

	AllPlacements.emplace_back(new G4PVPlacement(0,Sign*(G4ThreeVector(0. ,0. , HypHI_FiberTracker3_posZ + Systematic_shift)-transMFLD_new), FiberD3_MothVol_log, "FiberDetector3", MFLD_log,false,0)); 
	AllPlacements.emplace_back(new G4PVPlacement(0,Sign*(G4ThreeVector(0. ,0. , posZshift[0])), FiberD3_MothVol_log_x, "FiberD3_x", FiberD3_MothVol_log,false,0)); 
	AllPlacements.emplace_back(new G4PVPlacement(0,Sign*(G4ThreeVector(0. ,0. , posZshift[1])), FiberD3_MothVol_log_u, "FiberD3_u", FiberD3_MothVol_log,false,0)); 
	AllPlacements.emplace_back(new G4PVPlacement(0,Sign*(G4ThreeVector(0. ,0. , posZshift[2])), FiberD3_MothVol_log_v, "FiberD3_v", FiberD3_MothVol_log,false,0)); 	

	std::vector<G4double> FD3_startX1 = {-105.325*mm,-74.525*mm ,-39.325*mm, -4.125*mm,  31.075*mm, 66.275*mm};
	std::vector<G4double> FD3_startX2 = {-105.05*mm, -65.45*mm,  -30.25*mm,   4.95*mm,   40.15*mm,  75.35*mm};
	std::vector<G4double> numFibersTopLayer3 = {56, 64, 64, 64, 64, 72};


	G4VSolid* FiberD3_Core_solid_x = new G4Tubs("FiberD3_Core_solid_x",0, 0.24*mm, (1.5*10.5)*cm, 0.*deg, 360.*deg);
	G4VSolid* FiberD3_Core_solid_u = new G4Tubs("FiberD3_Core_solid_u",0, 0.24*mm, (1.5*10.5/sqrt(3)*2 )*cm, 0.*deg, 360.*deg);
	G4VSolid* FiberD3_Core_solid_v = new G4Tubs("FiberD3_Core_solid_v",0, 0.24*mm, (1.5*10.5/sqrt(3)*2 )*cm, 0.*deg, 360.*deg);
	G4LogicalVolume* FiberD3_Core_log_x = new G4LogicalVolume(FiberD3_Core_solid_x,FiberCoreScinti, "FiberD3_Core_log_x", 0, 0, 0);
	G4LogicalVolume* FiberD3_Core_log_u = new G4LogicalVolume(FiberD3_Core_solid_u,FiberCoreScinti, "FiberD3_Core_log_u", 0, 0, 0);	
	G4LogicalVolume* FiberD3_Core_log_v = new G4LogicalVolume(FiberD3_Core_solid_v,FiberCoreScinti, "FiberD3_Core_log_v", 0, 0, 0);	

	G4VSolid* FiberD3_Cladding_solid_x = new G4Tubs("FiberD3_Cladding_solid_x",0.24*mm, 0.25*mm, (1.5*10.5)*cm, 0.*deg, 360.*deg);
	G4VSolid* FiberD3_Cladding_solid_u = new G4Tubs("FiberD3_Cladding_solid_u",0.24*mm, 0.25*mm, (1.5*10.5/sqrt(3)*2 )*cm, 0.*deg, 360.*deg);
	G4VSolid* FiberD3_Cladding_solid_v = new G4Tubs("FiberD3_Cladding_solid_v",0.24*mm, 0.25*mm, (1.5*10.5/sqrt(3)*2 )*cm, 0.*deg, 360.*deg);
	G4LogicalVolume* FiberD3_Cladding_log_x = new G4LogicalVolume(FiberD3_Cladding_solid_x, Scinti, "FiberD3_Cladding_log_x", 0, 0, 0);	
	G4LogicalVolume* FiberD3_Cladding_log_u = new G4LogicalVolume(FiberD3_Cladding_solid_u, Scinti, "FiberD3_Cladding_log_u", 0, 0, 0);
	G4LogicalVolume* FiberD3_Cladding_log_v = new G4LogicalVolume(FiberD3_Cladding_solid_v, Scinti, "FiberD3_Cladding_log_v", 0, 0, 0);

	for(G4int IdReadout = 0 ; IdReadout < 6; ++IdReadout){
		for(G4int IdFib = 0 ; IdFib < 128; ++IdFib){
			if(IdFib<numFibersTopLayer3[IdReadout]){
				posFib = G4ThreeVector(FD3_startX1[IdReadout] + IdFib*spacingX ,0. , -startZ1);
			}else{
				posFib = G4ThreeVector(FD3_startX2[IdReadout] + (IdFib-numFibersTopLayer3[IdReadout])*spacingX ,0. , -startZ2);
			}

			AllPlacements.emplace_back(new G4PVPlacement(rotFib1,Sign*posFib, FiberD3_Cladding_log_x, "FiberD3_Cladding_x", FiberD3_MothVol_log_x,false,IdReadout*128+IdFib)); 
			AllPlacements.emplace_back(new G4PVPlacement(rotFib1,Sign*posFib, FiberD3_Core_log_x, "FiberD3_Core_x", FiberD3_MothVol_log_x,false,IdReadout*128+IdFib));
			AllPlacements.emplace_back(new G4PVPlacement(rotFib2,Sign*posFib, FiberD3_Cladding_log_u, "FiberD3_Cladding_u", FiberD3_MothVol_log_u,false,IdReadout*128+IdFib));
			AllPlacements.emplace_back(new G4PVPlacement(rotFib2,Sign*posFib, FiberD3_Core_log_u, "FiberD3_Core_u", FiberD3_MothVol_log_u,false,IdReadout*128+IdFib));
			AllPlacements.emplace_back(new G4PVPlacement(rotFib3,Sign*posFib, FiberD3_Cladding_log_v, "FiberD3_Cladding_v", FiberD3_MothVol_log_v,false,IdReadout*128+IdFib));
			AllPlacements.emplace_back(new G4PVPlacement(rotFib3,Sign*posFib, FiberD3_Core_log_v, "FiberD3_Core_v", FiberD3_MothVol_log_v,false,IdReadout*128+IdFib));
		}
	}
	NameDetectorsSD.push_back(FiberD3_Core_log_x->GetName());
	NameDetectorsSD.push_back(FiberD3_Core_log_u->GetName());
	NameDetectorsSD.push_back(FiberD3_Core_log_v->GetName());

	FiberD3_Core_log_x->SetVisAttributes(visAttributes_x);
	FiberD3_Core_log_u->SetVisAttributes(visAttributes_u);
	FiberD3_Core_log_v->SetVisAttributes(visAttributes_v);
	FiberD3_Cladding_log_x->SetVisAttributes(visAttributes_x);
	FiberD3_Cladding_log_u->SetVisAttributes(visAttributes_u);
	FiberD3_Cladding_log_v->SetVisAttributes(visAttributes_v);
	FiberD3_MothVol_log->SetVisAttributes(G4VisAttributes::Invisible);
	FiberD3_MothVol_log_x->SetVisAttributes(G4VisAttributes::Invisible);
	FiberD3_MothVol_log_u->SetVisAttributes(G4VisAttributes::Invisible);
	FiberD3_MothVol_log_v->SetVisAttributes(G4VisAttributes::Invisible);
// ----------------------------- @@ -----------------------------
// ----------------------------- @@ -----------------------------




  G4VSolid* EndFMF2_box = new G4Box("EndFMF2_box",25.*cm,25.*cm,1.*mm);
  G4LogicalVolume* EndFMF2_log = new G4LogicalVolume(EndFMF2_box, Scinti, "FMF2_log", 0, 0, 0);
  
  NameDetectorsSD.push_back(EndFMF2_log->GetName());

  
  G4RotationMatrix* rotFMF2 = new G4RotationMatrix;
  // if(WasaSide==1)
  //   rotFMF2->rotateY(-180*degree);
  double FMF2_posZ = Par.Get<double>("FRS_FMF2_posZ");
  
  AllPlacements.emplace_back(new G4PVPlacement(rotFMF2,Sign*(G4ThreeVector(0., 0., FMF2_posZ+Systematic_shift)-transMFLD_new), EndFMF2_log, "FMF2_phys", MFLD_log,false,0));  
  AllPlacements.emplace_back(new G4PVPlacement(rotFMF2,Sign*(G4ThreeVector(0., 0., FMF2_posZ+Systematic_shift+1.*cm)-transMFLD_new), EndFMF2_log, "FMF2_phys1", MFLD_log,false,1));
  AllPlacements.emplace_back(new G4PVPlacement(rotFMF2,Sign*(G4ThreeVector(0., 0., FMF2_posZ+Systematic_shift+2.*cm)-transMFLD_new), EndFMF2_log, "FMF2_phys2", MFLD_log,false,2));

  G4VisAttributes *FMF2_att = new G4VisAttributes(Red);
  FMF2_att->SetForceWireframe(false);
  EndFMF2_log->SetVisAttributes(FMF2_att);
  
  
  double HypHI_EndCap_rmax = Par.Get<double>("HypHI_EndCap_maxR");
  double HypHI_EndCap_PosZ = Par.Get<double>("HypHI_EndCap_posZ");

  G4VSolid* HypHI_Endcap = new G4Tubs("HypHI_Endcap",0, HypHI_EndCap_rmax, 2*cm, 0,CLHEP::twopi);
  G4LogicalVolume* HypHI_Endcap_log  = new G4LogicalVolume(HypHI_Endcap, Air, "HypHI_Endcap_log",0,0,0);// CDCFieldMgr,0,0);

  G4RotationMatrix* rotEndCap = new G4RotationMatrix;
  // if(WasaSide==1)
  //   rotEndCap->rotateY(90*degree);

  AllPlacements.emplace_back(new G4PVPlacement(rotEndCap, Sign*(G4ThreeVector(0., 0., HypHI_EndCap_PosZ+Systematic_shift)-transMFLD_new), HypHI_Endcap_log, "HypHI_Endcap",
					       MFLD_log, false,0));
  
  //--- Visualization ---//
  HypHI_Endcap_log->SetVisAttributes(G4VisAttributes::Invisible);
  
  G4VSolid* HypHI_TrackerFwd = new G4Tubs("HypHI_TrackerFwd",BeamHoleSize*0.5, HypHI_EndCap_rmax, 2*mm, 0,CLHEP::twopi);
  G4LogicalVolume* HypHI_TrackerFwd_log = new G4LogicalVolume(HypHI_TrackerFwd,Air,"HypHI_TrackFwd_log",0,0,0);
  
  NameDetectorsSD.push_back(HypHI_TrackerFwd_log->GetName());
  
  AllPlacements.emplace_back(new G4PVPlacement(0,G4ThreeVector(0, 0, -1.5*cm),
					       HypHI_TrackerFwd_log, "HypHI_TrackerFwd0", HypHI_Endcap_log, false,0));
      
  AllPlacements.emplace_back(new G4PVPlacement(0,G4ThreeVector(0, 0, -1.*cm),
					       HypHI_TrackerFwd_log, "HypHI_TrackerFwd1", HypHI_Endcap_log, false,1));

  AllPlacements.emplace_back(new G4PVPlacement(0,G4ThreeVector(0, 0, -0.5*cm),
					       HypHI_TrackerFwd_log, "HypHI_TrackerFwd2", HypHI_Endcap_log, false,2));
            
  G4VSolid* HypHI_RPC_l = new G4Tubs("HypHI_RPC_segment_L",BeamHoleSize*0.5, HypHI_EndCap_rmax/2., 0.5*cm, -CLHEP::twopi/16.,2.*CLHEP::twopi/16.);
  G4LogicalVolume* HypHI_RPC_l_log = new G4LogicalVolume(HypHI_RPC_l, Air, "HypHI_RPC_l_log",0,0,0);
  NameDetectorsSD.push_back(HypHI_RPC_l_log->GetName());
  G4VSolid* HypHI_RPC_h = new G4Tubs("HypHI_RPC_segment_H",HypHI_EndCap_rmax/2., HypHI_EndCap_rmax, 0.5*cm, -CLHEP::twopi/16.,2.*CLHEP::twopi/16.);
  G4LogicalVolume* HypHI_RPC_h_log = new G4LogicalVolume(HypHI_RPC_h, Air, "HypHI_RPC_h_log",0,0,0);
  NameDetectorsSD.push_back(HypHI_RPC_h_log->GetName());
      
  for(int idRPC = 0; idRPC < 8 ;++idRPC)
    {
      G4RotationMatrix* rotRPC = new G4RotationMatrix;
      double rotAngle = CLHEP::twopi/8.*static_cast<double>(idRPC);
      rotRPC->rotateZ(rotAngle);
      std::string nameRPC ("HypHI_RPC_l");
      nameRPC+=std::to_string(idRPC);
      AllPlacements.emplace_back(new G4PVPlacement(rotRPC, G4ThreeVector(0, 0, 0.5*cm),
						   HypHI_RPC_l_log, nameRPC, HypHI_Endcap_log, false, idRPC));

      std::string nameRPC2 ("HypHI_RPC_h");
      nameRPC2+=std::to_string(idRPC);
      AllPlacements.emplace_back(new G4PVPlacement(rotRPC, G4ThreeVector(0, 0, 0.5*cm),
						   HypHI_RPC_h_log, nameRPC2, HypHI_Endcap_log, false, idRPC));
	  

    }

  //--- Visualization ---//
  G4VisAttributes *HypHI_RPC_att = new G4VisAttributes(Orange);
  HypHI_RPC_att->SetForceWireframe(false);
  HypHI_RPC_l_log->SetVisAttributes(HypHI_RPC_att);
  HypHI_RPC_h_log->SetVisAttributes(HypHI_RPC_att);
  G4VisAttributes *HypHI_Tracker_att = new G4VisAttributes(LightPurple);
  HypHI_Tracker_att->SetForceWireframe(false);
  HypHI_TrackerFwd_log->SetVisAttributes(HypHI_Tracker_att);



  
  

  experimentalHall_logOutRoot = world->GetLogicalVolume();
  experimentalHall_physOutRoot = world;



  
  return world;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WasaDetectorConstruction::DefineMaterials()
{ 
// Dummy, as materials are imported via VGM
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* WasaDetectorConstruction::DefineVolumes()
{
// Dummy, as geometry is imported via VGM

  return 0;
}

  //
  // end VGM demo
  //---------------------------------------------------------------------------


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WasaDetectorConstruction::ConstructSDandField()
{ 
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.

  // G4ThreeVector fieldValue = G4ThreeVector();
  // fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  // fMagFieldMessenger->SetVerboseLevel(1);

  G4SDManager *SDman = G4SDManager::GetSDMpointer();

  for(auto& CurrentName : NameDetectorsSD)
    {
      G4LogicalVolume* Det = FindVolume(CurrentName);
      G4Sol_SD_Det* SD = new G4Sol_SD_Det(CurrentName);
      //SD->Init();
      SDman->AddNewDetector(SD);
      Det->SetSensitiveDetector(SD);
    }  

  std::cout<<" Sensitive Detectors :"<<std::endl;
  for(auto NameD : NameDetectorsSD)
    std::cout<<NameD<<std::endl;

  G4Region* aDetectorRegion = new G4Region("DetectorRegion");
  
  for(auto& CurrentName : NameDetectorsSD)
    {
      G4LogicalVolume* Det = FindVolume(CurrentName);
      Det->SetRegion(aDetectorRegion);
      aDetectorRegion->AddRootLogicalVolume(Det);
    }
  std::vector<double> cutsDet (4,Par.Get<double>("DetectorRegionCut"));
  aDetectorRegion->SetProductionCuts(new G4ProductionCuts());
  aDetectorRegion->GetProductionCuts()->SetProductionCuts(cutsDet);

  // G4Region* aTargetRegion = new G4Region("TargetRegion");
  // HypHI_Target_log->SetRegion(aTargetRegion);
  // aTargetRegion->AddRootLogicalVolume(HypHI_Target_log);
  // std::vector<double> cutsTarget (4,Par.Get<double>("TargetRegionCut"));
  // aTargetRegion->SetProductionCuts(new G4ProductionCuts());
  // aTargetRegion->GetProductionCuts()->SetProductionCuts(cutsTarget);

  experimentalHall_log->SetUserLimits( new G4UserLimits(DBL_MAX,2*m,10*s,0.,0.) );  
  

  G4double fCDField    = 0.   *tesla;  

  if(Par.IsAvailable("Field_CDS_Bz"))
    fCDField = Par.Get<double>("Field_CDS_Bz");

  G4ThreeVector fCDC   (0.0, 0.0, fCDField);  

  fMagneticField = new G4SolSimpleMagneticField();
  fMagneticField->SetField(fCDC);
  
  fEquation = new G4Mag_UsualEqRhs(fMagneticField);
  //fStepper = new G4ClassicalRK4( fEquation );
  fStepper = new G4NystromRK4(fEquation);

  
  fFieldMgr = new G4FieldManager();
  //fFieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fFieldMgr->SetDetectorField(fMagneticField);
  //fFieldMgr->CreateChordFinder(fMagneticField);
  fChordFinder = new G4ChordFinder(fMagneticField, 1.e-2, fStepper);
  fFieldMgr->SetChordFinder(fChordFinder);
  
  G4bool forceToAllDaughters = true;
  INNER_log->SetFieldManager(fFieldMgr,forceToAllDaughters);
  //CDS_endcap_log->SetFieldManager(fFieldMgr,forceToAllDaughters);
  // if(DoModHypHI)
  //   HypHI_InTracker_log->SetFieldManager(fFieldMgr,forceToAllDaughters);
  
  G4AutoDelete::Register(fMagneticField);
  //G4AutoDelete::Register(fFieldMgr);
  G4AutoDelete::Register(fEquation);
  G4AutoDelete::Register(fStepper);
  G4AutoDelete::Register(fChordFinder);
  
  
  
  // G4FieldManager* fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();

  // SksFieldinRoot = new WasaMagneticField("./field/SksPlusField_1000A.root");
  
  // double SKS_field1_value = Par.Get_Geometry_SKSField1();
  // double SKS_field2_value = Par.Get_Geometry_SKSField2();

  // std::cout<<" Sks field :"<<SKS_field1_value<<" "<<SKS_field2_value<<std::endl;

  // SksFieldinRoot->InitField(-SKS_field1_value/0.8016,true,SKS_field2_value);
  // fieldMgr->SetDetectorField(SksFieldinRoot);

  // double pointTest[3] = {0.,0.,90.*cm};
  // double fieldTest[3] = {0.,0.,0.};
  // SksFieldinRoot->GetFieldValue(pointTest,fieldTest);
  // std::cout<<" Field 1 :"<<fieldTest[1]/tesla<<" T | ";
  // pointTest[2] = 300.*cm;
  // SksFieldinRoot->GetFieldValue(pointTest,fieldTest);
  // std::cout<<" Field 2 :"<<fieldTest[1]/tesla<<" T"<<std::endl;
  
  
  // G4Mag_UsualEqRhs* pMagFldEquation = new G4Mag_UsualEqRhs(SksFieldinRoot);
  // //G4MagIntegratorStepper* fStepper = new G4NystromRK4( pMagFldEquation );

  // G4MagIntegratorStepper* fStepper = new G4SimpleHeum(pMagFldEquation,8);

  // G4ChordFinder *pChordFinder = new G4ChordFinder(SksFieldinRoot,1.e-2*mm, fStepper);
  // pChordFinder->SetDeltaChord(1.0e-1*mm);
  
  // fieldMgr->SetChordFinder(pChordFinder);

  // // G4LogicalVolume* SetupLV = FindVolume("Setup");
  // // SetupLV->SetFieldManager(fieldMgr,true);

  // // Register the field messenger for deleting
  // G4AutoDelete::Register(SksFieldinRoot);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* WasaDetectorConstruction::FindVolPhys(const G4String& name)
{
  G4PhysicalVolumeStore* pvStore = G4PhysicalVolumeStore::GetInstance();
  return pvStore->GetVolume(name);
}
