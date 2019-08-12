// --------------------------------------------------------------------
// Implementation of the WasaSimpleDetectorConstruction class
// Created by C.Rappold (c.rappold@gsi.de)
//---------------------------------------------------------------------

// ====================================================================
//    WasaSimpleDetectorConstruction.cc
//
//
// ====================================================================

#include "WasaSimpleDetectorConstruction.hh"
//#include "KnuclRegionInformation.hh"
#include "G4Box.hh"
#include "G4Color.hh"
#include "G4Element.hh"
#include "G4GenericMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4ProductionCuts.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4TwistedTubs.hh"
#include "G4UserLimits.hh"
#include "G4VSensitiveDetector.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"
//#include "KnuclMaterialManager.hh"
//#include "KnuclFieldSetup.hh"
//#include "KnuclCommon.h"
#include "G4SolSensitiveD.hh"
#include "G4SolSimpleMagneticField.hh"

//#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
//#include "G4SolidStore.hh"

#include "G4AutoDelete.hh"
#include "G4FieldManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4TransportationManager.hh"

#include <iostream>
#include <set>

G4ThreadLocal G4SolSimpleMagneticField* WasaSimpleDetectorConstruction::fMagneticField = 0;
G4ThreadLocal G4FieldManager* WasaSimpleDetectorConstruction::fFieldMgr                = 0;

// ====================================================================
//
// class description
//
// ====================================================================

//////////////////////////////////////////////////////
WasaSimpleDetectorConstruction::WasaSimpleDetectorConstruction(const G4SolConfig& conf)
    : G4SolVDetectorConstruction(conf), experimentalHall_box(0), experimentalHall_log(0), experimentalHall_phys(0),
      CD_tube(0), CD_log(0), CD_phys(0)

//////////////////////////////////////////////////////
{

  TargetMaterialID    = 1;
  ChamberGasSelection = "ArCO2Methan89-10-1";

  FieldInCD = 0.5;

  TargetLength = Par.Get<double>("Target_Size"); // 1.0*cm;

  DoCD       = true;
  DoPSB      = true;
  DoEndFMF2  = true;
  DoModHypHI = true;

  DoOnlySense     = 0;
  ReductionRadius = 1.;

  // if(Par.IsAvailable("ReductionFactor"))
  //   {
  //     ReductionRadius = Par.Get<double>("ReductionFactor");
  //     if(std::abs(ReductionRadius -1.)>1e-4)
  // 	{
  // 	  CD_Rmin *= ReductionRadius;
  // 	  CD_Rmax *= ReductionRadius;
  // 	  for(auto& layerR : CDC_RADIUS)
  // 	    layerR*= ReductionRadius;
  // 	}
  //   }
  if(Par.IsAvailable("ConvertRoot"))
    DoForRoot = true;
  else
    DoForRoot = false;

  DefineCommands();
}

//////////////////////////////////////////////////////
WasaSimpleDetectorConstruction::~WasaSimpleDetectorConstruction()
//////////////////////////////////////////////////////
{
  // if (fEmFieldSetupKURAMA) delete fEmFieldSetupKURAMA;
  // if (fEmFieldSetupCDC)    delete fEmFieldSetupCDC;
}

void WasaSimpleDetectorConstruction::ConstructMaterials()
{
  G4NistManager* materialMgr = G4NistManager::Instance();
  // G4Material* He4 = materialMgr->FindOrBuildMaterial("G4_He");
  G4Material* Argon            = materialMgr->FindOrBuildMaterial("G4_Ar");
  double density               = 2.67 * mg / cm3;
  std::vector<G4String> nameEl = {"C", "H"};
  std::vector<G4int> atoms     = {4, 10};
  G4Material* isobutane        = materialMgr->ConstructNewMaterial("isoC4H10", nameEl, atoms, density, true, kStateGas);

  density     = 1.74514 * mg / cm3;
  ArIsoButane = std::unique_ptr<G4Material>(
      new G4Material("ArIsoButane", density, 2)); // std::make_unique<G4Material>("ArIsoButane", density,2);
  ArIsoButane->AddMaterial(Argon, 0.90);
  ArIsoButane->AddMaterial(isobutane, 0.10);

  // G4Material* HeEthan_50_50 = materialMgr-> GetMaterial("He_Ethan_50_50");
  // G4Material* ArCH4_90_10  = materialMgr-> GetMaterial("ArCH4_90_10");
  // G4Material* ArEthan_50_50 = materialMgr-> GetMaterial("ArEthan_50_50");
  // G4Material* HeIsobutane_80_20 = materialMgr-> GetMaterial("He_Isobutane_80_20");
  // G4Material* ArEthan_80_20 = materialMgr-> GetMaterial("ArEthan_80_20");
  // G4Material* ArCO2Methan_89_10_1 = materialMgr-> GetMaterial("ArCO2Methan_89_10_1");
  // G4Material* CFRP = materialMgr-> GetMaterial("CFRP");
  // G4Material* LHe3 = materialMgr->FindOrBuildMaterial("LHelium-3");
  // G4Material* D2Gas = materialMgr->FindOrBuildMaterial("DeuteriumGas");
  // G4Material* SUS = materialMgr->FindOrBuildMaterial("SUS");
}

//////////////////////////////////////////////////////
G4VPhysicalVolume* WasaSimpleDetectorConstruction::Construct()
//////////////////////////////////////////////////////
{

  // ==============================================================
  // materials
  // ==============================================================

  // KnuclMaterialManager* materialMgr = new KnuclMaterialManager();
  ConstructMaterials();
  G4NistManager* materialMgr = G4NistManager::Instance();

  G4Material* Air = materialMgr->FindOrBuildMaterial("G4_AIR");
  // G4Material* Vacuum = materialMgr->FindOrBuildMaterial("G4_Galactic");
  // G4Material* HeGas = materialMgr->FindOrBuildMaterial("G4_He");

  G4cout << "!!! " << ChamberGasSelection << " was selected as MDC " << G4endl;

  G4cout << "Detectors :" << DoCD << " " << DoPSB << "\n";

  // ==============================================================
  // experimental hall (world volume)
  //   --- beam line along z axis ---
  // ==============================================================

  G4double expHall_x   = 10.0 * m;
  G4double expHall_y   = 2.0 * m;
  G4double expHall_z   = 20.0 * m;
  experimentalHall_box = new G4Box("expHall_box", expHall_x, expHall_y, expHall_z);
  experimentalHall_log = new G4LogicalVolume(experimentalHall_box, Air, "expHall_log", 0, 0, 0);
  if(DoForRoot)
    experimentalHall_logOutRoot = new G4LogicalVolume(experimentalHall_box, Air, "expHall_logR", 0, 0, 0);

  experimentalHall_phys = new G4PVPlacement(0, G4ThreeVector(), experimentalHall_log, "expHall", 0, false, 0);
  if(DoForRoot)
    experimentalHall_physOutRoot =
        new G4PVPlacement(0, G4ThreeVector(), experimentalHall_logOutRoot, "expHallR", 0, false, 0);

  experimentalHall_log->SetVisAttributes(G4VisAttributes::Invisible);

  // ==============================================================
  // CDS
  // ==============================================================

  const G4double TargetPosX = Par.Get<double>("Target_PosX");
  const G4double TargetPosY = Par.Get<double>("Target_PosY");
  const G4double TargetPosZ = Par.Get<double>("Target_PosZ");
  // const G4double RelativePos_CDSTarget = Par.Get<double>("CDS_RelativePosTarget");
  const G4double AbsPos_CDTarget = Par.Get<double>("CD_AbsPosTarget");
  const G4double cdsPos_x        = TargetPosX;                   // 0.0*m;
  const G4double cdsPos_y        = TargetPosY;                   // 0.0*m;
  const G4double cdsPos_z        = TargetPosZ + AbsPos_CDTarget; //-RelativePos_CDSTarget* CD_Z;//0.0*m;

  if(DoCD == true)
    ConstructCD(CD_Rmax, CD_Z, cdsPos_x, cdsPos_y, cdsPos_z);

  if(DoModHypHI == true)
    ConstructDownTracker(CD_Z, AbsPos_CDTarget, TargetPosX, TargetPosY, TargetPosZ);

  if(DoPSB == true)
    ConstructPSB(); // cdsPos_x, cdsPos_y, cdsPos_z);

  if(DoEndFMF2 == true)
    ConstructEndFMF2();

  G4cout << "WasaSimpleDetectorConstruction completed" << G4endl;

  return experimentalHall_phys;
}

void WasaSimpleDetectorConstruction::ConstructSDandField()
{

  // ==============================================================
  // Definition of sensitive detectors
  // ==============================================================

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  for(auto& CurrentName : NameDetectorsSD)
    {
      G4LogicalVolume* Det = FindVolume(CurrentName);
      G4Sol_SD_Det* SD     = new G4Sol_SD_Det(CurrentName);
      // SD->Init();
      SDman->AddNewDetector(SD);
      Det->SetSensitiveDetector(SD);
    }

  std::cout << " Sensitive Detectors :" << std::endl;
  for(auto NameD : NameDetectorsSD)
    std::cout << NameD << std::endl;

  G4Region* aDetectorRegion = new G4Region("DetectorRegion");

  for(auto& CurrentName : NameDetectorsSD)
    {
      G4LogicalVolume* Det = FindVolume(CurrentName);
      Det->SetRegion(aDetectorRegion);
      aDetectorRegion->AddRootLogicalVolume(Det);
    }
  std::vector<double> cutsDet(4, Par.Get<double>("DetectorRegionCut"));
  aDetectorRegion->SetProductionCuts(new G4ProductionCuts());
  aDetectorRegion->GetProductionCuts()->SetProductionCuts(cutsDet);

  G4Region* aTargetRegion = new G4Region("TargetRegion");
  HypHI_Target_log->SetRegion(aTargetRegion);
  aTargetRegion->AddRootLogicalVolume(HypHI_Target_log);
  std::vector<double> cutsTarget(4, Par.Get<double>("TargetRegionCut"));
  aTargetRegion->SetProductionCuts(new G4ProductionCuts());
  aTargetRegion->GetProductionCuts()->SetProductionCuts(cutsTarget);

  experimentalHall_log->SetUserLimits(new G4UserLimits(DBL_MAX, 2 * m, 10 * s, 0., 0.));

  // ==============================================================
  // magnetic field
  // ==============================================================

  G4double fCDField = -FieldInCD * tesla;

  if(Par.IsAvailable("Field_CDS_Bz"))
    fCDField = Par.Get<double>("Field_CDS_Bz");

  G4ThreeVector fCDC(0.0, 0.0, fCDField);

  fMagneticField = new G4SolSimpleMagneticField();
  fFieldMgr      = new G4FieldManager();
  fFieldMgr->SetDetectorField(fMagneticField);
  fFieldMgr->CreateChordFinder(fMagneticField);

  G4bool forceToAllDaughters = true;
  CD_log->SetFieldManager(fFieldMgr, forceToAllDaughters);
  // CDS_endcap_log->SetFieldManager(fFieldMgr,forceToAllDaughters);
  // if(DoModHypHI)
  //   HypHI_InTracker_log->SetFieldManager(fFieldMgr,forceToAllDaughters);

  G4AutoDelete::Register(fMagneticField);
  G4AutoDelete::Register(fFieldMgr);
}

void WasaSimpleDetectorConstruction::ConstructCD(G4double cd_rmax, G4double cd_z, G4double cdsPos_x, G4double cdsPos_y,
                                                 G4double cdsPos_z)
{
  G4NistManager* materialMgr = G4NistManager::Instance();

  G4Material* Air = materialMgr->FindOrBuildMaterial("G4_AIR");
  // G4Material* Vacuum = materialMgr->FindOrBuildMaterial("G4_Galactic");
  // G4Material* Fe = materialMgr->FindOrBuildMaterial("G4_Fe");
  // G4Material* Mylar = materialMgr->FindOrBuildMaterial("G4_MYLAR");
  // G4Material* W = materialMgr->FindOrBuildMaterial("G4_W");
  G4Material* Al = materialMgr->FindOrBuildMaterial("G4_Al");

  G4Material* ChamberGas = ArIsoButane.get();

  //*************************//
  //*** CD virtual holder ***//
  //*************************//

  // G4FieldManager* CDCFieldMgr = fEmFieldSetupCDC->GetFieldManager();
  // double cd_rmin = CD_Rmin*0.95*mm;
  CD_tube = new G4Tubs("CD_tube", 0., cd_rmax, 0.5 * cd_z, 0.0, CLHEP::twopi);
  CD_log  = new G4LogicalVolume(CD_tube, Air, "CD_log", 0, 0, 0); // CDCFieldMgr,0,0);
  CD_phys =
      new G4PVPlacement(0, G4ThreeVector(cdsPos_x, cdsPos_y, cdsPos_z), CD_log, "CD", experimentalHall_log, false, 0);
  //--- Visualization ---//
  CD_log->SetVisAttributes(G4VisAttributes::Invisible);
  if(DoForRoot)
    {
      CD_logOutRoot  = new G4LogicalVolume(CD_tube, Air, "CD_logR", 0, 0, 0); // CDCFieldMgr,0,0);
      CD_physOutRoot = new G4PVPlacement(0, G4ThreeVector(cdsPos_x, cdsPos_y, cdsPos_z), CD_logOutRoot, "CDR",
                                         experimentalHall_logOutRoot, false, 0);
      CD_logOutRoot->SetVisAttributes(G4VisAttributes::Invisible);
    }

    //****************//
    //*** CD Yoke ***//
    //****************//

#ifdef NOYET_YOKE
  G4double cdsyoke_rmax = CDSYOKE_RMAX * mm;
  G4double cdsyoke_rmin = CDSYOKE_RMIN * mm;

  G4Tubs* CDSYoke_tube         = new G4Tubs("CDS_tube", cdsyoke_rmin, cdsyoke_rmax, cd_z, 0.0, CLHEP::twopi);
  G4LogicalVolume* CDSYoke_log = new G4LogicalVolume(CDSYoke_tube, Fe, "CDSYoke_log", 0, 0, 0);

  // G4PVPlacement* CDSYoke_phys =
  AllPlacements.emplace_back(new G4PVPlacement(0, G4ThreeVector(cdsPos_x, cdsPos_y, cdsPos_z), CDSYoke_log, "CDS_Yoke",
                                               experimentalHall_log, false, 0));
  if(DoForRoot)
    AllPlacements.emplace_back(new G4PVPlacement(0, G4ThreeVector(cdsPos_x, cdsPos_y, cdsPos_z), CDSYoke_log,
                                                 "CDS_YokeR", experimentalHall_logOutRoot, false, 0));

  //--- Visualization ---//
  G4VisAttributes* CDSYoke_att = new G4VisAttributes(Gray);
  CDSYoke_att->SetForceWireframe(false);
  CDSYoke_log->SetVisAttributes(CDSYoke_att);

  //**********************//
  //*** CDS endcap ***//
  //**********************//
  G4double cds_endcap_rmax = CDSYOKE_RMAX * mm;
  G4double cds_endcap_rmin = CDS_ENDCAP_RMIN * mm;
  G4double cds_endcap_z    = CDS_ENDCAP_Z / 2.0 * mm;

  G4VSolid* CDS_endcap_tube =
      new G4Tubs("CDS_endcap_tube", cds_endcap_rmin, cds_endcap_rmax, cds_endcap_z, 0.0, CLHEP::twopi);
  CDS_endcap_log = new G4LogicalVolume(CDS_endcap_tube, Fe, "CDS_endcap_log", 0, 0, 0); //, CDCFieldMgr,0,0);

  // G4PVPlacement* CDS_endcap_phys[2];
  // CDS_endcap_phys[0] =
  AllPlacements.emplace_back(new G4PVPlacement(0, G4ThreeVector(cdsPos_x, cdsPos_y, cdsPos_z - (cd_z + cds_endcap_z)),
                                               CDS_endcap_log, "CDS_endcap_up", experimentalHall_log, false, 0));
  if(DoForRoot)
    AllPlacements.emplace_back(new G4PVPlacement(0, G4ThreeVector(cdsPos_x, cdsPos_y, cdsPos_z - (cd_z + cds_endcap_z)),
                                                 CDS_endcap_log, "CDS_endcap_upR", experimentalHall_logOutRoot, false,
                                                 0));
    // CDS_endcap_phys[1] =
#endif

#ifdef NOYET_ENDCAP
  if(DoModHypHI == false)
    {
      AllPlacements.emplace_back(new G4PVPlacement(0,
                                                   G4ThreeVector(cdsPos_x, cdsPos_y, cdsPos_z + (cd_z + cds_endcap_z)),
                                                   CDS_endcap_log, "CDS_endcap_down", experimentalHall_log, false, 0));
    }
  else
    {

      double HypHI_rmax = cd_rmax * 1.2;

      HypHI_Endcap      = new G4Tubs("HypHI_Endcap", 0, HypHI_rmax, 10 * cm, 0, CLHEP::twopi);
      HypHI_Endcap_log  = new G4LogicalVolume(HypHI_Endcap, Air, "HypHI_Endcap_log", 0, 0, 0); // CDCFieldMgr,0,0);
      HypHI_Endcap_phys = new G4PVPlacement(0, G4ThreeVector(cdsPos_x, cdsPos_y, cdsPos_z + (cd_z + 12. * cm)),
                                            HypHI_Endcap_log, "HypHI_Endcap", experimentalHall_log, false, 0);
      if(DoForRoot)
        AllPlacements.emplace_back(new G4PVPlacement(0, G4ThreeVector(cdsPos_x, cdsPos_y, cdsPos_z + (cd_z + 12. * cm)),
                                                     HypHI_Endcap_log, "HypHI_EndcapR", experimentalHall_logOutRoot,
                                                     false, 0));

      //--- Visualization ---//
      HypHI_Endcap_log->SetVisAttributes(G4VisAttributes::Invisible);

      G4VSolid* HypHI_TrackerFwd = new G4Tubs("HypHI_TrackerFwd", cds_endcap_rmin, HypHI_rmax, 1 * cm, 0, CLHEP::twopi);
      HypHI_TrackerFwd_log       = new G4LogicalVolume(HypHI_TrackerFwd, Air, "HypHI_TrackFwd_log", 0, 0, 0);

      NameDetectorsSD.push_back(HypHI_TrackerFwd_log->GetName());

      AllPlacements.emplace_back(new G4PVPlacement(0, G4ThreeVector(0, 0, -9 * cm), HypHI_TrackerFwd_log,
                                                   "HypHI_TrackerFwd0", HypHI_Endcap_log, false, 0));

      AllPlacements.emplace_back(new G4PVPlacement(0, G4ThreeVector(0, 0, -6 * cm), HypHI_TrackerFwd_log,
                                                   "HypHI_TrackerFwd1", HypHI_Endcap_log, false, 1));

      AllPlacements.emplace_back(new G4PVPlacement(0, G4ThreeVector(0, 0, -3 * cm), HypHI_TrackerFwd_log,
                                                   "HypHI_TrackerFwd2", HypHI_Endcap_log, false, 2));

      G4VSolid* HypHI_RPC_l = new G4Tubs("HypHI_RPC_segment_L", cds_endcap_rmin, HypHI_rmax / 2., 5 * cm,
                                         -CLHEP::twopi / 16., 2. * CLHEP::twopi / 16.);
      HypHI_RPC_l_log       = new G4LogicalVolume(HypHI_RPC_l, Air, "HypHI_RPC_l_log", 0, 0, 0);
      NameDetectorsSD.push_back(HypHI_RPC_l_log->GetName());
      G4VSolid* HypHI_RPC_h = new G4Tubs("HypHI_RPC_segment_H", HypHI_rmax / 2., HypHI_rmax, 5 * cm,
                                         -CLHEP::twopi / 16., 2. * CLHEP::twopi / 16.);
      HypHI_RPC_h_log       = new G4LogicalVolume(HypHI_RPC_h, Air, "HypHI_RPC_h_log", 0, 0, 0);
      NameDetectorsSD.push_back(HypHI_RPC_h_log->GetName());

      for(int idRPC = 0; idRPC < 8; ++idRPC)
        {
          G4RotationMatrix* rotRPC = new G4RotationMatrix;
          double rotAngle          = CLHEP::twopi / 8. * static_cast<double>(idRPC);
          rotRPC->rotateZ(rotAngle);
          std::string nameRPC("HypHI_RPC_l");
          nameRPC += std::to_string(idRPC);
          AllPlacements.emplace_back(new G4PVPlacement(rotRPC, G4ThreeVector(0, 0, 4. * cm), HypHI_RPC_l_log, nameRPC,
                                                       HypHI_Endcap_log, false, idRPC));

          std::string nameRPC2("HypHI_RPC_h");
          nameRPC2 += std::to_string(idRPC);
          AllPlacements.emplace_back(new G4PVPlacement(rotRPC, G4ThreeVector(0, 0, 4. * cm), HypHI_RPC_h_log, nameRPC2,
                                                       HypHI_Endcap_log, false, idRPC));
        }

      //--- Visualization ---//
      G4VisAttributes* HypHI_RPC_att = new G4VisAttributes(Orange);
      HypHI_RPC_att->SetForceWireframe(false);
      HypHI_RPC_l_log->SetVisAttributes(HypHI_RPC_att);
      HypHI_RPC_h_log->SetVisAttributes(HypHI_RPC_att);
      G4VisAttributes* HypHI_Tracker_att = new G4VisAttributes(LightPurple);
      HypHI_Tracker_att->SetForceWireframe(false);
      HypHI_TrackerFwd_log->SetVisAttributes(HypHI_Tracker_att);
    }

  //--- Visualization ---//
  G4VisAttributes* CDS_endcap_att = new G4VisAttributes(Gray);
  CDS_endcap_att->SetForceWireframe(false);
  CDS_endcap_log->SetVisAttributes(CDS_endcap_att);

#endif

  //****************//
  //*** Solenoid ***//
  //****************//

  Solenoid_tube =
      new G4Tubs("Solenoid_tube", Solenoid_RadiusMin, Solenoid_RadiusMax, Solenoid_Length * 0.5, 0.0, CLHEP::twopi);
  Solenoid_log  = new G4LogicalVolume(Solenoid_tube, Al, "Solenoid_log", 0, 0, 0); // CDCFieldMgr,0,0);
  Solenoid_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), Solenoid_log, "Solenoid_phys", CD_log, false, 0);

  //***********//
  //*** MDC ***//
  //***********//
  G4double mdc_rmin = MDC_RadiusMin * 0.95 * mm;
  G4double mdc_rmax = MDC_RadiusMax * 1.01 * mm;
  G4double mdc_z    = Solenoid_Length * .5 * mm;

  G4VSolid* MDC_body_tube       = new G4Tubs("MDC_body_tube", mdc_rmin, mdc_rmax, mdc_z, 0.0, CLHEP::twopi);
  G4LogicalVolume* MDC_body_log = new G4LogicalVolume(MDC_body_tube, ChamberGas, "MDC_body_log", 0, 0, 0);

  AllPlacements.emplace_back(
      new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), MDC_body_log, "MDC_body_phys", CD_log, false, 0));
  G4LogicalVolume* MDC_body_logOutRoot = nullptr;
  if(DoForRoot)
    {
      MDC_body_logOutRoot = new G4LogicalVolume(MDC_body_tube, ChamberGas, "MDC_body_logR", 0, 0, 0);
      AllPlacements.emplace_back(new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), MDC_body_logOutRoot, "MDC_body_physR",
                                                   CD_logOutRoot, false, 0));
    }

  G4VSolid* MDC_Set[17];
  G4VSolid* MDC_tube[17];
  // G4VSolid* CDC_twist_tube[17];

  G4LogicalVolume* MDC_log[17];
  G4LogicalVolume* MDC_Setlog[17];
  G4LogicalVolume* MDC_SetlogOutRoot[17];

  std::vector<std::vector<G4PVPlacement*> > MDC_phys(MDC_NbLayer); //[3][j];

  G4int n = 0;

  G4double rmin = 0.0;
  G4double rmax = MDC_RadiusMin - MDC_Cell_Thinkness[0] / 2.;

  // Layer 4mm
  for(G4int i = 0; i < MDC_NbLayer; i++)
    {
      rmin = rmax + MDC_Cell_Thinkness[i] / 2. + MDC_LayerSpacing;
      if(i == 0)
        rmin -= MDC_LayerSpacing;
      rmax = rmin + MDC_Cell_Thinkness[i] / 2.;

      double Nb_MDC_cell = std::floor(CLHEP::pi * (rmin + rmax) / MDC_Cell_Thinkness[i]);
      int Nb_cell        = static_cast<int>(Nb_MDC_cell);
      MDC_phys[i].resize(Nb_cell, nullptr);

      std::string name_sol("CDC_sol");
      name_sol += std::to_string(i);
      std::string name_log("MDC_log");
      name_log += std::to_string(i);

      std::string nameSet("MDC_Set");
      nameSet += std::to_string(i);
      std::string nameSetLog("MDC_SetLog");
      nameSetLog += std::to_string(i);
      std::string nameSetPhys("MDC_SetPhys");
      nameSetPhys += std::to_string(i);

      MDC_Set[i]    = new G4Tubs(nameSet, rmin - 0.95 * mm, rmax + 0.95 * mm, mdc_z, 0, CLHEP::twopi);
      MDC_Setlog[i] = new G4LogicalVolume(MDC_Set[i], ChamberGas, nameSetLog, 0, 0, 0);
      // MDC_Setlog[i]->SetVisAttributes(G4VisAttributes::Invisible);
      AllPlacements.emplace_back(
          new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), MDC_Setlog[i], nameSetPhys, MDC_body_log, false, i));

      MDC_tube[i] =
          new G4Tubs(name_sol, rmin, rmax, mdc_z, -0.5 * CLHEP::twopi / Nb_MDC_cell, CLHEP::twopi / Nb_MDC_cell);
      MDC_log[i] = new G4LogicalVolume(MDC_tube[i], ChamberGas, name_log, 0, 0, 0);

      // n=0;
      for(G4int j = 0; j < Nb_cell; j++)
        {
          // 	  //sprintf(name_phy,"CDC_phys%2.2d%03d",i+1,j+1);
          std::string name_phy("CDC_phys");
          name_phy += std::to_string(i + 1);
          name_phy += "_";
          name_phy += std::to_string(j + 1);

          G4ThreeVector xyzCounter(0., 0., 0.);
          G4RotationMatrix* rotCounterMDCcell = new G4RotationMatrix;
          rotCounterMDCcell->rotateZ(CLHEP::twopi / Nb_MDC_cell * (j + 0.5 * MDC_Cell_XXP[i]));
          G4Transform3D posCounterMDCcell(*rotCounterMDCcell, xyzCounter);

          MDC_phys[i][j] = new G4PVPlacement(posCounterMDCcell, MDC_log[i], name_phy, MDC_Setlog[i], false, n);
          // 	  n++;
        }
    }

  for(auto Vlog : MDC_log)
    NameDetectorsSD.push_back(Vlog->GetName());

  if(DoForRoot)
    {
      for(int i = 0; i < MDC_NbLayer; ++i)
        {
          std::string nameSetLog("MDC_SetLogR_");
          std::string nameSetPhys("MDC_SetPhysR_");
          nameSetLog += std::to_string(i);
          nameSetPhys += std::to_string(i);
          MDC_SetlogOutRoot[i] = new G4LogicalVolume(MDC_Set[i], ChamberGas, nameSetLog, 0, 0, 0);
          // MDC_Setlog[i]->SetVisAttributes(G4VisAttributes::Invisible);
          AllPlacements.emplace_back(new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), MDC_SetlogOutRoot[i], nameSetPhys,
                                                       MDC_body_logOutRoot, false, i));
        }
    }

  //--- Visualization ---//
  G4VisAttributes* MDC_body_att = new G4VisAttributes(LightBlue);
  MDC_body_att->SetForceWireframe(false);
  MDC_body_log->SetVisAttributes(MDC_body_att);
  // MDC_att->SetForceSolid(true);

  for(G4int i = 0; i < MDC_NbLayer; i++)
    {
      MDC_Setlog[i]->SetVisAttributes(G4VisAttributes::Invisible);

      G4VisAttributes* MDC_att = new G4VisAttributes(ColorCDC[i]);
      MDC_att->SetForceWireframe(false);
      // MDC_log[i]->SetSensitiveDetector(chamberSD);
      // MDC_log[i]-> SetUserLimits( new G4UserLimits(1.0*mm) );
      // MDC_log[i]->SetVisAttributes(G4VisAttributes::Invisible);
      MDC_log[i]->SetVisAttributes(MDC_att);
    }
}

void WasaSimpleDetectorConstruction::ConstructDownTracker(G4double cd_z, G4double AbsPos, G4double TargetPos_x,
                                                          G4double TargetPos_y, G4double TargetPos_z)
{
  G4NistManager* materialMgr = G4NistManager::Instance();

  // G4Material* Air = materialMgr->FindOrBuildMaterial("G4_AIR");
  G4Material* Vacuum = materialMgr->FindOrBuildMaterial("G4_Galactic");
  G4Material* Si     = materialMgr->FindOrBuildMaterial("G4_Si"); //"Plastic");

  G4Material* Scinti = materialMgr->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"); // G4_POLYETHYLENE");//"Plastic");

  G4Material* Carbon = materialMgr->FindOrBuildMaterial("G4_C"); //"Plastic");

  G4VisAttributes* Si_att = new G4VisAttributes(Pink);

  const double lengthDownstream = std::min(Par.Get<double>("HypHI_Downstream_Size"), AbsPos - cd_z * 0.5); // 0.8*m;
  const double shift_log =
      -0.45 * lengthDownstream; // Par.Get<double>("HypHI_Downstream_Shift");//-0.45*m; //RelativePos*cd_z*0.5

  std::cout << " DownTracker :" << lengthDownstream << " " << shift_log << "\n";
  // double mdc_rmin = 0;
  HypHI_InTracker      = new G4Tubs("HypHI_InTracker", 0, 1. * m, .5 * lengthDownstream, 0, CLHEP::twopi);
  HypHI_InTracker_log  = new G4LogicalVolume(HypHI_InTracker, Vacuum, "HypHI_InTracker_log", 0, 0, 0);
  HypHI_InTracker_phys = new G4PVPlacement(0, G4ThreeVector(TargetPos_x, TargetPos_y, TargetPos_z - shift_log),
                                           HypHI_InTracker_log, "HypHI_InTracker_Phys", experimentalHall_log, false, 0);
  if(DoForRoot)
    AllPlacements.emplace_back(new G4PVPlacement(0, G4ThreeVector(TargetPos_x, TargetPos_y, TargetPos_z - shift_log),
                                                 HypHI_InTracker_log, "HypHI_InTracker_PhysR",
                                                 experimentalHall_logOutRoot, false, 0));

  //--- Visualization ---//
  if(!Par.IsAvailable("HypHI_InnerTrackerBox_Visible"))
    HypHI_InTracker_log->SetVisAttributes(G4VisAttributes::Invisible);

  HypHI_Target      = new G4Box("HypHI_Target", TargetLength, TargetLength, TargetLength);
  HypHI_Target_log  = new G4LogicalVolume(HypHI_Target, Carbon, "HypHI_Target_log", 0, 0, 0);
  HypHI_Target_phys = new G4PVPlacement(0, G4ThreeVector(0, 0, TargetPos_z + shift_log), HypHI_Target_log,
                                        "HypHI_Target_Phys", HypHI_InTracker_log, false, 0);

  G4int nb_panel = 16;
  G4VSolid* HypHI_SiliciumSeg =
      new G4Tubs("HypHI_SiSeg", 0. * cm, 4. * cm, 3. * mm, -CLHEP::twopi / static_cast<double>(2 * nb_panel),
                 2. * CLHEP::twopi / static_cast<double>(2 * nb_panel));

  std::vector<double> posZ = {3. * cm, 4. * cm, 5. * cm, 6. * cm};
  if(Par.IsAvailable("HypHI_InnerTracker_Spec"))
    {
      double PosZInTracker = Par.Get<double>("HypHI_InnerTracker_PosZ");
      int NbInTracker      = Par.Get<int>("HypHI_InnerTracker_Nb");
      double DistInTracker = Par.Get<double>("HypHI_InnerTracker_Spacing");
      posZ.resize(NbInTracker);
      for(size_t idL = 0; idL < posZ.size(); ++idL)
        posZ[idL] = TargetLength + TargetPos_z + shift_log + PosZInTracker + static_cast<double>(idL) * DistInTracker;
    }

  for(size_t idLayer = 0; idLayer < posZ.size(); ++idLayer)
    {
      std::string name_Si("HypHI_InSi_log");
      name_Si += std::to_string(idLayer);
      G4LogicalVolume* HypHI_InSi_log = new G4LogicalVolume(HypHI_SiliciumSeg, Si, name_Si, 0, 0, 0);
      NameDetectorsSD.push_back(HypHI_InSi_log->GetName());

      for(G4int IdSi = 0; IdSi < nb_panel; ++IdSi)
        {
          G4RotationMatrix* rotSi = new G4RotationMatrix;
          double rotAngle         = CLHEP::twopi / static_cast<double>(nb_panel) * static_cast<double>(IdSi);
          rotSi->rotateZ(rotAngle);
          std::string nameSi("HypHI_Layer");
          nameSi += std::to_string(idLayer);
          nameSi += "_SiSeg_";
          nameSi += std::to_string(IdSi);
          AllPlacements.emplace_back(new G4PVPlacement(rotSi, G4ThreeVector(0, 0, posZ[idLayer]), HypHI_InSi_log,
                                                       nameSi, HypHI_InTracker_log, false, IdSi));
        }

      Si_att->SetForceWireframe(false);
      HypHI_InSi_log->SetVisAttributes(Si_att);
    }

  G4VSolid* FiberTrk_box = new G4Box("FiberTrk_box", 15. * cm, 15. * cm, 1. * mm);

  std::vector<double> posZ_Fiber = {40. * cm, 55. * cm, 70. * cm};
  double additional_shift        = 0.;
  std::cout << " IF :" << lengthDownstream - (0.5 * lengthDownstream + shift_log) << " " << 70. * cm << "\n";
  if(lengthDownstream - (0.5 * lengthDownstream + shift_log) < 70. * cm)
    additional_shift =
        10 * cm; //(70*cm+TargetLength + TargetPos_z + shift_log) - (lengthDownstream-(0.5*lengthDownstream+shift_log));
  std::cout << " Add:" << additional_shift << "\n";
  for(auto& pos : posZ_Fiber)
    {
      std::cout << " pos:" << pos;
      pos += TargetLength + TargetPos_z + shift_log - additional_shift;
      std::cout << " " << pos << "\n";
    }
  for(size_t idLayer = 0; idLayer < posZ_Fiber.size(); ++idLayer)
    {
      std::string name_Fiber("HypHI_Fiber_log");
      name_Fiber += std::to_string(idLayer);
      G4LogicalVolume* FiberTrk_log = new G4LogicalVolume(FiberTrk_box, Scinti, name_Fiber, 0, 0, 0);
      NameDetectorsSD.push_back(FiberTrk_log->GetName());

      std::string nameFiber("HypHI_Fiber");
      nameFiber += std::to_string(idLayer);

      AllPlacements.emplace_back(new G4PVPlacement(0, G4ThreeVector(0., 0., posZ_Fiber[idLayer]), FiberTrk_log,
                                                   nameFiber, HypHI_InTracker_log, false, idLayer));
    }
}

void WasaSimpleDetectorConstruction::ConstructPSB()
{

  G4NistManager* materialMgr = G4NistManager::Instance();

  // G4Material* HeEthan_50_50 = materialMgr-> GetMaterial("He_Ethan_50_50");
  // G4Material* ArCH4_90_10  = materialMgr-> GetMaterial("ArCH4_90_10");
  // G4Material* ArEthan_50_50 = materialMgr-> GetMaterial("ArEthan_50_50");
  // G4Material* HeIsobutane_80_20 = materialMgr-> GetMaterial("He_Isobutane_80_20");
  // G4Material* ArEthan_80_20 = materialMgr-> GetMaterial("ArEthan_80_20");
  // G4Material* ArCO2Methan_89_10_1 = materialMgr-> GetMaterial("ArCO2Methan_89_10_1");

  // G4Material* Air = materialMgr->FindOrBuildMaterial("G4_AIR");
  G4Material* Vacuum = materialMgr->FindOrBuildMaterial("G4_Galactic");
  // G4Material* HeGas = materialMgr->FindOrBuildMaterial("G4_He");
  // G4Material* Fe = materialMgr->FindOrBuildMaterial("G4_Fe");
  G4Material* Scinti = materialMgr->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"); // G4_POLYETHYLENE");//"Plastic");
  // G4Material* CFRP = materialMgr-> GetMaterial("CFRP");
  // G4Material* Mylar = materialMgr->FindOrBuildMaterial("G4_MYLAR");
  // G4Material* W = materialMgr->FindOrBuildMaterial("G4_W");
  // G4Material* Al = materialMgr->FindOrBuildMaterial("G4_Al");
  // G4Material* LHe3 = materialMgr->FindOrBuildMaterial("LHelium-3");
  // G4Material* CH2 = materialMgr->FindOrBuildMaterial("CH2");
  // G4Material* Be = materialMgr->FindOrBuildMaterial("G4_Be");
  // G4Material* Kapton = materialMgr->FindOrBuildMaterial("G4_KAPTON");
  // G4Material* Cu = materialMgr->FindOrBuildMaterial("G4_Cu");
  // G4Material* D2Gas = materialMgr->FindOrBuildMaterial("DeuteriumGas");
  // G4Material* SUS = materialMgr->FindOrBuildMaterial("SUS");
  // G4Material* Ar = materialMgr-> GetMaterial("ArgonGas");

  // G4Material* ChamberGas=Vacuum;

  //***********//
  //*** PSB ***//
  //***********//
  G4double psb_rmin = MDC_RadiusMax * mm + 10 * mm;
  G4double psb_rmax = psb_rmin + PSB_Thickness;
  G4double psb_z    = PSB_Length / 2.0 * mm;

  double N_bars = std::floor(CLHEP::pi * (psb_rmin + psb_rmax) / PSB_Width);
  int Nb_bars   = static_cast<int>(N_bars);
  // char phys_name[100];

  G4Tubs* PSB_Set             = new G4Tubs("PSB_Set", psb_rmin, psb_rmax, psb_z, 0, CLHEP::twopi);
  G4LogicalVolume* PSB_Setlog = new G4LogicalVolume(PSB_Set, Vacuum, "PSB_Setlog", 0, 0, 0);
  AllPlacements.emplace_back(
      new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), PSB_Setlog, "PSB_SetPhys", CD_log, false, 0));
  if(DoForRoot)
    AllPlacements.emplace_back(
        new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), PSB_Setlog, "PSB_SetPhysR", CD_logOutRoot, false, 0));

  PSB_Setlog->SetVisAttributes(G4VisAttributes::Invisible);

  G4Tubs* PSB_tube         = new G4Tubs("PSB_tube", psb_rmin, psb_rmax, psb_z, 0, CLHEP::twopi / N_bars);
  G4LogicalVolume* PSB_log = new G4LogicalVolume(PSB_tube, Scinti, "PSB_log", 0, 0, 0);
  NameDetectorsSD.push_back(PSB_log->GetName());

  std::vector<G4PVPlacement*> PSB_phys(Nb_bars, nullptr);
  for(G4int i = 0; i < Nb_bars; i++)
    {
      // sprintf(phys_name,"PSB_phys%02d", i+1);
      std::string phys_name("PSB_phys");
      phys_name += std::to_string(i + 1);
      G4ThreeVector xyzCounter(0., 0., 0.);
      G4RotationMatrix* rotCounterPSB = new G4RotationMatrix;
      rotCounterPSB->rotateZ(CLHEP::twopi / N_bars * (i));

      G4Transform3D posCounterPSB(*rotCounterPSB, xyzCounter);

      PSB_phys[i] = new G4PVPlacement(posCounterPSB, PSB_log, phys_name, PSB_Setlog, false, i);
    }

  // PSB_log->SetSensitiveDetector(counterSD);

  //--- Visualization ---//
  G4VisAttributes* PSB_att = new G4VisAttributes(Red);
  PSB_att->SetForceWireframe(false);
  // PSB_att->SetForceSolid(true);
  PSB_log->SetVisAttributes(PSB_att);
}

void WasaSimpleDetectorConstruction::ConstructEndFMF2()
{
  //***********//
  //*** TOF ***//
  //***********//
  G4NistManager* materialMgr = G4NistManager::Instance();

  G4Material* Scinti = materialMgr->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"); // G4_POLYETHYLENE");//"Plastic");

  G4VSolid* EndFMF2_box        = new G4Box("EndFMF2_box", 25. * cm, 25. * cm, 1. * mm);
  G4LogicalVolume* EndFMF2_log = new G4LogicalVolume(EndFMF2_box, Scinti, "FMF2_log", 0, 0, 0);

  NameDetectorsSD.push_back(EndFMF2_log->GetName());

  // G4PVPlacement* KDV_phys  =
  AllPlacements.emplace_back(
      new G4PVPlacement(0, G4ThreeVector(0., 0., 2. * m), EndFMF2_log, "FMF2_phys", experimentalHall_log, false, 0));
  AllPlacements.emplace_back(new G4PVPlacement(0, G4ThreeVector(0., 0., 2. * m + 1. * cm), EndFMF2_log, "FMF2_phys1",
                                               experimentalHall_log, false, 1));
  AllPlacements.emplace_back(new G4PVPlacement(0, G4ThreeVector(0., 0., 2. * m + 2. * cm), EndFMF2_log, "FMF2_phys2",
                                               experimentalHall_log, false, 2));

  if(DoForRoot)
    {
      AllPlacements.emplace_back(new G4PVPlacement(0, G4ThreeVector(0., 0., 2. * m), EndFMF2_log, "FMF2_physR",
                                                   experimentalHall_logOutRoot, false, 0));
      AllPlacements.emplace_back(new G4PVPlacement(0, G4ThreeVector(0., 0., 2. * m + 1. * cm), EndFMF2_log,
                                                   "FMF2_physR1", experimentalHall_logOutRoot, false, 1));
      AllPlacements.emplace_back(new G4PVPlacement(0, G4ThreeVector(0., 0., 2. * m + 2. * cm), EndFMF2_log,
                                                   "FMF2_physR2", experimentalHall_logOutRoot, false, 2));
    }
  //--- Visualization ---//
  G4VisAttributes* FMF2_att = new G4VisAttributes(Red);
  FMF2_att->SetForceWireframe(false);
  // KDV_att->SetForceSolid(true);
  EndFMF2_log->SetVisAttributes(FMF2_att);
}

void WasaSimpleDetectorConstruction::DefineCommands()
{

  // Define /B5/detector command directory using generic messenger class
  fMessenger = new G4GenericMessenger(this, "/G4SolSimple/detector/", "Detector control");

  G4GenericMessenger::Command& CDSCmd = fMessenger->DeclareProperty("SetCD", DoCD);
  G4String guidance4                  = "Boolean flag for setting CD.\n";
  CDSCmd.SetGuidance(guidance4);
  CDSCmd.SetParameterName("SetCD", true);
  CDSCmd.SetDefaultValue("true");

  G4GenericMessenger::Command& PSBCmd = fMessenger->DeclareProperty("SetPSB", DoPSB);
  G4String guidance5                  = "Boolean flag for setting PSB.\n";
  PSBCmd.SetGuidance(guidance5);
  PSBCmd.SetParameterName("SetPSB", false);
  PSBCmd.SetDefaultValue("false");
}
