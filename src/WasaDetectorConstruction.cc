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

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
//#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4Colour.hh"
#include "G4GeometryManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4NistManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4VisAttributes.hh"

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

#include "G4AutoDelete.hh"
#include "G4FieldManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ProductionCuts.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4SDManager.hh"
#include "G4SolSensitiveD.hh"
#include "G4SolSimpleMagneticField.hh"
#include "G4SolWASAMapMagneticField.hh"
#include "G4TransportationManager.hh"
#include "G4UserLimits.hh"

#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4MagneticField* WasaDetectorConstruction::fMagneticField = 0;
G4ThreadLocal G4FieldManager* WasaDetectorConstruction::fFieldMgr       = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WasaDetectorConstruction::WasaDetectorConstruction(G4SolConfig& _par)
    : G4SolVDetectorConstruction(_par), experimentalHall_log(nullptr), experimentalHall_phys(nullptr),
      MFLD_log(nullptr), MFLD_phys(nullptr), fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// ---------------------------------------------------------------------------
// VGM demo
//

WasaDetectorConstruction::~WasaDetectorConstruction() {}

G4VPhysicalVolume* WasaDetectorConstruction::Construct()
{

  // Input calibration parameters ////////////

  // Fiber
  double fiber_mft1_pos_x = 0.;  double fiber_mft1_pos_y = 0.;
  double fiber_mft2_pos_x = 0.;  double fiber_mft2_pos_y = 0.;

  double fiber_mft1_off_x1 = 0.; double fiber_mft1_off_u1 = 0.; double fiber_mft1_off_v1 = 0.;
  double fiber_mft1_off_x2 = 0.; double fiber_mft1_off_u2 = 0.; double fiber_mft1_off_v2 = 0.;
  double fiber_mft2_off_x1 = 0.; double fiber_mft2_off_u1 = 0.; double fiber_mft2_off_v1 = 0.;
  double fiber_mft2_off_x2 = 0.; double fiber_mft2_off_u2 = 0.; double fiber_mft2_off_v2 = 0.;

  double fiber_uft1_pos_x = 0.;  double fiber_uft1_pos_y = 0.;
  double fiber_uft2_pos_x = 0.;  double fiber_uft2_pos_y = 0.;
  double fiber_uft3_pos_x = 0.;  double fiber_uft3_pos_y = 0.;

  double fiber_uft1_off_x = 0.;  double fiber_uft1_off_u = 0.;  double fiber_uft1_off_v = 0.;
  double fiber_uft2_off_x = 0.;  double fiber_uft2_off_u = 0.;  double fiber_uft2_off_v = 0.;
  double fiber_uft3_off_x = 0.;  double fiber_uft3_off_u = 0.;  double fiber_uft3_off_v = 0.;

  double fiber_dft1_pos_x = 0.;  double fiber_dft1_pos_y = 0.;
  double fiber_dft2_pos_x = 0.;  double fiber_dft2_pos_y = 0.;

  double fiber_dft1_off_x = 0.;  double fiber_dft1_off_u = 0.;  double fiber_dft1_off_v = 0.;
  double fiber_dft2_off_x = 0.;  double fiber_dft2_off_u = 0.;  double fiber_dft2_off_v = 0.;

  std::array<std::array<std::array<double, 384>, 3>, 7> fiber_offset;
  for(int i=0; i<7; ++i)
    for(int j=0; j<3; ++j)
      for(int k=0; k<384; ++k)
        fiber_offset[i][j][k] = 0.;
  std::array<std::array<std::array<double, 2>, 3>, 7> fiber_angle_offset;
  for(int i=0; i<7; ++i)
    for(int j=0; j<3; ++j)
      for(int k=0; k<2; ++k)
        fiber_angle_offset[i][j][k] = 0.;
  std::array<std::array<std::array<std::array<double, 3>, 2>, 3>, 2> fiber_mft_cor_par;
  for(int i=0; i<2; ++i)
    for(int j=0; j<3; ++j)
      for(int k=0; k<2; ++k)
        for(int l=0; l<3; ++l)
          fiber_mft_cor_par[i][j][k][l] = 0.;

  if(Par.IsAvailable("Calib_Fiber_ON") && Par.Get<int>("Calib_Fiber_ON")==1){

    fiber_mft1_pos_x = 2.1*mm;  fiber_mft1_pos_y = 2.5*mm;
    fiber_mft2_pos_x = 2.1*mm;  fiber_mft2_pos_y = 2.5*mm;

    fiber_mft1_off_x1 = ( 0.93 + 0.7 )*mm;  fiber_mft1_off_u1 = (1.49 - 0.45       )*mm;  fiber_mft1_off_v1 = ( 0.76 + 0.05 + 0.66)*mm;
    fiber_mft1_off_x2 = (-0.22 - 0.05)*mm;  fiber_mft1_off_u2 = (0.25 - 0.05 - 0.17)*mm;  fiber_mft1_off_v2 = (-0.52 - 0.29       )*mm;
    fiber_mft2_off_x1 = ( 0.70 + 0.93)*mm;  fiber_mft2_off_v1 = (0.61 + 0.8        )*mm;  fiber_mft2_off_u1 = ( 1.62 - 0.55       )*mm;
    fiber_mft2_off_x2 = (-0.13 - 0.15)*mm;  fiber_mft2_off_v2 = (0.10 - 0.50       )*mm;  fiber_mft2_off_u2 = (-0.08 - 0.20 - 0.08)*mm;


    fiber_uft1_pos_x = 0.*mm;  fiber_uft1_pos_y = 0.*mm;
    fiber_uft2_pos_x = 0.*mm;  fiber_uft2_pos_y = 0.*mm;
    fiber_uft3_pos_x = 0.*mm;  fiber_uft3_pos_y = 0.*mm;

    fiber_uft1_off_x =  0.294*mm; fiber_uft1_off_u =  0.102*mm; fiber_uft1_off_v =  0.273*mm;
    fiber_uft2_off_x = -0.258*mm; fiber_uft2_off_u = -0.100*mm; fiber_uft2_off_v = -0.247*mm;
    fiber_uft3_off_x = -1.194*mm; fiber_uft3_off_u = -0.422*mm; fiber_uft3_off_v = -0.887*mm;


    fiber_dft1_pos_x = 0.*mm;  fiber_dft1_pos_y =  12.*mm;
    fiber_dft2_pos_x = 0.*mm;  fiber_dft2_pos_y = -12.*mm;

    fiber_dft1_off_x = 0.040*mm; fiber_dft1_off_u  = -0.287*mm; fiber_dft1_off_v  =  0.892*mm;
    fiber_dft2_off_x = 0.183*mm; fiber_dft2_off_u  =  0.436*mm; fiber_dft2_off_v  = -0.377*mm;


    std::string fiber_name_offset = Par.Get<std::string>("CalibFile_Fiber_Offset");
    std::ifstream ifs_fiber ( fiber_name_offset );
    if(ifs_fiber.is_open())
    {
      const std::string CommentSymbol("#");

      std::string temp_line;
      while(std::getline(ifs_fiber,temp_line))
      {
        std::stringstream stream(temp_line);
        std::string testComment(stream.str());
        std::size_t it_comment = testComment.find(CommentSymbol);
        if(it_comment!=std::string::npos)
        {
          //std::cout<<"!> Skip comment"<<temp_line<<std::endl;
          continue;
        }

        int det_id, lay_id, fib_id;
        double offset;

        stream >> det_id >> lay_id >> fib_id >> offset;
        //printf("%d %d %d : %.2f\n", det_id, lay_id, fib_id, offset );

        fiber_offset[det_id][lay_id][fib_id] = offset*mm;
      }
      //std::cout << "done " << mdc_name_map << std::endl;
      printf("fiber offset loaded : %s\n", fiber_name_offset.c_str());
    }
    else
    {
      //std::cout << " ! fail to open " << mdc_name_map << std::endl;
      printf(" ! fail to open  : %s\n", fiber_name_offset.c_str());
      exit(-1);
    }


    std::string fiber_name_angleoffset = Par.Get<std::string>("CalibFile_Fiber_Angleoffset");
    std::ifstream ifs_anglefiber ( fiber_name_angleoffset );
    if(ifs_anglefiber.is_open())
    {
      const std::string CommentSymbol("#");

      std::string temp_line;
      while(std::getline(ifs_anglefiber,temp_line))
      {
        std::stringstream stream(temp_line);
        std::string testComment(stream.str());
        std::size_t it_comment = testComment.find(CommentSymbol);
        if(it_comment!=std::string::npos)
        {
          //std::cout<<"!> Skip comment"<<temp_line<<std::endl;
          continue;
        }

        int det_id, lay_id, seg_id;
        double angle;

        stream >> det_id >> lay_id >> seg_id >> angle;
        //printf("%d %d %d : %.2f\n", det_id, lay_id, fib_id, offset );

        fiber_angle_offset[det_id][lay_id][seg_id] = angle;
      }
      //std::cout << "done " << mdc_name_map << std::endl;
      printf("fiber angle offset loaded : %s\n", fiber_name_angleoffset.c_str());
    }
    else
    {
      //std::cout << " ! fail to open " << mdc_name_map << std::endl;
      printf(" ! fail to open  : %s\n", fiber_name_angleoffset.c_str());
      exit(-1);
    }

    std::string fiber_name_mftcor = Par.Get<std::string>("CalibFile_Fiber_MFTcor");
    std::ifstream ifs_mftcorfiber ( fiber_name_mftcor );
    if(ifs_mftcorfiber.is_open())
    {
      const std::string CommentSymbol("#");

      std::string temp_line;
      while(std::getline(ifs_mftcorfiber,temp_line))
      {
        std::stringstream stream(temp_line);
        std::string testComment(stream.str());
        std::size_t it_comment = testComment.find(CommentSymbol);
        if(it_comment!=std::string::npos)
        {
          //std::cout<<"!> Skip comment"<<temp_line<<std::endl;
          continue;
        }

        int det, lay, seg;
        double p0, p1, p2;

        stream >> det >> lay >> seg >> p0 >> p1 >> p2;
        //printf("%d %d %d : %.2f\n", det_id, lay_id, fib_id, offset );

        fiber_mft_cor_par[det][lay][seg][0] = p0;
        fiber_mft_cor_par[det][lay][seg][1] = p1;
        fiber_mft_cor_par[det][lay][seg][2] = p2;
      }
      //std::cout << "done " << mdc_name_map << std::endl;
      printf("fiber mft corrections loaded : %s\n", fiber_name_mftcor.c_str());
    }
    else
    {
      //std::cout << " ! fail to open " << mdc_name_map << std::endl;
      printf(" ! fail to open  : %s\n", fiber_name_mftcor.c_str());
      exit(-1);
    }

  }

  // MDC
  double mdc_pos_x = 0.;
  double mdc_pos_y = 0.;
  double mdc_pos_z = 0.;

  double mdc_rot_x = 0.;
  double mdc_rot_y = 0.;
  double mdc_rot_z = 0.;

  if(Par.IsAvailable("Calib_MDC_ON") && Par.Get<int>("Calib_MDC_ON")==1)
    {
      mdc_pos_x = 1.5;
      mdc_pos_y = 5.5;
      mdc_pos_z = -4.;

      mdc_rot_x = -0.32;
      mdc_rot_y = -0.1;
      mdc_rot_z = -0.4;
    }

  // PSB
  double psb_pos_x = 0.;
  double psb_pos_y = 0.;
  double psb_pos_z = 0.;

  double psb_rot_z = 0.;
  if(Par.IsAvailable("Calib_PSB_ON") && Par.Get<int>("Calib_PSB_ON")==1)
    {
      psb_pos_x = 0.5;
      psb_pos_y = 5.5;
      psb_pos_z = 0.;

      psb_rot_z = -0.4;

    }


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

  experimentalHall_log  = world->GetLogicalVolume();
  experimentalHall_phys = world;

  MFLD_log  = FindVolume("MFLD");
  MFLD_phys = FindVolPhys("MFLD");

  std::cout << "MFLD Phys:" << MFLD_phys << " " << MFLD_phys->GetInstanceID() << "\n";
  transMFLD = MFLD_phys->GetObjectTranslation();
  rotMFLD   = MFLD_phys->GetObjectRotationValue();
  std::cout << "Trans:" << transMFLD << "\n";
  std::cout << "Rot:" << rotMFLD << "\n";

  const G4double Wasa_Zshift      = Par.Get<double>("Wasa_ShiftZ");
  const G4double Systematic_shift = Par.Get<double>("Systematic_Shift");
  const G4double TargetPosX       = Par.Get<double>("Target_PosX");
  const G4double TargetPosY       = Par.Get<double>("Target_PosY");
  const G4double TargetPosZ       = Par.Get<double>("Target_PosZ");

  const int WasaSide = Par.Get<int>("Wasa_Side");
  double Sign = WasaSide == 0 ? 1. : -1.;
  G4ThreeVector transMFLD_move(0., 0., Wasa_Zshift + Systematic_shift);
  //G4ThreeVector transMFLD_new = transMFLD + transMFLD_move;
  G4ThreeVector transMFLD_new = transMFLD*Sign + transMFLD_move;

  if(WasaSide == 1)
    rotMFLD.rotateY(180 * degree);

  std::cout << "Trans:" << transMFLD_new << "\n";
  std::cout << "RotAfter:" << rotMFLD << "\n";

  MFLD_phys->SetTranslation(transMFLD_new);

  if(WasaSide == 1)
    MFLD_phys->SetRotation(&rotMFLD);

  // Get volumes from logical volume store by name
  // G4LogicalVolume* calorLV = FindVolume("Calorimeter");

  G4LogicalVolume* WasaLV  = FindVolume("WASA");
  G4LogicalVolume* InnerLV = FindVolume("INNER");
  if(!Par.IsAvailable("HypHI_InnerTrackerBox_Visible"))
    InnerLV->SetVisAttributes(G4VisAttributes::GetInvisible());

  INNER_log  = FindVolume("INNER");
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

  // G4Region* aDetectorRegion = new G4Region("DetectorRegion");

  bool MiniVis = Par.IsAvailable("Gui_MiniVis");

  std::vector<G4String> NameCDC1_Invisible = {"ME01", "ME02", "ME03", "ME04", "ME05", "ME06", "ME07", "ME08", "ME09",
                                              "ME10", "ME11", "ME12", "ME13", "ME14", "ME15", "ME16", "ME17"};
  std::vector<G4String> NameCDC2_Invisible = {"MD01", "MD02", "MD03", "MD04", "MD05", "MD06", "MD07", "MD08", "MD09",
                                              "MD10", "MD11", "MD12", "MD13", "MD14", "MD15", "MD16", "MD17"};
  for(auto& CurrentName : NameCDC1_Invisible)
    {
      G4LogicalVolume* UTracker = FindVolume(CurrentName);
      UTracker->SetVisAttributes(G4VisAttributes::GetInvisible());
    }
  int iColorT = 0;
  for( auto& CurrentName : NameCDC2_Invisible)
    {
      G4LogicalVolume* UTracker = FindVolume(CurrentName);
      if(MiniVis)
        {
          G4VisAttributes MG_Color(ColorCDC[iColorT]);
          UTracker->SetVisAttributes(MG_Color);
          ++iColorT;
        }
      else
        UTracker->SetVisAttributes(G4VisAttributes::GetInvisible());
    }

  FindVolume("SOL")->SetVisAttributes(G4VisAttributes::GetInvisible());

  // int iColorT = 0;

  double PosZShiftMDC = 0.;
  if(Par.IsAvailable("PosZ_ShiftMDC"))
    {
      std::vector<G4String> Name_MDC_Shift = {"MD01", "MD02", "MD03", "MD04", "MD05", "MD06", "MD07", "MD08",
                                              "MD09", "MD10", "MD11", "MD12", "MD13", "MD14", "MD15", "MD16",
                                              "MD17", "MDO",  "MDB",  "MDF",  "PTB0", "PTB1", "PTB2", "PTB3"};

      PosZShiftMDC = Par.Get<double>("PosZ_ShiftMDC");

      if(PosZShiftMDC < 0.)
        Name_MDC_Shift.emplace_back("PSB");
      else
        Name_MDC_Shift.emplace_back("PSF");

      if(TMath::Abs(PosZShiftMDC) > 1e-9)
        for(auto& CurrentName : Name_MDC_Shift)
          {
            G4VPhysicalVolume* MDC_temp = FindVolPhys(CurrentName);
            if(MDC_temp == nullptr)
              {
                std::cout << "E> no volume :" << CurrentName << "\n";
              }
            G4ThreeVector transMDC = MDC_temp->GetObjectTranslation();
            G4ThreeVector transMDC_move(0., 0., PosZShiftMDC);
            G4ThreeVector transMDC_new = transMDC + transMDC_move;
            MDC_temp->SetTranslation(transMDC_new);
            if(CurrentName == "MDB" && PosZShiftMDC < -7.5 * cm)
              {
                std::cout << "MDB -> outside : removed \n";
                // MDC_temp->SetMotherLogical(MFLD_log);
                INNER_log->RemoveDaughter(MDC_temp);
              }
            if(CurrentName == "MDF" && PosZShiftMDC > 2.95 * cm)
              {
                std::cout << "MDF -> outside : removed \n";
                // MDC_temp->SetMotherLogical(MFLD_log);
                INNER_log->RemoveDaughter(MDC_temp);
              }
          }
    }

  if(Par.IsAvailable("Calib_MDC_ON") && Par.Get<int>("Calib_MDC_ON")==1)
    {


      // std::vector<G4String> Name_MDC = {"MD01", "MD02", "MD03", "MD04", "MD05", "MD06", "MD07", "MD08",
      // 					"MD09", "MD10", "MD11", "MD12", "MD13", "MD14", "MD15", "MD16",
      // 					"MD17", "MDO",  "MDB",  "MDF",  "PTB0", "PTB1", "PTB2", "PTB3"};

      // for(size_t iName = 0; iName< Name_MDC.size() ; ++iName)
      // 	{
      // 	  auto& CurrentName = Name_MDC[iName];
      // 	  G4VPhysicalVolume* MDC_temp = FindVolPhys(CurrentName);
      // 	  if(MDC_temp == nullptr)
      // 	  std::cout << "E> no volume :" << CurrentName << "\n";

      // 	transMDC[iName] = MDC_temp->GetObjectTranslation();
      // 	G4ThreeVector transMDC_move(mdc_pos_x*mm * Sign, mdc_pos_y*mm, mdc_pos_z*mm * Sign);
      // 	G4ThreeVector transMDC_new = transMDC[iName] + transMDC_move;
      // 	rotMDC[iName] = MDC_temp->GetObjectRotationValue();
      // 	rotMDC[iName].rotateZ(mdc_rot_z*deg);// * Sign);
      // 	rotMDC[iName].rotateX(mdc_rot_x*deg);// * Sign);
      // 	rotMDC[iName].rotateY(mdc_rot_y*deg);

      // 	MDC_temp->SetTranslation(transMDC_new);
      // 	std::cout<<" volume : "<<CurrentName<<" "<<rotMDC.size()<<" "<<rotMDC[iName]<<"\n";

      // 	MDC_temp->SetRotation(&(rotMDC[iName]));
      // }

      G4VPhysicalVolume* MDC_temp = FindVolPhys("MDC");

      transMDC = MDC_temp->GetObjectTranslation();
      G4ThreeVector transMDC_move(mdc_pos_x*mm * Sign, mdc_pos_y*mm, mdc_pos_z*mm * Sign);
      G4ThreeVector transMDC_new = transMDC + transMDC_move;

      rotMDC = MDC_temp->GetObjectRotationValue();
      //std::cout << "rotMDC before : " << rotMDC << std::endl;
      rotMDC.rotateZ(mdc_rot_z * degree);
      rotMDC.rotateX(mdc_rot_x * degree);
      rotMDC.rotateY(mdc_rot_y * degree);
      //std::cout << "rotMDC after : " << rotMDC << std::endl;

      MDC_temp->SetTranslation(transMDC_new);
      MDC_temp->SetRotation(&rotMDC);
      //G4RotationMatrix rotMDC_buf = MDC_temp->GetObjectRotationValue();
      //std::cout << "rotMDC buf : " << rotMDC_buf << std::endl;

    }

  if(Par.IsAvailable("Calib_PSB_ON") && Par.Get<int>("Calib_PSB_ON")==1)
    {

      G4VPhysicalVolume* PSCE_temp = FindVolPhys("PSCEall");

      transPSCE = PSCE_temp->GetObjectTranslation();
      G4ThreeVector transPSCE_move(psb_pos_x*mm * Sign, psb_pos_y*mm, psb_pos_z*mm * Sign);
      G4ThreeVector transPSCE_new = transPSCE + transPSCE_move;

      rotPSCE = PSCE_temp->GetObjectRotationValue();
      //std::cout << "rotPSB before : " << rotPSB << std::endl;
      rotPSCE.rotateZ(psb_rot_z * degree);
      //std::cout << "rotPSB after : " << rotPSB << std::endl;

      PSCE_temp->SetTranslation(transPSCE_new);
      PSCE_temp->SetRotation(&rotPSCE);

      //G4RotationMatrix rotPSB_buf = PSB_temp->GetObjectRotationValue();
      //std::cout << "rotPSB buf : " << rotPSB_buf << std::endl;

      // for(auto& CurrentName : Name_PSCE){
      //   G4VPhysicalVolume* PSCE_temp = FindVolPhys(CurrentName);
      //   if(PSCE_temp == nullptr) std::cout << "E> no volume :" << CurrentName << "\n";
      //   transPSCE = PSCE_temp->GetObjectTranslation();
      //   G4ThreeVector transPSCE_move(psb_pos_x*mm * Sign, psb_pos_y*mm, psb_pos_z*mm * Sign);
      //   G4ThreeVector transPSCE_new = transPSCE + transPSCE_move;
      //   rotPSCE = PSCE_temp->GetObjectRotationValue();
      //   rotPSCE.rotateZ(psb_rot_z*deg * Sign);
      //   PSCE_temp->SetTranslation(transPSCE_new);
      // }

    }


  std::vector<G4String> NameSD = {"MG01", "MG02", "MG03", "MG04", "MG05", "MG06", "MG07", "MG08", "MG09", "MG10",
                                  "MG11", "MG12", "MG13", "MG14", "MG15", "MG16", "MG17", "PSCE", "PSBE", "PSFE"};

  NameDetectorsSD = NameSD;

  std::vector<G4String> NameSD_Color = {"MG01", "MG02", "MG03", "MG04", "MG05", "MG06", "MG07", "MG08", "MG09",
                                        "MG10", "MG11", "MG12", "MG13", "MG14", "MG15", "MG16", "MG17"};

  FindVolume("PSCE")->SetVisAttributes(G4VisAttributes(Blue));
  FindVolume("PSBE")->SetVisAttributes(G4VisAttributes(Blue));
  FindVolume("PSFE")->SetVisAttributes(G4VisAttributes(LightBlue));

  //     UTracker->SetRegion(aDetectorRegion);
  //     aDetectorRegion->AddRootLogicalVolume(UTracker);

  // UTracker->SetVisAttributes(VisDetectorSD);

  for(size_t iColor = 0; iColor < NameSD_Color.size(); ++iColor)
    {
      G4LogicalVolume* UTracker = FindVolume(NameSD_Color[iColor]);
      if(!MiniVis)
        {
          G4VisAttributes MG_Color(ColorCDC[iColor]);
          UTracker->SetVisAttributes(MG_Color);
        }
      else
        UTracker->SetVisAttributes(G4VisAttributes::GetInvisible());
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
  // worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

  // worldLV->SetUserLimits( new G4UserLimits(DBL_MAX,2*m,10*s,0.,0.) );

  if(WasaLV)
    WasaLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  if(MFLD_log)
    if(!Par.IsAvailable("HypHI_InnerTrackerBox_Visible"))
      MFLD_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  // if(VolGapLV)
  //   VolGapLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  // if(SecondMagFieldLV)
  //   SecondMagFieldLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  // if(SecondMagLV)
  //   SecondMagLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  // if(SetupLV)
  //   SetupLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  // if(IronQuadLV)
  //   IronQuadLV->SetVisAttributes(G4VisAttributes::GetInvisible());

  // G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  // simpleBoxVisAtt->SetVisibility(true);
  // if (calorLV) calorLV->SetVisAttributes(simpleBoxVisAtt);

  G4NistManager* materialMgr = G4NistManager::Instance();

  // G4Material* Air = materialMgr->FindOrBuildMaterial("G4_AIR");
  G4Material* Air    = materialMgr->FindOrBuildMaterial("G4_AIR");
  G4Material* Vacuum = materialMgr->FindOrBuildMaterial("G4_Galactic");
  G4Material* Si     = materialMgr->FindOrBuildMaterial("G4_Si"); //"Plastic");
  G4Material* Scinti =
      materialMgr->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"); // G4_POLYETHYLENE");//"Plastic");
  G4Material* FiberCoreScinti = materialMgr->FindOrBuildMaterial("G4_POLYSTYRENE");
  G4Material* Carbon          = materialMgr->FindOrBuildMaterial("G4_C"); //"Plastic");
  materialMgr->BuildMaterialWithNewDensity("G4_Diamond","G4_C",3.*g/cm3);

  G4Material* Diamond         = materialMgr->FindOrBuildMaterial("G4_Diamond"); //"Plastic");

  G4VisAttributes* Si_att = new G4VisAttributes(Pink);

  G4double TargetLengthX = 0., TargetLengthY = 0., TargetLengthZ = 0.;
  if(Par.IsAvailable("Target_Size"))
    {
      TargetLengthX = Par.Get<double>("Target_Size"); // 1.0*cm;
      TargetLengthY = Par.Get<double>("Target_Size"); // 1.0*cm;
      TargetLengthZ = Par.Get<double>("Target_Size"); // 1.0*cm;
    }
  if(Par.IsAvailable("Target_SizeX"))
    TargetLengthX = Par.Get<double>("Target_SizeX"); // 1.0*cm;
  if(Par.IsAvailable("Target_SizeY"))
    TargetLengthY = Par.Get<double>("Target_SizeY"); // 1.0*cm;
  if(Par.IsAvailable("Target_SizeZ"))
    TargetLengthZ = Par.Get<double>("Target_SizeZ"); // 1.0*cm;

  G4double BeamHoleSize = 0. * cm;
  if(Par.IsAvailable("HypHI_BeamHole"))
    BeamHoleSize = Par.Get<double>("HypHI_BeamHole");

  G4VSolid* HypHI_Target = new G4Box("HypHI_Target", 0.5*TargetLengthX, 0.5*TargetLengthY, 0.5*TargetLengthZ);
  // G4LogicalVolume* HypHI_Target_log = new G4LogicalVolume(HypHI_Target, Carbon,"HypHI_Target_log", 0,0,0);
  G4Material* MatTarget             = Par.IsAvailable("Target_Diamond") ? Diamond : Par.IsAvailable("Target_Carbon") ? Carbon : Vacuum;
  G4LogicalVolume* HypHI_Target_log = new G4LogicalVolume(HypHI_Target, MatTarget, "HypHI_Target_log", 0, 0, 0);
  // G4PVPlacement*   HypHI_Target_phys =
  G4ThreeVector TargetTrans   = G4ThreeVector(TargetPosX, TargetPosY, TargetPosZ) - transMFLD_new;
  G4RotationMatrix* TargetRot = new G4RotationMatrix;
  // if(WasaSide==1)
  //   TargetRot->rotateY(-180*degree);
  G4Transform3D posTarget(*TargetRot, Sign * TargetTrans);
  AllPlacements.emplace_back(new G4PVPlacement(posTarget, HypHI_Target_log, "HypHI_Target_Phys", MFLD_log, false, 0));

  if(Par.IsAvailable("HypHI_InnerTracker_On"))
    {
      G4int nb_panel         = 16;
      G4double HypHI_Si_minR = Par.Get<double>("HypHI_Si_minR");
      G4double HypHI_Si_maxR = Par.Get<double>("HypHI_Si_maxR");

      G4VSolid* HypHI_SiliciumSeg = new G4Tubs("HypHI_SiSeg", HypHI_Si_minR, HypHI_Si_maxR, 3 * mm,
                                               -CLHEP::twopi / static_cast<double>(2 * nb_panel),
                                               2. * CLHEP::twopi / static_cast<double>(2 * nb_panel));

      std::vector<double> posZ = {20. * cm, 25. * cm, 30. * cm, 40. * cm};
      if(Par.IsAvailable("HypHI_InnerTracker_Spec"))
        {
          double PosZInTracker = Par.Get<double>("HypHI_InnerTracker_PosZ");
          int NbInTracker      = Par.Get<int>("HypHI_InnerTracker_Nb");
          double DistInTracker = Par.Get<double>("HypHI_InnerTracker_Spacing");
          posZ.resize(NbInTracker);
          for(size_t idL = 0; idL < posZ.size(); ++idL)
            posZ[idL] = TargetLengthZ + TargetPosZ + PosZInTracker + static_cast<double>(idL) * DistInTracker;
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
              G4ThreeVector SiTrans = G4ThreeVector(0., 0., posZ[idLayer]) - transMFLD_new;
              // if(WasaSide==1)
              //   rotSi->rotateY(-180*degree);
              G4Transform3D posSi(*rotSi, Sign * SiTrans);

              AllPlacements.emplace_back(new G4PVPlacement(posSi, HypHI_InSi_log, nameSi, MFLD_log, false, IdSi));
            }
        }
    }

  // ----------------------------- @@ -----------------------------
  //         Silicon Detectors
  // ----------------------------- @@ -----------------------------

  if(Par.IsAvailable("HypHI_Si1_On"))
    {
      G4double HypHI_Si1_length     = Par.Get<double>("HypHI_Si1_length");
      G4double HypHI_Si1_thickness  = Par.Get<double>("HypHI_Si1_thickness");
      G4double HypHI_Si1_stripwidth = Par.Get<double>("HypHI_Si1_stripwidth");
      G4double HypHI_Si1_posZ       = Par.Get<double>("HypHI_Si1_posZ");
      int HypHI_Si1_SingleSided     = Par.IsAvailable("HypHI_Si1_SingleSided") ? Par.Get<int>("HypHI_Si1_SingleSided") : 0;

      double scaling_thickness = HypHI_Si1_SingleSided == 1 ? 2 : 1 ;

      G4VSolid* Si1_MothVol_solid =
	new G4Box("Si1_solid", HypHI_Si1_length / 2., HypHI_Si1_length / 2., (scaling_thickness * HypHI_Si1_thickness + (0.25*mm)*HypHI_Si1_SingleSided)/ 2.);
      G4VSolid* Si1_layerX_MothVol_solid =
	new G4Box("Si1_layerX_solid", HypHI_Si1_length / 2., HypHI_Si1_length / 2., scaling_thickness * HypHI_Si1_thickness / 4.);
      G4VSolid* Si1_layerY_MothVol_solid =
          new G4Box("Si1_layerY_solid", HypHI_Si1_length / 2., HypHI_Si1_length / 2., scaling_thickness * HypHI_Si1_thickness / 4.);
      G4LogicalVolume* Si1_MothVol_log   = new G4LogicalVolume(Si1_MothVol_solid, Air, "Si1_log", 0, 0, 0);
      G4LogicalVolume* Si1_MothVol_log_x = new G4LogicalVolume(Si1_layerX_MothVol_solid, Air, "Si1_log_x", 0, 0, 0);
      G4LogicalVolume* Si1_MothVol_log_y = new G4LogicalVolume(Si1_layerY_MothVol_solid, Air, "Si1_log_y", 0, 0, 0);

      AllPlacements.emplace_back(
          new G4PVPlacement(0, Sign * (G4ThreeVector(0., 0., HypHI_Si1_posZ + Systematic_shift) - transMFLD_new),
                            Si1_MothVol_log, "Silicon1", MFLD_log, false, 0));
      AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(0., 0., -scaling_thickness * HypHI_Si1_thickness/ 4. - (0.125*mm)*HypHI_Si1_SingleSided)),
                                                   Si1_MothVol_log_x, "Si1_x", Si1_MothVol_log, false, 0));
      AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(0., 0., +scaling_thickness * HypHI_Si1_thickness / 4. + (0.125*mm)*HypHI_Si1_SingleSided)),
                                                   Si1_MothVol_log_y, "Si1_y", Si1_MothVol_log, false, 0));

      G4VSolid* Si1_Strip_solid_x =
          new G4Box("Si1_Strip_solid_x", HypHI_Si1_stripwidth / 2., HypHI_Si1_length / 2., scaling_thickness * HypHI_Si1_thickness / 4.);
      G4VSolid* Si1_Strip_solid_y =
          new G4Box("Si1_Strip_solid_y", HypHI_Si1_length / 2., HypHI_Si1_stripwidth / 2., scaling_thickness * HypHI_Si1_thickness / 4.);
      G4LogicalVolume* Si1_Strip_log_x = new G4LogicalVolume(Si1_Strip_solid_x, Si, "Si1_Strip_log_x", 0, 0, 0);
      G4LogicalVolume* Si1_Strip_log_y = new G4LogicalVolume(Si1_Strip_solid_y, Si, "Si1_Strip_log_y", 0, 0, 0);

      G4int HypHI_Si1_Nstrips = (int)(HypHI_Si1_length / HypHI_Si1_stripwidth);

      for(G4int idStrip = 0; idStrip < HypHI_Si1_Nstrips; ++idStrip)
        {
          G4double posStrip = -HypHI_Si1_length / 2. + HypHI_Si1_stripwidth * (0.5 + idStrip);

          AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(posStrip, 0., 0.)), Si1_Strip_log_x,
                                                       "Si1_Strip_x", Si1_MothVol_log_x, false, idStrip));
          AllPlacements.emplace_back(new G4PVPlacement(0, G4ThreeVector(0., posStrip, 0.), Si1_Strip_log_y,
                                                       "Si1_Strip_y", Si1_MothVol_log_y, false, idStrip));
        }

      NameDetectorsSD.push_back(Si1_Strip_log_x->GetName());
      NameDetectorsSD.push_back(Si1_Strip_log_y->GetName());

      if(!MiniVis)
        {
          Si1_Strip_log_x->SetVisAttributes(Si_att);
          Si1_Strip_log_y->SetVisAttributes(Si_att);
          Si1_MothVol_log->SetVisAttributes(G4VisAttributes::GetInvisible());
          Si1_MothVol_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
          Si1_MothVol_log_y->SetVisAttributes(G4VisAttributes::GetInvisible());
        }
      else
        {
          Si1_Strip_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
          Si1_Strip_log_y->SetVisAttributes(G4VisAttributes::GetInvisible());
          Si1_MothVol_log->SetVisAttributes(Si_att);
          Si1_MothVol_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
          Si1_MothVol_log_y->SetVisAttributes(G4VisAttributes::GetInvisible());
        }
    }

  if(Par.IsAvailable("HypHI_Si2_On"))
    {
      G4double HypHI_Si2_length     = Par.Get<double>("HypHI_Si2_length");
      G4double HypHI_Si2_thickness  = Par.Get<double>("HypHI_Si2_thickness");
      G4double HypHI_Si2_stripwidth = Par.Get<double>("HypHI_Si2_stripwidth");
      G4double HypHI_Si2_posZ       = Par.Get<double>("HypHI_Si2_posZ");
      int HypHI_Si2_SingleSided     = Par.IsAvailable("HypHI_Si2_SingleSided") ? Par.Get<int>("HypHI_Si2_SingleSided") : 0;

      double scaling_thickness = HypHI_Si2_SingleSided == 1 ? 2 : 1 ;

      G4VSolid* Si2_MothVol_solid =
	new G4Box("Si2_solid", HypHI_Si2_length / 2., HypHI_Si2_length / 2., (scaling_thickness * HypHI_Si2_thickness + (0.25*mm)*HypHI_Si2_SingleSided)/ 2. );
      G4VSolid* Si2_layerX_MothVol_solid =
          new G4Box("Si2_layerX_solid", HypHI_Si2_length / 2., HypHI_Si2_length / 2., scaling_thickness * HypHI_Si2_thickness / 4.);
      G4VSolid* Si2_layerY_MothVol_solid =
          new G4Box("Si2_layerY_solid", HypHI_Si2_length / 2., HypHI_Si2_length / 2., scaling_thickness * HypHI_Si2_thickness / 4.);
      G4LogicalVolume* Si2_MothVol_log   = new G4LogicalVolume(Si2_MothVol_solid, Air, "Si2_log", 0, 0, 0);
      G4LogicalVolume* Si2_MothVol_log_x = new G4LogicalVolume(Si2_layerX_MothVol_solid, Air, "Si2_log_x", 0, 0, 0);
      G4LogicalVolume* Si2_MothVol_log_y = new G4LogicalVolume(Si2_layerY_MothVol_solid, Air, "Si2_log_y", 0, 0, 0);

      AllPlacements.emplace_back(
          new G4PVPlacement(0, Sign * (G4ThreeVector(0., 0., HypHI_Si2_posZ + Systematic_shift) - transMFLD_new),
                            Si2_MothVol_log, "Silicon1", MFLD_log, false, 0));
      AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(0., 0., -scaling_thickness * HypHI_Si2_thickness / 4. - (0.125*mm)*HypHI_Si2_SingleSided)),
                                                   Si2_MothVol_log_x, "Si2_x", Si2_MothVol_log, false, 0));
      AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(0., 0., +scaling_thickness * HypHI_Si2_thickness / 4. + (0.125*mm)*HypHI_Si2_SingleSided)),
                                                   Si2_MothVol_log_y, "Si2_y", Si2_MothVol_log, false, 0));

      G4VSolid* Si2_Strip_solid_x =
          new G4Box("Si2_Strip_solid_x", HypHI_Si2_stripwidth / 2., HypHI_Si2_length / 2., scaling_thickness * HypHI_Si2_thickness / 4.);
      G4VSolid* Si2_Strip_solid_y =
          new G4Box("Si2_Strip_solid_y", HypHI_Si2_length / 2., HypHI_Si2_stripwidth / 2., scaling_thickness * HypHI_Si2_thickness / 4.);
      G4LogicalVolume* Si2_Strip_log_x = new G4LogicalVolume(Si2_Strip_solid_x, Si, "Si2_Strip_log_x", 0, 0, 0);
      G4LogicalVolume* Si2_Strip_log_y = new G4LogicalVolume(Si2_Strip_solid_y, Si, "Si2_Strip_log_y", 0, 0, 0);

      G4int HypHI_Si2_Nstrips = (int)(HypHI_Si2_length / HypHI_Si2_stripwidth);

      for(G4int idStrip = 0; idStrip < HypHI_Si2_Nstrips; ++idStrip)
        {
          G4double posStrip = -HypHI_Si2_length / 2. + HypHI_Si2_stripwidth * (0.5 + idStrip);

          AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(posStrip, 0., 0.)), Si2_Strip_log_x,
                                                       "Si2_Strip_x", Si2_MothVol_log_x, false, idStrip));
          AllPlacements.emplace_back(new G4PVPlacement(0, G4ThreeVector(0., posStrip, 0.), Si2_Strip_log_y,
                                                       "Si2_Strip_y", Si2_MothVol_log_y, false, idStrip));
        }

      NameDetectorsSD.push_back(Si2_Strip_log_x->GetName());
      NameDetectorsSD.push_back(Si2_Strip_log_y->GetName());
      if(!MiniVis)
        {
          Si2_Strip_log_x->SetVisAttributes(Si_att);
          Si2_Strip_log_y->SetVisAttributes(Si_att);
          Si2_MothVol_log->SetVisAttributes(G4VisAttributes::GetInvisible());
          Si2_MothVol_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
          Si2_MothVol_log_y->SetVisAttributes(G4VisAttributes::GetInvisible());
        }
      else
        {
          Si2_Strip_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
          Si2_Strip_log_y->SetVisAttributes(G4VisAttributes::GetInvisible());
          Si2_MothVol_log->SetVisAttributes(Si_att);
          Si2_MothVol_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
          Si2_MothVol_log_y->SetVisAttributes(G4VisAttributes::GetInvisible());
        }
    }

  if(Par.IsAvailable("HypHI_SD1_On"))
    {
      G4double HypHI_SD1_length      = 98.77 * mm;
      G4double HypHI_SD1_pcb         = 250 * mm;
      G4double HypHI_SD1_thickness   = 0.320 * mm;
      G4double HypHI_SD1_stripwidth  = 0.190 * mm;
      G4double HypHI_SD1_striplength = 48.365 * mm;
      G4double HypHI_SD1_gap         = Par.IsAvailable("HypHI_SD_NoGap") ? 0.0 * mm : 0.5 * mm;
      G4double HypHI_SD1_padding     = (0.235 + 0.510) * mm;
      G4int HypHI_SD1_Nch            = Par.Get<int>("HypHI_SD1_Nch");
      G4int HypHI_SD1_stripPerCh = Par.IsAvailable("HypHI_SD1_stripPerCh") ? Par.Get<int>("HypHI_SD1_stripPerCh") : 1;
      if(HypHI_SD1_Nch * HypHI_SD1_stripPerCh != 512)
        {
          std::cout << "E> SD1 number of channel does not correspond to the design ! " << HypHI_SD1_Nch << " "
                    << HypHI_SD1_stripPerCh << "\n";
          exit(-1);
        }
      if(HypHI_SD1_stripPerCh != 1)
        HypHI_SD1_stripwidth *= static_cast<double>(HypHI_SD1_stripPerCh);
      G4double HypHI_SD1_posZ = Par.Get<double>("HypHI_SD1_posZ");

      {
        G4VSolid* SD1u_MothVol_solid =
            new G4Box("SD1u_solid", HypHI_SD1_length / 2., HypHI_SD1_length / 2., HypHI_SD1_thickness / 2.);
        G4VSolid* SD1u_PCBext_Vol_solid =
            new G4Box("SD1u_PCBext_solid", HypHI_SD1_pcb / 2., HypHI_SD1_pcb / 2., HypHI_SD1_thickness / 2.);
        G4VSolid* SD1u_PCB_Vol_solid =
            new G4SubtractionSolid("SD1u_PCB_solid", SD1u_PCBext_Vol_solid, SD1u_MothVol_solid);
        G4VSolid* SD1u_layer1_MothVol_solid =
            new G4Box("SD1u_layer1_solid", HypHI_SD1_striplength / 2., HypHI_SD1_Nch * (HypHI_SD1_stripwidth / 2.),
                      HypHI_SD1_thickness / 2.);
        G4VSolid* SD1u_layer2_MothVol_solid =
            new G4Box("SD1u_layer2_solid", HypHI_SD1_striplength / 2., HypHI_SD1_Nch * (HypHI_SD1_stripwidth / 2.),
                      HypHI_SD1_thickness / 2.);
        G4LogicalVolume* SD1u_MothVol_log = new G4LogicalVolume(SD1u_MothVol_solid, Air, "SD1u_log", 0, 0, 0);
        G4LogicalVolume* SD1u_PCBVol_log  = new G4LogicalVolume(SD1u_PCB_Vol_solid, Air, "SD1u_PCB_log", 0, 0, 0);
        G4LogicalVolume* SD1u_MothVol_log_1 =
            new G4LogicalVolume(SD1u_layer1_MothVol_solid, Air, "SD1u_log_t", 0, 0, 0);
        G4LogicalVolume* SD1u_MothVol_log_2 =
            new G4LogicalVolume(SD1u_layer2_MothVol_solid, Air, "SD1u_log_b", 0, 0, 0);

        G4RotationMatrix* rotStereoSD1_u = nullptr;
        if(Par.IsAvailable("HypHI_SD1_stereoAngle"))
          {
            G4double HypHI_SD1_stereoAngle = Par.Get<double>("HypHI_SD1_stereoAngle");
            rotStereoSD1_u                 = new G4RotationMatrix;
            rotStereoSD1_u->rotateZ(HypHI_SD1_stereoAngle - 90 * deg);
          }
        G4ThreeVector TransSD1u = G4ThreeVector(0., 0., HypHI_SD1_posZ + Systematic_shift) - transMFLD_new;

        AllPlacements.emplace_back(
            new G4PVPlacement(rotStereoSD1_u, Sign * TransSD1u, SD1u_MothVol_log, "Silicon1u", MFLD_log, false, 0));

        if(Par.IsAvailable("HypHI_SD1_PCB"))
          AllPlacements.emplace_back(new G4PVPlacement(rotStereoSD1_u, Sign * TransSD1u, SD1u_PCBVol_log,
                                                       "Silicon1u_PCB", MFLD_log, false, 0));

        AllPlacements.emplace_back(
            new G4PVPlacement(0, Sign * (G4ThreeVector(-0.5 * (HypHI_SD1_striplength + HypHI_SD1_gap), 0., 0.)),
                              SD1u_MothVol_log_1, "SD1_u1", SD1u_MothVol_log, false, 0));
        AllPlacements.emplace_back(
            new G4PVPlacement(0, Sign * (G4ThreeVector(0.5 * (HypHI_SD1_striplength + HypHI_SD1_gap), 0., 0.)),
                              SD1u_MothVol_log_2, "SD1_u2", SD1u_MothVol_log, false, 0));

        G4VSolid* SD1_Strip_solid_u      = new G4Box("SD1_Strip_solid_u", HypHI_SD1_striplength / 2.,
                                                HypHI_SD1_stripwidth / 2., HypHI_SD1_thickness / 2.);
        G4LogicalVolume* SD1_Strip_log_u = new G4LogicalVolume(SD1_Strip_solid_u, Si, "SD1_Strip_log_u", 0, 0, 0);

        // G4LogicalVolume* SD1_Strip_log_y = new G4LogicalVolume(SD1_Strip_solid_y, Si, "SD1_Strip_log_y", 0, 0, 0);

        for(G4int idStrip = 0; idStrip < HypHI_SD1_Nch; ++idStrip)
          {
            G4double posStrip = -HypHI_SD1_Nch * (HypHI_SD1_stripwidth / 2.) + HypHI_SD1_stripwidth * (0.5 + idStrip);

            AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(0., posStrip, 0.)), SD1_Strip_log_u,
                                                         "SD1_Strip_u", SD1u_MothVol_log_1, false, idStrip));
            AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(0., posStrip, 0.)), SD1_Strip_log_u,
                                                         "SD1_Strip_u", SD1u_MothVol_log_2, false,
                                                         idStrip + 2 * HypHI_SD1_Nch));
          }

        NameDetectorsSD.push_back(SD1_Strip_log_u->GetName());
        // NameDetectorsSD.push_back(Si1_Strip_log_y->GetName());
        if(!MiniVis)
          {
            SD1_Strip_log_u->SetVisAttributes(Si_att);
            SD1u_MothVol_log->SetVisAttributes(G4VisAttributes::GetInvisible());
            SD1u_MothVol_log_1->SetVisAttributes(G4VisAttributes::GetInvisible());
            SD1u_MothVol_log_2->SetVisAttributes(G4VisAttributes::GetInvisible());
          }
        else
          {
            SD1_Strip_log_u->SetVisAttributes(G4VisAttributes::GetInvisible());
            SD1u_MothVol_log->SetVisAttributes(Si_att);
            SD1u_MothVol_log_1->SetVisAttributes(G4VisAttributes::GetInvisible());
            SD1u_MothVol_log_2->SetVisAttributes(G4VisAttributes::GetInvisible());
          }
      }
      {
        G4VSolid* SD1v_MothVol_solid =
            new G4Box("SD1v_solid", HypHI_SD1_length / 2., HypHI_SD1_length / 2., HypHI_SD1_thickness / 2.);
        G4VSolid* SD1v_PCBext_Vol_solid =
            new G4Box("SD1v_PCBext_solid", HypHI_SD1_pcb / 2., HypHI_SD1_pcb / 2., HypHI_SD1_thickness / 2.);
        G4VSolid* SD1v_PCB_Vol_solid =
            new G4SubtractionSolid("SD1v_PCB_solid", SD1v_PCBext_Vol_solid, SD1v_MothVol_solid);
        G4VSolid* SD1v_layer1_MothVol_solid =
            new G4Box("SD1v_layer1_solid", HypHI_SD1_Nch * (HypHI_SD1_stripwidth / 2.), HypHI_SD1_striplength / 2.,
                      HypHI_SD1_thickness / 2.);
        G4VSolid* SD1v_layer2_MothVol_solid =
            new G4Box("SD1v_layer2_solid", HypHI_SD1_Nch * (HypHI_SD1_stripwidth / 2.), HypHI_SD1_striplength / 2.,
                      HypHI_SD1_thickness / 2.);
        G4LogicalVolume* SD1v_MothVol_log = new G4LogicalVolume(SD1v_MothVol_solid, Air, "SD1v_log", 0, 0, 0);
        G4LogicalVolume* SD1v_PCBVol_log  = new G4LogicalVolume(SD1v_PCB_Vol_solid, Air, "SD1v_PCB_log", 0, 0, 0);
        G4LogicalVolume* SD1v_MothVol_log_1 =
            new G4LogicalVolume(SD1v_layer1_MothVol_solid, Air, "SD1v_log_t", 0, 0, 0);
        G4LogicalVolume* SD1v_MothVol_log_2 =
            new G4LogicalVolume(SD1v_layer2_MothVol_solid, Air, "SD1v_log_b", 0, 0, 0);

        G4RotationMatrix* rotStereoSD1_v = nullptr;
        if(Par.IsAvailable("HypHI_SD1_stereoAngle"))
          {
            G4double HypHI_SD1_stereoAngle = Par.Get<double>("HypHI_SD1_stereoAngle");
            rotStereoSD1_v                 = new G4RotationMatrix;
            rotStereoSD1_v->rotateZ(-HypHI_SD1_stereoAngle);
          }
        G4ThreeVector TransSD1v = G4ThreeVector(0., 0., HypHI_SD1_posZ + 0.5 * mm + Systematic_shift) - transMFLD_new;

        AllPlacements.emplace_back(
            new G4PVPlacement(rotStereoSD1_v, Sign * TransSD1v, SD1v_MothVol_log, "Silicon1v", MFLD_log, false, 0));

        if(Par.IsAvailable("HypHI_SD1_PCB"))
          AllPlacements.emplace_back(new G4PVPlacement(rotStereoSD1_v, Sign * TransSD1v, SD1v_PCBVol_log,
                                                       "Silicon1v_PCB", MFLD_log, false, 0));

        AllPlacements.emplace_back(
            new G4PVPlacement(0, Sign * (G4ThreeVector(0., -0.5 * (HypHI_SD1_striplength + HypHI_SD1_gap), 0.)),
                              SD1v_MothVol_log_1, "SD1_v1", SD1v_MothVol_log, false, 0));
        AllPlacements.emplace_back(
            new G4PVPlacement(0, Sign * (G4ThreeVector(0., 0.5 * (HypHI_SD1_striplength + HypHI_SD1_gap), 0.)),
                              SD1v_MothVol_log_2, "SD1_v2", SD1v_MothVol_log, false, 0));

        G4VSolid* SD1_Strip_solid_v      = new G4Box("SD1_Strip_solid_v", HypHI_SD1_stripwidth / 2.,
                                                HypHI_SD1_striplength / 2., HypHI_SD1_thickness / 2.);
        G4LogicalVolume* SD1_Strip_log_v = new G4LogicalVolume(SD1_Strip_solid_v, Si, "SD1_Strip_log_v", 0, 0, 0);

        // G4LogicalVolume* SD1_Strip_log_y = new G4LogicalVolume(SD1_Strip_solid_y, Si, "SD1_Strip_log_y", 0, 0, 0);

        for(G4int idStrip = 0; idStrip < HypHI_SD1_Nch; ++idStrip)
          {
            G4double posStrip = -HypHI_SD1_Nch * (HypHI_SD1_stripwidth / 2.) + HypHI_SD1_stripwidth * (0.5 + idStrip);

            AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(posStrip, 0., 0.)), SD1_Strip_log_v,
                                                         "SD1_Strip_v", SD1v_MothVol_log_1, false, idStrip));
            AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(posStrip, 0., 0.)), SD1_Strip_log_v,
                                                         "SD1_Strip_v", SD1v_MothVol_log_2, false,
                                                         idStrip + 2 * HypHI_SD1_Nch));
          }

        NameDetectorsSD.push_back(SD1_Strip_log_v->GetName());
        // NameDetectorsSD.push_back(Si1_Strip_log_y->GetName());
        if(!MiniVis)
          {
            SD1_Strip_log_v->SetVisAttributes(Si_att);
            SD1v_MothVol_log->SetVisAttributes(G4VisAttributes::GetInvisible());
            SD1v_MothVol_log_1->SetVisAttributes(G4VisAttributes::GetInvisible());
            SD1v_MothVol_log_2->SetVisAttributes(G4VisAttributes::GetInvisible());
          }
        else
          {
            SD1_Strip_log_v->SetVisAttributes(G4VisAttributes::GetInvisible());
            SD1v_MothVol_log->SetVisAttributes(Si_att);
            SD1v_MothVol_log_1->SetVisAttributes(G4VisAttributes::GetInvisible());
            SD1v_MothVol_log_2->SetVisAttributes(G4VisAttributes::GetInvisible());
          }
      }
    }
  if(Par.IsAvailable("HypHI_SD2_On"))
    {
      G4double HypHI_SD2_length      = 98.77 * mm;
      G4double HypHI_SD2_pcb         = 250 * mm;
      G4double HypHI_SD2_thickness   = 0.320 * mm;
      G4double HypHI_SD2_stripwidth  = 0.190 * mm;
      G4double HypHI_SD2_striplength = 48.365 * mm;
      G4double HypHI_SD2_gap         = Par.IsAvailable("HypHI_SD_NoGap") ? 0.0 * mm : 0.5 * mm;
      G4double HypHI_SD2_padding     = (0.235 + 0.510) * mm;
      G4int HypHI_SD2_Nch            = Par.Get<int>("HypHI_SD2_Nch");
      G4int HypHI_SD2_stripPerCh = Par.IsAvailable("HypHI_SD2_stripPerCh") ? Par.Get<int>("HypHI_SD2_stripPerCh") : 1;
      if(HypHI_SD2_Nch * HypHI_SD2_stripPerCh != 512)
        {
          std::cout << "E> SD2 number of channel does not correspond to the design ! " << HypHI_SD2_Nch << " "
                    << HypHI_SD2_stripPerCh << "\n";
          exit(-1);
        }
      if(HypHI_SD2_stripPerCh != 1)
        HypHI_SD2_stripwidth *= static_cast<double>(HypHI_SD2_stripPerCh);

      G4double HypHI_SD2_posZ = Par.Get<double>("HypHI_SD2_posZ");

      {
        G4VSolid* SD2u_MothVol_solid =
            new G4Box("SD2u_solid", HypHI_SD2_length / 2., HypHI_SD2_length / 2., HypHI_SD2_thickness / 2.);
        G4VSolid* SD2u_PCBext_Vol_solid =
            new G4Box("SD2u_PCBext_solid", HypHI_SD2_pcb / 2., HypHI_SD2_pcb / 2., HypHI_SD2_thickness / 2.);
        G4VSolid* SD2u_PCB_Vol_solid =
            new G4SubtractionSolid("SD2u_PCB_solid", SD2u_PCBext_Vol_solid, SD2u_MothVol_solid);
        G4VSolid* SD2u_layer1_MothVol_solid =
            new G4Box("SD2u_layer1_solid", HypHI_SD2_striplength / 2., HypHI_SD2_Nch * (HypHI_SD2_stripwidth / 2.),
                      HypHI_SD2_thickness / 2.);
        G4VSolid* SD2u_layer2_MothVol_solid =
            new G4Box("SD2u_layer2_solid", HypHI_SD2_striplength / 2., HypHI_SD2_Nch * (HypHI_SD2_stripwidth / 2.),
                      HypHI_SD2_thickness / 2.);
        G4LogicalVolume* SD2u_MothVol_log = new G4LogicalVolume(SD2u_MothVol_solid, Air, "SD2u_log", 0, 0, 0);
        G4LogicalVolume* SD2u_PCBVol_log  = new G4LogicalVolume(SD2u_PCB_Vol_solid, Air, "SD2u_PCB_log", 0, 0, 0);
        G4LogicalVolume* SD2u_MothVol_log_1 =
            new G4LogicalVolume(SD2u_layer1_MothVol_solid, Air, "SD2u_log_t", 0, 0, 0);
        G4LogicalVolume* SD2u_MothVol_log_2 =
            new G4LogicalVolume(SD2u_layer2_MothVol_solid, Air, "SD2u_log_b", 0, 0, 0);

        G4RotationMatrix* rotStereoSD2_u = nullptr;
        if(Par.IsAvailable("HypHI_SD2_stereoAngle"))
          {
            G4double HypHI_SD2_stereoAngle = Par.Get<double>("HypHI_SD2_stereoAngle");
            rotStereoSD2_u                 = new G4RotationMatrix;
            rotStereoSD2_u->rotateZ(HypHI_SD2_stereoAngle - 90 * deg);
          }
        G4ThreeVector TransSD2u = G4ThreeVector(0., 0., HypHI_SD2_posZ + Systematic_shift) - transMFLD_new;

        AllPlacements.emplace_back(
            new G4PVPlacement(rotStereoSD2_u, Sign * TransSD2u, SD2u_MothVol_log, "Silicon2u", MFLD_log, false, 0));

        if(Par.IsAvailable("HypHI_SD2_PCB"))
          AllPlacements.emplace_back(new G4PVPlacement(rotStereoSD2_u, Sign * TransSD2u, SD2u_PCBVol_log,
                                                       "Silicon2u_PCB", MFLD_log, false, 0));

        AllPlacements.emplace_back(
            new G4PVPlacement(0, Sign * (G4ThreeVector(-0.5 * (HypHI_SD2_striplength + HypHI_SD2_gap), 0., 0.)),
                              SD2u_MothVol_log_1, "SD2_u1", SD2u_MothVol_log, false, 0));
        AllPlacements.emplace_back(
            new G4PVPlacement(0, Sign * (G4ThreeVector(0.5 * (HypHI_SD2_striplength + HypHI_SD2_gap), 0., 0.)),
                              SD2u_MothVol_log_2, "SD2_u2", SD2u_MothVol_log, false, 0));

        G4VSolid* SD2_Strip_solid_u      = new G4Box("SD2_Strip_solid_u", HypHI_SD2_striplength / 2.,
                                                HypHI_SD2_stripwidth / 2., HypHI_SD2_thickness / 2.);
        G4LogicalVolume* SD2_Strip_log_u = new G4LogicalVolume(SD2_Strip_solid_u, Si, "SD2_Strip_log_u", 0, 0, 0);

        // G4LogicalVolume* SD2_Strip_log_y = new G4LogicalVolume(SD2_Strip_solid_y, Si, "SD2_Strip_log_y", 0, 0,
        // 0);

        for(G4int idStrip = 0; idStrip < HypHI_SD2_Nch; ++idStrip)
          {
            G4double posStrip = -HypHI_SD2_Nch * (HypHI_SD2_stripwidth / 2.) + HypHI_SD2_stripwidth * (0.5 + idStrip);

            AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(0., posStrip, 0.)), SD2_Strip_log_u,
                                                         "SD2_Strip_u", SD2u_MothVol_log_1, false, idStrip));
            AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(0., posStrip, 0.)), SD2_Strip_log_u,
                                                         "SD2_Strip_u", SD2u_MothVol_log_2, false,
                                                         idStrip + 2 * HypHI_SD2_Nch));
          }

        NameDetectorsSD.push_back(SD2_Strip_log_u->GetName());
        // NameDetectorsSD.push_back(Si1_Strip_log_y->GetName());
        if(!MiniVis)
          {
            SD2_Strip_log_u->SetVisAttributes(Si_att);
            SD2u_MothVol_log->SetVisAttributes(G4VisAttributes::GetInvisible());
            SD2u_MothVol_log_1->SetVisAttributes(G4VisAttributes::GetInvisible());
            SD2u_MothVol_log_2->SetVisAttributes(G4VisAttributes::GetInvisible());
          }
        else
          {
            SD2_Strip_log_u->SetVisAttributes(G4VisAttributes::GetInvisible());
            SD2u_MothVol_log->SetVisAttributes(Si_att);
            SD2u_MothVol_log_1->SetVisAttributes(G4VisAttributes::GetInvisible());
            SD2u_MothVol_log_2->SetVisAttributes(G4VisAttributes::GetInvisible());
          }
      }
      {
        G4VSolid* SD2v_MothVol_solid =
            new G4Box("SD2v_solid", HypHI_SD2_length / 2., HypHI_SD2_length / 2., HypHI_SD2_thickness / 2.);
        G4VSolid* SD2v_PCBext_Vol_solid =
            new G4Box("SD2u_PCBext_solid", HypHI_SD2_pcb / 2., HypHI_SD2_pcb / 2., HypHI_SD2_thickness / 2.);
        G4VSolid* SD2v_PCB_Vol_solid =
            new G4SubtractionSolid("SD2u_PCB_solid", SD2v_PCBext_Vol_solid, SD2v_MothVol_solid);
        G4VSolid* SD2v_layer1_MothVol_solid =
            new G4Box("SD2v_layer1_solid", HypHI_SD2_Nch * (HypHI_SD2_stripwidth / 2.), HypHI_SD2_striplength / 2.,
                      HypHI_SD2_thickness / 2.);
        G4VSolid* SD2v_layer2_MothVol_solid =
            new G4Box("SD2v_layer2_solid", HypHI_SD2_Nch * (HypHI_SD2_stripwidth / 2.), HypHI_SD2_striplength / 2.,
                      HypHI_SD2_thickness / 2.);
        G4LogicalVolume* SD2v_MothVol_log = new G4LogicalVolume(SD2v_MothVol_solid, Air, "SD2v_log", 0, 0, 0);
        G4LogicalVolume* SD2v_PCBVol_log  = new G4LogicalVolume(SD2v_PCB_Vol_solid, Air, "SD2v_PCB_log", 0, 0, 0);
        G4LogicalVolume* SD2v_MothVol_log_1 =
            new G4LogicalVolume(SD2v_layer1_MothVol_solid, Air, "SD2v_log_t", 0, 0, 0);
        G4LogicalVolume* SD2v_MothVol_log_2 =
            new G4LogicalVolume(SD2v_layer2_MothVol_solid, Air, "SD2v_log_b", 0, 0, 0);

        G4RotationMatrix* rotStereoSD2_v = nullptr;
        if(Par.IsAvailable("HypHI_SD2_stereoAngle"))
          {
            G4double HypHI_SD2_stereoAngle = Par.Get<double>("HypHI_SD2_stereoAngle");
            rotStereoSD2_v                 = new G4RotationMatrix;
            rotStereoSD2_v->rotateZ(-HypHI_SD2_stereoAngle);
          }
        G4ThreeVector TransSD2v = G4ThreeVector(0., 0., HypHI_SD2_posZ + 0.5 * mm + Systematic_shift) - transMFLD_new;

        AllPlacements.emplace_back(
            new G4PVPlacement(rotStereoSD2_v, Sign * TransSD2v, SD2v_MothVol_log, "Silicon2v", MFLD_log, false, 0));

        if(Par.IsAvailable("HypHI_SD2_PCB"))
          AllPlacements.emplace_back(new G4PVPlacement(rotStereoSD2_v, Sign * TransSD2v, SD2v_PCBVol_log,
                                                       "Silicon2v_PCB", MFLD_log, false, 0));

        AllPlacements.emplace_back(
            new G4PVPlacement(0, Sign * (G4ThreeVector(0., -0.5 * (HypHI_SD2_striplength + HypHI_SD2_gap), 0.)),
                              SD2v_MothVol_log_1, "SD2_v1", SD2v_MothVol_log, false, 0));
        AllPlacements.emplace_back(
            new G4PVPlacement(0, Sign * (G4ThreeVector(0., 0.5 * (HypHI_SD2_striplength + HypHI_SD2_gap), 0.)),
                              SD2v_MothVol_log_2, "SD2_v2", SD2v_MothVol_log, false, 0));

        G4VSolid* SD2_Strip_solid_v      = new G4Box("SD2_Strip_solid_v", HypHI_SD2_stripwidth / 2.,
                                                HypHI_SD2_striplength / 2., HypHI_SD2_thickness / 2.);
        G4LogicalVolume* SD2_Strip_log_v = new G4LogicalVolume(SD2_Strip_solid_v, Si, "SD2_Strip_log_v", 0, 0, 0);

        // G4LogicalVolume* SD2_Strip_log_y = new G4LogicalVolume(SD2_Strip_solid_y, Si, "SD2_Strip_log_y", 0, 0,
        // 0);

        for(G4int idStrip = 0; idStrip < HypHI_SD2_Nch; ++idStrip)
          {
            G4double posStrip = -HypHI_SD2_Nch * (HypHI_SD2_stripwidth / 2.) + HypHI_SD2_stripwidth * (0.5 + idStrip);

            AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(posStrip, 0., 0.)), SD2_Strip_log_v,
                                                         "SD2_Strip_v", SD2v_MothVol_log_1, false, idStrip));
            AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(posStrip, 0., 0.)), SD2_Strip_log_v,
                                                         "SD2_Strip_v", SD2v_MothVol_log_2, false,
                                                         idStrip + 2 * HypHI_SD2_Nch));
          }

        NameDetectorsSD.push_back(SD2_Strip_log_v->GetName());
        // NameDetectorsSD.push_back(Si1_Strip_log_y->GetName());
        if(!MiniVis)
          {
            SD2_Strip_log_v->SetVisAttributes(Si_att);
            SD2v_MothVol_log->SetVisAttributes(G4VisAttributes::GetInvisible());
            SD2v_MothVol_log_1->SetVisAttributes(G4VisAttributes::GetInvisible());
            SD2v_MothVol_log_2->SetVisAttributes(G4VisAttributes::GetInvisible());
          }
        else
          {
            SD2_Strip_log_v->SetVisAttributes(G4VisAttributes::GetInvisible());
            SD2v_MothVol_log->SetVisAttributes(Si_att);
            SD2v_MothVol_log_1->SetVisAttributes(G4VisAttributes::GetInvisible());
            SD2v_MothVol_log_2->SetVisAttributes(G4VisAttributes::GetInvisible());
          }
      }
    }

  if(Par.IsAvailable("HypHI_SDpad1_On"))
    {
      G4double HypHI_SD_length = 20480 * um;
      // G4double HypHI_SD_pcb         = 1650 * um;
      G4double HypHI_SD_thickness   = 285 * um;
      G4double HypHI_SD_stripwidth  = 80 * um;
      G4double HypHI_SD_striplength = 20480 * um; // 20192 * um;
      G4double HypHI_SD_gap =
          Par.IsAvailable("HypHI_SDpad_NoGap")
              ? 0.0 * mm
              : Par.IsAvailable("HypHI_SDpad1_gapSize") ? Par.Get<double>("HypHI_SDpad1_gapSize") : 0.5 * mm;

      G4int HypHI_SD_Nch = Par.Get<int>("HypHI_SDpad1_Nch");
      G4int HypHI_SD_stripPerCh =
          Par.IsAvailable("HypHI_SDpad1_stripPerCh") ? Par.Get<int>("HypHI_SDpad1_stripPerCh") : 1;
      if(HypHI_SD_Nch * HypHI_SD_stripPerCh != 256)
        {
          std::cout << "E> SDpad1 number of channel does not correspond to the design ! " << HypHI_SD_Nch << " "
                    << HypHI_SD_stripPerCh << "\n";
          exit(-1);
        }
      if(HypHI_SD_stripPerCh != 1)
        HypHI_SD_stripwidth *= static_cast<double>(HypHI_SD_stripPerCh);

      G4double HypHI_SD1_posZ = Par.Get<double>("HypHI_SD1_posZ");

      G4int HypHI_SD1_padX = Par.Get<int>("HypHI_SD1_padX");
      G4int HypHI_SD1_padY = Par.Get<int>("HypHI_SD1_padY");

      {
        G4VSolid* SDu_MothVol_solid =
            new G4Box("SD1u_solid", ((HypHI_SD_length + HypHI_SD_gap) * HypHI_SD1_padX) / 2.,
                      ((HypHI_SD_length + HypHI_SD_gap) * HypHI_SD1_padY) / 2., HypHI_SD_thickness / 2.);

        G4LogicalVolume* SDu_MothVol_log = new G4LogicalVolume(SDu_MothVol_solid, Air, "SD1u_log", 0, 0, 0);

        G4RotationMatrix* rotSD1_u = nullptr;
        if(Par.IsAvailable("HypHI_SDpad1_rotate"))
          {
            G4double HypHI_SD1_rotate = Par.Get<double>("HypHI_SDpad1_rotate");
            rotSD1_u                  = new G4RotationMatrix;
            double signRot            = WasaSide == 1 ? -1 : 1;
            rotSD1_u->rotateZ(signRot * HypHI_SD1_rotate);
          }

        G4ThreeVector TransSDu = G4ThreeVector(0., 0., HypHI_SD1_posZ + Systematic_shift) - transMFLD_new;

        AllPlacements.emplace_back(
            new G4PVPlacement(rotSD1_u, Sign * TransSDu, SDu_MothVol_log, "Silicon1u", MFLD_log, false, 0));

        G4VSolid* SD_Strip_solid_u = new G4Box("SD1_Strip_solid_u", HypHI_SD_striplength / 2., HypHI_SD_stripwidth / 2.,
                                               HypHI_SD_thickness / 2.);

        G4LogicalVolume* SD_Strip_log_u = new G4LogicalVolume(SD_Strip_solid_u, Si, "SD1pad_Strip_log_u", 0, 0, 0);

        for(int iX = 0; iX < HypHI_SD1_padX; ++iX)
          {
            std::string nameX("pad_X");
            nameX += std::to_string(iX);
            for(int iY = 0; iY < HypHI_SD1_padY; ++iY)
              {
                std::string nameY("_Y");
                nameY += std::to_string(iY);

                std::string nameSolid("SD1u_solid_");
                nameSolid += nameX;
                nameSolid += nameY;
                G4VSolid* SDpad_solid =
                    new G4Box(nameSolid, HypHI_SD_length / 2., HypHI_SD_striplength / 2., HypHI_SD_thickness / 2.);

                std::string nameLogic("SD1u_logic_");
                nameLogic += nameX;
                nameLogic += nameY;
                G4LogicalVolume* SDu_pad_log_1 = new G4LogicalVolume(SDpad_solid, Air, nameLogic, 0, 0, 0);

                double minDistX            = HypHI_SD_length + 2 * HypHI_SD_gap;
                double minDistY            = HypHI_SD_length + 2 * HypHI_SD_gap;
                G4ThreeVector TransSDu_pad = G4ThreeVector(iX * minDistX - 0.5 * minDistX * (HypHI_SD1_padX - 1),
                                                           iY * minDistY - 0.5 * minDistY * (HypHI_SD1_padY - 1), 0.);

                std::string namePlace("SD1u_pad_");
                namePlace += nameX;
                namePlace += nameY;
                AllPlacements.emplace_back(
                    new G4PVPlacement(0, Sign * (TransSDu_pad), SDu_pad_log_1, namePlace, SDu_MothVol_log, false, 0));

                for(G4int idStrip = 0; idStrip < HypHI_SD_Nch; ++idStrip)
                  {
                    G4double posStrip =
                        -HypHI_SD_Nch * (HypHI_SD_stripwidth / 2.) + HypHI_SD_stripwidth * (0.5 + idStrip);

                    AllPlacements.emplace_back(new G4PVPlacement(
                        0, Sign * (G4ThreeVector(0., posStrip, 0.)), SD_Strip_log_u, "SD1_Strip_u", SDu_pad_log_1,
                        false, idStrip + 2 * (iX + HypHI_SD1_padX * iY) * HypHI_SD_Nch));
                  }
              }
          }

        NameDetectorsSD.push_back(SD_Strip_log_u->GetName());
        if(!MiniVis)
          {
            SD_Strip_log_u->SetVisAttributes(Si_att);
            SDu_MothVol_log->SetVisAttributes(G4VisAttributes::GetInvisible());
          }
        else
          {
            SD_Strip_log_u->SetVisAttributes(G4VisAttributes::GetInvisible());
            SDu_MothVol_log->SetVisAttributes(Si_att);
          }
      }

      {
        G4VSolid* SDv_MothVol_solid =
            new G4Box("SD1v_solid", ((HypHI_SD_length + HypHI_SD_gap) * HypHI_SD1_padX) / 2.,
                      ((HypHI_SD_length + HypHI_SD_gap) * HypHI_SD1_padY) / 2., HypHI_SD_thickness / 2.);

        G4LogicalVolume* SDv_MothVol_log = new G4LogicalVolume(SDv_MothVol_solid, Air, "SD1v_log", 0, 0, 0);

        G4RotationMatrix* rotSD1_v = nullptr;
        if(Par.IsAvailable("HypHI_SDpad1_rotate"))
          {
            G4double HypHI_SD1_rotate = Par.Get<double>("HypHI_SDpad1_rotate");
            rotSD1_v                  = new G4RotationMatrix;
            double signRot            = WasaSide == 1 ? -1 : 1;
            rotSD1_v->rotateZ(signRot * HypHI_SD1_rotate);
          }

        G4ThreeVector TransSDv = G4ThreeVector(0., 0., HypHI_SD1_posZ + 0.5 * mm + Systematic_shift) - transMFLD_new;

        AllPlacements.emplace_back(
            new G4PVPlacement(rotSD1_v, Sign * TransSDv, SDv_MothVol_log, "Silicon1v", MFLD_log, false, 0));

        G4VSolid* SD_Strip_solid_v = new G4Box("SD1_Strip_solid_v", HypHI_SD_stripwidth / 2., HypHI_SD_striplength / 2.,
                                               HypHI_SD_thickness / 2.);

        G4LogicalVolume* SD_Strip_log_v = new G4LogicalVolume(SD_Strip_solid_v, Si, "SD1pad_Strip_log_v", 0, 0, 0);

        for(int iX = 0; iX < HypHI_SD1_padX; ++iX)
          {
            std::string nameX("pad_X");
            nameX += std::to_string(iX);
            for(int iY = 0; iY < HypHI_SD1_padY; ++iY)
              {
                std::string nameY("_Y");
                nameY += std::to_string(iY);

                std::string nameSolid("SD1v_solid_");
                nameSolid += nameX;
                nameSolid += nameY;
                G4VSolid* SDpad_solid =
                    new G4Box(nameSolid, HypHI_SD_length / 2., HypHI_SD_striplength / 2., HypHI_SD_thickness / 2.);

                std::string nameLogic("SD1v_logic_");
                nameLogic += nameX;
                nameLogic += nameY;
                G4LogicalVolume* SDv_pad_log_1 = new G4LogicalVolume(SDpad_solid, Air, nameLogic, 0, 0, 0);

                double minDistX            = HypHI_SD_length + 2 * HypHI_SD_gap;
                double minDistY            = HypHI_SD_length + 2 * HypHI_SD_gap;
                G4ThreeVector TransSDv_pad = G4ThreeVector(iX * minDistX - 0.5 * minDistX * (HypHI_SD1_padX - 1),
                                                           iY * minDistY - 0.5 * minDistY * (HypHI_SD1_padY - 1), 0.);

                std::string namePlace("SD1v_pad_");
                namePlace += nameX;
                namePlace += nameY;
                AllPlacements.emplace_back(
                    new G4PVPlacement(0, Sign * (TransSDv_pad), SDv_pad_log_1, namePlace, SDv_MothVol_log, false, 0));

                for(G4int idStrip = 0; idStrip < HypHI_SD_Nch; ++idStrip)
                  {
                    G4double posStrip =
                        -HypHI_SD_Nch * (HypHI_SD_stripwidth / 2.) + HypHI_SD_stripwidth * (0.5 + idStrip);

                    AllPlacements.emplace_back(new G4PVPlacement(
                        0, Sign * (G4ThreeVector(posStrip, 0., 0.)), SD_Strip_log_v, "SD1_Strip_v", SDv_pad_log_1,
                        false, idStrip + 2 * (iX + HypHI_SD1_padX * iY) * HypHI_SD_Nch));
                  }
              }
          }

        NameDetectorsSD.push_back(SD_Strip_log_v->GetName());
        if(!MiniVis)
          {
            SD_Strip_log_v->SetVisAttributes(Si_att);
            SDv_MothVol_log->SetVisAttributes(G4VisAttributes::GetInvisible());
          }
        else
          {
            SD_Strip_log_v->SetVisAttributes(G4VisAttributes::GetInvisible());
            SDv_MothVol_log->SetVisAttributes(Si_att);
          }
      }
    }

  if(Par.IsAvailable("HypHI_SDpad2_On"))
    {
      G4double HypHI_SD_length = 20480 * um;
      // G4double HypHI_SD_pcb         = 1650 * um;
      G4double HypHI_SD_thickness   = 285 * um;
      G4double HypHI_SD_stripwidth  = 80 * um;
      G4double HypHI_SD_striplength = 20480 * um; // 20192 * um;
      G4double HypHI_SD_gap =
          Par.IsAvailable("HypHI_SDpad_NoGap")
              ? 0.0 * mm
              : Par.IsAvailable("HypHI_SDpad2_gapSize") ? Par.Get<double>("HypHI_SDpad2_gapSize") : 0.5 * mm;

      G4int HypHI_SD_Nch = Par.Get<int>("HypHI_SDpad2_Nch");
      G4int HypHI_SD_stripPerCh =
          Par.IsAvailable("HypHI_SDpad2_stripPerCh") ? Par.Get<int>("HypHI_SDpad2_stripPerCh") : 1;
      if(HypHI_SD_Nch * HypHI_SD_stripPerCh != 256)
        {
          std::cout << "E> SDpad2 number of channel does not correspond to the design ! " << HypHI_SD_Nch << " "
                    << HypHI_SD_stripPerCh << "\n";
          exit(-1);
        }
      if(HypHI_SD_stripPerCh != 1)
        HypHI_SD_stripwidth *= static_cast<double>(HypHI_SD_stripPerCh);

      G4double HypHI_SD2_posZ = Par.Get<double>("HypHI_SD2_posZ");

      G4int HypHI_SD2_padX = Par.Get<int>("HypHI_SD2_padX");
      G4int HypHI_SD2_padY = Par.Get<int>("HypHI_SD2_padY");

      {
        G4VSolid* SDu_MothVol_solid =
            new G4Box("SD2u_solid", ((HypHI_SD_length + HypHI_SD_gap) * HypHI_SD2_padX) / 2.,
                      ((HypHI_SD_length + HypHI_SD_gap) * HypHI_SD2_padY) / 2., HypHI_SD_thickness / 2.);

        G4LogicalVolume* SDu_MothVol_log = new G4LogicalVolume(SDu_MothVol_solid, Air, "SD2u_log", 0, 0, 0);

        G4RotationMatrix* rotSD2_u = nullptr;
        if(Par.IsAvailable("HypHI_SDpad2_rotate"))
          {
            G4double HypHI_SD2_rotate = Par.Get<double>("HypHI_SDpad2_rotate");
            rotSD2_u                  = new G4RotationMatrix;
            double signRot            = WasaSide == 1 ? -1 : 1;
            rotSD2_u->rotateZ(signRot * HypHI_SD2_rotate);
          }

        G4ThreeVector TransSDu = G4ThreeVector(0., 0., HypHI_SD2_posZ + Systematic_shift) - transMFLD_new;

        AllPlacements.emplace_back(
            new G4PVPlacement(rotSD2_u, Sign * TransSDu, SDu_MothVol_log, "Silicon2u", MFLD_log, false, 0));

        G4VSolid* SD_Strip_solid_u = new G4Box("SD2_Strip_solid_u", HypHI_SD_striplength / 2., HypHI_SD_stripwidth / 2.,
                                               HypHI_SD_thickness / 2.);

        G4LogicalVolume* SD_Strip_log_u = new G4LogicalVolume(SD_Strip_solid_u, Si, "SD2pad_Strip_log_u", 0, 0, 0);

        for(int iX = 0; iX < HypHI_SD2_padX; ++iX)
          {
            std::string nameX("pad_X");
            nameX += std::to_string(iX);
            for(int iY = 0; iY < HypHI_SD2_padY; ++iY)
              {
                std::string nameY("_Y");
                nameY += std::to_string(iY);

                std::string nameSolid("SD2u_solid_");
                nameSolid += nameX;
                nameSolid += nameY;
                G4VSolid* SDpad_solid =
                    new G4Box(nameSolid, HypHI_SD_length / 2., HypHI_SD_striplength / 2., HypHI_SD_thickness / 2.);

                std::string nameLogic("SD2u_logic_");
                nameLogic += nameX;
                nameLogic += nameY;
                G4LogicalVolume* SDu_pad_log_1 = new G4LogicalVolume(SDpad_solid, Air, nameLogic, 0, 0, 0);

                double minDistX            = HypHI_SD_length + 2 * HypHI_SD_gap;
                double minDistY            = HypHI_SD_length + 2 * HypHI_SD_gap;
                G4ThreeVector TransSDu_pad = G4ThreeVector(iX * minDistX - 0.5 * minDistX * (HypHI_SD2_padX - 1),
                                                           iY * minDistY - 0.5 * minDistY * (HypHI_SD2_padY - 1), 0.);

                std::string namePlace("SD2u_pad_");
                namePlace += nameX;
                namePlace += nameY;
                AllPlacements.emplace_back(
                    new G4PVPlacement(0, Sign * (TransSDu_pad), SDu_pad_log_1, namePlace, SDu_MothVol_log, false, 0));

                for(G4int idStrip = 0; idStrip < HypHI_SD_Nch; ++idStrip)
                  {
                    G4double posStrip =
                        -HypHI_SD_Nch * (HypHI_SD_stripwidth / 2.) + HypHI_SD_stripwidth * (0.5 + idStrip);

                    AllPlacements.emplace_back(new G4PVPlacement(
                        0, Sign * (G4ThreeVector(0., posStrip, 0.)), SD_Strip_log_u, "SD2_Strip_u", SDu_pad_log_1,
                        false, idStrip + 2 * (iX + HypHI_SD2_padX * iY) * HypHI_SD_Nch));
                  }
              }
          }

        NameDetectorsSD.push_back(SD_Strip_log_u->GetName());
        if(!MiniVis)
          {
            SD_Strip_log_u->SetVisAttributes(Si_att);
            SDu_MothVol_log->SetVisAttributes(G4VisAttributes::GetInvisible());
          }
        else
          {
            SD_Strip_log_u->SetVisAttributes(G4VisAttributes::GetInvisible());
            SDu_MothVol_log->SetVisAttributes(Si_att);
          }
      }

      {
        G4VSolid* SDv_MothVol_solid =
            new G4Box("SD2v_solid", ((HypHI_SD_length + HypHI_SD_gap) * HypHI_SD2_padX) / 2.,
                      ((HypHI_SD_length + HypHI_SD_gap) * HypHI_SD2_padY) / 2., HypHI_SD_thickness / 2.);

        G4LogicalVolume* SDv_MothVol_log = new G4LogicalVolume(SDv_MothVol_solid, Air, "SD2v_log", 0, 0, 0);

        G4RotationMatrix* rotSD2_v = nullptr;
        if(Par.IsAvailable("HypHI_SDpad2_rotate"))
          {
            G4double HypHI_SD2_rotate = Par.Get<double>("HypHI_SDpad2_rotate");
            rotSD2_v                  = new G4RotationMatrix;
            double signRot            = WasaSide == 1 ? -1 : 1;
            rotSD2_v->rotateZ(signRot * HypHI_SD2_rotate);
          }

        G4ThreeVector TransSDv = G4ThreeVector(0., 0., HypHI_SD2_posZ + 0.5 * mm + Systematic_shift) - transMFLD_new;

        AllPlacements.emplace_back(
            new G4PVPlacement(rotSD2_v, Sign * TransSDv, SDv_MothVol_log, "Silicon2v", MFLD_log, false, 0));

        G4VSolid* SD_Strip_solid_v = new G4Box("SD2_Strip_solid_v", HypHI_SD_stripwidth / 2., HypHI_SD_striplength / 2.,
                                               HypHI_SD_thickness / 2.);

        G4LogicalVolume* SD_Strip_log_v = new G4LogicalVolume(SD_Strip_solid_v, Si, "SD2pad_Strip_log_v", 0, 0, 0);

        for(int iX = 0; iX < HypHI_SD2_padX; ++iX)
          {
            std::string nameX("pad_X");
            nameX += std::to_string(iX);
            for(int iY = 0; iY < HypHI_SD2_padY; ++iY)
              {
                std::string nameY("_Y");
                nameY += std::to_string(iY);

                std::string nameSolid("SD2v_solid_");
                nameSolid += nameX;
                nameSolid += nameY;
                G4VSolid* SDpad_solid =
                    new G4Box(nameSolid, HypHI_SD_length / 2., HypHI_SD_striplength / 2., HypHI_SD_thickness / 2.);

                std::string nameLogic("SD2v_logic_");
                nameLogic += nameX;
                nameLogic += nameY;
                G4LogicalVolume* SDv_pad_log_1 = new G4LogicalVolume(SDpad_solid, Air, nameLogic, 0, 0, 0);

                double minDistX            = HypHI_SD_length + 2 * HypHI_SD_gap;
                double minDistY            = HypHI_SD_length + 2 * HypHI_SD_gap;
                G4ThreeVector TransSDv_pad = G4ThreeVector(iX * minDistX - 0.5 * minDistX * (HypHI_SD2_padX - 1),
                                                           iY * minDistY - 0.5 * minDistY * (HypHI_SD2_padY - 1), 0.);

                std::string namePlace("SD2v_pad_");
                namePlace += nameX;
                namePlace += nameY;
                AllPlacements.emplace_back(
                    new G4PVPlacement(0, Sign * (TransSDv_pad), SDv_pad_log_1, namePlace, SDv_MothVol_log, false, 0));

                for(G4int idStrip = 0; idStrip < HypHI_SD_Nch; ++idStrip)
                  {
                    G4double posStrip =
                        -HypHI_SD_Nch * (HypHI_SD_stripwidth / 2.) + HypHI_SD_stripwidth * (0.5 + idStrip);

                    AllPlacements.emplace_back(new G4PVPlacement(
                        0, Sign * (G4ThreeVector(posStrip, 0., 0.)), SD_Strip_log_v, "SD2_Strip_v", SDv_pad_log_1,
                        false, idStrip + 2 * (iX + HypHI_SD2_padX * iY) * HypHI_SD_Nch));
                  }
              }
          }
        NameDetectorsSD.push_back(SD_Strip_log_v->GetName());
        if(!MiniVis)
          {
            SD_Strip_log_v->SetVisAttributes(Si_att);
            SDv_MothVol_log->SetVisAttributes(G4VisAttributes::GetInvisible());
          }
        else
          {
            SD_Strip_log_v->SetVisAttributes(G4VisAttributes::GetInvisible());
            SDv_MothVol_log->SetVisAttributes(Si_att);
          }
      }
    }

  // ----------------------------- @@ -----------------------------
  //         Virtual Detectors
  // ----------------------------- @@ -----------------------------

  if(Par.IsAvailable("HypHI_VirtualTR_On"))
    {
      const double TR1_posZ = Par.Get<double>("HypHI_TR1_posZ");
      G4VSolid* TR1_box     = nullptr;
      if(Par.IsAvailable("HypHI_BeamHole"))
        {
          G4VSolid* TR1_box_init = new G4Box("TR1_box_init", 15. * cm, 15. * cm, 1. * mm);
          G4VSolid* TR1_box_hole = new G4Box("TR1_box", BeamHoleSize * 0.5, BeamHoleSize * 0.5, 1. * mm);
          TR1_box                = new G4SubtractionSolid("TR1_box", TR1_box_init, TR1_box_hole);
        }
      else
        TR1_box = new G4Box("TR1_box", 15. * cm, 15. * cm, 1. * mm);

      G4LogicalVolume* TR1_log = new G4LogicalVolume(TR1_box, Scinti, "TR1_log", 0, 0, 0);
      AllPlacements.emplace_back(
          new G4PVPlacement(0, Sign * (G4ThreeVector(0., 0., TR1_posZ + Systematic_shift) - transMFLD_new), TR1_log,
                            "TR1_phys", MFLD_log, false, 0));
      NameDetectorsSD.push_back(TR1_log->GetName());

      const double TR2_posZ = Par.Get<double>("HypHI_TR2_posZ");

      G4VSolid* TR2_box = nullptr;
      if(Par.IsAvailable("HypHI_BeamHole"))
        {
          G4VSolid* TR2_box_init = new G4Box("TR2_box_init", 15. * cm, 15. * cm, 1. * mm);
          G4VSolid* TR2_box_hole = new G4Box("TR2_box_hole", BeamHoleSize * 0.5, BeamHoleSize * 0.5, 1. * mm);
          TR2_box                = new G4SubtractionSolid("TR2_box", TR2_box_init, TR2_box_hole);
        }
      else
        TR2_box = new G4Box("TR2_box", 15. * cm, 15. * cm, 1. * mm);

      G4LogicalVolume* TR2_log = new G4LogicalVolume(TR2_box, Scinti, "TR2_log", 0, 0, 0);
      AllPlacements.emplace_back(
          new G4PVPlacement(0, Sign * (G4ThreeVector(0., 0., TR2_posZ + Systematic_shift) - transMFLD_new), TR2_log,
                            "TR2_phys", MFLD_log, false, 0));
      NameDetectorsSD.push_back(TR2_log->GetName());
    }

  // ----------------------------- @@ -----------------------------
  // 		     Fiber Detectors
  // ----------------------------- @@ -----------------------------
  if(Par.IsAvailable("HypHI_FiberTR_On"))
    {

      const double HypHI_FiberTracker1_posZ = Par.Get<double>("HypHI_FiberTracker1_posZ");
      const double HypHI_FiberTracker2_posZ = Par.Get<double>("HypHI_FiberTracker2_posZ");
      const double HypHI_FiberTracker3_posZ = Par.Get<double>("HypHI_FiberTracker3_posZ");
      const double HypHI_FiberTracker4_posZ = Par.Get<double>("HypHI_FiberTracker4_posZ");
      const double HypHI_FiberTracker5_posZ = Par.Get<double>("HypHI_FiberTracker5_posZ");

      G4RotationMatrix* rotFib1 = new G4RotationMatrix;
      rotFib1->rotateX(90. * deg);
      G4RotationMatrix* rotFib2 = new G4RotationMatrix;
      rotFib2->rotateZ(30. * deg);
      rotFib2->rotateX(90. * deg);
      G4RotationMatrix* rotFib3 = new G4RotationMatrix;
      rotFib3->rotateZ(-30. * deg);
      rotFib3->rotateX(90. * deg);

      const G4double spacingX               = 0.55 * mm;
      const G4double startZ1                = 0. * mm;
      const G4double startZ2                = 0.47631397 * mm;
      const std::vector<G4double> posZshift = {-4. * mm, 0. * mm, 4. * mm};

      G4ThreeVector posFibX;
      G4ThreeVector posFibUV;

      G4ThreeVector shiftFibX  = G4ThreeVector(spacingX, 0., 0.);
      G4ThreeVector shiftFibUV = G4ThreeVector(spacingX / cos(30. * deg), 0., 0.);

      G4VisAttributes* visAttributes_x = new G4VisAttributes(color_x);
      G4VisAttributes* visAttributes_u = new G4VisAttributes(color_u);
      G4VisAttributes* visAttributes_v = new G4VisAttributes(color_v);



      // -------------------------- First Fiber Detector --------------------------
      G4VSolid* FiberD1_MothVol_solid        = new G4Box("FiberDetector1", 15. * cm, 20. * cm, 6. * mm);
      G4VSolid* FiberD1_layerX_MothVol_solid = new G4Box("FiberD1_layerX_solid", 15. * cm, 20. * cm, 2. * mm);
      G4VSolid* FiberD1_layerU_MothVol_solid = new G4Box("FiberD1_layerU_solid", 15. * cm, 20. * cm, 2. * mm);
      G4VSolid* FiberD1_layerV_MothVol_solid = new G4Box("FiberD1_layerV_solid", 15. * cm, 20. * cm, 2. * mm);
      G4LogicalVolume* FiberD1_MothVol_log   = new G4LogicalVolume(FiberD1_MothVol_solid, Air, "FiberD1_log", 0, 0, 0);
      G4LogicalVolume* FiberD1_MothVol_log_x =
          new G4LogicalVolume(FiberD1_layerX_MothVol_solid, Air, "FiberD1_log_x", 0, 0, 0);
      G4LogicalVolume* FiberD1_MothVol_log_u =
          new G4LogicalVolume(FiberD1_layerU_MothVol_solid, Air, "FiberD1_log_u", 0, 0, 0);
      G4LogicalVolume* FiberD1_MothVol_log_v =
          new G4LogicalVolume(FiberD1_layerV_MothVol_solid, Air, "FiberD1_log_v", 0, 0, 0);

      AllPlacements.emplace_back(new G4PVPlacement(
          0, Sign * (G4ThreeVector(fiber_uft1_pos_x, fiber_uft1_pos_y*Sign, HypHI_FiberTracker1_posZ + Systematic_shift) - transMFLD_new),
          FiberD1_MothVol_log, "FiberDetector1", MFLD_log, false, 0));
      AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(0., 0., posZshift[0])),
                                                   FiberD1_MothVol_log_x, "FiberD1_x", FiberD1_MothVol_log, false, 0));
      AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(0., 0., posZshift[1])),
                                                   FiberD1_MothVol_log_u, "FiberD1_u", FiberD1_MothVol_log, false, 0));
      AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(0., 0., posZshift[2])),
                                                   FiberD1_MothVol_log_v, "FiberD1_v", FiberD1_MothVol_log, false, 0));

      const G4double FD1_startX1 = -69.9875 * mm;
      const G4double FD1_startX2 = -70.2625 * mm;

      G4VSolid* FiberD1_Core_solid_x =
          new G4Tubs("FiberD1_Core_solid_x", 0, 0.24 * mm, 10.5 * cm, 0. * deg, 360. * deg);
      G4VSolid* FiberD1_Core_solid_u =
          new G4Tubs("FiberD1_Core_solid_u", 0, 0.24 * mm, (10.5 / sqrt(3) * 2) * cm, 0. * deg, 360. * deg);
      G4VSolid* FiberD1_Core_solid_v =
          new G4Tubs("FiberD1_Core_solid_v", 0, 0.24 * mm, (10.5 / sqrt(3) * 2) * cm, 0. * deg, 360. * deg);
      G4LogicalVolume* FiberD1_Core_log_x =
          new G4LogicalVolume(FiberD1_Core_solid_x, FiberCoreScinti, "FiberD1_Core_log_x", 0, 0, 0);
      G4LogicalVolume* FiberD1_Core_log_u =
          new G4LogicalVolume(FiberD1_Core_solid_u, FiberCoreScinti, "FiberD1_Core_log_u", 0, 0, 0);
      G4LogicalVolume* FiberD1_Core_log_v =
          new G4LogicalVolume(FiberD1_Core_solid_v, FiberCoreScinti, "FiberD1_Core_log_v", 0, 0, 0);

      G4VSolid* FiberD1_Cladding_solid_x =
          new G4Tubs("FiberD1_Cladding_solid_x", 0.24 * mm, 0.25 * mm, 10.5 * cm, 0. * deg, 360. * deg);
      G4VSolid* FiberD1_Cladding_solid_u =
          new G4Tubs("FiberD1_Cladding_solid_u", 0.24 * mm, 0.25 * mm, (10.5 / sqrt(3) * 2) * cm, 0. * deg, 360. * deg);
      G4VSolid* FiberD1_Cladding_solid_v =
          new G4Tubs("FiberD1_Cladding_solid_v", 0.24 * mm, 0.25 * mm, (10.5 / sqrt(3) * 2) * cm, 0. * deg, 360. * deg);
      G4LogicalVolume* FiberD1_Cladding_log_x =
          new G4LogicalVolume(FiberD1_Cladding_solid_x, Scinti, "FiberD1_Cladding_log_x", 0, 0, 0);
      G4LogicalVolume* FiberD1_Cladding_log_u =
          new G4LogicalVolume(FiberD1_Cladding_solid_u, Scinti, "FiberD1_Cladding_log_u", 0, 0, 0);
      G4LogicalVolume* FiberD1_Cladding_log_v =
          new G4LogicalVolume(FiberD1_Cladding_solid_v, Scinti, "FiberD1_Cladding_log_v", 0, 0, 0);

      for(G4int IdUD = 0; IdUD < 2; ++IdUD)
        for(G4int IdFib = 0; IdFib < 128; ++IdFib)
          {
            if(IdUD == 0)
              {
                posFibX  = G4ThreeVector( FD1_startX1 + IdFib * 2.*spacingX                  , 0., startZ1);
                posFibUV = G4ThreeVector((FD1_startX1 + IdFib * 2.*spacingX) / cos(30. * deg), 0., startZ1);
              }
            else
              {
                posFibX  = G4ThreeVector( FD1_startX2 + IdFib * 2.* spacingX                  , 0., startZ2);
                posFibUV = G4ThreeVector((FD1_startX2 + IdFib * 2.* spacingX) / cos(30. * deg), 0., startZ2);
              }

            int fiberID = IdFib*4 + (1-IdUD)*2;

            G4ThreeVector posFib_x1;
            G4ThreeVector posFib_u1;
            G4ThreeVector posFib_v1;
            G4ThreeVector posFib_x2;
            G4ThreeVector posFib_u2;
            G4ThreeVector posFib_v2;

            double off_x_buf = fiber_uft1_off_x + fiber_offset[0][0][fiberID/2];
            double off_u_buf = fiber_uft1_off_u + fiber_offset[0][1][fiberID/2];
            double off_v_buf = fiber_uft1_off_v + fiber_offset[0][2][fiberID/2];
            off_u_buf /= cos(30. * deg);
            off_v_buf /= cos(30. * deg);
            double pos_x1_buf = posFibX.x()  + off_x_buf;
            double pos_x2_buf = posFibX.x()  + off_x_buf + spacingX;
            double pos_u1_buf = posFibUV.x() + off_u_buf;
            double pos_u2_buf = posFibUV.x() + off_u_buf + spacingX / cos(30. * deg);
            double pos_v1_buf = posFibUV.x() + off_v_buf;
            double pos_v2_buf = posFibUV.x() + off_v_buf + spacingX / cos(30. * deg);
            posFib_x1 = G4ThreeVector(pos_x1_buf, 0., posFibX.z());
            posFib_x2 = G4ThreeVector(pos_x2_buf, 0., posFibX.z());
            posFib_u1 = G4ThreeVector(pos_u1_buf, 0., posFibUV.z());
            posFib_u2 = G4ThreeVector(pos_u2_buf, 0., posFibUV.z());
            posFib_v1 = G4ThreeVector(pos_v1_buf, 0., posFibUV.z());
            posFib_v2 = G4ThreeVector(pos_v2_buf, 0., posFibUV.z());

            AllPlacements.emplace_back(new G4PVPlacement(rotFib1, Sign * posFib_x1, FiberD1_Cladding_log_x, "FiberD1_Cladding_x",
                                                         FiberD1_MothVol_log_x, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib1, Sign * posFib_x1, FiberD1_Core_log_x    , "FiberD1_Core_x",
                                                         FiberD1_MothVol_log_x, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib2, Sign * posFib_u1, FiberD1_Cladding_log_u, "FiberD1_Cladding_u",
                                                         FiberD1_MothVol_log_u, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib2, Sign * posFib_u1, FiberD1_Core_log_u    , "FiberD1_Core_u",
                                                         FiberD1_MothVol_log_u, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib3, Sign * posFib_v1, FiberD1_Cladding_log_v, "FiberD1_Cladding_v",
                                                         FiberD1_MothVol_log_v, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib3, Sign * posFib_v1, FiberD1_Core_log_v    , "FiberD1_Core_v",
                                                         FiberD1_MothVol_log_v, false, fiberID));

            //posFibX  += shiftFibX;
            //posFibUV += shiftFibUV;
            AllPlacements.emplace_back(new G4PVPlacement(rotFib1, Sign * posFib_x2, FiberD1_Cladding_log_x, "FiberD1_Cladding_x",
                                                         FiberD1_MothVol_log_x, false, fiberID + 1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib1, Sign * posFib_x2, FiberD1_Core_log_x    , "FiberD1_Core_x",
                                                         FiberD1_MothVol_log_x, false, fiberID + 1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib2, Sign * posFib_u2, FiberD1_Cladding_log_u, "FiberD1_Cladding_u",
                                                         FiberD1_MothVol_log_u, false, fiberID + 1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib2, Sign * posFib_u2, FiberD1_Core_log_u    , "FiberD1_Core_u",
                                                         FiberD1_MothVol_log_u, false, fiberID + 1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib3, Sign * posFib_v2, FiberD1_Cladding_log_v, "FiberD1_Cladding_v",
                                                         FiberD1_MothVol_log_v, false, fiberID + 1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib3, Sign * posFib_v2, FiberD1_Core_log_v    , "FiberD1_Core_v",
                                                         FiberD1_MothVol_log_v, false, fiberID + 1));
          }

      NameDetectorsSD.push_back(FiberD1_Core_log_x->GetName());
      NameDetectorsSD.push_back(FiberD1_Core_log_u->GetName());
      NameDetectorsSD.push_back(FiberD1_Core_log_v->GetName());
      if(!MiniVis)
        {
          FiberD1_Core_log_x->SetVisAttributes(visAttributes_x);
          FiberD1_Core_log_u->SetVisAttributes(visAttributes_u);
          FiberD1_Core_log_v->SetVisAttributes(visAttributes_v);
          FiberD1_Cladding_log_x->SetVisAttributes(visAttributes_x);
          FiberD1_Cladding_log_u->SetVisAttributes(visAttributes_u);
          FiberD1_Cladding_log_v->SetVisAttributes(visAttributes_v);
          FiberD1_MothVol_log->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD1_MothVol_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD1_MothVol_log_u->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD1_MothVol_log_v->SetVisAttributes(G4VisAttributes::GetInvisible());
        }
      else
        {
          FiberD1_Core_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD1_Core_log_u->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD1_Core_log_v->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD1_Cladding_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD1_Cladding_log_u->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD1_Cladding_log_v->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD1_MothVol_log->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD1_MothVol_log_x->SetVisAttributes(visAttributes_x);
          FiberD1_MothVol_log_u->SetVisAttributes(visAttributes_u);
          FiberD1_MothVol_log_v->SetVisAttributes(visAttributes_v);
        }



      // -------------------------- Second Fiber Detector --------------------------
      G4VSolid* FiberD2_MothVol_solid        = new G4Box("FiberDetector2", 15. * cm, 20. * cm, 6. * mm);
      G4VSolid* FiberD2_layerX_MothVol_solid = new G4Box("FiberD2_layerX_solid", 15. * cm, 20. * cm, 2. * mm);
      G4VSolid* FiberD2_layerU_MothVol_solid = new G4Box("FiberD2_layerU_solid", 15. * cm, 20. * cm, 2. * mm);
      G4VSolid* FiberD2_layerV_MothVol_solid = new G4Box("FiberD2_layerV_solid", 15. * cm, 20. * cm, 2. * mm);
      G4LogicalVolume* FiberD2_MothVol_log   = new G4LogicalVolume(FiberD2_MothVol_solid, Air, "FiberD2_log", 0, 0, 0);
      G4LogicalVolume* FiberD2_MothVol_log_x =
          new G4LogicalVolume(FiberD2_layerX_MothVol_solid, Air, "FiberD2_log_x", 0, 0, 0);
      G4LogicalVolume* FiberD2_MothVol_log_u =
          new G4LogicalVolume(FiberD2_layerU_MothVol_solid, Air, "FiberD2_log_u", 0, 0, 0);
      G4LogicalVolume* FiberD2_MothVol_log_v =
          new G4LogicalVolume(FiberD2_layerV_MothVol_solid, Air, "FiberD2_log_v", 0, 0, 0);

      AllPlacements.emplace_back(new G4PVPlacement(
          0, Sign * (G4ThreeVector(fiber_uft2_pos_x, fiber_uft2_pos_y*Sign, HypHI_FiberTracker2_posZ + Systematic_shift) - transMFLD_new),
          FiberD2_MothVol_log, "FiberDetector2", MFLD_log, false, 0));
      AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(0., 0., posZshift[0])),
                                                   FiberD2_MothVol_log_x, "FiberD2_x", FiberD2_MothVol_log, false, 0));
      AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(0., 0., posZshift[1])),
                                                   FiberD2_MothVol_log_u, "FiberD2_u", FiberD2_MothVol_log, false, 0));
      AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(0., 0., posZshift[2])),
                                                   FiberD2_MothVol_log_v, "FiberD2_v", FiberD2_MothVol_log, false, 0));

      const G4double FD2_startX1 = -69.9875 * mm;
      const G4double FD2_startX2 = -70.2625 * mm;

      G4VSolid* FiberD2_Core_solid_x =
          new G4Tubs("FiberD2_Core_solid_x", 0, 0.24 * mm, 10.5 * cm, 0. * deg, 360. * deg);
      G4VSolid* FiberD2_Core_solid_u =
          new G4Tubs("FiberD2_Core_solid_u", 0, 0.24 * mm, (10.5 / sqrt(3) * 2) * cm, 0. * deg, 360. * deg);
      G4VSolid* FiberD2_Core_solid_v =
          new G4Tubs("FiberD2_Core_solid_v", 0, 0.24 * mm, (10.5 / sqrt(3) * 2) * cm, 0. * deg, 360. * deg);
      G4LogicalVolume* FiberD2_Core_log_x =
          new G4LogicalVolume(FiberD2_Core_solid_x, FiberCoreScinti, "FiberD2_Core_log_x", 0, 0, 0);
      G4LogicalVolume* FiberD2_Core_log_u =
          new G4LogicalVolume(FiberD2_Core_solid_u, FiberCoreScinti, "FiberD2_Core_log_u", 0, 0, 0);
      G4LogicalVolume* FiberD2_Core_log_v =
          new G4LogicalVolume(FiberD2_Core_solid_v, FiberCoreScinti, "FiberD2_Core_log_v", 0, 0, 0);

      G4VSolid* FiberD2_Cladding_solid_x =
          new G4Tubs("FiberD2_Cladding_solid_x", 0.24 * mm, 0.25 * mm, 10.5 * cm, 0. * deg, 360. * deg);
      G4VSolid* FiberD2_Cladding_solid_u =
          new G4Tubs("FiberD2_Cladding_solid_u", 0.24 * mm, 0.25 * mm, (10.5 / sqrt(3) * 2) * cm, 0. * deg, 360. * deg);
      G4VSolid* FiberD2_Cladding_solid_v =
          new G4Tubs("FiberD2_Cladding_solid_v", 0.24 * mm, 0.25 * mm, (10.5 / sqrt(3) * 2) * cm, 0. * deg, 360. * deg);
      G4LogicalVolume* FiberD2_Cladding_log_x =
          new G4LogicalVolume(FiberD2_Cladding_solid_x, Scinti, "FiberD2_Cladding_log_x", 0, 0, 0);
      G4LogicalVolume* FiberD2_Cladding_log_u =
          new G4LogicalVolume(FiberD2_Cladding_solid_u, Scinti, "FiberD2_Cladding_log_u", 0, 0, 0);
      G4LogicalVolume* FiberD2_Cladding_log_v =
          new G4LogicalVolume(FiberD2_Cladding_solid_v, Scinti, "FiberD2_Cladding_log_v", 0, 0, 0);

      for(G4int IdUD = 0; IdUD < 2; ++IdUD)
        for(G4int IdFib = 0; IdFib < 128; ++IdFib)
          {
            if(IdUD == 0)
              {
                posFibX  = G4ThreeVector( FD2_startX1 + IdFib * 2.* spacingX                  , 0., startZ1);
                posFibUV = G4ThreeVector((FD2_startX1 + IdFib * 2.* spacingX) / cos(30. * deg), 0., startZ1);
              }
            else
              {
                posFibX  = G4ThreeVector( FD2_startX2 + IdFib * 2.* spacingX                  , 0., startZ2);
                posFibUV = G4ThreeVector((FD2_startX2 + IdFib * 2.* spacingX) / cos(30. * deg), 0., startZ2);
              }

            int fiberID = IdFib*4 + (1-IdUD)*2;

            G4ThreeVector posFib_x1;
            G4ThreeVector posFib_u1;
            G4ThreeVector posFib_v1;
            G4ThreeVector posFib_x2;
            G4ThreeVector posFib_u2;
            G4ThreeVector posFib_v2;

            double off_x_buf = fiber_uft2_off_x + fiber_offset[1][0][fiberID/2];
            double off_u_buf = fiber_uft2_off_u + fiber_offset[1][1][fiberID/2];
            double off_v_buf = fiber_uft2_off_v + fiber_offset[1][2][fiberID/2];
            off_u_buf /= cos(30. * deg);
            off_v_buf /= cos(30. * deg);
            double pos_x1_buf = posFibX.x()  + off_x_buf;
            double pos_x2_buf = posFibX.x()  + off_x_buf + spacingX;
            double pos_u1_buf = posFibUV.x() + off_u_buf;
            double pos_u2_buf = posFibUV.x() + off_u_buf + spacingX / cos(30. * deg);
            double pos_v1_buf = posFibUV.x() + off_v_buf;
            double pos_v2_buf = posFibUV.x() + off_v_buf + spacingX / cos(30. * deg);
            posFib_x1 = G4ThreeVector(pos_x1_buf, 0., posFibX.z());
            posFib_x2 = G4ThreeVector(pos_x2_buf, 0., posFibX.z());
            posFib_u1 = G4ThreeVector(pos_u1_buf, 0., posFibUV.z());
            posFib_u2 = G4ThreeVector(pos_u2_buf, 0., posFibUV.z());
            posFib_v1 = G4ThreeVector(pos_v1_buf, 0., posFibUV.z());
            posFib_v2 = G4ThreeVector(pos_v2_buf, 0., posFibUV.z());


            AllPlacements.emplace_back(new G4PVPlacement(rotFib1, Sign * posFib_x1, FiberD2_Cladding_log_x, "FiberD2_Cladding_x",
                                                         FiberD2_MothVol_log_x, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib1, Sign * posFib_x1, FiberD2_Core_log_x    , "FiberD2_Core_x",
                                                         FiberD2_MothVol_log_x, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib2, Sign * posFib_u1, FiberD2_Cladding_log_u, "FiberD2_Cladding_u",
                                                         FiberD2_MothVol_log_u, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib2, Sign * posFib_u1, FiberD2_Core_log_u    , "FiberD2_Core_u",
                                                         FiberD2_MothVol_log_u, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib3, Sign * posFib_v1, FiberD2_Cladding_log_v, "FiberD2_Cladding_v",
                                                         FiberD2_MothVol_log_v, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib3, Sign * posFib_v1, FiberD2_Core_log_v    , "FiberD2_Core_v",
                                                         FiberD2_MothVol_log_v, false, fiberID));

            //posFibX  += shiftFibX;
            //posFibUV += shiftFibUV;
            AllPlacements.emplace_back(new G4PVPlacement(rotFib1, Sign * posFib_x2, FiberD2_Cladding_log_x, "FiberD2_Cladding_x",
                                                         FiberD2_MothVol_log_x, false, fiberID + 1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib1, Sign * posFib_x2, FiberD2_Core_log_x    , "FiberD2_Core_x",
                                                         FiberD2_MothVol_log_x, false, fiberID + 1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib2, Sign * posFib_u2, FiberD2_Cladding_log_u, "FiberD2_Cladding_u",
                                                         FiberD2_MothVol_log_u, false, fiberID + 1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib2, Sign * posFib_u2, FiberD2_Core_log_u    , "FiberD2_Core_u",
                                                         FiberD2_MothVol_log_u, false, fiberID + 1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib3, Sign * posFib_v2, FiberD2_Cladding_log_v, "FiberD2_Cladding_v",
                                                         FiberD2_MothVol_log_v, false, fiberID + 1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib3, Sign * posFib_v2, FiberD2_Core_log_v    , "FiberD2_Core_v",
                                                         FiberD2_MothVol_log_v, false, fiberID + 1));
          }

      NameDetectorsSD.push_back(FiberD2_Core_log_x->GetName());
      NameDetectorsSD.push_back(FiberD2_Core_log_u->GetName());
      NameDetectorsSD.push_back(FiberD2_Core_log_v->GetName());

      if(!MiniVis)
        {
          FiberD2_Core_log_x->SetVisAttributes(visAttributes_x);
          FiberD2_Core_log_u->SetVisAttributes(visAttributes_u);
          FiberD2_Core_log_v->SetVisAttributes(visAttributes_v);
          FiberD2_Cladding_log_x->SetVisAttributes(visAttributes_x);
          FiberD2_Cladding_log_u->SetVisAttributes(visAttributes_u);
          FiberD2_Cladding_log_v->SetVisAttributes(visAttributes_v);
          FiberD2_MothVol_log->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD2_MothVol_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD2_MothVol_log_u->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD2_MothVol_log_v->SetVisAttributes(G4VisAttributes::GetInvisible());
        }
      else
        {
          FiberD2_Core_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD2_Core_log_u->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD2_Core_log_v->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD2_Cladding_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD2_Cladding_log_u->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD2_Cladding_log_v->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD2_MothVol_log->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD2_MothVol_log_x->SetVisAttributes(visAttributes_x);
          FiberD2_MothVol_log_u->SetVisAttributes(visAttributes_u);
          FiberD2_MothVol_log_v->SetVisAttributes(visAttributes_v);
        }

      // -------------------------- Third Fiber Detector --------------------------
      G4VSolid* FiberD3_MothVol_solid        = new G4Box("FiberDetector3", 30. * cm, 30. * cm, 6. * mm);
      G4VSolid* FiberD3_layerX_MothVol_solid = new G4Box("FiberD3_layerX_solid", 25. * cm, 20. * cm, 2. * mm);
      G4VSolid* FiberD3_layerU_MothVol_solid = new G4Box("FiberD3_layerU_solid", 25. * cm, 20. * cm, 2. * mm);
      G4VSolid* FiberD3_layerV_MothVol_solid = new G4Box("FiberD3_layerV_solid", 25. * cm, 20. * cm, 2. * mm);
      G4LogicalVolume* FiberD3_MothVol_log   = new G4LogicalVolume(FiberD3_MothVol_solid, Air, "FiberD3_log", 0, 0, 0);
      G4LogicalVolume* FiberD3_MothVol_log_x =
          new G4LogicalVolume(FiberD3_layerX_MothVol_solid, Air, "FiberD3_log_x", 0, 0, 0);
      G4LogicalVolume* FiberD3_MothVol_log_u =
          new G4LogicalVolume(FiberD3_layerU_MothVol_solid, Air, "FiberD3_log_u", 0, 0, 0);
      G4LogicalVolume* FiberD3_MothVol_log_v =
          new G4LogicalVolume(FiberD3_layerV_MothVol_solid, Air, "FiberD3_log_v", 0, 0, 0);

      AllPlacements.emplace_back(new G4PVPlacement(
          0, Sign * (G4ThreeVector(fiber_uft3_pos_x, fiber_uft3_pos_y*Sign, HypHI_FiberTracker3_posZ + Systematic_shift) - transMFLD_new),
          FiberD3_MothVol_log, "FiberDetector3", MFLD_log, false, 0));
      AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(0., 0., posZshift[0])),
                                                   FiberD3_MothVol_log_x, "FiberD3_x", FiberD3_MothVol_log, false, 0));
      AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(0., 0., posZshift[1])),
                                                   FiberD3_MothVol_log_u, "FiberD3_u", FiberD3_MothVol_log, false, 0));
      AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(0., 0., posZshift[2])),
                                                   FiberD3_MothVol_log_v, "FiberD3_v", FiberD3_MothVol_log, false, 0));

      const G4double FD3_startX1  = 105.4625 * mm;
      const G4double FD3_startX2  = 105.1875 * mm;
      const G4double FD3_startUV1 = 105.1875 * mm;
      const G4double FD3_startUV2 = 105.4625 * mm;

      G4VSolid* FiberD3_Core_solid_x =
          new G4Tubs("FiberD3_Core_solid_x", 0, 0.24 * mm, (1.5 * 10.5) * cm, 0. * deg, 360. * deg);
      G4VSolid* FiberD3_Core_solid_u =
          new G4Tubs("FiberD3_Core_solid_u", 0, 0.24 * mm, (1.5 * 10.5 / sqrt(3) * 2) * cm, 0. * deg, 360. * deg);
      G4VSolid* FiberD3_Core_solid_v =
          new G4Tubs("FiberD3_Core_solid_v", 0, 0.24 * mm, (1.5 * 10.5 / sqrt(3) * 2) * cm, 0. * deg, 360. * deg);
      G4LogicalVolume* FiberD3_Core_log_x =
          new G4LogicalVolume(FiberD3_Core_solid_x, FiberCoreScinti, "FiberD3_Core_log_x", 0, 0, 0);
      G4LogicalVolume* FiberD3_Core_log_u =
          new G4LogicalVolume(FiberD3_Core_solid_u, FiberCoreScinti, "FiberD3_Core_log_u", 0, 0, 0);
      G4LogicalVolume* FiberD3_Core_log_v =
          new G4LogicalVolume(FiberD3_Core_solid_v, FiberCoreScinti, "FiberD3_Core_log_v", 0, 0, 0);

      G4VSolid* FiberD3_Cladding_solid_x =
          new G4Tubs("FiberD3_Cladding_solid_x", 0.24 * mm, 0.25 * mm, (1.5 * 10.5) * cm, 0. * deg, 360. * deg);
      G4VSolid* FiberD3_Cladding_solid_u = new G4Tubs("FiberD3_Cladding_solid_u", 0.24 * mm, 0.25 * mm,
                                                      (1.5 * 10.5 / sqrt(3) * 2) * cm, 0. * deg, 360. * deg);
      G4VSolid* FiberD3_Cladding_solid_v = new G4Tubs("FiberD3_Cladding_solid_v", 0.24 * mm, 0.25 * mm,
                                                      (1.5 * 10.5 / sqrt(3) * 2) * cm, 0. * deg, 360. * deg);
      G4LogicalVolume* FiberD3_Cladding_log_x =
          new G4LogicalVolume(FiberD3_Cladding_solid_x, Scinti, "FiberD3_Cladding_log_x", 0, 0, 0);
      G4LogicalVolume* FiberD3_Cladding_log_u =
          new G4LogicalVolume(FiberD3_Cladding_solid_u, Scinti, "FiberD3_Cladding_log_u", 0, 0, 0);
      G4LogicalVolume* FiberD3_Cladding_log_v =
          new G4LogicalVolume(FiberD3_Cladding_solid_v, Scinti, "FiberD3_Cladding_log_v", 0, 0, 0);


      for(G4int IdUD = 0; IdUD < 2; ++IdUD)
        for(G4int IdFib = 0; IdFib < 192; ++IdFib)
          {
            if(IdUD == 0)
              {
                posFibX  = G4ThreeVector( FD3_startX1  - IdFib * 2.*spacingX                  , 0., startZ1);
                posFibUV = G4ThreeVector((FD3_startUV1 - IdFib * 2.*spacingX) / cos(30. * deg), 0., startZ1);
              }
            else
              {
                posFibX  = G4ThreeVector( FD3_startX2  - IdFib * 2.* spacingX                  , 0., startZ2);
                posFibUV = G4ThreeVector((FD3_startUV2 - IdFib * 2.* spacingX) / cos(30. * deg), 0., startZ2);
              }

            int fiberID_X  = IdFib*4 + IdUD*2;
            int fiberID_UV = IdFib*4 + (1-IdUD)*2;

            G4ThreeVector posFib_x1;
            G4ThreeVector posFib_u1;
            G4ThreeVector posFib_v1;
            G4ThreeVector posFib_x2;
            G4ThreeVector posFib_u2;
            G4ThreeVector posFib_v2;

            double off_x_buf = fiber_uft3_off_x + fiber_offset[2][0][fiberID_X /2];
            double off_u_buf = fiber_uft3_off_u + fiber_offset[2][1][fiberID_UV/2];
            double off_v_buf = fiber_uft3_off_v + fiber_offset[2][2][fiberID_UV/2];
            off_u_buf /= cos(30. * deg);
            off_v_buf /= cos(30. * deg);
            double pos_x1_buf = posFibX.x()  + off_x_buf;
            double pos_x2_buf = posFibX.x()  + off_x_buf - spacingX;
            double pos_u1_buf = posFibUV.x() + off_u_buf;
            double pos_u2_buf = posFibUV.x() + off_u_buf - spacingX / cos(30. * deg);
            double pos_v1_buf = posFibUV.x() + off_v_buf;
            double pos_v2_buf = posFibUV.x() + off_v_buf - spacingX / cos(30. * deg);
            posFib_x1 = G4ThreeVector(pos_x1_buf, 0., posFibX.z());
            posFib_x2 = G4ThreeVector(pos_x2_buf, 0., posFibX.z());
            posFib_u1 = G4ThreeVector(pos_u1_buf, 0., posFibUV.z());
            posFib_u2 = G4ThreeVector(pos_u2_buf, 0., posFibUV.z());
            posFib_v1 = G4ThreeVector(pos_v1_buf, 0., posFibUV.z());
            posFib_v2 = G4ThreeVector(pos_v2_buf, 0., posFibUV.z());

            AllPlacements.emplace_back(new G4PVPlacement(rotFib1, Sign * posFib_x1, FiberD3_Cladding_log_x, "FiberD3_Cladding_x",
                                                         FiberD3_MothVol_log_x, false, fiberID_X));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib1, Sign * posFib_x1, FiberD3_Core_log_x    , "FiberD3_Core_x",
                                                         FiberD3_MothVol_log_x, false, fiberID_X));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib2, Sign * posFib_u1, FiberD3_Cladding_log_u, "FiberD3_Cladding_u",
                                                         FiberD3_MothVol_log_u, false, fiberID_UV));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib2, Sign * posFib_u1, FiberD3_Core_log_u    , "FiberD3_Core_u",
                                                         FiberD3_MothVol_log_u, false, fiberID_UV));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib3, Sign * posFib_u1, FiberD3_Cladding_log_v, "FiberD3_Cladding_v",
                                                         FiberD3_MothVol_log_v, false, fiberID_UV));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib3, Sign * posFib_u1, FiberD3_Core_log_v    , "FiberD3_Core_v",
                                                         FiberD3_MothVol_log_v, false, fiberID_UV));
            //posFibX  -= shiftFibX;
            //posFibUV -= shiftFibUV;
            AllPlacements.emplace_back(new G4PVPlacement(rotFib1, Sign * posFib_x2, FiberD3_Cladding_log_x, "FiberD3_Cladding_x",
                                                         FiberD3_MothVol_log_x, false, fiberID_X + 1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib1, Sign * posFib_x2, FiberD3_Core_log_x    , "FiberD3_Core_x",
                                                         FiberD3_MothVol_log_x, false, fiberID_X + 1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib2, Sign * posFib_u2, FiberD3_Cladding_log_u, "FiberD3_Cladding_u",
                                                         FiberD3_MothVol_log_u, false, fiberID_UV + 1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib2, Sign * posFib_u2, FiberD3_Core_log_u    , "FiberD3_Core_u",
                                                         FiberD3_MothVol_log_u, false, fiberID_UV + 1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib3, Sign * posFib_u2, FiberD3_Cladding_log_v, "FiberD3_Cladding_v",
                                                         FiberD3_MothVol_log_v, false, fiberID_UV + 1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib3, Sign * posFib_u2, FiberD3_Core_log_v    , "FiberD3_Core_v",
                                                         FiberD3_MothVol_log_v, false, fiberID_UV + 1));
          }

      NameDetectorsSD.push_back(FiberD3_Core_log_x->GetName());
      NameDetectorsSD.push_back(FiberD3_Core_log_u->GetName());
      NameDetectorsSD.push_back(FiberD3_Core_log_v->GetName());

      if(!MiniVis)
        {
          FiberD3_Core_log_x->SetVisAttributes(visAttributes_x);
          FiberD3_Core_log_u->SetVisAttributes(visAttributes_u);
          FiberD3_Core_log_v->SetVisAttributes(visAttributes_v);
          FiberD3_Cladding_log_x->SetVisAttributes(visAttributes_x);
          FiberD3_Cladding_log_u->SetVisAttributes(visAttributes_u);
          FiberD3_Cladding_log_v->SetVisAttributes(visAttributes_v);
          FiberD3_MothVol_log->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD3_MothVol_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD3_MothVol_log_u->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD3_MothVol_log_v->SetVisAttributes(G4VisAttributes::GetInvisible());
        }
      else
        {
          FiberD3_Core_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD3_Core_log_u->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD3_Core_log_v->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD3_Cladding_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD3_Cladding_log_u->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD3_Cladding_log_v->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD3_MothVol_log->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD3_MothVol_log_x->SetVisAttributes(visAttributes_x);
          FiberD3_MothVol_log_u->SetVisAttributes(visAttributes_u);
          FiberD3_MothVol_log_v->SetVisAttributes(visAttributes_v);
        }

      // -------------------------- Fourth Fiber Detector --------------------------
      G4VSolid* FiberD4_MothVol_solid        = new G4Box("FiberDetector4", 15. * cm, 20. * cm, 6. * mm);
      G4VSolid* FiberD4_layerV_MothVol_solid = new G4Box("FiberD4_layerV_solid", 15. * cm, 20. * cm, 2. * mm);
      G4VSolid* FiberD4_layerU_MothVol_solid = new G4Box("FiberD4_layerU_solid", 15. * cm, 20. * cm, 2. * mm);
      G4VSolid* FiberD4_layerX_MothVol_solid = new G4Box("FiberD4_layerX_solid", 15. * cm, 20. * cm, 2. * mm);
      G4LogicalVolume* FiberD4_MothVol_log   = new G4LogicalVolume(FiberD4_MothVol_solid, Air, "FiberD4_log", 0, 0, 0);
      G4LogicalVolume* FiberD4_MothVol_log_v =
          new G4LogicalVolume(FiberD4_layerV_MothVol_solid, Air, "FiberD4_log_v", 0, 0, 0);
      G4LogicalVolume* FiberD4_MothVol_log_u =
          new G4LogicalVolume(FiberD4_layerU_MothVol_solid, Air, "FiberD4_log_u", 0, 0, 0);
      G4LogicalVolume* FiberD4_MothVol_log_x =
          new G4LogicalVolume(FiberD4_layerX_MothVol_solid, Air, "FiberD4_log_x", 0, 0, 0);

      AllPlacements.emplace_back(new G4PVPlacement(
          0, Sign * (G4ThreeVector(fiber_dft1_pos_x, fiber_dft1_pos_y*Sign, HypHI_FiberTracker4_posZ + Systematic_shift) - transMFLD_new),
          FiberD4_MothVol_log, "FiberDetector4", MFLD_log, false, 0));
      AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(0., 0., posZshift[0])),
                                                   FiberD4_MothVol_log_v, "FiberD4_v", FiberD4_MothVol_log, false, 0));
      AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(0., 0., posZshift[1])),
                                                   FiberD4_MothVol_log_u, "FiberD4_u", FiberD4_MothVol_log, false, 0));
      AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(0., 0., posZshift[2])),
                                                   FiberD4_MothVol_log_x, "FiberD4_x", FiberD4_MothVol_log, false, 0));

      const G4double FD4_startX1 = 70.2625 * mm;
      const G4double FD4_startX2 = 69.9875 * mm;

      G4VSolid* FiberD4_Core_solid_v =
          new G4Tubs("FiberD4_Core_solid_v", 0, 0.24 * mm, (10.5 / sqrt(3) * 2) * cm, 0. * deg, 360. * deg);
      G4VSolid* FiberD4_Core_solid_u =
          new G4Tubs("FiberD4_Core_solid_u", 0, 0.24 * mm, (10.5 / sqrt(3) * 2) * cm, 0. * deg, 360. * deg);
      G4VSolid* FiberD4_Core_solid_x =
          new G4Tubs("FiberD4_Core_solid_x", 0, 0.24 * mm, 10.5 * cm, 0. * deg, 360. * deg);
      G4LogicalVolume* FiberD4_Core_log_v =
          new G4LogicalVolume(FiberD4_Core_solid_v, FiberCoreScinti, "FiberD4_Core_log_v", 0, 0, 0);
      G4LogicalVolume* FiberD4_Core_log_u =
          new G4LogicalVolume(FiberD4_Core_solid_u, FiberCoreScinti, "FiberD4_Core_log_u", 0, 0, 0);
      G4LogicalVolume* FiberD4_Core_log_x =
          new G4LogicalVolume(FiberD4_Core_solid_x, FiberCoreScinti, "FiberD4_Core_log_x", 0, 0, 0);

      G4VSolid* FiberD4_Cladding_solid_v =
          new G4Tubs("FiberD4_Cladding_solid_v", 0.24 * mm, 0.25 * mm, (10.5 / sqrt(3) * 2) * cm, 0. * deg, 360. * deg);
      G4VSolid* FiberD4_Cladding_solid_u =
          new G4Tubs("FiberD4_Cladding_solid_u", 0.24 * mm, 0.25 * mm, (10.5 / sqrt(3) * 2) * cm, 0. * deg, 360. * deg);
      G4VSolid* FiberD4_Cladding_solid_x =
          new G4Tubs("FiberD4_Cladding_solid_x", 0.24 * mm, 0.25 * mm, 10.5 * cm, 0. * deg, 360. * deg);
      G4LogicalVolume* FiberD4_Cladding_log_v =
          new G4LogicalVolume(FiberD4_Cladding_solid_v, Scinti, "FiberD4_Cladding_log_v", 0, 0, 0);
      G4LogicalVolume* FiberD4_Cladding_log_u =
          new G4LogicalVolume(FiberD4_Cladding_solid_u, Scinti, "FiberD4_Cladding_log_u", 0, 0, 0);
      G4LogicalVolume* FiberD4_Cladding_log_x =
          new G4LogicalVolume(FiberD4_Cladding_solid_x, Scinti, "FiberD4_Cladding_log_x", 0, 0, 0);

      for(G4int IdUD = 0; IdUD < 2; ++IdUD)
        for(G4int IdFib = 0; IdFib < 128; ++IdFib)
          {
            if(IdUD == 0)
              {
                posFibX  = G4ThreeVector( FD4_startX1 - IdFib * 2.*spacingX                  , 0., startZ1);
                posFibUV = G4ThreeVector((FD4_startX1 - IdFib * 2.*spacingX) / cos(30. * deg), 0., startZ1);
              }
            else
              {
                posFibX  = G4ThreeVector( FD4_startX2 - IdFib * 2.* spacingX                  , 0., startZ2);
                posFibUV = G4ThreeVector((FD4_startX2 - IdFib * 2.* spacingX) / cos(30. * deg), 0., startZ2);
              }

            int fiberID = IdFib*4 + IdUD*2;

            G4ThreeVector posFib_x1;
            G4ThreeVector posFib_u1;
            G4ThreeVector posFib_v1;
            G4ThreeVector posFib_x2;
            G4ThreeVector posFib_u2;
            G4ThreeVector posFib_v2;

            double off_x_buf = fiber_dft1_off_x + fiber_offset[5][0][fiberID/2];
            double off_u_buf = fiber_dft1_off_u + fiber_offset[5][1][fiberID/2];
            double off_v_buf = fiber_dft1_off_v + fiber_offset[5][2][fiberID/2];
            off_u_buf /= cos(30. * deg);
            off_v_buf /= cos(30. * deg);
            double pos_x1_buf = posFibX.x()  + off_x_buf;
            double pos_x2_buf = posFibX.x()  + off_x_buf - spacingX;
            double pos_u1_buf = posFibUV.x() + off_u_buf;
            double pos_u2_buf = posFibUV.x() + off_u_buf - spacingX / cos(30. * deg);
            double pos_v1_buf = posFibUV.x() + off_v_buf;
            double pos_v2_buf = posFibUV.x() + off_v_buf - spacingX / cos(30. * deg);
            posFib_x1 = G4ThreeVector(pos_x1_buf, 0., posFibX.z());
            posFib_x2 = G4ThreeVector(pos_x2_buf, 0., posFibX.z());
            posFib_u1 = G4ThreeVector(pos_u1_buf, 0., posFibUV.z());
            posFib_u2 = G4ThreeVector(pos_u2_buf, 0., posFibUV.z());
            posFib_v1 = G4ThreeVector(pos_v1_buf, 0., posFibUV.z());
            posFib_v2 = G4ThreeVector(pos_v2_buf, 0., posFibUV.z());

            AllPlacements.emplace_back(new G4PVPlacement(rotFib1, Sign * posFib_x1, FiberD4_Cladding_log_x, "FiberD4_Cladding_x",
                                                         FiberD4_MothVol_log_x, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib1, Sign * posFib_x1, FiberD4_Core_log_x    , "FiberD4_Core_x",
                                                         FiberD4_MothVol_log_x, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib3, Sign * posFib_u1, FiberD4_Cladding_log_u, "FiberD4_Cladding_u",
                                                         FiberD4_MothVol_log_u, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib3, Sign * posFib_u1, FiberD4_Core_log_u    , "FiberD4_Core_u",
                                                         FiberD4_MothVol_log_u, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib2, Sign * posFib_v1, FiberD4_Cladding_log_v, "FiberD4_Cladding_v",
                                                         FiberD4_MothVol_log_v, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib2, Sign * posFib_v1, FiberD4_Core_log_v    , "FiberD4_Core_v",
                                                         FiberD4_MothVol_log_v, false, fiberID));
            //posFibX  -= shiftFibX;
            //posFibUV -= shiftFibUV;
            AllPlacements.emplace_back(new G4PVPlacement(rotFib1, Sign * posFib_x2, FiberD4_Cladding_log_x, "FiberD4_Cladding_x",
                                                         FiberD4_MothVol_log_x, false, fiberID + 1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib1, Sign * posFib_x2, FiberD4_Core_log_x    , "FiberD4_Core_x",
                                                         FiberD4_MothVol_log_x, false, fiberID + 1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib3, Sign * posFib_u2, FiberD4_Cladding_log_u, "FiberD4_Cladding_u",
                                                         FiberD4_MothVol_log_u, false, fiberID + 1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib3, Sign * posFib_u2, FiberD4_Core_log_u    , "FiberD4_Core_u",
                                                         FiberD4_MothVol_log_u, false, fiberID + 1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib2, Sign * posFib_v2, FiberD4_Cladding_log_v, "FiberD4_Cladding_v",
                                                         FiberD4_MothVol_log_v, false, fiberID + 1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib2, Sign * posFib_v2, FiberD4_Core_log_v    , "FiberD4_Core_v",
                                                         FiberD4_MothVol_log_v, false, fiberID + 1));
          }

      NameDetectorsSD.push_back(FiberD4_Core_log_v->GetName());
      NameDetectorsSD.push_back(FiberD4_Core_log_u->GetName());
      NameDetectorsSD.push_back(FiberD4_Core_log_x->GetName());

      if(!MiniVis)
        {
          FiberD4_Core_log_v->SetVisAttributes(visAttributes_v);
          FiberD4_Core_log_u->SetVisAttributes(visAttributes_u);
          FiberD4_Core_log_x->SetVisAttributes(visAttributes_x);
          FiberD4_Cladding_log_v->SetVisAttributes(visAttributes_v);
          FiberD4_Cladding_log_u->SetVisAttributes(visAttributes_u);
          FiberD4_Cladding_log_x->SetVisAttributes(visAttributes_x);
          FiberD4_MothVol_log->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD4_MothVol_log_v->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD4_MothVol_log_u->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD4_MothVol_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
        }
      else
        {
          FiberD4_Core_log_v->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD4_Core_log_u->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD4_Core_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD4_Cladding_log_v->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD4_Cladding_log_u->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD4_Cladding_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD4_MothVol_log->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD4_MothVol_log_v->SetVisAttributes(visAttributes_v);
          FiberD4_MothVol_log_u->SetVisAttributes(visAttributes_u);
          FiberD4_MothVol_log_x->SetVisAttributes(visAttributes_x);
        }

      // -------------------------- Fifth Fiber Detector --------------------------
      G4VSolid* FiberD5_MothVol_solid        = new G4Box("FiberDetector5", 15. * cm, 20. * cm, 6. * mm);
      G4VSolid* FiberD5_layerX_MothVol_solid = new G4Box("FiberD5_layerX_solid", 15. * cm, 20. * cm, 2. * mm);
      G4VSolid* FiberD5_layerU_MothVol_solid = new G4Box("FiberD5_layerU_solid", 15. * cm, 20. * cm, 2. * mm);
      G4VSolid* FiberD5_layerV_MothVol_solid = new G4Box("FiberD5_layerV_solid", 15. * cm, 20. * cm, 2. * mm);
      G4LogicalVolume* FiberD5_MothVol_log   = new G4LogicalVolume(FiberD5_MothVol_solid, Air, "FiberD5_log", 0, 0, 0);
      G4LogicalVolume* FiberD5_MothVol_log_x =
          new G4LogicalVolume(FiberD5_layerX_MothVol_solid, Air, "FiberD5_log_x", 0, 0, 0);
      G4LogicalVolume* FiberD5_MothVol_log_u =
          new G4LogicalVolume(FiberD5_layerU_MothVol_solid, Air, "FiberD5_log_u", 0, 0, 0);
      G4LogicalVolume* FiberD5_MothVol_log_v =
          new G4LogicalVolume(FiberD5_layerV_MothVol_solid, Air, "FiberD5_log_v", 0, 0, 0);

      AllPlacements.emplace_back(new G4PVPlacement(
          0, Sign * (G4ThreeVector(fiber_dft2_pos_x, fiber_dft2_pos_y*Sign, HypHI_FiberTracker5_posZ + Systematic_shift) - transMFLD_new),
          FiberD5_MothVol_log, "FiberDetector5", MFLD_log, false, 0));
      AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(0., 0., posZshift[0])),
                                                   FiberD5_MothVol_log_x, "FiberD5_x", FiberD5_MothVol_log, false, 0));
      AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(0., 0., posZshift[1])),
                                                   FiberD5_MothVol_log_u, "FiberD5_u", FiberD5_MothVol_log, false, 0));
      AllPlacements.emplace_back(new G4PVPlacement(0, Sign * (G4ThreeVector(0., 0., posZshift[2])),
                                                   FiberD5_MothVol_log_v, "FiberD5_v", FiberD5_MothVol_log, false, 0));

      const G4double FD5_startX1 = 69.9875 * mm;
      const G4double FD5_startX2 = 70.2625 * mm;

      G4VSolid* FiberD5_Core_solid_x =
          new G4Tubs("FiberD5_Core_solid_x", 0, 0.24 * mm, 10.5 * cm, 0. * deg, 360. * deg);
      G4VSolid* FiberD5_Core_solid_u =
          new G4Tubs("FiberD5_Core_solid_u", 0, 0.24 * mm, (10.5 / sqrt(3) * 2) * cm, 0. * deg, 360. * deg);
      G4VSolid* FiberD5_Core_solid_v =
          new G4Tubs("FiberD5_Core_solid_v", 0, 0.24 * mm, (10.5 / sqrt(3) * 2) * cm, 0. * deg, 360. * deg);
      G4LogicalVolume* FiberD5_Core_log_x =
          new G4LogicalVolume(FiberD5_Core_solid_x, FiberCoreScinti, "FiberD5_Core_log_x", 0, 0, 0);
      G4LogicalVolume* FiberD5_Core_log_u =
          new G4LogicalVolume(FiberD5_Core_solid_u, FiberCoreScinti, "FiberD5_Core_log_u", 0, 0, 0);
      G4LogicalVolume* FiberD5_Core_log_v =
          new G4LogicalVolume(FiberD5_Core_solid_v, FiberCoreScinti, "FiberD5_Core_log_v", 0, 0, 0);

      G4VSolid* FiberD5_Cladding_solid_x =
          new G4Tubs("FiberD5_Cladding_solid_x", 0.24 * mm, 0.25 * mm, 10.5 * cm, 0. * deg, 360. * deg);
      G4VSolid* FiberD5_Cladding_solid_u =
          new G4Tubs("FiberD5_Cladding_solid_u", 0.24 * mm, 0.25 * mm, (10.5 / sqrt(3) * 2) * cm, 0. * deg, 360. * deg);
      G4VSolid* FiberD5_Cladding_solid_v =
          new G4Tubs("FiberD5_Cladding_solid_v", 0.24 * mm, 0.25 * mm, (10.5 / sqrt(3) * 2) * cm, 0. * deg, 360. * deg);
      G4LogicalVolume* FiberD5_Cladding_log_x =
          new G4LogicalVolume(FiberD5_Cladding_solid_x, Scinti, "FiberD5_Cladding_log_x", 0, 0, 0);
      G4LogicalVolume* FiberD5_Cladding_log_u =
          new G4LogicalVolume(FiberD5_Cladding_solid_u, Scinti, "FiberD5_Cladding_log_u", 0, 0, 0);
      G4LogicalVolume* FiberD5_Cladding_log_v =
          new G4LogicalVolume(FiberD5_Cladding_solid_v, Scinti, "FiberD5_Cladding_log_v", 0, 0, 0);

      for(G4int IdUD = 0; IdUD < 2; ++IdUD)
        for(G4int IdFib = 0; IdFib < 128; ++IdFib)
          {
            if(IdUD == 0)
              {
                posFibX  = G4ThreeVector( FD5_startX1 - IdFib * 2.*spacingX                  , 0., startZ1);
                posFibUV = G4ThreeVector((FD5_startX1 - IdFib * 2.*spacingX) / cos(30. * deg), 0., startZ1);
              }
            else
              {
                posFibX  = G4ThreeVector( FD5_startX2 - IdFib * 2.* spacingX                  , 0., startZ2);
                posFibUV = G4ThreeVector((FD5_startX2 - IdFib * 2.* spacingX) / cos(30. * deg), 0., startZ2);
              }

            int fiberID = IdFib*4 + (1-IdUD)*2;

            G4ThreeVector posFib_x1;
            G4ThreeVector posFib_u1;
            G4ThreeVector posFib_v1;
            G4ThreeVector posFib_x2;
            G4ThreeVector posFib_u2;
            G4ThreeVector posFib_v2;

            double off_x_buf = fiber_dft2_off_x + fiber_offset[6][0][fiberID/2];
            double off_u_buf = fiber_dft2_off_u + fiber_offset[6][1][fiberID/2];
            double off_v_buf = fiber_dft2_off_v + fiber_offset[6][2][fiberID/2];
            off_u_buf /= cos(30. * deg);
            off_v_buf /= cos(30. * deg);
            double pos_x1_buf = posFibX.x()  + off_x_buf;
            double pos_x2_buf = posFibX.x()  + off_x_buf - spacingX;
            double pos_u1_buf = posFibUV.x() + off_u_buf;
            double pos_u2_buf = posFibUV.x() + off_u_buf - spacingX / cos(30. * deg);
            double pos_v1_buf = posFibUV.x() + off_v_buf;
            double pos_v2_buf = posFibUV.x() + off_v_buf - spacingX / cos(30. * deg);
            posFib_x1 = G4ThreeVector(pos_x1_buf, 0., posFibX.z());
            posFib_x2 = G4ThreeVector(pos_x2_buf, 0., posFibX.z());
            posFib_u1 = G4ThreeVector(pos_u1_buf, 0., posFibUV.z());
            posFib_u2 = G4ThreeVector(pos_u2_buf, 0., posFibUV.z());
            posFib_v1 = G4ThreeVector(pos_v1_buf, 0., posFibUV.z());
            posFib_v2 = G4ThreeVector(pos_v2_buf, 0., posFibUV.z());

            AllPlacements.emplace_back(new G4PVPlacement(rotFib1, Sign * posFib_x1, FiberD5_Cladding_log_x, "FiberD5_Cladding_x",
                                                         FiberD5_MothVol_log_x, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib1, Sign * posFib_x1, FiberD5_Core_log_x    , "FiberD5_Core_x",
                                                         FiberD5_MothVol_log_x, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib2, Sign * posFib_u1, FiberD5_Cladding_log_u, "FiberD5_Cladding_u",
                                                         FiberD5_MothVol_log_u, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib2, Sign * posFib_u1, FiberD5_Core_log_u    , "FiberD5_Core_u",
                                                         FiberD5_MothVol_log_u, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib3, Sign * posFib_v1, FiberD5_Cladding_log_v, "FiberD5_Cladding_v",
                                                         FiberD5_MothVol_log_v, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib3, Sign * posFib_v1, FiberD5_Core_log_v    , "FiberD5_Core_v",
                                                         FiberD5_MothVol_log_v, false, fiberID));

            //posFibX  -= shiftFibX;
            //posFibUV -= shiftFibUV;
            AllPlacements.emplace_back(new G4PVPlacement(rotFib1, Sign * posFib_x2, FiberD5_Cladding_log_x, "FiberD5_Cladding_x",
                                                         FiberD5_MothVol_log_x, false, fiberID + 1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib1, Sign * posFib_x2, FiberD5_Core_log_x    , "FiberD5_Core_x",
                                                         FiberD5_MothVol_log_x, false, fiberID + 1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib2, Sign * posFib_u2, FiberD5_Cladding_log_u, "FiberD5_Cladding_u",
                                                         FiberD5_MothVol_log_u, false, fiberID + 1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib2, Sign * posFib_u2, FiberD5_Core_log_u    , "FiberD5_Core_u",
                                                         FiberD5_MothVol_log_u, false, fiberID + 1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib3, Sign * posFib_v2, FiberD5_Cladding_log_v, "FiberD5_Cladding_v",
                                                         FiberD5_MothVol_log_v, false, fiberID + 1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFib3, Sign * posFib_v2, FiberD5_Core_log_v    , "FiberD5_Core_v",
                                                         FiberD5_MothVol_log_v, false, fiberID + 1));
          }

      NameDetectorsSD.push_back(FiberD5_Core_log_x->GetName());
      NameDetectorsSD.push_back(FiberD5_Core_log_u->GetName());
      NameDetectorsSD.push_back(FiberD5_Core_log_v->GetName());

      if(!MiniVis)
        {
          FiberD5_Core_log_x->SetVisAttributes(visAttributes_x);
          FiberD5_Core_log_u->SetVisAttributes(visAttributes_u);
          FiberD5_Core_log_v->SetVisAttributes(visAttributes_v);
          FiberD5_Cladding_log_x->SetVisAttributes(visAttributes_x);
          FiberD5_Cladding_log_u->SetVisAttributes(visAttributes_u);
          FiberD5_Cladding_log_v->SetVisAttributes(visAttributes_v);
          FiberD5_MothVol_log->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD5_MothVol_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD5_MothVol_log_u->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD5_MothVol_log_v->SetVisAttributes(G4VisAttributes::GetInvisible());
        }
      else
        {
          FiberD5_Core_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD5_Core_log_u->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD5_Core_log_v->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD5_Cladding_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD5_Cladding_log_u->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD5_Cladding_log_v->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD5_MothVol_log->SetVisAttributes(G4VisAttributes::GetInvisible());
          FiberD5_MothVol_log_x->SetVisAttributes(visAttributes_x);
          FiberD5_MothVol_log_u->SetVisAttributes(visAttributes_u);
          FiberD5_MothVol_log_v->SetVisAttributes(visAttributes_v);
        }
    }
  // ----------------------------- @@ -----------------------------
  // ----------------------------- @@ -----------------------------

  // ----------------------------- @@ -----------------------------
  // Mini Fiber Detector
  // ----------------------------- @@ -----------------------------

  if(Par.IsAvailable("HypHI_MiniFiberTR_On"))
    {

      const double HypHI_MiniFiberTracker1_posX = fiber_mft1_pos_x;
      const double HypHI_MiniFiberTracker1_posY = fiber_mft1_pos_y;
      const double HypHI_MiniFiberTracker1_posZ = Par.IsAvailable("HypHI_MiniFiberTracker1_posZ") ? Par.Get<double>("HypHI_MiniFiberTracker1_posZ") : 2269.3*mm;

      const double HypHI_MiniFiberTracker2_posX = fiber_mft2_pos_x;
      const double HypHI_MiniFiberTracker2_posY = fiber_mft2_pos_y;
      const double HypHI_MiniFiberTracker2_posZ = Par.IsAvailable("HypHI_MiniFiberTracker2_posZ") ? Par.Get<double>("HypHI_MiniFiberTracker2_posZ") : 2309.3*mm;

      G4RotationMatrix* rotFibM1 = new G4RotationMatrix;
      rotFibM1->rotateZ(0. * deg);
      G4RotationMatrix* rotFibM2 = new G4RotationMatrix;
      rotFibM2->rotateZ(60. * deg);
      G4RotationMatrix* rotFibM3 = new G4RotationMatrix;
      rotFibM3->rotateZ(-60. * deg);

      G4RotationMatrix* rotFibM1x1 = new G4RotationMatrix;
      G4RotationMatrix* rotFibM1u1 = new G4RotationMatrix;
      G4RotationMatrix* rotFibM1v1 = new G4RotationMatrix;
      G4RotationMatrix* rotFibM1x2 = new G4RotationMatrix;
      G4RotationMatrix* rotFibM1u2 = new G4RotationMatrix;
      G4RotationMatrix* rotFibM1v2 = new G4RotationMatrix;

      G4RotationMatrix* rotFibM2x1 = new G4RotationMatrix;
      G4RotationMatrix* rotFibM2u1 = new G4RotationMatrix;
      G4RotationMatrix* rotFibM2v1 = new G4RotationMatrix;
      G4RotationMatrix* rotFibM2x2 = new G4RotationMatrix;
      G4RotationMatrix* rotFibM2u2 = new G4RotationMatrix;
      G4RotationMatrix* rotFibM2v2 = new G4RotationMatrix;

      rotFibM1x1->rotateZ( fiber_angle_offset[3][0][0] * deg);
      rotFibM1u1->rotateZ( fiber_angle_offset[3][1][0] * deg);
      rotFibM1v1->rotateZ( fiber_angle_offset[3][2][0] * deg);
      rotFibM1x2->rotateZ( fiber_angle_offset[3][0][1] * deg);
      rotFibM1u2->rotateZ( fiber_angle_offset[3][1][1] * deg);
      rotFibM1v2->rotateZ( fiber_angle_offset[3][2][1] * deg);

      rotFibM2x1->rotateZ( fiber_angle_offset[4][0][0] * deg);
      rotFibM2u1->rotateZ( fiber_angle_offset[4][2][0] * deg);
      rotFibM2v1->rotateZ( fiber_angle_offset[4][1][0] * deg);
      rotFibM2x2->rotateZ( fiber_angle_offset[4][0][1] * deg);
      rotFibM2u2->rotateZ( fiber_angle_offset[4][2][1] * deg);
      rotFibM2v2->rotateZ( fiber_angle_offset[4][1][1] * deg);

      rotFibM1x1->rotateX( 90. * deg);
      rotFibM1u1->rotateX( 90. * deg);
      rotFibM1v1->rotateX( 90. * deg);
      rotFibM1x2->rotateX( 90. * deg);
      rotFibM1u2->rotateX( 90. * deg);
      rotFibM1v2->rotateX( 90. * deg);

      rotFibM2x1->rotateX( 90. * deg);
      rotFibM2u1->rotateX( 90. * deg);
      rotFibM2v1->rotateX( 90. * deg);
      rotFibM2x2->rotateX( 90. * deg);
      rotFibM2u2->rotateX( 90. * deg);
      rotFibM2v2->rotateX( 90. * deg);

      G4RotationMatrix* rotFib1 = new G4RotationMatrix;
      rotFib1->rotateX(90. * deg);

      const G4double spacingX               = 0.55 * mm;
      const G4double startZ1                = 0. * mm;
      const G4double startZ2                = 0.47631397 * mm;
      const std::vector<G4double> posZshift = {-4. * mm, 0. * mm, 4. * mm};

      G4ThreeVector posFib;
      G4ThreeVector shiftFib  = G4ThreeVector(spacingX, 0., 0.);

      G4VisAttributes* visAttributes_x = new G4VisAttributes(color_x);
      G4VisAttributes* visAttributes_u = new G4VisAttributes(color_u);
      G4VisAttributes* visAttributes_v = new G4VisAttributes(color_v);

      // -------------------------- First MiniFiber Detector --------------------------
      G4VSolid* MiniFiberD1_MothVol_solid         = new G4Box("MiniFiberDetector1", 30. * cm, 30. * cm, 1.5 * cm);
      G4VSolid* MiniFiberD1_layerX_MothVol_solid = new G4Box("MiniFiberD1_layerX_solid", 20. * cm, 20. * cm, 2. * mm);
      G4VSolid* MiniFiberD1_layerU_MothVol_solid = new G4Box("MiniFiberD1_layerU_solid", 20. * cm, 20. * cm, 2. * mm);
      G4VSolid* MiniFiberD1_layerV_MothVol_solid = new G4Box("MiniFiberD1_layerV_solid", 20. * cm, 20. * cm, 2. * mm);
      G4LogicalVolume* MiniFiberD1_MothVol_log =
          new G4LogicalVolume(MiniFiberD1_MothVol_solid, Air, "MiniFiberD1_log", 0, 0, 0);
      G4LogicalVolume* MiniFiberD1_MothVol_log_x =
          new G4LogicalVolume(MiniFiberD1_layerX_MothVol_solid, Air, "MiniFiberD1_log_x", 0, 0, 0);
      G4LogicalVolume* MiniFiberD1_MothVol_log_u =
          new G4LogicalVolume(MiniFiberD1_layerU_MothVol_solid, Air, "MiniFiberD1_log_u", 0, 0, 0);
      G4LogicalVolume* MiniFiberD1_MothVol_log_v =
          new G4LogicalVolume(MiniFiberD1_layerV_MothVol_solid, Air, "MiniFiberD1_log_v", 0, 0, 0);

      AllPlacements.emplace_back(new G4PVPlacement(
            0, Sign * (G4ThreeVector( HypHI_MiniFiberTracker1_posX, HypHI_MiniFiberTracker1_posY*Sign, HypHI_MiniFiberTracker1_posZ + Systematic_shift) - transMFLD_new),
          MiniFiberD1_MothVol_log, "MiniFiberDetector1", MFLD_log, false, 0));

      AllPlacements.emplace_back(new G4PVPlacement(rotFibM1, Sign * (G4ThreeVector(0., 0., posZshift[0])),
                                                   MiniFiberD1_MothVol_log_x, "MiniFiberD1_x",
                                                   MiniFiberD1_MothVol_log, false, 0));
      AllPlacements.emplace_back(new G4PVPlacement(rotFibM3, Sign * (G4ThreeVector(0., 0., posZshift[1])),
                                                   MiniFiberD1_MothVol_log_u, "MiniFiberD1_u",
                                                   MiniFiberD1_MothVol_log, false, 0));
      AllPlacements.emplace_back(new G4PVPlacement(rotFibM2, Sign * (G4ThreeVector(0., 0., posZshift[2])),
                                                   MiniFiberD1_MothVol_log_v, "MiniFiberD1_v",
                                                   MiniFiberD1_MothVol_log, false, 0));

      const std::vector<G4double> MFD1_startX1 = { -82.4   * mm, 12.55 * mm};
      const std::vector<G4double> MFD1_startX2 = { -82.675 * mm, 12.825 * mm};

      G4VSolid* MiniFiberD1_Core_solid_x =
          new G4Tubs("MiniFiberD1_Core_solid_x", 0, 0.24 * mm, 10.0 * cm, 0. * deg, 360. * deg);
      G4VSolid* MiniFiberD1_Core_solid_u =
          new G4Tubs("MiniFiberD1_Core_solid_u", 0, 0.24 * mm, 10.0 * cm, 0. * deg, 360. * deg);
      G4VSolid* MiniFiberD1_Core_solid_v =
          new G4Tubs("MiniFiberD1_Core_solid_v", 0, 0.24 * mm, 10.0 * cm, 0. * deg, 360. * deg);
      G4LogicalVolume* MiniFiberD1_Core_log_x =
          new G4LogicalVolume(MiniFiberD1_Core_solid_x, FiberCoreScinti, "MiniFiberD1_Core_log_x", 0, 0, 0);
      G4LogicalVolume* MiniFiberD1_Core_log_u =
          new G4LogicalVolume(MiniFiberD1_Core_solid_u, FiberCoreScinti, "MiniFiberD1_Core_log_u", 0, 0, 0);
      G4LogicalVolume* MiniFiberD1_Core_log_v =
          new G4LogicalVolume(MiniFiberD1_Core_solid_v, FiberCoreScinti, "MiniFiberD1_Core_log_v", 0, 0, 0);

      G4VSolid* MiniFiberD1_Cladding_solid_x =
          new G4Tubs("MiniFiberD1_Cladding_solid_x", 0.24 * mm, 0.25 * mm, 10.0 * cm, 0. * deg, 360. * deg);
      G4VSolid* MiniFiberD1_Cladding_solid_u =
          new G4Tubs("MiniFiberD1_Cladding_solid_u", 0.24 * mm, 0.25 * mm, 10.0 * cm, 0. * deg, 360. * deg);
      G4VSolid* MiniFiberD1_Cladding_solid_v =
          new G4Tubs("MiniFiberD1_Cladding_solid_v", 0.24 * mm, 0.25 * mm, 10.0 * cm, 0. * deg, 360. * deg);
      G4LogicalVolume* MiniFiberD1_Cladding_log_x =
          new G4LogicalVolume(MiniFiberD1_Cladding_solid_x, Scinti, "MiniFiberD1_Cladding_log_x", 0, 0, 0);
      G4LogicalVolume* MiniFiberD1_Cladding_log_u =
          new G4LogicalVolume(MiniFiberD1_Cladding_solid_u, Scinti, "MiniFiberD1_Cladding_log_u", 0, 0, 0);
      G4LogicalVolume* MiniFiberD1_Cladding_log_v =
          new G4LogicalVolume(MiniFiberD1_Cladding_solid_v, Scinti, "MiniFiberD1_Cladding_log_v", 0, 0, 0);


      for(G4int IdUD = 0; IdUD < 2; ++IdUD)
        for(G4int IdFib = 0; IdFib < 128; ++IdFib)
          {
            int leftright=IdFib>63;
            if(IdUD == 0)
              {
                posFib  = G4ThreeVector(MFD1_startX1[leftright] + (IdFib - leftright*64) * 2.*spacingX, 0., startZ1);
              }
            else
              {
                posFib = G4ThreeVector(MFD1_startX2[leftright] + (IdFib - leftright*64) * 2.* spacingX, 0., startZ2);
              }

            int fiberID = IdFib*4 + (1-IdUD)*2;
            if(leftright) fiberID = IdFib*4 + IdUD*2;

            G4ThreeVector posFib_x1;
            G4ThreeVector posFib_u1;
            G4ThreeVector posFib_v1;
            G4ThreeVector posFib_x2;
            G4ThreeVector posFib_u2;
            G4ThreeVector posFib_v2;

            G4RotationMatrix* rotFibMx;
            G4RotationMatrix* rotFibMu;
            G4RotationMatrix* rotFibMv;
            if(leftright==0){
              double off_x_buf = fiber_mft1_off_x1 + fiber_offset[3][0][fiberID/2];
              double off_u_buf = fiber_mft1_off_u1 + fiber_offset[3][1][fiberID/2];
              double off_v_buf = fiber_mft1_off_v1 + fiber_offset[3][2][fiberID/2];
              double pos_x1_buf = posFib.x() + off_x_buf;
              double pos_x2_buf = posFib.x() + off_x_buf + spacingX;
              double pos_u1_buf = posFib.x() + off_u_buf;
              double pos_u2_buf = posFib.x() + off_u_buf + spacingX;
              double pos_v1_buf = posFib.x() + off_v_buf;
              double pos_v2_buf = posFib.x() + off_v_buf + spacingX;
              pos_x1_buf /= cos(fiber_angle_offset[3][0][0]*deg);
              pos_x2_buf /= cos(fiber_angle_offset[3][0][0]*deg);
              pos_u1_buf /= cos(fiber_angle_offset[3][1][0]*deg);
              pos_u2_buf /= cos(fiber_angle_offset[3][1][0]*deg);
              pos_v1_buf /= cos(fiber_angle_offset[3][2][0]*deg);
              pos_v2_buf /= cos(fiber_angle_offset[3][2][0]*deg);
              posFib_x1 = G4ThreeVector(pos_x1_buf, 0., posFib.z());
              posFib_x2 = G4ThreeVector(pos_x2_buf, 0., posFib.z());
              posFib_u1 = G4ThreeVector(pos_u1_buf, 0., posFib.z());
              posFib_u2 = G4ThreeVector(pos_u2_buf, 0., posFib.z());
              posFib_v1 = G4ThreeVector(pos_v1_buf, 0., posFib.z());
              posFib_v2 = G4ThreeVector(pos_v2_buf, 0., posFib.z());
              rotFibMx = rotFibM1x1;
              rotFibMu = rotFibM1u1;
              rotFibMv = rotFibM1v1;
            }
            else{
              double off_x_buf = fiber_mft1_off_x2 + fiber_offset[3][0][fiberID/2];
              double off_u_buf = fiber_mft1_off_u2 + fiber_offset[3][1][fiberID/2];
              double off_v_buf = fiber_mft1_off_v2 + fiber_offset[3][2][fiberID/2];
              double pos_x1_buf = posFib.x() + off_x_buf;
              double pos_x2_buf = posFib.x() + off_x_buf + spacingX;
              double pos_u1_buf = posFib.x() + off_u_buf;
              double pos_u2_buf = posFib.x() + off_u_buf + spacingX;
              double pos_v1_buf = posFib.x() + off_v_buf;
              double pos_v2_buf = posFib.x() + off_v_buf + spacingX;
              pos_x1_buf /= cos(fiber_angle_offset[3][0][1]*deg);
              pos_x2_buf /= cos(fiber_angle_offset[3][0][1]*deg);
              pos_u1_buf /= cos(fiber_angle_offset[3][1][1]*deg);
              pos_u2_buf /= cos(fiber_angle_offset[3][1][1]*deg);
              pos_v1_buf /= cos(fiber_angle_offset[3][2][1]*deg);
              pos_v2_buf /= cos(fiber_angle_offset[3][2][1]*deg);
              posFib_x1 = G4ThreeVector(pos_x1_buf, 0., posFib.z());
              posFib_x2 = G4ThreeVector(pos_x2_buf, 0., posFib.z());
              posFib_u1 = G4ThreeVector(pos_u1_buf, 0., posFib.z());
              posFib_u2 = G4ThreeVector(pos_u2_buf, 0., posFib.z());
              posFib_v1 = G4ThreeVector(pos_v1_buf, 0., posFib.z());
              posFib_v2 = G4ThreeVector(pos_v2_buf, 0., posFib.z());
              rotFibMx = rotFibM1x2;
              rotFibMu = rotFibM1u2;
              rotFibMv = rotFibM1v2;
            }

            AllPlacements.emplace_back(new G4PVPlacement(rotFibMx, Sign * posFib_x1, MiniFiberD1_Cladding_log_x,
                                                         "MiniFiberD1_Cladding_x", MiniFiberD1_MothVol_log_x, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFibMx, Sign * posFib_x1, MiniFiberD1_Core_log_x, "MiniFiberD1_Core_x",
                                                         MiniFiberD1_MothVol_log_x, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFibMu, Sign * posFib_u1, MiniFiberD1_Cladding_log_u,
                                                         "MiniFiberD1_Cladding_u", MiniFiberD1_MothVol_log_u, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFibMu, Sign * posFib_u1, MiniFiberD1_Core_log_u, "MiniFiberD1_Core_u",
                                                         MiniFiberD1_MothVol_log_u, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFibMv, Sign * posFib_v1, MiniFiberD1_Cladding_log_v,
                                                         "MiniFiberD1_Cladding_v", MiniFiberD1_MothVol_log_v, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFibMv, Sign * posFib_v1, MiniFiberD1_Core_log_v, "MiniFiberD1_Core_v",
                                                         MiniFiberD1_MothVol_log_v, false, fiberID));

            //posFib += shiftFib;
            AllPlacements.emplace_back(new G4PVPlacement(rotFibMx, Sign * posFib_x2, MiniFiberD1_Cladding_log_x,
                                                         "MiniFiberD1_Cladding_x", MiniFiberD1_MothVol_log_x, false, fiberID+1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFibMx, Sign * posFib_x2, MiniFiberD1_Core_log_x, "MiniFiberD1_Core_x",
                                                         MiniFiberD1_MothVol_log_x, false, fiberID+1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFibMu, Sign * posFib_u2, MiniFiberD1_Cladding_log_u,
                                                         "MiniFiberD1_Cladding_u", MiniFiberD1_MothVol_log_u, false, fiberID+1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFibMu, Sign * posFib_u2, MiniFiberD1_Core_log_u, "MiniFiberD1_Core_u",
                                                         MiniFiberD1_MothVol_log_u, false, fiberID+1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFibMv, Sign * posFib_v2, MiniFiberD1_Cladding_log_v,
                                                         "MiniFiberD1_Cladding_v", MiniFiberD1_MothVol_log_v, false, fiberID+1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFibMv, Sign * posFib_v2, MiniFiberD1_Core_log_v, "MiniFiberD1_Core_v",
                                                         MiniFiberD1_MothVol_log_v, false, fiberID+1));
          }

      NameDetectorsSD.push_back(MiniFiberD1_Core_log_x->GetName());
      NameDetectorsSD.push_back(MiniFiberD1_Core_log_u->GetName());
      NameDetectorsSD.push_back(MiniFiberD1_Core_log_v->GetName());


      if(!MiniVis)
        {
          MiniFiberD1_Core_log_x->SetVisAttributes(visAttributes_x);
          MiniFiberD1_Core_log_u->SetVisAttributes(visAttributes_u);
          MiniFiberD1_Core_log_v->SetVisAttributes(visAttributes_v);

          MiniFiberD1_Cladding_log_x->SetVisAttributes(visAttributes_x);
          MiniFiberD1_Cladding_log_u->SetVisAttributes(visAttributes_u);
          MiniFiberD1_Cladding_log_v->SetVisAttributes(visAttributes_v);

          MiniFiberD1_MothVol_log->SetVisAttributes(G4VisAttributes::GetInvisible());
          MiniFiberD1_MothVol_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
          MiniFiberD1_MothVol_log_u->SetVisAttributes(G4VisAttributes::GetInvisible());
          MiniFiberD1_MothVol_log_v->SetVisAttributes(G4VisAttributes::GetInvisible());

        }
      else
        {
          MiniFiberD1_Core_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
          MiniFiberD1_Core_log_u->SetVisAttributes(G4VisAttributes::GetInvisible());
          MiniFiberD1_Core_log_v->SetVisAttributes(G4VisAttributes::GetInvisible());

          MiniFiberD1_Cladding_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
          MiniFiberD1_Cladding_log_u->SetVisAttributes(G4VisAttributes::GetInvisible());
          MiniFiberD1_Cladding_log_v->SetVisAttributes(G4VisAttributes::GetInvisible());


          MiniFiberD1_MothVol_log->SetVisAttributes(G4VisAttributes::GetInvisible());
          MiniFiberD1_MothVol_log_x->SetVisAttributes(visAttributes_x);
          MiniFiberD1_MothVol_log_u->SetVisAttributes(visAttributes_u);
          MiniFiberD1_MothVol_log_v->SetVisAttributes(visAttributes_v);

        }


      // -------------------------- Second MiniFiber Detector --------------------------
      G4VSolid* MiniFiberD2_MothVol_solid         = new G4Box("MiniFiberDetector2", 30. * cm, 30. * cm, 1.5 * cm);
      G4VSolid* MiniFiberD2_layerX_MothVol_solid = new G4Box("MiniFiberD2_layerX_solid", 20. * cm, 20. * cm, 2. * mm);
      G4VSolid* MiniFiberD2_layerV_MothVol_solid = new G4Box("MiniFiberD2_layerV_solid", 20. * cm, 20. * cm, 2. * mm);
      G4VSolid* MiniFiberD2_layerU_MothVol_solid = new G4Box("MiniFiberD2_layerU_solid", 20. * cm, 20. * cm, 2. * mm);
      G4LogicalVolume* MiniFiberD2_MothVol_log =
          new G4LogicalVolume(MiniFiberD2_MothVol_solid, Air, "MiniFiberD2_log", 0, 0, 0);
      G4LogicalVolume* MiniFiberD2_MothVol_log_x =
          new G4LogicalVolume(MiniFiberD2_layerX_MothVol_solid, Air, "MiniFiberD2_log_x", 0, 0, 0);
      G4LogicalVolume* MiniFiberD2_MothVol_log_v =
          new G4LogicalVolume(MiniFiberD2_layerV_MothVol_solid, Air, "MiniFiberD2_log_v", 0, 0, 0);
      G4LogicalVolume* MiniFiberD2_MothVol_log_u =
          new G4LogicalVolume(MiniFiberD2_layerU_MothVol_solid, Air, "MiniFiberD2_log_u", 0, 0, 0);

      AllPlacements.emplace_back(new G4PVPlacement(
            0, Sign * (G4ThreeVector(HypHI_MiniFiberTracker2_posX, HypHI_MiniFiberTracker2_posY*Sign, HypHI_MiniFiberTracker2_posZ + Systematic_shift) - transMFLD_new),
               MiniFiberD2_MothVol_log, "MiniFiberDetector2", MFLD_log, false, 0));

      AllPlacements.emplace_back(new G4PVPlacement(rotFibM1, Sign * (G4ThreeVector(0., 0., posZshift[0])),
                                                   MiniFiberD2_MothVol_log_x, "MiniFiberD2_x",
                                                   MiniFiberD2_MothVol_log, false, 0));
      AllPlacements.emplace_back(new G4PVPlacement(rotFibM2, Sign * (G4ThreeVector(0., 0., posZshift[1])),
                                                   MiniFiberD2_MothVol_log_v, "MiniFiberD2_v",
                                                   MiniFiberD2_MothVol_log, false, 0));
      AllPlacements.emplace_back(new G4PVPlacement(rotFibM3, Sign * (G4ThreeVector(0., 0., posZshift[2])),
                                                   MiniFiberD2_MothVol_log_u, "MiniFiberD2_u",
                                                   MiniFiberD2_MothVol_log, false, 0));

      const std::vector<G4double> MFD2_startX1 = { -82.675 * mm, 12.825 * mm};
      const std::vector<G4double> MFD2_startX2 = { -82.4   * mm, 12.55 * mm};

      G4VSolid* MiniFiberD2_Core_solid_x =
          new G4Tubs("MiniFiberD2_Core_solid_x", 0, 0.24 * mm, 10.0 * cm, 0. * deg, 360. * deg);
      G4VSolid* MiniFiberD2_Core_solid_v =
          new G4Tubs("MiniFiberD2_Core_solid_v", 0, 0.24 * mm, 10.0 * cm, 0. * deg, 360. * deg);
      G4VSolid* MiniFiberD2_Core_solid_u =
          new G4Tubs("MiniFiberD2_Core_solid_u", 0, 0.24 * mm, 10.0 * cm, 0. * deg, 360. * deg);
      G4LogicalVolume* MiniFiberD2_Core_log_x =
          new G4LogicalVolume(MiniFiberD2_Core_solid_x, FiberCoreScinti, "MiniFiberD2_Core_log_x", 0, 0, 0);
      G4LogicalVolume* MiniFiberD2_Core_log_v =
          new G4LogicalVolume(MiniFiberD2_Core_solid_v, FiberCoreScinti, "MiniFiberD2_Core_log_v", 0, 0, 0);
      G4LogicalVolume* MiniFiberD2_Core_log_u =
          new G4LogicalVolume(MiniFiberD2_Core_solid_u, FiberCoreScinti, "MiniFiberD2_Core_log_u", 0, 0, 0);

      G4VSolid* MiniFiberD2_Cladding_solid_x =
          new G4Tubs("MiniFiberD2_Cladding_solid_x", 0.24 * mm, 0.25 * mm, 10.0 * cm, 0. * deg, 360. * deg);
      G4VSolid* MiniFiberD2_Cladding_solid_v =
          new G4Tubs("MiniFiberD2_Cladding_solid_v", 0.24 * mm, 0.25 * mm, 10.0 * cm, 0. * deg, 360. * deg);
      G4VSolid* MiniFiberD2_Cladding_solid_u =
          new G4Tubs("MiniFiberD2_Cladding_solid_u", 0.24 * mm, 0.25 * mm, 10.0 * cm, 0. * deg, 360. * deg);
      G4LogicalVolume* MiniFiberD2_Cladding_log_x =
          new G4LogicalVolume(MiniFiberD2_Cladding_solid_x, Scinti, "MiniFiberD2_Cladding_log_x", 0, 0, 0);
      G4LogicalVolume* MiniFiberD2_Cladding_log_v =
          new G4LogicalVolume(MiniFiberD2_Cladding_solid_v, Scinti, "MiniFiberD2_Cladding_log_v", 0, 0, 0);
      G4LogicalVolume* MiniFiberD2_Cladding_log_u =
          new G4LogicalVolume(MiniFiberD2_Cladding_solid_u, Scinti, "MiniFiberD2_Cladding_log_u", 0, 0, 0);


      for(G4int IdUD = 0; IdUD < 2; ++IdUD)
        for(G4int IdFib = 0; IdFib < 128; ++IdFib)
          {
            int leftright=IdFib>63;
            if(IdUD == 0)
              {
                posFib  = G4ThreeVector(MFD2_startX1[leftright] + (IdFib - leftright*64) * 2.*spacingX, 0., startZ1);
              }
            else
              {
                posFib = G4ThreeVector(MFD2_startX2[leftright] + (IdFib - leftright*64) * 2.* spacingX, 0., startZ2);
              }

            int fiberID = IdFib*4 + IdUD*2;
            if(leftright) fiberID = IdFib*4 + (1-IdUD)*2;

            G4ThreeVector posFib_x1;
            G4ThreeVector posFib_u1;
            G4ThreeVector posFib_v1;
            G4ThreeVector posFib_x2;
            G4ThreeVector posFib_u2;
            G4ThreeVector posFib_v2;

            G4RotationMatrix* rotFibMx;
            G4RotationMatrix* rotFibMu;
            G4RotationMatrix* rotFibMv;
            if(leftright==0){
              double off_x_buf = fiber_mft2_off_x1 + fiber_offset[4][0][fiberID/2];
              double off_u_buf = fiber_mft2_off_u1 + fiber_offset[4][2][fiberID/2];
              double off_v_buf = fiber_mft2_off_v1 + fiber_offset[4][1][fiberID/2];
              double pos_x1_buf = posFib.x() + off_x_buf;
              double pos_x2_buf = posFib.x() + off_x_buf + spacingX;
              double pos_u1_buf = posFib.x() + off_u_buf;
              double pos_u2_buf = posFib.x() + off_u_buf + spacingX;
              double pos_v1_buf = posFib.x() + off_v_buf;
              double pos_v2_buf = posFib.x() + off_v_buf + spacingX;
              pos_x1_buf /= cos(fiber_angle_offset[4][0][0]*deg);
              pos_x2_buf /= cos(fiber_angle_offset[4][0][0]*deg);
              pos_u1_buf /= cos(fiber_angle_offset[4][2][0]*deg);
              pos_u2_buf /= cos(fiber_angle_offset[4][2][0]*deg);
              pos_v1_buf /= cos(fiber_angle_offset[4][1][0]*deg);
              pos_v2_buf /= cos(fiber_angle_offset[4][1][0]*deg);
              posFib_x1 = G4ThreeVector(pos_x1_buf, 0., posFib.z());
              posFib_x2 = G4ThreeVector(pos_x2_buf, 0., posFib.z());
              posFib_u1 = G4ThreeVector(pos_u1_buf, 0., posFib.z());
              posFib_u2 = G4ThreeVector(pos_u2_buf, 0., posFib.z());
              posFib_v1 = G4ThreeVector(pos_v1_buf, 0., posFib.z());
              posFib_v2 = G4ThreeVector(pos_v2_buf, 0., posFib.z());
              rotFibMx = rotFibM2x1;
              rotFibMu = rotFibM2u1;
              rotFibMv = rotFibM2v1;
            }
            else{
              double off_x_buf = fiber_mft2_off_x2 + fiber_offset[4][0][fiberID/2];
              double off_u_buf = fiber_mft2_off_u2 + fiber_offset[4][2][fiberID/2];
              double off_v_buf = fiber_mft2_off_v2 + fiber_offset[4][1][fiberID/2];
              double pos_x1_buf = posFib.x() + off_x_buf;
              double pos_x2_buf = posFib.x() + off_x_buf + spacingX;
              double pos_u1_buf = posFib.x() + off_u_buf;
              double pos_u2_buf = posFib.x() + off_u_buf + spacingX;
              double pos_v1_buf = posFib.x() + off_v_buf;
              double pos_v2_buf = posFib.x() + off_v_buf + spacingX;
              pos_x1_buf /= cos(fiber_angle_offset[4][0][1]*deg);
              pos_x2_buf /= cos(fiber_angle_offset[4][0][1]*deg);
              pos_u1_buf /= cos(fiber_angle_offset[4][2][1]*deg);
              pos_u2_buf /= cos(fiber_angle_offset[4][2][1]*deg);
              pos_v1_buf /= cos(fiber_angle_offset[4][1][1]*deg);
              pos_v2_buf /= cos(fiber_angle_offset[4][1][1]*deg);
              posFib_x1 = G4ThreeVector(pos_x1_buf, 0., posFib.z());
              posFib_x2 = G4ThreeVector(pos_x2_buf, 0., posFib.z());
              posFib_u1 = G4ThreeVector(pos_u1_buf, 0., posFib.z());
              posFib_u2 = G4ThreeVector(pos_u2_buf, 0., posFib.z());
              posFib_v1 = G4ThreeVector(pos_v1_buf, 0., posFib.z());
              posFib_v2 = G4ThreeVector(pos_v2_buf, 0., posFib.z());
              rotFibMx = rotFibM2x2;
              rotFibMu = rotFibM2u2;
              rotFibMv = rotFibM2v2;
            }

            AllPlacements.emplace_back(new G4PVPlacement(rotFibMx, Sign * posFib_x1, MiniFiberD2_Cladding_log_x,
                                                         "MiniFiberD2_Cladding_x", MiniFiberD2_MothVol_log_x, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFibMx, Sign * posFib_x1, MiniFiberD2_Core_log_x, "MiniFiberD2_Core_x",
                                                         MiniFiberD2_MothVol_log_x, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFibMu, Sign * posFib_u1, MiniFiberD2_Cladding_log_u,
                                                         "MiniFiberD2_Cladding_u", MiniFiberD2_MothVol_log_u, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFibMu, Sign * posFib_u1, MiniFiberD2_Core_log_u, "MiniFiberD2_Core_u",
                                                         MiniFiberD2_MothVol_log_u, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFibMv, Sign * posFib_v1, MiniFiberD2_Cladding_log_v,
                                                         "MiniFiberD2_Cladding_v", MiniFiberD2_MothVol_log_v, false, fiberID));
            AllPlacements.emplace_back(new G4PVPlacement(rotFibMv, Sign * posFib_v1, MiniFiberD2_Core_log_v, "MiniFiberD2_Core_v",
                                                         MiniFiberD2_MothVol_log_v, false, fiberID));

            //posFib  += shiftFib;
            AllPlacements.emplace_back(new G4PVPlacement(rotFibMx, Sign * posFib_x2, MiniFiberD2_Cladding_log_x,
                                                         "MiniFiberD2_Cladding_x", MiniFiberD2_MothVol_log_x, false, fiberID+1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFibMx, Sign * posFib_x2, MiniFiberD2_Core_log_x, "MiniFiberD2_Core_x",
                                                         MiniFiberD2_MothVol_log_x, false, fiberID+1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFibMu, Sign * posFib_u2, MiniFiberD2_Cladding_log_u,
                                                         "MiniFiberD2_Cladding_u", MiniFiberD2_MothVol_log_u, false, fiberID+1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFibMu, Sign * posFib_u2, MiniFiberD2_Core_log_u, "MiniFiberD2_Core_u",
                                                         MiniFiberD2_MothVol_log_u, false, fiberID+1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFibMv, Sign * posFib_v2, MiniFiberD2_Cladding_log_v,
                                                         "MiniFiberD2_Cladding_v", MiniFiberD2_MothVol_log_v, false, fiberID+1));
            AllPlacements.emplace_back(new G4PVPlacement(rotFibMv, Sign * posFib_v2, MiniFiberD2_Core_log_v, "MiniFiberD2_Core_v",
                                                         MiniFiberD2_MothVol_log_v, false, fiberID+1));
          }

      NameDetectorsSD.push_back(MiniFiberD2_Core_log_x->GetName());
      NameDetectorsSD.push_back(MiniFiberD2_Core_log_u->GetName());
      NameDetectorsSD.push_back(MiniFiberD2_Core_log_v->GetName());


      if(!MiniVis)
        {
          MiniFiberD2_Core_log_x->SetVisAttributes(visAttributes_x);
          MiniFiberD2_Core_log_u->SetVisAttributes(visAttributes_u);
          MiniFiberD2_Core_log_v->SetVisAttributes(visAttributes_v);

          MiniFiberD2_Cladding_log_x->SetVisAttributes(visAttributes_x);
          MiniFiberD2_Cladding_log_u->SetVisAttributes(visAttributes_u);
          MiniFiberD2_Cladding_log_v->SetVisAttributes(visAttributes_v);

          MiniFiberD2_MothVol_log->SetVisAttributes(G4VisAttributes::GetInvisible());
          MiniFiberD2_MothVol_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
          MiniFiberD2_MothVol_log_u->SetVisAttributes(G4VisAttributes::GetInvisible());
          MiniFiberD2_MothVol_log_v->SetVisAttributes(G4VisAttributes::GetInvisible());

        }
      else
        {
          MiniFiberD2_Core_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
          MiniFiberD2_Core_log_u->SetVisAttributes(G4VisAttributes::GetInvisible());
          MiniFiberD2_Core_log_v->SetVisAttributes(G4VisAttributes::GetInvisible());

          MiniFiberD2_Cladding_log_x->SetVisAttributes(G4VisAttributes::GetInvisible());
          MiniFiberD2_Cladding_log_u->SetVisAttributes(G4VisAttributes::GetInvisible());
          MiniFiberD2_Cladding_log_v->SetVisAttributes(G4VisAttributes::GetInvisible());


          MiniFiberD2_MothVol_log->SetVisAttributes(G4VisAttributes::GetInvisible());
          MiniFiberD2_MothVol_log_x->SetVisAttributes(visAttributes_x);
          MiniFiberD2_MothVol_log_u->SetVisAttributes(visAttributes_u);
          MiniFiberD2_MothVol_log_v->SetVisAttributes(visAttributes_v);

        }
    }

  // ----------------------------- @@ -----------------------------
  // ----------------------------- @@ -----------------------------
  if(Par.IsAvailable("FMF2_On"))
    {
      G4VSolid* EndFMF2_box        = new G4Box("EndFMF2_box", 25. * cm, 25. * cm, 1. * mm);
      G4LogicalVolume* EndFMF2_log = new G4LogicalVolume(EndFMF2_box, Scinti, "FMF2_log", 0, 0, 0);

      NameDetectorsSD.push_back(EndFMF2_log->GetName());

      G4RotationMatrix* rotFMF2 = new G4RotationMatrix;
      // if(WasaSide==1)
      //   rotFMF2->rotateY(-180*degree);
      double FMF2_posZ = Par.Get<double>("FRS_FMF2_posZ");

      AllPlacements.emplace_back(
				 new G4PVPlacement(rotFMF2, Sign * (G4ThreeVector(0., 0., FMF2_posZ + Systematic_shift) - transMFLD_new),
						   EndFMF2_log, "FMF2_phys", MFLD_log, false, 0));
      AllPlacements.emplace_back(
				 new G4PVPlacement(rotFMF2, Sign * (G4ThreeVector(0., 0., FMF2_posZ + Systematic_shift + 1. * cm) - transMFLD_new),
						   EndFMF2_log, "FMF2_phys1", MFLD_log, false, 1));
      AllPlacements.emplace_back(
				 new G4PVPlacement(rotFMF2, Sign * (G4ThreeVector(0., 0., FMF2_posZ + Systematic_shift + 2. * cm) - transMFLD_new),
						   EndFMF2_log, "FMF2_phys2", MFLD_log, false, 2));

      G4VisAttributes* FMF2_att = new G4VisAttributes(Red);
      FMF2_att->SetForceWireframe(false);
      EndFMF2_log->SetVisAttributes(FMF2_att);
    }

  if(Par.IsAvailable("WASA2_FrontCap"))
    {
      double WASA2_FrontCap_rmin = Par.Get<double>("WASA2_FrontCap_minR");
      double WASA2_FrontCap_rmax = Par.Get<double>("WASA2_FrontCap_maxR");
      double WASA2_FrontCap_PosZ = Par.Get<double>("WASA2_FrontCap_posZ");
      double WASA2_FrontCap_Thickness = 0.3*mm;

      G4VSolid* WASA2_FrontCap = new G4Tubs("WASA2_FrontCap", WASA2_FrontCap_rmin, WASA2_FrontCap_rmax, WASA2_FrontCap_Thickness, 0, CLHEP::twopi);
      G4LogicalVolume* WASA2_FrontCap_log = new G4LogicalVolume(WASA2_FrontCap, Air, "WASA2_FrontCap_log", 0, 0, 0); // CDCFieldMgr,0,0);

      G4RotationMatrix* rotEndCap = new G4RotationMatrix;
      AllPlacements.emplace_back(
                                 new G4PVPlacement(rotEndCap, Sign * (G4ThreeVector(0., 0., WASA2_FrontCap_PosZ + Systematic_shift) - transMFLD_new),
                                                   WASA2_FrontCap_log, "WASA2_FrontCap", MFLD_log, false, 0));

      NameDetectorsSD.push_back(WASA2_FrontCap_log->GetName());
    }
  
  if(Par.IsAvailable("HypHI_EndCap"))
    {
      double HypHI_EndCap_rmax = Par.Get<double>("HypHI_EndCap_maxR");
      double HypHI_EndCap_PosZ = Par.Get<double>("HypHI_EndCap_posZ");

      G4VSolid* HypHI_Endcap = new G4Tubs("HypHI_Endcap", 0, HypHI_EndCap_rmax, 2 * cm, 0, CLHEP::twopi);
      G4LogicalVolume* HypHI_Endcap_log =
	new G4LogicalVolume(HypHI_Endcap, Air, "HypHI_Endcap_log", 0, 0, 0); // CDCFieldMgr,0,0);

      G4RotationMatrix* rotEndCap = new G4RotationMatrix;
      // if(WasaSide==1)
      //   rotEndCap->rotateY(90*degree);

      AllPlacements.emplace_back(
				 new G4PVPlacement(rotEndCap, Sign * (G4ThreeVector(0., 0., HypHI_EndCap_PosZ + Systematic_shift) - transMFLD_new),
						   HypHI_Endcap_log, "HypHI_Endcap", MFLD_log, false, 0));

      G4VSolid* HypHI_TrackerFwd =
	new G4Tubs("HypHI_TrackerFwd", BeamHoleSize * 0.5, HypHI_EndCap_rmax, 2 * mm, 0, CLHEP::twopi);
      G4LogicalVolume* HypHI_TrackerFwd_log = new G4LogicalVolume(HypHI_TrackerFwd, Air, "HypHI_TrackFwd_log", 0, 0, 0);

      NameDetectorsSD.push_back(HypHI_TrackerFwd_log->GetName());

      AllPlacements.emplace_back(new G4PVPlacement(0, G4ThreeVector(0, 0, -1.5 * cm), HypHI_TrackerFwd_log,
						   "HypHI_TrackerFwd0", HypHI_Endcap_log, false, 0));

      AllPlacements.emplace_back(new G4PVPlacement(0, G4ThreeVector(0, 0, -1. * cm), HypHI_TrackerFwd_log,
						   "HypHI_TrackerFwd1", HypHI_Endcap_log, false, 1));

      AllPlacements.emplace_back(new G4PVPlacement(0, G4ThreeVector(0, 0, -0.5 * cm), HypHI_TrackerFwd_log,
						   "HypHI_TrackerFwd2", HypHI_Endcap_log, false, 2));

      G4VSolid* HypHI_RPC_l = new G4Tubs("HypHI_RPC_segment_L", BeamHoleSize * 0.5, HypHI_EndCap_rmax / 2., 0.5 * cm,
					 -CLHEP::twopi / 16., 2. * CLHEP::twopi / 16.);
      G4LogicalVolume* HypHI_RPC_l_log = new G4LogicalVolume(HypHI_RPC_l, Air, "HypHI_RPC_l_log", 0, 0, 0);
      NameDetectorsSD.push_back(HypHI_RPC_l_log->GetName());
      G4VSolid* HypHI_RPC_h = new G4Tubs("HypHI_RPC_segment_H", HypHI_EndCap_rmax / 2., HypHI_EndCap_rmax, 0.5 * cm,
					 -CLHEP::twopi / 16., 2. * CLHEP::twopi / 16.);
      G4LogicalVolume* HypHI_RPC_h_log = new G4LogicalVolume(HypHI_RPC_h, Air, "HypHI_RPC_h_log", 0, 0, 0);
      NameDetectorsSD.push_back(HypHI_RPC_h_log->GetName());

      for(int idRPC = 0; idRPC < 8; ++idRPC)
	{
	  G4RotationMatrix* rotRPC = new G4RotationMatrix;
	  double rotAngle          = CLHEP::twopi / 8. * static_cast<double>(idRPC);
	  rotRPC->rotateZ(rotAngle);
	  std::string nameRPC("HypHI_RPC_l");
	  nameRPC += std::to_string(idRPC);
	  AllPlacements.emplace_back(new G4PVPlacement(rotRPC, G4ThreeVector(0, 0, 0.5 * cm), HypHI_RPC_l_log, nameRPC,
						       HypHI_Endcap_log, false, idRPC));

	  std::string nameRPC2("HypHI_RPC_h");
	  nameRPC2 += std::to_string(idRPC);
	  AllPlacements.emplace_back(new G4PVPlacement(rotRPC, G4ThreeVector(0, 0, 0.5 * cm), HypHI_RPC_h_log, nameRPC2,
						       HypHI_Endcap_log, false, idRPC));
	}

      //--- Visualization ---//
      if(!MiniVis)
	{
	  HypHI_Endcap_log->SetVisAttributes(G4VisAttributes::GetInvisible());
	  G4VisAttributes* HypHI_RPC_att = new G4VisAttributes(Orange);
	  HypHI_RPC_att->SetForceWireframe(false);
	  HypHI_RPC_l_log->SetVisAttributes(HypHI_RPC_att);
	  HypHI_RPC_h_log->SetVisAttributes(HypHI_RPC_att);
	  G4VisAttributes* HypHI_Tracker_att = new G4VisAttributes(LightPurple);
	  HypHI_Tracker_att->SetForceWireframe(false);
	  HypHI_TrackerFwd_log->SetVisAttributes(HypHI_Tracker_att);
	}
      else
	{
	  G4VisAttributes* HypHI_RPC_att = new G4VisAttributes(Orange);
	  HypHI_Endcap_log->SetVisAttributes(HypHI_RPC_att);
	  HypHI_RPC_l_log->SetVisAttributes(G4VisAttributes::GetInvisible());
	  HypHI_RPC_h_log->SetVisAttributes(G4VisAttributes::GetInvisible());
	  HypHI_TrackerFwd_log->SetVisAttributes(G4VisAttributes::GetInvisible());
	}
    }
  // Region

  G4Region* aDetectorRegion =
      G4RegionStore::GetInstance()->FindOrCreateRegion("DetectorRegion"); // new G4Region("DetectorRegion");

  for(auto& CurrentName : NameDetectorsSD)
    {
      G4LogicalVolume* Det = FindVolume(CurrentName);
      Det->SetRegion(aDetectorRegion);
      aDetectorRegion->AddRootLogicalVolume(Det);
    }
  std::vector<double> cutsDet(4, Par.Get<double>("DetectorRegionCut"));
  aDetectorRegion->SetProductionCuts(new G4ProductionCuts());

  aDetectorRegion->GetProductionCuts()->SetProductionCuts(cutsDet);

  G4Region* aTargetRegion = G4RegionStore::GetInstance()->FindOrCreateRegion("TargetRegion");
  HypHI_Target_log->SetRegion(aTargetRegion);
  aTargetRegion->AddRootLogicalVolume(HypHI_Target_log);
  std::vector<double> cutsTarget(4, Par.Get<double>("TargetRegionCut"));
  aTargetRegion->SetProductionCuts(new G4ProductionCuts());
  aTargetRegion->GetProductionCuts()->SetProductionCuts(cutsTarget);

  experimentalHall_logOutRoot  = world->GetLogicalVolume();
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

  experimentalHall_log->SetUserLimits(new G4UserLimits(DBL_MAX, 2 * m, 10 * s, 0., 0.));

  bool isFieldMap = Par.IsAvailable("Field_WASAMap");

  if(isFieldMap)
    {
      const std::string field_name = Par.Get<std::string>("Field_WASAMap");
      double max_valueField =
          Par.IsAvailable("Field_WASAMapMaxField") ? Par.Get<double>("Field_WASAMapMaxField") : 1. * tesla;

      std::cout << "Field origin MFLD Phys:" << MFLD_phys << " " << MFLD_phys->GetInstanceID() << "\n";
      transMFLD = MFLD_phys->GetObjectTranslation();
      rotMFLD   = MFLD_phys->GetObjectRotationValue();

      std::cout << "Trans:" << transMFLD << "\n";
      std::cout << "Rot:" << rotMFLD << "\n";

      double signDir = rotMFLD.xx()< 0 ? 1. : -1. ;
      double originF[3] = {transMFLD.x(), transMFLD.y(), transMFLD.z()};

      fMagneticField = new G4SolWASAMapMagneticField(field_name);
      dynamic_cast<G4SolWASAMapMagneticField*>(fMagneticField)->SetOriginField(originF,signDir);
      dynamic_cast<G4SolWASAMapMagneticField*>(fMagneticField)->SetMaxField(max_valueField);
      dynamic_cast<G4SolWASAMapMagneticField*>(fMagneticField)->InitField();
    }
  else
    {
      G4double fCDField = 0. * tesla;
      if(Par.IsAvailable("Field_CDS_Bz"))
        fCDField = Par.Get<double>("Field_CDS_Bz");

      G4ThreeVector fCDC(0.0, 0.0, fCDField);

      fMagneticField = new G4SolSimpleMagneticField();
      dynamic_cast<G4SolSimpleMagneticField*>(fMagneticField)->SetField(fCDC);
    }

  fEquation = new G4Mag_UsualEqRhs(fMagneticField);
  // fStepper = new G4ClassicalRK4( fEquation );
  fStepper = new G4NystromRK4(fEquation);

  fFieldMgr = new G4FieldManager();
  // fFieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  fFieldMgr->SetDetectorField(fMagneticField);
  // fFieldMgr->CreateChordFinder(fMagneticField);
  fChordFinder = new G4ChordFinder(fMagneticField, 1.e-2, fStepper);
  fFieldMgr->SetChordFinder(fChordFinder);

  G4bool forceToAllDaughters = true;
  if(isFieldMap)
    MFLD_log->SetFieldManager(fFieldMgr, forceToAllDaughters);
  else
    INNER_log->SetFieldManager(fFieldMgr, forceToAllDaughters);
  // CDS_endcap_log->SetFieldManager(fFieldMgr,forceToAllDaughters);
  // if(DoModHypHI)
  //   HypHI_InTracker_log->SetFieldManager(fFieldMgr,forceToAllDaughters);

  G4AutoDelete::Register(fMagneticField);
  // G4AutoDelete::Register(fFieldMgr);
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
