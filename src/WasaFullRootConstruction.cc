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
// $Id: WasaFullRootConstruction.cc 77601 2013-11-26 17:08:44Z gcosmo $
//
/// \file WasaFullRootConstruction.cc
/// \brief Implementation of the WasaFullRootConstruction class

#include "WasaFullRootConstruction.hh"

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

#include "TFile.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal G4MagneticField* WasaFullRootConstruction::fMagneticField = 0;
G4ThreadLocal G4FieldManager* WasaFullRootConstruction::fFieldMgr       = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

WasaFullRootConstruction::WasaFullRootConstruction(G4SolConfig& _par)
  : G4SolVDetectorConstruction(_par), experimentalHall_log(nullptr), experimentalHall_phys(nullptr)
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// ---------------------------------------------------------------------------
// VGM demo
//

WasaFullRootConstruction::~WasaFullRootConstruction() {}

G4VPhysicalVolume* WasaFullRootConstruction::Construct()
{


  const G4Colour Blue        = {0.0, 0.0, 1.0};
  const G4Colour Gray        = {0.5, 0.5, 0.5};
  const G4Colour Red         = {1.0, 0.0, 0.0};
  const G4Colour LightBlue   = {0.0, 0.0, 1.0, 0.7};
  const G4Colour Yellow      = {1.0, 1.0, 0.0};
  const G4Colour Purple      = {.5, 0.0, .5};
  const G4Colour LightPurple = {106. / 255., 27. / 255., 189. / 255.};
  const G4Colour Green       = {0.0, 1.0, 0.0};
  const G4Colour Orange      = {1.0, 0.647, 0};
  const G4Colour Pink        = {1.0, 0.753, 0.796};

  const G4Colour  color_x     = {1., 1., 0.8};
  const G4Colour  color_u     = {1., 0.8, 1.};
  const G4Colour  color_v     = {0.8, 1., 1.};

  // G4Colour  white   (1.0, 1.0, 1.0) ;
  // G4Colour  grey    (0.5, 0.5, 0.5) ;
  // G4Colour  lgrey   (.85, .85, .85) ;
  // G4Colour  red     (1.0, 0.0, 0.0) ;
  // G4Colour  blue    (0.0, 0.0, 1.0) ;
  // G4Colour  cyan    (0.0, 1.0, 1.0) ;
  // G4Colour  magenta (1.0, 0.0, 1.0) ;
  // G4Colour  yellow  (1.0, 1.0, 0.0) ;
  // G4Colour  orange  (.75, .55, 0.0) ;
  // G4Colour  lblue   (0.0, 0.0, .75) ;
  // G4Colour  lgreen  (0.0, .75, 0.0) ;
  // G4Colour  green   (0.0, 1.0, 0.0) ;
  // G4Colour  brown   (0.7, 0.4, 0.1) ;

  const G4Colour PaletteOrange[5] = {{254 / 255., 237 / 255., 222 / 255.},
                                     {253 / 255., 190 / 255., 133 / 255.},
                                     {253 / 255., 141 / 255., 060 / 255.},
                                     {230 / 255., 85 / 255., 13 / 255.},
                                     {166 / 255., 54 / 255., 3 / 255.}};
  const G4Colour PaletteGreen[5]  = {{237 / 255., 248 / 255., 233 / 255.},
				     {186 / 255., 228 / 255., 179 / 255.},
				     {116 / 255., 196 / 255., 118 / 255.},
				     {49 / 255., 163 / 255., 84 / 255.},
				     {0, 109 / 255., 44 / 255.}};
  const G4Colour PaletteRed[5]    = {{254 / 255., 229 / 255., 217 / 255.},
				     {252 / 255., 174 / 255., 145 / 255.},
				     {251 / 255., 106 / 255., 74 / 255.},
				     {222 / 255., 45 / 255., 38 / 255.},
				     {165 / 255., 15 / 255., 21 / 255.}};

  const G4Colour ColorCDC[17] = {
    PaletteOrange[0], PaletteOrange[0], PaletteOrange[1], PaletteOrange[1], PaletteOrange[1], PaletteOrange[2],
    PaletteOrange[2], PaletteOrange[2], PaletteOrange[3], PaletteOrange[3], PaletteOrange[3], PaletteOrange[3],
    PaletteOrange[4], PaletteOrange[4], PaletteOrange[4], PaletteOrange[4], PaletteOrange[4]};


  G4VisAttributes* visAttributes_x = new G4VisAttributes(color_x);
  G4VisAttributes* visAttributes_u = new G4VisAttributes(color_u);
  G4VisAttributes* visAttributes_v = new G4VisAttributes(color_v);

  G4VisAttributes* Si_att = new G4VisAttributes(Pink);
  G4VisAttributes* FMF2_att = new G4VisAttributes(Red);
  G4VisAttributes* HypHI_RPC_att = new G4VisAttributes(Orange);
  G4VisAttributes* HypHI_Tracker_att = new G4VisAttributes(LightPurple);


  //
  // Import geometry from Root
  //

  // Import geometry from the root file
  new TGeoManager("Wasa", "Geant4 basic example Wasa");

  const std::string nameGeometry = Par.Get<std::string>("Geometry_Namefile");
  gGeoManager->Import(nameGeometry.c_str());

  TFile* f_geo = TFile::Open(nameGeometry.c_str());
  std::vector<std::string>* nameDet = (std::vector<std::string>*)(f_geo->Get("nameDet"));

  if(nameDet != nullptr)
    for(size_t i=0;i<nameDet->size();++i)
      NameDetectorsSD.push_back(nameDet->at(i));

  std::map<std::string, double>* simParameters = (std::map<std::string, double>*)(f_geo->Get("simParameters"));

  if(simParameters!= nullptr)
    Par.ParsePreviousParams(simParameters);

  f_geo->Close();
  f_geo->Delete();
  f_geo = nullptr;

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

  // Visibility

  G4LogicalVolume* InnerLV = FindVolume("INNER");
  if(!Par.IsAvailable("HypHI_InnerTrackerBox_Visible"))
    InnerLV->SetVisAttributes(G4VisAttributes::GetInvisible());

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
  FindVolume("PSCE")->SetVisAttributes(G4VisAttributes(Blue));
  FindVolume("PSBE")->SetVisAttributes(G4VisAttributes(Blue));
  FindVolume("PSFE")->SetVisAttributes(G4VisAttributes(LightBlue));

  std::vector<G4String> NameSD_Color = {"MG01", "MG02", "MG03", "MG04", "MG05", "MG06", "MG07", "MG08", "MG09",
                                        "MG10", "MG11", "MG12", "MG13", "MG14", "MG15", "MG16", "MG17"};

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

  auto WasaLV  = FindVolume("WASA");
  if(WasaLV)
    FindVolume("WASA")->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto MFLD_log  = FindVolume("MFLD");
  if(MFLD_log)
    if(!Par.IsAvailable("HypHI_InnerTrackerBox_Visible"))
      FindVolume("MFLD")->SetVisAttributes(G4VisAttributes::GetInvisible());

  if(!MiniVis)
    {
      auto tempVol = FindVolume("Si1_Strip_log_x");
      if(tempVol != nullptr)
	{
	  FindVolume("Si1_Strip_log_x")->SetVisAttributes(Si_att);
	  FindVolume("Si1_Strip_log_y")->SetVisAttributes(Si_att);
	  FindVolume("Si1_log")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("Si1_log_x")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("Si1_log_y")->SetVisAttributes(G4VisAttributes::GetInvisible());
	}
      tempVol = FindVolume("Si2_Strip_log_x");
      if(tempVol != nullptr)
	{
	  FindVolume("Si2_Strip_log_x")->SetVisAttributes(Si_att);
	  FindVolume("Si2_Strip_log_y")->SetVisAttributes(Si_att);
	  FindVolume("Si2_log")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("Si2_log_x")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("Si2_log_y")->SetVisAttributes(G4VisAttributes::GetInvisible());
	}
      tempVol = FindVolume("SD1_Strip_log_u");
      if(tempVol != nullptr)
	{
	  FindVolume("SD1_Strip_log_u")->SetVisAttributes(Si_att);
	  FindVolume("SD1u_log")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD1u_log_t")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD1u_log_b")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD1_Strip_log_v")->SetVisAttributes(Si_att);
	  FindVolume("SD1v_log")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD1v_log_t")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD1v_log_b")->SetVisAttributes(G4VisAttributes::GetInvisible());
	}
      tempVol = FindVolume("SD2_Strip_log_u");
      if(tempVol != nullptr)
	{
	  FindVolume("SD2_Strip_log_u")->SetVisAttributes(Si_att);
	  FindVolume("SD2u_log")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD2u_log_t")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD2u_log_b")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD2_Strip_log_v")->SetVisAttributes(Si_att);
	  FindVolume("SD2v_log")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD2v_log_t")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD2v_log_b")->SetVisAttributes(G4VisAttributes::GetInvisible());
	}
      tempVol = FindVolume("SD1pad_Strip_log_u");
      if(tempVol != nullptr)
	{
	  FindVolume("SD1pad_Strip_log_u")->SetVisAttributes(Si_att);
	  FindVolume("SD1u_log")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD1pad_Strip_log_v")->SetVisAttributes(Si_att);
	  FindVolume("SD1v_log")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD2pad_Strip_log_u")->SetVisAttributes(Si_att);
	  FindVolume("SD2u_log")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD2pad_Strip_log_v")->SetVisAttributes(Si_att);
	  FindVolume("SD2v_log")->SetVisAttributes(G4VisAttributes::GetInvisible());
	}

      tempVol = FindVolume("FiberD1_Core_log_x");
      if(tempVol != nullptr)
	{
	  FindVolume("FiberD1_Core_log_x")->SetVisAttributes(visAttributes_x);
	  FindVolume("FiberD1_Core_log_u")->SetVisAttributes(visAttributes_u);
	  FindVolume("FiberD1_Core_log_v")->SetVisAttributes(visAttributes_v);
	  FindVolume("FiberD1_Cladding_log_x")->SetVisAttributes(visAttributes_x);
	  FindVolume("FiberD1_Cladding_log_u")->SetVisAttributes(visAttributes_u);
	  FindVolume("FiberD1_Cladding_log_v")->SetVisAttributes(visAttributes_v);
	  FindVolume("FiberD1_log")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD1_log_x")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD1_log_u")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD1_log_v")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD2_Core_log_x")->SetVisAttributes(visAttributes_x);
	  FindVolume("FiberD2_Core_log_u")->SetVisAttributes(visAttributes_u);
	  FindVolume("FiberD2_Core_log_v")->SetVisAttributes(visAttributes_v);
	  FindVolume("FiberD2_Cladding_log_x")->SetVisAttributes(visAttributes_x);
	  FindVolume("FiberD2_Cladding_log_u")->SetVisAttributes(visAttributes_u);
	  FindVolume("FiberD2_Cladding_log_v")->SetVisAttributes(visAttributes_v);
	  FindVolume("FiberD2_log")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD2_log_x")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD2_log_u")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD2_log_v")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD3_Core_log_x")->SetVisAttributes(visAttributes_x);
	  FindVolume("FiberD3_Core_log_u")->SetVisAttributes(visAttributes_u);
	  FindVolume("FiberD3_Core_log_v")->SetVisAttributes(visAttributes_v);
	  FindVolume("FiberD3_Cladding_log_x")->SetVisAttributes(visAttributes_x);
	  FindVolume("FiberD3_Cladding_log_u")->SetVisAttributes(visAttributes_u);
	  FindVolume("FiberD3_Cladding_log_v")->SetVisAttributes(visAttributes_v);
	  FindVolume("FiberD3_log")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD3_log_x")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD3_log_u")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD3_log_v")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD4_Core_log_v")->SetVisAttributes(visAttributes_v);
	  FindVolume("FiberD4_Core_log_u")->SetVisAttributes(visAttributes_u);
	  FindVolume("FiberD4_Core_log_x")->SetVisAttributes(visAttributes_x);
	  FindVolume("FiberD4_Cladding_log_v")->SetVisAttributes(visAttributes_v);
	  FindVolume("FiberD4_Cladding_log_u")->SetVisAttributes(visAttributes_u);
	  FindVolume("FiberD4_Cladding_log_x")->SetVisAttributes(visAttributes_x);
	  FindVolume("FiberD4_log")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD4_log_v")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD4_log_u")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD4_log_x")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD5_Core_log_x")->SetVisAttributes(visAttributes_x);
	  FindVolume("FiberD5_Core_log_u")->SetVisAttributes(visAttributes_u);
	  FindVolume("FiberD5_Core_log_v")->SetVisAttributes(visAttributes_v);
	  FindVolume("FiberD5_Cladding_log_x")->SetVisAttributes(visAttributes_x);
	  FindVolume("FiberD5_Cladding_log_u")->SetVisAttributes(visAttributes_u);
	  FindVolume("FiberD5_Cladding_log_v")->SetVisAttributes(visAttributes_v);
	  FindVolume("FiberD5_log")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD5_log_x")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD5_log_u")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD5_log_v")->SetVisAttributes(G4VisAttributes::GetInvisible());
	}

      tempVol = FindVolume("MiniFiberD1_Core_log_x1");
      if(tempVol != nullptr)
	{
	  FindVolume("MiniFiberD1_Core_log_x1")->SetVisAttributes(visAttributes_x);
	  FindVolume("MiniFiberD1_Core_log_u1")->SetVisAttributes(visAttributes_u);
	  FindVolume("MiniFiberD1_Core_log_v1")->SetVisAttributes(visAttributes_v);
	  FindVolume("MiniFiberD1_Core_log_x2")->SetVisAttributes(visAttributes_x);
	  FindVolume("MiniFiberD1_Core_log_v2")->SetVisAttributes(visAttributes_v);
	  FindVolume("MiniFiberD1_Core_log_u2")->SetVisAttributes(visAttributes_u);
	  FindVolume("MiniFiberD1_Cladding_log_x1")->SetVisAttributes(visAttributes_x);
	  FindVolume("MiniFiberD1_Cladding_log_u1")->SetVisAttributes(visAttributes_u);
	  FindVolume("MiniFiberD1_Cladding_log_v1")->SetVisAttributes(visAttributes_v);
	  FindVolume("MiniFiberD1_Cladding_log_x2")->SetVisAttributes(visAttributes_x);
	  FindVolume("MiniFiberD1_Cladding_log_v2")->SetVisAttributes(visAttributes_v);
	  FindVolume("MiniFiberD1_Cladding_log_u2")->SetVisAttributes(visAttributes_u);
	  FindVolume("MiniFiberD1_log")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("MiniFiberD1_log_x1")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("MiniFiberD1_log_u1")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("MiniFiberD1_log_v1")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("MiniFiberD1_log_x2")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("MiniFiberD1_log_v2")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("MiniFiberD1_log_u2")->SetVisAttributes(G4VisAttributes::GetInvisible());
	}
      tempVol = FindVolume("FMF2_log");
      if(tempVol != nullptr)
	{
	  FMF2_att->SetForceWireframe(false);
	  FindVolume("FMF2_log")->SetVisAttributes(FMF2_att);
	}
      tempVol = FindVolume("HypHI_Endcap_log");
      if(tempVol != nullptr)
	{
	  FindVolume("HypHI_Endcap_log")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  HypHI_RPC_att->SetForceWireframe(false);
	  FindVolume("HypHI_RPC_l_log")->SetVisAttributes(HypHI_RPC_att);
	  FindVolume("HypHI_RPC_h_log")->SetVisAttributes(HypHI_RPC_att);
	  HypHI_Tracker_att->SetForceWireframe(false);
	  FindVolume("HypHI_TrackFwd_log")->SetVisAttributes(HypHI_Tracker_att);
	}
    }
  else
    {
      auto tempVol = FindVolume("Si1_Strip_log_x");
      if(tempVol != nullptr)
	{
	  FindVolume("Si1_Strip_log_x")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("Si1_Strip_log_y")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("Si1_log")->SetVisAttributes(Si_att);
	  FindVolume("Si1_log_x")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("Si1_log_y")->SetVisAttributes(G4VisAttributes::GetInvisible());
	}
      tempVol = FindVolume("Si2_Strip_log_x");
      if(tempVol != nullptr)
	{
	  FindVolume("Si2_Strip_log_x")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("Si2_Strip_log_y")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("Si2_log")->SetVisAttributes(Si_att);
	  FindVolume("Si2_log_x")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("Si2_log_y")->SetVisAttributes(G4VisAttributes::GetInvisible());
	}
      tempVol = FindVolume("SD1_Strip_log_u");
      if(tempVol != nullptr)
	{
	  FindVolume("SD1_Strip_log_u")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD1u_log")->SetVisAttributes(Si_att);
	  FindVolume("SD1u_log_t")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD1u_log_b")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD1_Strip_log_v")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD1v_log")->SetVisAttributes(Si_att);
	  FindVolume("SD1v_log_t")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD1v_log_b")->SetVisAttributes(G4VisAttributes::GetInvisible());
	}
      tempVol = FindVolume("SD2_Strip_log_u");
      if(tempVol != nullptr)
	{
	  FindVolume("SD2_Strip_log_u")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD2u_log")->SetVisAttributes(Si_att);
	  FindVolume("SD2u_log_t")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD2u_log_b")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD2_Strip_log_v")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD2v_log")->SetVisAttributes(Si_att);
	  FindVolume("SD2v_log_t")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD2v_log_b")->SetVisAttributes(G4VisAttributes::GetInvisible());
	}
      tempVol = FindVolume("SD1pad_Strip_log_u");
      if(tempVol != nullptr)
	{
	  FindVolume("SD1pad_Strip_log_u")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD1u_log")->SetVisAttributes(Si_att);
	  FindVolume("SD1pad_Strip_log_v")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD1v_log")->SetVisAttributes(Si_att);
	  FindVolume("SD2pad_Strip_log_u")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD2u_log")->SetVisAttributes(Si_att);
	  FindVolume("SD2pad_Strip_log_v")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("SD2v_log")->SetVisAttributes(Si_att);
	}
      tempVol = FindVolume("FiberD1_Core_log_x");
      if(tempVol != nullptr)
	{

	  FindVolume("FiberD1_Core_log_x")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD1_Core_log_u")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD1_Core_log_v")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD1_Cladding_log_x")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD1_Cladding_log_u")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD1_Cladding_log_v")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD1_log")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD1_log_x")->SetVisAttributes(visAttributes_x);
	  FindVolume("FiberD1_log_u")->SetVisAttributes(visAttributes_u);
	  FindVolume("FiberD1_log_v")->SetVisAttributes(visAttributes_v);
	  FindVolume("FiberD2_Core_log_x")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD2_Core_log_u")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD2_Core_log_v")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD2_Cladding_log_x")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD2_Cladding_log_u")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD2_Cladding_log_v")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD2_log")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD2_log_x")->SetVisAttributes(visAttributes_x);
	  FindVolume("FiberD2_log_u")->SetVisAttributes(visAttributes_u);
	  FindVolume("FiberD2_log_v")->SetVisAttributes(visAttributes_v);
	  FindVolume("FiberD3_Core_log_x")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD3_Core_log_u")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD3_Core_log_v")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD3_Cladding_log_x")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD3_Cladding_log_u")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD3_Cladding_log_v")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD3_log")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD3_log_x")->SetVisAttributes(visAttributes_x);
	  FindVolume("FiberD3_log_u")->SetVisAttributes(visAttributes_u);
	  FindVolume("FiberD3_log_v")->SetVisAttributes(visAttributes_v);
	  FindVolume("FiberD4_Core_log_v")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD4_Core_log_u")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD4_Core_log_x")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD4_Cladding_log_v")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD4_Cladding_log_u")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD4_Cladding_log_x")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD4_log")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD4_log_v")->SetVisAttributes(visAttributes_v);
	  FindVolume("FiberD4_log_u")->SetVisAttributes(visAttributes_u);
	  FindVolume("FiberD4_log_x")->SetVisAttributes(visAttributes_x);
	  FindVolume("FiberD5_Core_log_x")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD5_Core_log_u")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD5_Core_log_v")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD5_Cladding_log_x")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD5_Cladding_log_u")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD5_Cladding_log_v")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD5_log")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("FiberD5_log_x")->SetVisAttributes(visAttributes_x);
	  FindVolume("FiberD5_log_u")->SetVisAttributes(visAttributes_u);
	  FindVolume("FiberD5_log_v")->SetVisAttributes(visAttributes_v);
	}
      tempVol = FindVolume("MiniFiberD1_Core_log_x1");
      if(tempVol != nullptr)
	{
	  FindVolume("MiniFiberD1_Core_log_x1")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("MiniFiberD1_Core_log_u1")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("MiniFiberD1_Core_log_v1")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("MiniFiberD1_Core_log_x2")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("MiniFiberD1_Core_log_v2")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("MiniFiberD1_Core_log_u2")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("MiniFiberD1_Cladding_log_x1")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("MiniFiberD1_Cladding_log_u1")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("MiniFiberD1_Cladding_log_v1")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("MiniFiberD1_Cladding_log_x2")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("MiniFiberD1_Cladding_log_v2")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("MiniFiberD1_Cladding_log_u2")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("MiniFiberD1_log")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("MiniFiberD1_log_x1")->SetVisAttributes(visAttributes_x);
	  FindVolume("MiniFiberD1_log_u1")->SetVisAttributes(visAttributes_u);
	  FindVolume("MiniFiberD1_log_v1")->SetVisAttributes(visAttributes_v);
	  FindVolume("MiniFiberD1_log_x2")->SetVisAttributes(visAttributes_x);
	  FindVolume("MiniFiberD1_log_v2")->SetVisAttributes(visAttributes_v);
	  FindVolume("MiniFiberD1_log_u2")->SetVisAttributes(visAttributes_u);
	}
      tempVol = FindVolume("HypHI_Endcap_log");
      if(tempVol != nullptr)
	{
	  FindVolume("HypHI_Endcap_log")->SetVisAttributes(HypHI_RPC_att);
	  FindVolume("HypHI_RPC_l_log")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("HypHI_RPC_h_log")->SetVisAttributes(G4VisAttributes::GetInvisible());
	  FindVolume("HypHI_TrackFwd_log")->SetVisAttributes(G4VisAttributes::GetInvisible());
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

  G4LogicalVolume* HypHI_Target_log = FindVolume("HypHI_Target_log");
  if(HypHI_Target_log!=nullptr)
    {
      HypHI_Target_log->SetRegion(aTargetRegion);
      aTargetRegion->AddRootLogicalVolume(HypHI_Target_log);
    }
  std::vector<double> cutsTarget(4, Par.Get<double>("TargetRegionCut"));
  aTargetRegion->SetProductionCuts(new G4ProductionCuts());
  aTargetRegion->GetProductionCuts()->SetProductionCuts(cutsTarget);

  experimentalHall_logOutRoot  = world->GetLogicalVolume();
  experimentalHall_physOutRoot = world;

  return world;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WasaFullRootConstruction::DefineMaterials()
{
  // Dummy, as materials are imported via VGM
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* WasaFullRootConstruction::DefineVolumes()
{
  // Dummy, as geometry is imported via VGM

  return 0;
}

//
// end VGM demo
//---------------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void WasaFullRootConstruction::ConstructSDandField()
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

  auto MFLD_log  = FindVolume("MFLD");
  auto MFLD_phys = FindVolPhys("MFLD");

  if(isFieldMap)
    {
      const std::string field_name = Par.Get<std::string>("Field_WASAMap");
      double max_valueField =
	Par.IsAvailable("Field_WASAMapMaxField") ? Par.Get<double>("Field_WASAMapMaxField") : 1. * tesla;

    
      std::cout << "Field origin MFLD Phys:" << MFLD_phys << " " << MFLD_phys->GetInstanceID() << "\n";
      auto transMFLD = MFLD_phys->GetObjectTranslation();
      auto rotMFLD   = MFLD_phys->GetObjectRotationValue();

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
    {
      auto INNER_log  = FindVolume("INNER");
      INNER_log->SetFieldManager(fFieldMgr, forceToAllDaughters);
    }
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

G4VPhysicalVolume* WasaFullRootConstruction::FindVolPhys(const G4String& name)
{
  G4PhysicalVolumeStore* pvStore = G4PhysicalVolumeStore::GetInstance();
  return pvStore->GetVolume(name);
}
