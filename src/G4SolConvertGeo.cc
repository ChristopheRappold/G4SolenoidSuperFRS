// ---------------------------------------------------------
// Implementation of the G4SolConvertGeo class
// Created by C.Rappold (c.rappold@gsi.de)
//----------------------------------------------------------

#include "G4SolConvertGeo.hh"

#include "Geant4GM/volumes/Factory.h"
#include "RootGM/volumes/Factory.h"
#include "TColor.h"
#include "TFile.h"
#include "TGeoManager.h"
#include "TROOT.h"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

G4SolConvertGeo::G4SolConvertGeo(G4VPhysicalVolume* w) : physiWorld(w) {}
G4SolConvertGeo::~G4SolConvertGeo() {}

int G4SolConvertGeo::Convert(const std::string& nameOut, const G4String& nameGeometry, const std::vector<G4String>& NameDetectorsSD, const G4SolConfig& Config)
{
  // Import Geant4 geometry to VGM

  if(gGeoManager != nullptr)
    {
      std::cout << "!> G4SolConvertGeo: gGeoManager not null !\n";
      gGeoManager->Print();
      TGeoManager* oldgeo1 = gGeoManager;
      gGeoManager          = nullptr;
      delete oldgeo1;
      oldgeo1 = nullptr;
    }

  Geant4GM::Factory g4Factory;
  g4Factory.Import(physiWorld);
  // where physiWorld is of G4VPhysicalVolume* type

  // Export VGM geometry to Root
  RootGM::Factory rtFactory;
  g4Factory.Export(&rtFactory);

  const int ColorCDC[15] = {0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4};

  auto set_color = [](const std::string& nameV, Color_t indexColor) {
    TGeoVolume* geoVol = nullptr;
    geoVol             = gGeoManager->GetVolume(nameV.c_str());
    if(geoVol != nullptr)
      geoVol->SetLineColor(indexColor);
    else
      std::cout << "E> set_color : no GeoVolume with this name :" << nameV << std::endl;
  };

  auto set_invisible = [](const std::string& nameV) {
    TGeoVolume* geoVol = nullptr;
    geoVol             = gGeoManager->GetVolume(nameV.c_str());
    if(geoVol != nullptr)
      geoVol->SetVisibility(false);
    else
      std::cout << "E> set_invisible : no GeoVolume with this name :" << nameV << std::endl;
  };

  if(nameGeometry == "CDS" && nameGeometry == "WasaSimple")
    {
      std::string tempNameGeo("CDC_SetLogR_");
      for(int i = 0; i < 15; ++i)
        {
          std::string nameTemp(tempNameGeo);
          nameTemp += std::to_string(i);
          set_color(nameTemp, kOrange + ColorCDC[i]);
        }
      set_color("CDC_win_out_log", kViolet - 6);
      set_color("CDC_win_log", kViolet - 6);

      set_color("CDH_log", kRed);

      set_color("CDSYoke_log", kBlue);
      set_color("CDS_endcap_log", kBlue);
    }
  else if(nameGeometry == "Wasa")
    {

      set_color("SOAL", kSpring - 8);
      set_color("COIL", kSpring - 8);

      set_color("IRC", kBlue);
      set_color("SOC0", kBlue);
      set_color("SOC1", kBlue);

      set_color("PSCE", kRed - 4);
      set_color("PSBE", kRed + 1);
      set_color("PSFE", kRed + 1);

      set_color("PSL0", kYellow + 1);
      set_color("PSL1", kYellow + 1);

      std::vector<std::string> namesDC = {"ME01", "ME02", "ME03", "ME04", "ME05", "ME06", "ME07", "ME08", "ME09",
                                          "ME10", "ME11", "ME12", "ME13", "ME14", "ME15", "ME16", "ME17"};

      int index_color = -8;

      for(auto nameDC : namesDC)
        {
          set_color(nameDC, kOrange + index_color);
          ++index_color;
        }
      set_color("PTB0", kGray + 1);
      set_color("PTB1", kGray + 1);
      set_color("PTB2", kGray + 1);
      set_color("PTB3", kGray + 1);

      set_color("MDO", kMagenta - 8);
      set_color("MDB", kMagenta - 8);
      set_color("MDF", kMagenta - 8);

      std::vector<std::string> namesSub1 = {"_Core_log_","_Cladding_log_"};
      std::vector<std::string> namesSub2 = {"x","u","v"};
      std::vector<std::string> namesFiber = {"FiberD1","FiberD2","FiberD3","FiberD4","FiberD5"};

      for(auto nameF1 : namesFiber)
	for(auto nameS2 : namesSub2)
	  set_color(nameF1+namesSub1[0]+nameS2, kViolet -1);

      for(auto nameF1 : namesFiber)
	for(auto nameS2 : namesSub2)
	  set_invisible(nameF1+namesSub1[1]+nameS2);


      std::vector<std::string> namesMiniFiber = {"MiniFiberD1_Core_log_x1","MiniFiberD1_Core_log_u1","MiniFiberD1_Core_log_v1","MiniFiberD1_Core_log_x2","MiniFiberD1_Core_log_u2","MiniFiberD1_Core_log_v2"};
      for(auto nameF : namesMiniFiber)
	set_color(nameF, kMagenta - 9);

      std::vector<std::string> namesMiniFiberInv = {"MiniFiberD1_Cladding_log_x1","MiniFiberD1_Cladding_log_u1","MiniFiberD1_Cladding_log_v1","MiniFiberD1_Cladding_log_x2","MiniFiberD1_Cladding_log_u2","MiniFiberD1_Cladding_log_v2"};
      for(auto nameF : namesMiniFiberInv)
	set_invisible(nameF);

    }
  else if(nameGeometry == "HIAFSimple")
    {
      std::string tempNameGeo("GEM_CD_log");
      for(int i = 0; i < 10; ++i)
        {
          std::string nameTemp(tempNameGeo);
          nameTemp += std::to_string(i);
          set_color(nameTemp, kOrange + ColorCDC[i]);
        }
      set_color("Solenoid_log", kGreen - 7);
      set_color("PSB_log", kRed + 1);
      set_color("Feyoke_up", kGreen - 7);
      set_color("Feyoke_down", kGreen - 7);
      set_color("Feyoke_left", kGreen - 7);
      set_color("Feyoke_right", kGreen - 7);
      set_color("LastTrackerPi_log", kMagenta + 1);
      set_color("LastTrackerFrag_log", kMagenta + 1);
      set_color("LastRPCPi_log", kOrange - 1);
      set_color("LastRPCFrag_log", kOrange - 1);
    }

  set_color("HypHI_TrackFwd_log", kMagenta + 2);

  set_color("HypHI_RPC_l_log", kOrange - 1);
  set_color("HypHI_RPC_h_log", kOrange - 1);

  set_color("FMF2_log", kRed);

  set_color("HypHI_Target_log", kWhite);

  std::string nameSi("HypHI_InSi_log");
  for(int i = 0; i < 4; ++i)
    {
      std::string nameTemp(nameSi);
      nameTemp += std::to_string(i);
      set_color(nameTemp, kPink + 1);
    }

  set_color("Si1_Strip_log_x",kPink+1);
  set_color("Si1_Strip_log_y",kPink+1);
  set_color("Si2_Strip_log_x",kPink+2);
  set_color("Si2_Strip_log_y",kPink+2);

  set_color("SD1_Strip_log_u",kPink+1);
  set_color("SD1_Strip_log_v",kPink+1);
  set_color("SD2_Strip_log_u",kPink+2);
  set_color("SD2_Strip_log_v",kPink+2);

  set_color("SD1pad_Strip_log_u",kPink+1);
  set_color("SD1pad_Strip_log_v",kPink+1);
  set_color("SD2pad_Strip_log_u",kPink+2);
  set_color("SD2pad_Strip_log_v",kPink+2);


  gGeoManager->CloseGeometry();

  TFile* f_outGeo = new TFile(nameOut.c_str(), "RECREATE");
  f_outGeo->cd();
  gGeoManager->Write();

  std::map<std::string, double> parameterToFile;
  parameterToFile.insert(std::make_pair("Target_Size", Config.Get<double>("Target_Size") / cm));
  parameterToFile.insert(std::make_pair("Target_PosX", Config.Get<double>("Target_PosX") / cm));
  parameterToFile.insert(std::make_pair("Target_PosY", Config.Get<double>("Target_PosY") / cm));
  parameterToFile.insert(std::make_pair("Target_PosZ", Config.Get<double>("Target_PosZ") / cm));
  parameterToFile.insert(std::make_pair("Wasa_ShiftZ", Config.Get<double>("Wasa_ShiftZ") / cm));
  parameterToFile.insert(std::make_pair("Wasa_Side", Config.Get<int>("Wasa_Side")));
  parameterToFile.insert(std::make_pair("Systematic_Shift", Config.Get<double>("Systematic_Shift") / cm));
  if(Config.IsAvailable("Field_CDS_Bz"))
    parameterToFile.insert(std::make_pair("Field_CDS_Bz", Config.Get<double>("Field_CDS_Bz") / tesla));
  else if(Config.IsAvailable("Field_WASAMapMaxField"))
    {
      parameterToFile.insert(std::make_pair("Field_CDS_Bz", Config.Get<double>("Field_WASAMapMaxField") / tesla));
      parameterToFile.insert(std::make_pair("Field_CDS_FieldMap", 1.0));
    }

  f_outGeo->WriteObjectAny(&NameDetectorsSD, "std::vector<std::string>", "nameDet");
  f_outGeo->WriteObjectAny(&parameterToFile, "std::map<std::string,double>", "simParameters");
  
  f_outGeo->Close();

  // rtFactory.World();
  // returns Root top node, of TGeoNode* type
  return 1;
}
