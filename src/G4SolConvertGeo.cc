// --------------------------------------------------------- 
// Implementation of the G4SolConvertGeo class
// Created by C.Rappold (c.rappold@gsi.de)
//----------------------------------------------------------

#include "G4SolConvertGeo.hh"

#include "Geant4GM/volumes/Factory.h" 
#include "RootGM/volumes/Factory.h" 

#include "TGeoManager.h"
#include "TFile.h"
#include "TROOT.h"                                                        
#include "TColor.h"

G4SolConvertGeo::G4SolConvertGeo(G4VPhysicalVolume* w):physiWorld(w)
{ }
G4SolConvertGeo::~G4SolConvertGeo()
{ }

int G4SolConvertGeo::Convert(const std::string& nameOut, const G4String& nameGeometry)
{
  // Import Geant4 geometry to VGM                       

  if(gGeoManager!=nullptr)
    {
      std::cout<<"!> G4SolConvertGeo: gGeoManager not null !\n";
      gGeoManager->Print();
      TGeoManager* oldgeo1 = gGeoManager;
      gGeoManager = nullptr;
      delete oldgeo1;
      oldgeo1 = nullptr;
    }
  
  Geant4GM::Factory g4Factory;
  g4Factory.Import(physiWorld);                          
  // where physiWorld is of G4VPhysicalVolume* type 
  
  // Export VGM geometry to Root                         
  RootGM::Factory rtFactory;
  g4Factory.Export(&rtFactory);
   
  const int ColorCDC[15] = {0,0,0,
			    1,1,1,1,
			    2,2,
			    3,3,3,3,
			    4,4};
		      

  auto set_color = [] (const std::string& nameV, Color_t indexColor) {
    TGeoVolume* geoVol = nullptr;
    geoVol = gGeoManager->GetVolume(nameV.c_str());
    if(geoVol!=nullptr)
      geoVol->SetLineColor(indexColor);
    else
      std::cout<<"E> set_color : no GeoVolume with this name :"<<nameV<<std::endl;
  };

  if(nameGeometry == "CDS" && nameGeometry == "WasaSimple")
    {
      std::string tempNameGeo ("CDC_SetLogR_");
      for(int i=0;i<15;++i)
	{
	  std::string nameTemp(tempNameGeo);
	  nameTemp+= std::to_string(i);
	  set_color(nameTemp,kOrange+ColorCDC[i]);
	}
      set_color("CDC_win_out_log",kViolet-6);
      set_color("CDC_win_log",kViolet-6);
      
      set_color("CDH_log",kRed);
      
      set_color("CDSYoke_log",kBlue);
      set_color("CDS_endcap_log",kBlue);      
		
    }
  else if(nameGeometry == "Wasa")
    {
      
      set_color("SOAL",kSpring-8);
      set_color("COIL",kSpring-8);

      set_color("IRC",kBlue);
      set_color("SOC0",kBlue);
      set_color("SOC1",kBlue);

      set_color("PSCE",kRed-4);
      set_color("PSBE",kRed+1);
      set_color("PSFE",kRed+1);
      
      set_color("PSL0",kYellow+1);
      set_color("PSL1",kYellow+1);

      std::vector<std::string> namesDC = {"MG01","MG02","MG03","MG04","MG05","MG06","MG07","MG08","MG09","MG10","MG11","MG12","MG13","MG14","MG15","MG16","MG17"};
      int index_color = -8;
      for(auto nameDC : namesDC)
	{
	  set_color(nameDC, kOrange+index_color);
	  ++index_color;
	}
      set_color("PTB0",kGray+1);
      set_color("PTB1",kGray+1);
      set_color("PTB2",kGray+1);
      set_color("PTB3",kGray+1);

      set_color("MDO",kMagenta-8);
      set_color("MDB",kMagenta-8);
      set_color("MDF",kMagenta-8);
    }
  else if(nameGeometry == "HIAFSimple")
    {
      std::string tempNameGeo ("GEM_CD_log");
      for(int i=0;i<10;++i)
	{
	  std::string nameTemp(tempNameGeo);
	  nameTemp+= std::to_string(i);
	  set_color(nameTemp,kOrange+ColorCDC[i]);
	}
      set_color("Solenoid_log",kGreen-7);
      set_color("PSB_log",kRed+1);
      set_color("Feyoke_up",kGreen-7);
      set_color("Feyoke_down",kGreen-7);
      set_color("Feyoke_left",kGreen-7);
      set_color("Feyoke_right",kGreen-7);
      set_color("LastTrackerPi_log",kMagenta+1);
      set_color("LastTrackerFrag_log",kMagenta+1);
      set_color("LastRPCPi_log",kOrange-1);
      set_color("LastRPCFrag_log",kOrange-1);
    }
  
  set_color("HypHI_TrackFwd_log",kMagenta+2);
  
  set_color("HypHI_RPC_l_log",kOrange-1);
  set_color("HypHI_RPC_h_log",kOrange-1);

  set_color("FMF2_log",kRed);
  
  set_color("HypHI_Target_log",kWhite);
  
  std::string nameSi("HypHI_InSi_log");
  for(int i=0;i<4;++i)
    {
      std::string nameTemp(nameSi);
      nameTemp+=std::to_string(i);
      set_color(nameTemp,kPink+1);
    }


  gGeoManager->CloseGeometry();

  TFile* f_outGeo = new TFile(nameOut.c_str(),"RECREATE");
  f_outGeo->cd();
  gGeoManager->Write();
  f_outGeo->Close();

  //rtFactory.World();                              
  // returns Root top node, of TGeoNode* type  
  return 1;
}
