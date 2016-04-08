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

int G4SolConvertGeo::Convert(const std::string& nameOut)
{
  // Import Geant4 geometry to VGM                       
  
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
  
  set_color("HypHI_TrackFwd_log",kMagenta+2);

  set_color("HypHI_RPC_l_log",kOrange-1);
  set_color("HypHI_RPC_h_log",kOrange-1);

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
