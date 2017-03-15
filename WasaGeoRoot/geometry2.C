#include "TFile.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoPcon.h"
#include "TGeoVolume.h"

#include "TEveGeoNode.h"
#include "TEveManager.h"
#include "TROOT.h"
//#include "TEveVSDStructs.h"
//#include "TEveTrackPropagator.h"
//#include "TEveTrack.h"
//#include "TEveViewer.h"
#include "TGLViewer.h"
#include "TMath.h"
#include "TSystem.h"
//#include "TEvePointSet.h"
#include "Riostream.h"
#include "TTree.h"

#include "TCanvas.h"

#include <vector>

#include "geometryVolEMC.C"
#include "geometryVolFW.C"
#include "geometryVolume.C"

void geometry()
{
  //
  //  This file has been generated automatically via the root
  //  utility g2root from an interactive version of GEANT
  //   (see ROOT class TGeoManager for an example of use)
  //
  gSystem->Load("libGeom");
  TGeoRotation* rot;
  TGeoNode *Node, *Node1;

  TGeoManager* oldgeo1 = gGeoManager;
  gGeoManager = 0;

  TEveManager* oldeve1 = gEve;
  gEve = 0;

  std::cout << " Change gGeoManager !" << oldgeo1 << " " << oldeve1 << std::endl;
  if(0 != oldgeo1)
    {
      delete oldgeo1;
      oldgeo1 = 0;
    }
  if(oldeve1 != 0)
    {
      delete oldeve1;
      oldeve1 = 0;
    }

  TGeoManager* geometry = new TGeoManager("geometry", "geometry.C");

  std::vector<TGeoMaterial*> list_mat;
  //-----------List of Materials and Mixtures--------------

  TGeoMaterial* mat1 = new TGeoMaterial(" VAC", 14.61, 7.3, 0.1205000E-05);
  mat1->SetUniqueID(1);
  list_mat.push_back(mat1);
  TGeoMaterial* mat2 = new TGeoMaterial(" LH2", 1.01, 1, 0.7080001E-01);
  mat2->SetUniqueID(2);
  list_mat.push_back(mat2);
  TGeoMaterial* mat6 = new TGeoMaterial("  BE", 9.01, 4, 1.848);
  mat6->SetUniqueID(6);
  list_mat.push_back(mat6);

  TGeoMaterial* mat7 = new TGeoMaterial("   C", 12.01, 6, 2.265);
  mat7->SetUniqueID(7);
  list_mat.push_back(mat7);

  TGeoMaterial* mat11 = new TGeoMaterial("  AL", 26.98, 13, 2.7);
  mat11->SetUniqueID(11);
  list_mat.push_back(mat11);

  TGeoMaterial* mat15 = new TGeoMaterial("  FE", 55.85, 26, 7.87);
  mat15->SetUniqueID(15);
  list_mat.push_back(mat15);

  TGeoMaterial* mat16 = new TGeoMaterial("  CU", 63.54, 29, 8.96);
  mat16->SetUniqueID(16);
  list_mat.push_back(mat16);

  TGeoMaterial* mat18 = new TGeoMaterial("  PB", 207.19, 82, 11.35);
  mat18->SetUniqueID(18);
  list_mat.push_back(mat18);

  TGeoMaterial* mat22 = new TGeoMaterial(" MYL", 12.876, 6.456, 1.39);
  mat22->SetUniqueID(22);
  list_mat.push_back(mat22);

  TGeoMixture* mat23 = new TGeoMixture("PSCI", 2, 1.03200);
  mat23->SetUniqueID(23);
  mat23->DefineElement(0, 1.01, 1, 0.8495516E-01);
  mat23->DefineElement(1, 12.01, 6, 0.9150448);
  list_mat.push_back(mat23);
  TGeoMaterial* mat26 = new TGeoMaterial(" CSI", 129.974, 54.023, 4.5);
  mat26->SetUniqueID(26);
  list_mat.push_back(mat26);
  TGeoMaterial* mat28 = new TGeoMaterial("PXGL", 11.157, 5.612, 1.2);
  mat28->SetUniqueID(28);
  list_mat.push_back(mat28);
  TGeoMaterial* mat29 = new TGeoMaterial(" H2O", 14.333, 7.222, 1);
  mat29->SetUniqueID(29);
  list_mat.push_back(mat29);
  TGeoMaterial* mat30 = new TGeoMaterial("SIO2", 21.649, 7.077, 2.203);
  mat30->SetUniqueID(30);
  list_mat.push_back(mat30);
  TGeoMaterial* mat31 = new TGeoMaterial(" AIR", 14.61, 7.3, 0.1205000E-02);
  mat31->SetUniqueID(31);
  list_mat.push_back(mat31);
  TGeoMaterial* mat39 = new TGeoMaterial(" CO2", 14.911, 7.454, 0.1842000E-02);
  mat39->SetUniqueID(39);
  list_mat.push_back(mat39);
  //-----------List of Tracking Media--------------

  std::vector<TGeoMedium*> list_med;

  TGeoMedium* med1 = new TGeoMedium("VAC", 1, 1, 0, 0, 0, 5, 0.1000000E+11, 0.2499637, 0.5000000E-02, 30.65114);
  TGeoMedium* med2 = new TGeoMedium("VAC+MFIELD", 2, 1, 0, 1, 11.69974, 5, 0.1000000E+11, 0.2499637, 0.5000000E-02, 30.65114);
  TGeoMedium* med3 = new TGeoMedium("AIR", 3, 31, 0, 0, 0, 5, 0.1000000E+11, 0.2488529, 0.5000000E-02, 0.9692742);
  TGeoMedium* med4 = new TGeoMedium("STEEL", 4, 15, 0, 0, 0, 5, 0.1000000E+11, 0.25, 0.5000000E-02, 0.2419303E-01);
  TGeoMedium* med5 = new TGeoMedium("STEEL+MFIELD", 5, 15, 0, 1, 11.69974, 5, 0.1000000E+11, 0.25, 0.5000000E-02, 0.2419303E-01);
  TGeoMedium* med6 = new TGeoMedium("BE", 6, 6, 0, 1, 11.69974, 5, 0.1000000E+11, 0.2163378, 0.5000000E-02, 0.1960746E-01);
  TGeoMedium* med7 = new TGeoMedium("AL", 7, 11, 0, 0, 0, 5, 0.1000000E+11, 0.1829598, 0.5000000E-02, 0.2794475E-01);
  TGeoMedium* med8 = new TGeoMedium("CU", 8, 16, 0, 0, 0, 5, 0.1000000E+11, 0.25, 0.5000000E-02, 0.2436345E-01);
  TGeoMedium* med9 = new TGeoMedium("MYLAR", 9, 22, 0, 1, 11.69974, 5, 0.1000000E+11, 0.2126021, 0.5000000E-02, 0.2694496E-01);
  TGeoMedium* med10 = new TGeoMedium("PXGLASS", 10, 28, 0, 0, 0, 5, 0.1000000E+11, 0.2159003, 0.5000000E-02, 0.2795223E-01);
  TGeoMedium* med11 = new TGeoMedium("STEEL", 11, 15, 1, 0, 0, 5, 0.1000000E+11, 0.1507557, 0.5000000E-02, 0.4559045E-01);
  TGeoMedium* med16 = new TGeoMedium("FD.FLANGES", 16, 15, 0, 0, 0, 5, 0.1000000E+11, 0.25, 0.5000000E-02, 0.2419303E-01);
  TGeoMedium* med21 = new TGeoMedium("CD.FLANGES", 21, 15, 0, 0, 0, 5, 0.1000000E+11, 0.25, 0.5000000E-02, 0.2419303E-01);
  TGeoMedium* med22 = new TGeoMedium("CD.FLANGES", 22, 15, 0, 0, 0, 5, 0.1000000E+11, 0.25, 0.5000000E-02, 0.2419303E-01);
  TGeoMedium* med23 = new TGeoMedium("CD.FLANGES", 23, 15, 1, 0, 0, 5, 0.1000000E+11, 0.1507557, 0.5000000E-02, 0.4559045E-01);
  TGeoMedium* med26 = new TGeoMedium("CD.BARS", 26, 15, 0, 0, 0, 5, 0.1000000E+11, 0.25, 0.5000000E-02, 0.2419303E-01);
  TGeoMedium* med27 = new TGeoMedium("CD.BARS", 27, 15, 0, 0, 0, 5, 0.1000000E+11, 0.25, 0.5000000E-02, 0.2419303E-01);
  TGeoMedium* med28 = new TGeoMedium("CD.BARS", 28, 15, 0, 0, 0, 5, 0.1000000E+11, 0.25, 0.5000000E-02, 0.2419303E-01);
  TGeoMedium* med29 = new TGeoMedium("CD.BARS", 29, 15, 0, 0, 0, 5, 0.1000000E+11, 0.25, 0.5000000E-02, 0.2419303E-01);
  TGeoMedium* med34 = new TGeoMedium("SOL.CHIMNEY", 34, 15, 0, 0, 0, 5, 0.1000000E+11, 0.25, 0.5000000E-02, 0.2419303E-01);
  TGeoMedium* med35 = new TGeoMedium("SOL.CHIMNEY", 35, 15, 0, 0, 0, 5, 0.1000000E+11, 0.25, 0.5000000E-02, 0.2419303E-01);
  TGeoMedium* med42 = new TGeoMedium("PEL.TUB", 42, 15, 0, 0, 0, 5, 0.1000000E+11, 0.25, 0.5000000E-02, 0.2419303E-01);
  TGeoMedium* med43 = new TGeoMedium("PEL.TUB", 43, 15, 0, 0, 0, 5, 0.1000000E+11, 0.25, 0.5000000E-02, 0.2419303E-01);
  TGeoMedium* med44 = new TGeoMedium("PEL.TUB", 44, 15, 0, 0, 0, 5, 0.1000000E+11, 0.25, 0.5000000E-02, 0.2419303E-01);
  TGeoMedium* med45 = new TGeoMedium("PEL.TUB", 45, 15, 0, 0, 0, 5, 0.1000000E+11, 0.25, 0.5000000E-02, 0.2419303E-01);
  TGeoMedium* med51 = new TGeoMedium("FHD.COVER", 51, 22, 0, 0, 0, 5, 0.1000000E+11, 0.2126021, 0.5000000E-02, 0.2694496E-01);
  TGeoMedium* med58 = new TGeoMedium("FD.ABSORBER", 58, 28, 0, 0, 0, 5, 0.1000000E+11, 0.2159003, 0.5000000E-02, 0.2795223E-01);
  TGeoMedium* med59 = new TGeoMedium("FD.ABSORBER", 59, 28, 0, 0, 0, 5, 0.1000000E+11, 0.2159003, 0.5000000E-02, 0.2795223E-01);
  TGeoMedium* med60 = new TGeoMedium("FD.ABSORBER", 60, 28, 0, 0, 0, 5, 0.1000000E+11, 0.2159003, 0.5000000E-02, 0.2795223E-01);
  TGeoMedium* med61 = new TGeoMedium("FD.ABSORBER", 61, 28, 0, 0, 0, 5, 0.1000000E+11, 0.2159003, 0.5000000E-02, 0.2795223E-01);
  TGeoMedium* med62 = new TGeoMedium("FD.ABSORBER", 62, 15, 0, 0, 0, 5, 0.1000000E+11, 0.25, 0.5000000E-02, 0.2419303E-01);
  TGeoMedium* med68 = new TGeoMedium("FD.MWPC", 68, 31, 1, 0, 0, 5, 0.1000000E+11, 0.1147079E-02, 0.5000000E-02, 0.1389792E-01);
  TGeoMedium* med69 = new TGeoMedium("CO2", 69, 39, 1, 1, 11.69974, 5, 0.1000000E+11, 0.1478443E-02, 0.5000000E-02, 0.1512986E-01);
  TGeoMedium* med71 = new TGeoMedium("FD.SCINT", 71, 23, 1, 0, 0, 5, 0.1000000E+11, 0.3081252E-01, 0.5000000E-02, 0.1028338E-01);
  TGeoMedium* med72 = new TGeoMedium("CD.SCINT", 72, 23, 1, 1, 11.69974, 5, 0.1000000E+11, 0.3081252E-01, 0.5000000E-02, 0.1028338E-01);
  TGeoMedium* med73 = new TGeoMedium("SEC", 73, 26, 1, 1, 11.69974, 5, 0.1000000E+11, 0.1466471, 0.5000000E-02, 0.9235660E-01);

  list_med.push_back(med1);
  list_med.push_back(med2);
  list_med.push_back(med3);
  list_med.push_back(med4);
  list_med.push_back(med5);
  list_med.push_back(med6);
  list_med.push_back(med7);
  list_med.push_back(med8);
  list_med.push_back(med9);
  list_med.push_back(med10);
  list_med.push_back(med11);
  list_med.push_back(med16);
  list_med.push_back(med21);
  list_med.push_back(med22);
  list_med.push_back(med23);
  list_med.push_back(med26);
  list_med.push_back(med27);
  list_med.push_back(med28);
  list_med.push_back(med29);
  list_med.push_back(med34);
  list_med.push_back(med35);
  list_med.push_back(med42);
  list_med.push_back(med43);
  list_med.push_back(med44);
  list_med.push_back(med45);
  list_med.push_back(med51);
  list_med.push_back(med58);
  list_med.push_back(med59);
  list_med.push_back(med60);
  list_med.push_back(med61);
  list_med.push_back(med62);
  list_med.push_back(med68);
  list_med.push_back(med69);
  list_med.push_back(med71);
  list_med.push_back(med72);
  list_med.push_back(med73);

  TGeoVolume* WASA = gGeoManager->MakeTube("WASA", med1, 0, 150, 300);

  // TGeoVolumeMulti *ESC_ = gGeoManager->MakeVolumeMulti("ESC", med11);
  // ESC_->AddVolume(gGeoManager->MakeTube("ESC",med11,0,150,0.5000000E-01));
  // ESC_->AddVolume(gGeoManager->MakeTube("ESC",med11,0,150,0.5000000E-01));
  // ESC_->AddVolume(gGeoManager->MakeTube("ESC",med11,149.9,150,300));

  TGeoVolume* MFLD = geometryVolumeCentral(WASA, list_mat, list_med);

  geometryVolEMC(MFLD, list_mat, list_med);
  geometryVolFW(WASA, list_mat, list_med);

  //-----------List of Nodes--------------

  gGeoManager->SetTopVolume(WASA);

  // WASA->AddNode(ESC_->GetVolume(0),1,new TGeoTranslation(0,0,-299.95));
  // WASA->AddNode(ESC_->GetVolume(1),3,new TGeoTranslation(0,0,299.95));
  // WASA->AddNode(ESC_->GetVolume(2),2,new TGeoTranslation(0,0,-2.5));

  gGeoManager->CloseGeometry();

  TEveManager::Create();
  TGeoNode* node = gGeoManager->GetTopNode();
  TEveGeoTopNode* en = new TEveGeoTopNode(gGeoManager, node);
  en->SetVisLevel(4);
  en->GetNode()->GetVolume()->SetVisibility(kFALSE);

  gEve->AddGlobalElement(en);
  gEve->Redraw3D(kTRUE);
}
