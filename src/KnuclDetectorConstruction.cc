// ====================================================================
//    KnuclDetectorConstruction.cc
//
//
// ====================================================================

#include "KnuclDetectorConstruction.hh"
//#include "KnuclRegionInformation.hh"
#include "globals.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4TwistedTubs.hh"

#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"

#include "G4VisAttributes.hh"
#include "G4Color.hh"

#include "G4UserLimits.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"

#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
//#include "KnuclMaterialManager.hh"
//#include "KnuclFieldSetup.hh"
//#include "KnuclCommon.h"
#include "G4SolSimpleMagneticField.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4AutoDelete.hh"


G4ThreadLocal G4SolSimpleMagneticField* KnuclDetectorConstruction::fMagneticField = 0;
G4ThreadLocal G4FieldManager* KnuclDetectorConstruction::fFieldMgr = 0;


// ====================================================================
//
// class description
//
// ====================================================================

//////////////////////////////////////////////////////
KnuclDetectorConstruction::KnuclDetectorConstruction()//KnuclAnaManager* ana)
  :  experimentalHall_box(0), experimentalHall_log(0), experimentalHall_phys(0),
     Kurama_box(0),           Kurama_log(0),           Kurama_phys(0),
     TOF_box(0),              TOF_log(0),              TOF_phys(0),
     CDS_tube(0),             CDS_log(0),              CDS_phys(0)
     //////////////////////////////////////////////////////
{
  TofRefPos           = 0.;//ana->GetTofRefPos(); 
  Ncdh                = 75;//ana->GetNdch();

  TargetMaterialID    = 1;//ana->GetTargetMaterialID();
  ChamberGasSelection = "ArCO2Methan89-10-1";//ana->GetChamberGas();
  CDCType             = "A";//ana->GetCDCType();
  
  FieldInCDC         = 0.5;
  FieldInKurama      = 1.0;

  ThicknessOfTgtCell = 0.000001;
  ThicknessOfTgtAl   = 0.000001;
  ThicknessOfTgtCFRP = 0.000001;
  ThicknessOfCDCwCFRP= 0.000001;

  TargetLength       = 1.0*cm;
  //BindingEnergy      = 100.0;
  //DecayWidth         =   0.0;

  CDSLength          = 1000.0*mm;
  ChamberSmear       = 1;
  CDCResolution      = 0.25*mm;
  ZVCLength          = 300.0*mm;
  ZVertexChamber     = 1;
  ZVertexChamberCell = 10.0*mm;



}

//////////////////////////////////////////////////////
KnuclDetectorConstruction::~KnuclDetectorConstruction()
//////////////////////////////////////////////////////
{
  //if (fEmFieldSetupKURAMA) delete fEmFieldSetupKURAMA;
  //if (fEmFieldSetupCDC)    delete fEmFieldSetupCDC;
}

void KnuclDetectorConstruction::ConstructMaterials()
{
  // G4NistManager* materialMgr = G4NistManager::Instance();
  // G4Material* He4 = materialMgr->FindOrBuildMaterial("G4_He");
  
  // double density = 
  // G4Material* HeEthan_50_50 = new G4Material(name="Water", density, ncomponents=2);
  
  // G4Material* HeEthan_50_50 = materialMgr-> GetMaterial("He_Ethan_50_50");
  // G4Material* ArCH4_90_10  = materialMgr-> GetMaterial("ArCH4_90_10");
  // G4Material* ArEthan_50_50 = materialMgr-> GetMaterial("ArEthan_50_50");
  // G4Material* HeIsobutane_80_20 = materialMgr-> GetMaterial("He_Isobutane_80_20");
  // G4Material* ArEthan_80_20 = materialMgr-> GetMaterial("ArEthan_80_20");
  // G4Material* ArCO2Methan_89_10_1 = materialMgr-> GetMaterial("ArCO2Methan_89_10_1");
  //G4Material* CFRP = materialMgr-> GetMaterial("CFRP");
  //G4Material* LHe3 = materialMgr->FindOrBuildMaterial("LHelium-3");
  //G4Material* D2Gas = materialMgr->FindOrBuildMaterial("DeuteriumGas");
  //G4Material* SUS = materialMgr->FindOrBuildMaterial("SUS");
  
}

//////////////////////////////////////////////////////
G4VPhysicalVolume* KnuclDetectorConstruction::Construct()
//////////////////////////////////////////////////////
{

  G4Colour  Blue    (0.0, 0.0, 1.0) ;
  G4Colour  Gray (0.5, 0.5, 0.5) ;
  G4Colour  Red     (1.0, 0.0, 0.0) ;
  G4Colour  LightBlue    (0.0, 0.0, 1.0,0.7) ;
  G4Colour  Yellow  (1.0, 1.0, 0.0) ; 
  G4Colour  Purple (.5, 0.0, .5) ;
  G4Colour  Green   (0.0, 1.0, 0.0);
  G4Colour  Orange (1.0,0.647,0);
  G4Colour  Pink (1.0,0.753,0.796);  
  // ==============================================================
  // materials
  // ==============================================================

  //KnuclMaterialManager* materialMgr = new KnuclMaterialManager();
  ConstructMaterials();
  G4NistManager* materialMgr = G4NistManager::Instance();
  
  // G4Material* HeEthan_50_50 = materialMgr-> GetMaterial("He_Ethan_50_50");
  // G4Material* ArCH4_90_10  = materialMgr-> GetMaterial("ArCH4_90_10");
  // G4Material* ArEthan_50_50 = materialMgr-> GetMaterial("ArEthan_50_50");
  // G4Material* HeIsobutane_80_20 = materialMgr-> GetMaterial("He_Isobutane_80_20");
  // G4Material* ArEthan_80_20 = materialMgr-> GetMaterial("ArEthan_80_20");
  // G4Material* ArCO2Methan_89_10_1 = materialMgr-> GetMaterial("ArCO2Methan_89_10_1");

  G4Material* Air = materialMgr->FindOrBuildMaterial("G4_AIR");
  G4Material* Vacuum = materialMgr->FindOrBuildMaterial("G4_Galactic");
  G4Material* HeGas = materialMgr->FindOrBuildMaterial("G4_He");
  G4Material* Fe = materialMgr->FindOrBuildMaterial("G4_Fe");
  G4Material* Scinti = materialMgr->FindOrBuildMaterial("G4_POLYETHYLENE");//"Plastic");
  //G4Material* CFRP = materialMgr-> GetMaterial("CFRP");
  G4Material* Mylar = materialMgr->FindOrBuildMaterial("G4_MYLAR");
  G4Material* W = materialMgr->FindOrBuildMaterial("G4_W");
  G4Material* Al = materialMgr->FindOrBuildMaterial("G4_Al");
  //G4Material* LHe3 = materialMgr->FindOrBuildMaterial("LHelium-3");
  G4Material* CH2 = materialMgr->FindOrBuildMaterial("CH2");
  G4Material* Be = materialMgr->FindOrBuildMaterial("G4_Be");
  G4Material* Kapton = materialMgr->FindOrBuildMaterial("G4_KAPTON");
  G4Material* Cu = materialMgr->FindOrBuildMaterial("G4_Cu");
  //G4Material* D2Gas = materialMgr->FindOrBuildMaterial("DeuteriumGas");
  //G4Material* SUS = materialMgr->FindOrBuildMaterial("SUS");
  //G4Material* Ar = materialMgr-> GetMaterial("ArgonGas");


  G4Material* ChamberGas=Vacuum;
  // if (ChamberGasSelection =="HeEthan50-50")
  //   ChamberGas = HeEthan_50_50;
  // else if (ChamberGasSelection =="P10")
  //     ChamberGas = ArCH4_90_10;
  // else if (ChamberGasSelection=="ArEthan50-50")
  //   ChamberGas = ArEthan_50_50;
  // else if (ChamberGasSelection=="HeIsobutane80-20")
  //   ChamberGas = HeIsobutane_80_20;
  // else if (ChamberGasSelection=="ArEthan80-20")
  //   ChamberGas = ArEthan_80_20;
  // else if (ChamberGasSelection=="ArCO2Methan89-10-1")
  //   ChamberGas = ArCO2Methan_89_10_1;
  // else if (ChamberGasSelection=="Vacuum")
  //   ChamberGas = Vacuum;
  // else
  //   {
  //     G4cout << "!!! Selected Chamber Gas " << ChamberGasSelection.Data() << " is not supported yet. " << G4endl; 
  //   }

  G4cout << "!!! " << ChamberGasSelection << " was selected as chamber gas for CDC " << G4endl;
  

  constexpr double twopi = 6.28318530717959;


  bool DoKurama = true;
  bool DoTOFn = true;
  bool DoTOFp = true;
  bool DoCDS = true;
  bool DoCDH = true;
  bool DoTargetChamber = true;
  bool DoAC = true;


  
  // ==============================================================
  // experimental hall (world volume)
  //   --- beam line along z axis ---
  // ==============================================================

  G4double expHall_x =10.0*m;
  G4double expHall_y = 2.0*m;
  G4double expHall_z =20.0*m;
  experimentalHall_box  = new G4Box("expHall_box",expHall_x,expHall_y,expHall_z);
  experimentalHall_log  = new G4LogicalVolume(experimentalHall_box,Air,"expHall_log",0,0,0);
  experimentalHall_phys = new G4PVPlacement(0,G4ThreeVector(),experimentalHall_log,"expHall",0,false,0);

  experimentalHall_log->SetVisAttributes(G4VisAttributes::Invisible);

#ifdef NOYET  
  G4Region* defaultRegion = (*(G4RegionStore::GetInstance()))[0];
  KnuclRegionInformation* defaultRInfo = new KnuclRegionInformation();
  defaultRInfo->SetWorld();
  defaultRegion->SetUserInformation(defaultRInfo);
#endif

  // ==============================================================
  // Kurama
  // ==============================================================

  //**************//
  //*** Kurama ***//
  //**************//
  if(DoKurama == true)
    {
      
      G4double kurama_x = 1.4/2.*m;
      G4double kurama_y = 1.4/2.*m;
      G4double kurama_z = 1.4/2.*m;

      G4double kuramaAperture_x = 0.8/2.*m;
      G4double kuramaAperture_y = 0.4/2.*m;
      G4double kuramaAperture_z = 0.65/2.*m;
      
      G4double kuramaPos_x = 0.0*m;
      G4double kuramaPos_y = 0.0*m;
      G4double kuramaPos_z = 2.8*m;
      
      Kurama_box  = new G4Box("Kurama_box",kurama_x, kurama_y,kurama_z);
      Kurama_log  = new G4LogicalVolume(Kurama_box,Vacuum,"Kurama_log",0,0,0);
      Kurama_phys = new G4PVPlacement(0,G4ThreeVector(kuramaPos_x,kuramaPos_y,kuramaPos_z),Kurama_log,"Kurama",experimentalHall_log,false,0);
      
      //************************//
      //*** Kurama Aperture ****//
      //************************//
      G4VSolid* KuramaAperture_box = new G4Box("KuramaAp_box",kuramaAperture_x,kuramaAperture_y,kuramaAperture_z  );
      
      //G4FieldManager* KuramaFieldMgr = fEmFieldSetupKURAMA->GetFieldManager(); 
      
      //  KuramaAperture_log = new G4LogicalVolume(KuramaAperture_box, HeGas, 
      KuramaAperture_log = new G4LogicalVolume(KuramaAperture_box, Air,"KuramaAperture_log",0,0,0);//, KuramaFieldMgr, 0, 0);
      G4PVPlacement* KuramaAperture_phys= new G4PVPlacement(0, G4ThreeVector(0.0,0.0,0.0),KuramaAperture_log,"KuramaAperture",Kurama_log,false,0);
      
      //*******************//
      //*** Kurama Yoke ***//
      //*******************//
      //    Kurama Yoke Up and Down side
      G4double kuramaYokeUD_x = 1.4/2.*m;
      G4double kuramaYokeUD_y = 0.5/2.*m;
      G4double kuramaYokeUD_z = 1.4/2.*m;
      G4VSolid* KuramaYokeUD_Box = new G4Box("KuramaYokeUD_box",kuramaYokeUD_x,kuramaYokeUD_y,kuramaYokeUD_z );
      //    Kurama Yoke Left and Right side
      G4double kuramaYokeLR_x = 0.3/2.*m;
      G4double kuramaYokeLR_y = 0.4/2.*m;
      G4double kuramaYokeLR_z = 1.4/2.*m;
      G4VSolid* KuramaYokeLR_Box = new G4Box("KuramaYokeLR_box",kuramaYokeLR_x,kuramaYokeLR_y,kuramaYokeLR_z );
      G4LogicalVolume* KuramaYoke_log[4];
      G4PVPlacement* KuramaYoke_phys[4];
      KuramaYoke_log[0] = new G4LogicalVolume(KuramaYokeUD_Box, Fe, "KuramaYoke0_log",0, 0, 0 );
      KuramaYoke_phys[0]= new G4PVPlacement(0, G4ThreeVector(0.0, 45.0*cm,0.0),KuramaYoke_log[0],"KuramaYoke0",Kurama_log,false,0);
      
      KuramaYoke_log[1] = new G4LogicalVolume(KuramaYokeUD_Box, Fe, "KuramaYoke1_log",0, 0, 0 );
      KuramaYoke_phys[1]= new G4PVPlacement(0, G4ThreeVector(0.0,-45.0*cm,0.0),KuramaYoke_log[1],"KuramaYoke1",Kurama_log,false,0);
      
      KuramaYoke_log[2] = new G4LogicalVolume(KuramaYokeLR_Box, Fe, "KuramaYoke2_log",0, 0, 0 );
      KuramaYoke_phys[2]= new G4PVPlacement(0, G4ThreeVector(55.0*cm,0.0,0.0),KuramaYoke_log[2],"KuramaYoke2",Kurama_log,false,0);
      
      KuramaYoke_log[3] = new G4LogicalVolume(KuramaYokeLR_Box, Fe, "KuramaYoke3_log",0, 0, 0 );
      KuramaYoke_phys[3]= new G4PVPlacement(0, G4ThreeVector(-55.0*cm,0.0,0.0),KuramaYoke_log[3],"KuramaYoke3",Kurama_log,false,0);
      
      //--- Visualization ---//
      Kurama_log->SetVisAttributes(G4VisAttributes::Invisible);
      KuramaAperture_log->SetVisAttributes(G4VisAttributes::Invisible);
      G4VisAttributes* KuramaYoke_att = new G4VisAttributes(Gray);
      KuramaYoke_att->SetForceWireframe(true);
      for (auto& KuramaLogV : KuramaYoke_log)
	KuramaLogV->SetVisAttributes(KuramaYoke_att);

            // ==============================================================
      // Drift Chamber for Proton (KURAMA-UP) by fujioka
      // ==============================================================
      const G4double SDC1_SDC1_space = 5.0*cm;

      //***********//
      //*** SDC1 ***//
      //***********//
      G4double SDC1_block_x = 24./2.*cm;
      G4double SDC1_block_y = 24./2.*cm;
      G4double SDC1_block_z = 5.0/2.*cm;
      G4VSolid* SDC1_box = new G4Box("SDC1_box",SDC1_block_x, SDC1_block_y,SDC1_block_z);
      G4LogicalVolume* SDC1_log[2] = {0,0};
      char log_name[100];
      for (G4int i=0; i<2; i++)
	{
	  sprintf(log_name, "SDC1_log%d", i+1);
	  SDC1_log[i] = new G4LogicalVolume(SDC1_box,Vacuum,log_name,0,0,0);
	}
      G4double SDC1_pos_x = 0.0*m;
      G4double SDC1_pos_y = 0.0*m;
      G4double SDC1_pos_z_base = kuramaPos_z-kurama_z-5.0*cm-SDC1_block_z;

      G4double SDC1_pos_z[2] = {SDC1_pos_z_base-SDC1_block_z*2-SDC1_SDC1_space,
				SDC1_pos_z_base};

      G4PVPlacement* SDC1_phys[2];
      char phys_name[100];
      for (G4int i=0; i<2; i++)
	{
	  sprintf(phys_name, "SDC1_phys%d", i);
	  SDC1_phys[i] = new G4PVPlacement(0, G4ThreeVector(SDC1_pos_x,SDC1_pos_y,SDC1_pos_z[i]),SDC1_log[i], phys_name, experimentalHall_log, false, 0);
	}
      const G4double SDC1_cell_size = 8.0*mm; // full width
      G4double SDC1_seg_x = SDC1_cell_size/2;
      G4double SDC1_seg_y = SDC1_block_y;
      G4double SDC1_seg_z = 0.000001/2*mm;
      G4VSolid* sSDC1_box = new G4Box("sSDC1_box",SDC1_seg_x,SDC1_seg_y,SDC1_seg_z);
      G4LogicalVolume* sSDC1_log = new G4LogicalVolume(sSDC1_box, ChamberGas, "sSDC1_log", 0,0,0);
      G4double SDC1_layer_pos[4]={-1.5*cm,-0.5*cm,0.5*cm,1.5*cm};
      G4int SDC1_layer_n[4]={30,29,30,29};
      G4double SDC1_layer_rot[4]={0*deg,0*deg,90*deg,90*deg};
      G4double SDC1_layer_sta[4]={-SDC1_block_x+SDC1_cell_size/2,-SDC1_block_x+SDC1_cell_size,
				  -SDC1_block_x+SDC1_cell_size/2,-SDC1_block_x+SDC1_cell_size};

    
      std::vector<std::vector< std::vector<G4PVPlacement*> > > sSDC1_phys(2, std::vector< std::vector<G4PVPlacement*> >(4) );
      for (G4int i=0; i<2; i++)
	{
	  for (G4int j=0; j<4; j++)
	    {
	      G4int n=0;
	      sSDC1_phys[i][j].resize(SDC1_layer_n[j],nullptr);
	      for (G4int k=0; k<SDC1_layer_n[j]; k++)
		{
		  sprintf(phys_name, "sSDC1_phys%d%d%02d", i+1, j+1, k+1);
		  G4ThreeVector xyzCounter;
		  if (j<2)
		    xyzCounter=G4ThreeVector(SDC1_layer_sta[j]+SDC1_cell_size*(k),0,SDC1_layer_pos[j]);
		  else
		    xyzCounter=G4ThreeVector(0,SDC1_layer_sta[j]+SDC1_cell_size*(k),SDC1_layer_pos[j]);
		  G4RotationMatrix* rotCounterSDC = new G4RotationMatrix;
		  rotCounterSDC->rotateZ(SDC1_layer_rot[j]);
		  G4Transform3D posCounterSDC(*rotCounterSDC, xyzCounter);
		  sSDC1_phys[i][j][k] = new G4PVPlacement(posCounterSDC, sSDC1_log, phys_name, SDC1_log[i], false, n);
		  n++;
		}
	    }
	}

	
      //--- Visualization ---//
      G4VisAttributes *SDC1_att = new G4VisAttributes(Green);
      SDC1_att->SetForceWireframe(true);
      for (G4int i=0; i<2; i++)
	{
	  //BLC_log[i-1]->SetVisAttributes(BLC_att);
	  SDC1_log[i]->SetVisAttributes(G4VisAttributes::Invisible);
	}
      G4VisAttributes *sSDC1_att = new G4VisAttributes(Yellow);
      sSDC1_att->SetForceWireframe(true);
      sSDC1_log->SetVisAttributes(sSDC1_att);      
  
      //***********//
      //*** SDC2 ***//
      //***********//
      G4double SDC2_block_x = 80./2.*cm;
      G4double SDC2_block_y = 80./2.*cm;
      G4double SDC2_block_z = 5.0/2.*cm;
      G4VSolid* SDC2_box = new G4Box("SDC2_box",SDC2_block_x, SDC2_block_y,SDC2_block_z);
      G4LogicalVolume* SDC2_log[1] = {0}; 
      for (G4int i=0; i<1; i++)
	{
	  sprintf(log_name, "SDC2_log%d", i);
	  SDC2_log[i] = new G4LogicalVolume(SDC2_box,Vacuum,log_name,0,0,0);
	}
      G4double SDC2_pos_x = 0.0*m;
      G4double SDC2_pos_y = 0.0*m;
      G4double SDC2_pos_z_base = kuramaPos_z+kurama_z+5.0*cm+SDC2_block_z;

      G4double SDC2_pos_z[1] = {SDC2_pos_z_base};
      G4PVPlacement* SDC2_phys[2];

      for (G4int i=0; i<1; i++)
	{
	  sprintf(phys_name, "SDC2_phys%d", i);
	  SDC2_phys[i] = new G4PVPlacement(0, G4ThreeVector(SDC2_pos_x,SDC2_pos_y,SDC2_pos_z[i]),
					   SDC2_log[i], phys_name, experimentalHall_log, false, 0);
	}
      const G4double SDC2_cell_size = 20.0*mm; // full width
      G4double SDC2_seg_x = SDC2_cell_size/2;
      G4double SDC2_seg_y = SDC2_block_y;
      G4double SDC2_seg_z = 0.000001/2*mm;
      G4VSolid* sSDC2_box = new G4Box("sSDC2_box",SDC2_seg_x,SDC2_seg_y,SDC2_seg_z);
      G4LogicalVolume* sSDC2_log = new G4LogicalVolume(sSDC2_box, ChamberGas, "sSDC2_log", 0,0,0);
      G4double SDC2_layer_pos[4]={-1.5*cm,-0.5*cm,0.5*cm,1.5*cm};
      G4int SDC2_layer_n[4]={40,39,40,39};
      G4double SDC2_layer_rot[4]={0*deg,0*deg,90*deg,90*deg};
      G4double SDC2_layer_sta[4]={-SDC2_block_x+SDC2_cell_size/2,-SDC2_block_x+SDC2_cell_size,
				  -SDC2_block_x+SDC2_cell_size/2,-SDC2_block_x+SDC2_cell_size};

      std::vector< std::vector< std::vector<G4PVPlacement*> > > sSDC2_phys(1,std::vector<std::vector<G4PVPlacement*> >(4) );
      for (G4int i=0; i<1; i++)
	for (G4int j=0; j<4; j++)
	  {
	    G4int n=0;
	    sSDC2_phys[i][j].resize(SDC2_layer_n[j],nullptr);
	    for (G4int k=0; k<SDC2_layer_n[j]; k++)
	      {
		sprintf(phys_name, "sSDC2_phys%d%d%02d", i+1, j+1, k+1);
		G4ThreeVector xyzCounter;
		if (j<2)
		  xyzCounter=G4ThreeVector(SDC2_layer_sta[j]+SDC2_cell_size*(k),0,SDC2_layer_pos[j]);
		else
		  xyzCounter=G4ThreeVector(0,SDC2_layer_sta[j]+SDC2_cell_size*(k),SDC2_layer_pos[j]);
		G4RotationMatrix* rotCounterSDC = new G4RotationMatrix;
		rotCounterSDC->rotateZ(SDC2_layer_rot[j]);
		G4Transform3D posCounterSDC(*rotCounterSDC, xyzCounter);
		sSDC2_phys[i][j][k] = new G4PVPlacement(posCounterSDC,sSDC2_log, phys_name, SDC2_log[i], false, n);
		n++;
	      }
	  }



      //--- Visualization ---//
      G4VisAttributes *SDC2_att = new G4VisAttributes(Green);
      SDC2_att->SetForceWireframe(true);
      for (G4int i=0; i<1; i++)
	{
	  //BLC_log[i]->SetVisAttributes(BLC_att);
	  SDC2_log[i]->SetVisAttributes(G4VisAttributes::Invisible);
	}
      G4VisAttributes *sSDC2_att = new G4VisAttributes(Yellow);
      sSDC2_att->SetForceWireframe(true);
      sSDC2_log->SetVisAttributes(sSDC2_att);      



    }
  //******************************//
  //*** VETO counter in Kurama ***//
  //******************************//
//  G4double VetoScintiUD_x = 1.00/2.0*m;
//  G4double VetoScintiUD_y = 0.01/2.0*m;
//  G4double VetoScintiUD_z = 0.80/2.0*m;
//  VetoScintiUD_box = new G4Box("VetoScintiUD_box",VetoScintiUD_x,
//                                                  VetoScintiUD_y,
//                                                  VetoScintiUD_z );
//
//  G4double VetoScintiLR_x = 0.01/2.0*m;
//  G4double VetoScintiLR_y = 0.48/2.0*m;
//  G4double VetoScintiLR_z = 0.80/2.0*m;
//  VetoScintiLR_box = new G4Box("VetoScintiLR_box",VetoScintiLR_x,
//                                                  VetoScintiLR_y,
//                                                  VetoScintiLR_z );
//
//  VetoScinti_log[0] =  new G4LogicalVolume(VetoScintiUD_box,Scinti,"VetoScinti0_log",0,0,0);
//  VetoScinti_log[1] =  new G4LogicalVolume(VetoScintiUD_box,Scinti,"VetoScinti1_log",0,0,0);
//  VetoScinti_log[2] =  new G4LogicalVolume(VetoScintiLR_box,Scinti,"VetoScinti2_log",0,0,0);
//  VetoScinti_log[3] =  new G4LogicalVolume(VetoScintiLR_box,Scinti,"VetoScinti3_log",0,0,0);
//  
//  VetoScinti_phys[0] = new G4PVPlacement(0, G4ThreeVector(  0.0,  24.5*cm,0.0),
//                           VetoScinti_log[0],"VetoScinti0",Kurama_log,false,0);
//  VetoScinti_phys[1] = new G4PVPlacement(0, G4ThreeVector(  0.0, -24.5*cm,0.0),
//                           VetoScinti_log[1],"VetoScinti1",Kurama_log,false,0);
//  VetoScinti_phys[2] = new G4PVPlacement(0, G4ThreeVector( 49.5*cm,0.0,0.0),
//                           VetoScinti_log[2],"VetoScinti2",Kurama_log,false,0);
//  VetoScinti_phys[3] = new G4PVPlacement(0, G4ThreeVector(-49.5*cm,0.0*cm,0.0),
//                           VetoScinti_log[3],"VetoScinti3",Kurama_log,false,0);
//
//  VetoScinti_log[0]->SetSensitiveDetector(counterSD);
//  VetoScinti_log[1]->SetSensitiveDetector(counterSD);
//  VetoScinti_log[2]->SetSensitiveDetector(counterSD);
//  VetoScinti_log[3]->SetSensitiveDetector(counterSD);
//
//  //--- Visualization ---//
//  G4VisAttributes *VetoScinti_att = new G4VisAttributes(Red);
//  VetoScinti_att->SetForceWireframe(true);
//  for (G4int i=0; i<=4; i++){
//    VetoScinti_log[i]->SetVisAttributes(VetoScinti_att);
//  }


  // ==============================================================
  // Neutron Counter
  // ==============================================================
  if(DoTOFn==true)
    {  
      //***********//
      //*** TOF ***//
      //***********//
      G4double TOF_block_x = 3.2/2.*m;
      G4double TOF_block_y = 1.5/2.*m;
      G4double TOF_block_z = 40.0/2.*cm;
      TOF_box = new G4Box("TOF_box",TOF_block_x, TOF_block_y,TOF_block_z);
      TOF_log = new G4LogicalVolume(TOF_box,Vacuum,"TOF_log",0,0,0);
      G4double TOF_blockPos_x = 0.0*m;
      G4double TOF_blockPos_y = 0.0*m;
      // to be modified
      G4double TOF_blockPos_z = TofRefPos*m+TOF_block_z+0.5*cm+1.*m;
      G4cout << " ToF length for NC        = " << TOF_blockPos_z << G4endl; 
      TOF_phys = new G4PVPlacement(0,G4ThreeVector(TOF_blockPos_x,TOF_blockPos_y,TOF_blockPos_z),TOF_log,"TOF",experimentalHall_log,false,0);

      //**********************************************//
      //*** sTOF ( counter element in TOF volume ) ***//
      //**********************************************//
      G4double TOF_seg_x = 20.0/2.*cm;
      G4double TOF_seg_y = 1.5/2.*m;
      G4double TOF_seg_z = 5.0/2.*cm;
      G4VSolid* sTOF_box = new G4Box("sTOF_box",TOF_seg_x,TOF_seg_y,TOF_seg_z);

      char log_name[100], phys_name[100];
      //std::string log_name,phys_name;
  
      G4LogicalVolume* sTOF_log = new G4LogicalVolume(sTOF_box, Scinti, log_name, 0,0,0);
      G4PVPlacement* sTOF_phys [7][16];
      G4int n=0;
      for (G4int i=0; i<7; ++i)
	{
	  n=0;
	  for (G4int j=0; j<16; j++)
	    {
	      sprintf(log_name, "sTOF_log%d%02d", i+1, j+1);
	      sprintf(phys_name, "sTOF_phys%d%02d", i+1, j+1);
	      G4ThreeVector sTOF_Pos((double)(j+1.-8.5)*TOF_seg_x*2.0, 0.0,(double)(i+1.-4.0)*TOF_seg_z*2.0);
	      sTOF_phys[i][j] = new G4PVPlacement(0, sTOF_Pos, sTOF_log, phys_name, TOF_log, false, n);
	      n++;
	    }
	}

      //sTOF_log->SetSensitiveDetector(counterSD);

      //--- Visualization ---//
      TOF_log->SetVisAttributes(G4VisAttributes::Invisible);
      G4VisAttributes *sTOF_att = new G4VisAttributes(Red);
      sTOF_att->SetForceWireframe(true);
      sTOF_log->SetVisAttributes(sTOF_att);

    }
  // ==============================================================
  // Proton Counter (central)
  // ==============================================================
  if(DoTOFp==true)
    {
      //***********//
      //*** TOF ***//
      //***********//
      G4double pTOF_block_x = 3.433/2.*m;
      G4double pTOF_block_y = 1.5/2.*m;
      G4double pTOF_block_z = 3.0/2.*cm;
      G4VSolid* pTOF_box = new G4Box("pTOF_box",pTOF_block_x, pTOF_block_y,pTOF_block_z);
      G4LogicalVolume* pTOF_log = new G4LogicalVolume(pTOF_box,Vacuum,"pTOF_log",0,0,0);
      //  G4double pTOF_blockPos_x = -1.6*m-pTOF_block_x;
      G4double pTOF_blockPos_x = 0.0*m;
      G4double pTOF_blockPos_y = 0.0*m;
      G4double pTOF_blockPos_z = TofRefPos*m+pTOF_block_z+0.5*cm-1.0*m;
      G4cout << " ToF length for PB        = " << pTOF_blockPos_z << G4endl; 
      G4PVPlacement* pTOF_phys = new G4PVPlacement(0,G4ThreeVector(pTOF_blockPos_x,pTOF_blockPos_y,pTOF_blockPos_z), pTOF_log,"pTOF",experimentalHall_log,false,0);

      //**********************************************//
      //*** sTOF ( counter element in TOF volume ) ***//
      //**********************************************//
      G4double pTOF_seg_x = 0.1/2.*m;
      G4double pTOF_seg_y = 1.5/2.*m;
      G4double pTOF_seg_z = 3.0/2.*cm;
      G4VSolid* spTOF_box = new G4Box("spTOF_box",pTOF_seg_x,pTOF_seg_y,pTOF_seg_z);

      char log_name[100], phys_name[100];
      
      G4LogicalVolume* spTOF_log = new G4LogicalVolume(spTOF_box, Scinti,log_name, 0,0,0);
      G4PVPlacement* spTOF_phys [1][61];
      G4int n=0;
      for (G4int i=0; i<1; i++)
	{
	  n=0;
	  for (G4int j=0; j<34; j++)
	    {
	      sprintf(log_name, "sptof_log%d%02d", i+1, j+1);
	      sprintf(phys_name, "sptof_phys%d%02d", i+1, j+1);
	      G4ThreeVector spTOF_Pos(-(double)(j+1-17.5)*(pTOF_seg_x+0.5*mm)*2.0, 0.0,(double)(i+1.-1.0)*pTOF_seg_z*2.0);
	      spTOF_phys[i][j] = new G4PVPlacement(0, spTOF_Pos, spTOF_log, phys_name, pTOF_log, false, n);
	      n++;
	    }
	}


      // ==============================================================
      // Proton Counter (side)
      // ==============================================================

      //***********//
      //*** TOF ***//
      //***********//
      G4double pTOF2_block_x = 2.726/2.*m;
      G4double pTOF2_block_y = 1.5/2.*m;
      //  G4double pTOF2_block_z = 3.0/2.*cm;
      G4double pTOF2_block_z = 4.316*cm;
      G4VSolid* pTOF2_box = new G4Box("pTOF2_box",pTOF2_block_x, pTOF2_block_y,pTOF2_block_z);
      G4LogicalVolume* pTOF2_log = new G4LogicalVolume(pTOF2_box,Vacuum,"pTOF2_log",0,0,0);
      //  G4double pTOF_blockPos_x = -1.6*m-pTOF_block_x;

      G4double pTOF2_tiltangle = 16.0*deg;
      G4double pTOF2_blockPos_x = -pTOF_block_x-pTOF2_block_x*cos(pTOF2_tiltangle);
      G4double pTOF2_blockPos_y = 0.0*m;
      G4double pTOF2_blockPos_z = TofRefPos*m+pTOF_block_z+0.5*cm-1.0*m-pTOF2_block_x*sin(pTOF2_tiltangle);

      G4ThreeVector pTOF2Counter(pTOF2_blockPos_x,pTOF2_blockPos_y,pTOF2_blockPos_z);
      G4RotationMatrix* rotCounterTof2 = new G4RotationMatrix;
      rotCounterTof2->rotateY(-pTOF2_tiltangle);
      G4Transform3D posCounterTof2(*rotCounterTof2, pTOF2Counter);
      G4PVPlacement* pTOF2_phys = new G4PVPlacement(posCounterTof2, pTOF2_log, "pTOF2", experimentalHall_log, false, 0);  


      //   pTOF2_phys = new G4PVPlacement(0,
      //              G4ThreeVector(pTOF2_blockPos_x,pTOF2_blockPos_y,pTOF2_blockPos_z),
      //              pTOF2_log,"pTOF2",experimentalHall_log,false,0);

      //**********************************************//
      //*** sTOF ( counter element in TOF volume ) ***//
      //**********************************************//

      for (G4int i=0; i<1; i++)
	{
	  n=0;
	  for (G4int j=34; j<61; j++)
	    {
	      sprintf(log_name, "sptof_log%d%02d", i+1, j+1);
	      sprintf(phys_name, "sptof_phys%d%02d", i+1, j+1);
	      G4ThreeVector spTOF_Pos(-(double)(j+1-48.0)*(pTOF_seg_x+0.5*mm)*2.0, 0.0,-pTOF2_block_z+pTOF_seg_z);
	      spTOF_phys[i][j] = new G4PVPlacement(0, spTOF_Pos, spTOF_log, phys_name, pTOF2_log, false, n);
	      n++;
	    }
	}
  
      //spTOF_log->SetSensitiveDetector(counterSD);
  
      //--- Visualization ---//
      pTOF_log->SetVisAttributes(G4VisAttributes::Invisible);
      G4VisAttributes *spTOF_att = new G4VisAttributes(Red);
      spTOF_att->SetForceWireframe(true);
      spTOF_log->SetVisAttributes(spTOF_att);
    }
  // ==============================================================
  // CDS
  // ==============================================================

  const G4double cds_rmax = CDS_RMAX*mm;
  const G4double cds_z    = CDS_Z/2.0*mm;
  
  const G4double cdsPos_x = 0.0*m;
  const G4double cdsPos_y = 0.0*m;
  const G4double cdsPos_z = 0.0*m;
      
  if(DoCDS == true)
    {
      //***********//
      //*** CDS ***//
      //***********//

      //G4FieldManager* CDCFieldMgr = fEmFieldSetupCDC->GetFieldManager();
 
      CDS_tube = new G4Tubs("CDS_tube", 0.0, cds_rmax, cds_z, 0.0, twopi);
      CDS_log  = new G4LogicalVolume(CDS_tube, Air, "CDS_log",0,0,0);// CDCFieldMgr,0,0);
      CDS_phys = new G4PVPlacement(0, G4ThreeVector(cdsPos_x, cdsPos_y, cdsPos_z), CDS_log, "CDS", experimentalHall_log, false,0);
      //--- Visualization ---//
      CDS_log->SetVisAttributes(G4VisAttributes::Invisible);
  
      //****************//
      //*** CDS Yoke ***//
      //****************//
      G4double cdsyoke_rmax = CDSYOKE_RMAX*mm;
      G4double cdsyoke_rmin = CDSYOKE_RMIN*mm;

      G4Tubs* CDSYoke_tube= new G4Tubs("CDS_tube", cdsyoke_rmin, cdsyoke_rmax, cds_z, 0.0, twopi);
      G4LogicalVolume* CDSYoke_log = new G4LogicalVolume(CDSYoke_tube, Fe, "CDSYoke_log", 0,0,0);

      G4PVPlacement* CDSYoke_phys = new G4PVPlacement(0, G4ThreeVector(cdsPos_x, cdsPos_y, cdsPos_z), CDSYoke_log, "CDS_Yoke", experimentalHall_log, false,0);

      //--- Visualization ---//
      G4VisAttributes *CDSYoke_att = new G4VisAttributes(Gray);
      CDSYoke_att->SetForceWireframe(true);
      CDSYoke_log->SetVisAttributes(CDSYoke_att);

      //**********************//
      //*** CDS endcap ***//
      //**********************//
      G4double cds_endcap_rmax = CDSYOKE_RMAX*mm;
      G4double cds_endcap_rmin = CDS_ENDCAP_RMIN*mm;
      G4double cds_endcap_z    = CDS_ENDCAP_Z/2.0*mm;

      G4VSolid* CDS_endcap_tube = new G4Tubs("CDS_endcap_tube", cds_endcap_rmin, cds_endcap_rmax, cds_endcap_z, 0.0, twopi);
      CDS_endcap_log  = new G4LogicalVolume(CDS_endcap_tube, Fe, "CDS_endcap_log",0,0,0);//, CDCFieldMgr,0,0);

      G4PVPlacement* CDS_endcap_phys[2];
      CDS_endcap_phys[0] = new G4PVPlacement(0, G4ThreeVector(cdsPos_x, cdsPos_y, cdsPos_z-(cds_z+cds_endcap_z)),
					     CDS_endcap_log, "CDS_endcap_up", experimentalHall_log, false,0);
      CDS_endcap_phys[1] = new G4PVPlacement(0,G4ThreeVector(cdsPos_x, cdsPos_y, cdsPos_z+(cds_z+cds_endcap_z)),
					     CDS_endcap_log, "CDS_endcap_down", experimentalHall_log, false,0);

      //--- Visualization ---//
      G4VisAttributes *CDS_endcap_att = new G4VisAttributes(Gray);
      CDS_endcap_att->SetForceWireframe(true);
      CDS_endcap_log->SetVisAttributes(CDS_endcap_att);

      //***********//
      //*** CDC ***//
      //***********//
      G4double cdc_rmin  = CDC_RMIN*mm;
      G4double cdc_rmax  = CDC_RMAX*mm;
      G4double cdc_z     = CDSLength/2.0*mm; 

      G4VSolid* CDC_body_tube = new G4Tubs("CDC_body_tube", cdc_rmin, cdc_rmax, cdc_z, 0.0, twopi);
      G4LogicalVolume* CDC_body_log  = new G4LogicalVolume(CDC_body_tube,ChamberGas,"CDC_body_log",0,0,0);
      G4PVPlacement* CDC_body_phys = new G4PVPlacement(0, G4ThreeVector(cdsPos_x, cdsPos_y, cdsPos_z), CDC_body_log, "CDC_body_phys", CDS_log, false, 0);  

      const G4double cdc_wire_dist = 0.45*cm;
      const G4double cdc_off = cdc_wire_dist;

      G4double rmin = 0.0;
      G4double rmax = 0.0;
      char name_sol[64];
      char name_log[64];
      char name_phy[64];
      //#if 1
      G4double cdc_cell_twist[15];
      G4double cdc_radius[15];
      for (G4int i=0; i<15; i++)
	{
	  cdc_cell_twist[i]  = twopi/(double)N_CDC_CELL[i]*CDC_CELL_OFFSET[i];
	  cdc_radius[i] = CDC_RADIUS[i]*mm;
	}

      G4VSolid* CDC_tube[15];
      G4VSolid* CDC_twist_tube[15];
      G4LogicalVolume* CDC_log[15];
      std::vector< std::vector<G4PVPlacement*> > CDC_phys(15);//[3][j];
      for(size_t i=0;i< CDC_phys.size();++i)
	CDC_phys[i].resize(N_CDC_CELL[i],nullptr);

      G4int n=0;

      //X, X', X
      for (G4int i=0; i<3; i++)
	{
	  rmin = cdc_radius[i]-cdc_off;
	  rmax = cdc_radius[i]+cdc_off;
      
	  sprintf(name_sol,"CDC_sol_%2.2d",i);
	  sprintf(name_log,"CDC_log_%2.2d",i);
      
	  CDC_tube[i] = new G4Tubs(name_sol, rmin, rmax, cdc_z, -0.5*twopi/(double)N_CDC_CELL[i], twopi/(double)N_CDC_CELL[i]);
	  CDC_log [i] = new G4LogicalVolume(CDC_tube[i],ChamberGas,name_log,0,0,0);
	  n=0;
	  for (G4int j=0; j<N_CDC_CELL[i]; j++)
	    {
	      sprintf(name_phy,"CDC_phys%2.2d%03d",i+1,j+1);

	      G4ThreeVector xyzCounter(cdsPos_x, cdsPos_y, cdsPos_z);
	      G4RotationMatrix* rotCounterCDCcell = new G4RotationMatrix;
	      rotCounterCDCcell->rotateZ(twopi/(double)N_CDC_CELL[i]*(j+0.5*CDC_CELL_XXP[i]));
	      G4Transform3D posCounterCDCcell(*rotCounterCDCcell, xyzCounter);

	      CDC_phys[i][j] = new G4PVPlacement(posCounterCDCcell, CDC_log[i], name_phy, CDC_body_log, false, n);  
	      n++;
	    }
	}

      //U, U'
      for (G4int i=0; i<2; i++)
	{
	  rmin = cdc_radius[i+3]-cdc_off;
	  rmax = cdc_radius[i+3]+cdc_off;
      
	  sprintf(name_sol,"CDC_sol_%2.2d",i+3);
	  sprintf(name_log,"CDC_log_%2.2d",i+3);

	  CDC_twist_tube[i+3] = new G4TwistedTubs(name_sol, cdc_cell_twist[i+3], rmin, rmax,cdc_z, twopi/(double)N_CDC_CELL[i+3]);
	  CDC_log [i+3] = new G4LogicalVolume(CDC_twist_tube[i+3],ChamberGas,name_log,0,0,0);
	  n=0;
	  for (G4int j=0; j<N_CDC_CELL[i+3]; j++)
	    {
	      sprintf(name_phy,"CDC_phys%2.2d%03d",i+4,j+1);

	      G4ThreeVector xyzCounter(cdsPos_x, cdsPos_y, cdsPos_z);
	      G4RotationMatrix* rotCounterCDCcell = new G4RotationMatrix;
	      rotCounterCDCcell->rotateZ(twopi/(double)N_CDC_CELL[i+3]*(j+0.5*CDC_CELL_XXP[i+3]));
	      G4Transform3D posCounterCDCcell(*rotCounterCDCcell, xyzCounter);

	      CDC_phys[i+3][j] = new G4PVPlacement(posCounterCDCcell, CDC_log[i+3], name_phy, CDC_body_log, false, n);
	      n++;
	    }
	}

      //V, V'
      for (G4int i=0; i<2; i++)
	{
	  rmin = cdc_radius[i+5]-cdc_off;
	  rmax = cdc_radius[i+5]+cdc_off;
      
	  sprintf(name_sol,"CDC_sol_%2.2d",i+5);
	  sprintf(name_log,"CDC_log_%2.2d",i+5);
      
	  CDC_twist_tube[i+5] = new G4TwistedTubs(name_sol, cdc_cell_twist[i+5], rmin, rmax, cdc_z, twopi/(double)N_CDC_CELL[i+5]);
	  CDC_log [i+5] = new G4LogicalVolume(CDC_twist_tube[i+5],ChamberGas,name_log,0,0,0);
	  n=0;
	  for (G4int j=0; j<N_CDC_CELL[i+5]; j++)
	    {
	      sprintf(name_phy,"CDC_phys%2.2d%03d",i+6,j+1);
	  
	      G4ThreeVector xyzCounter(cdsPos_x, cdsPos_y, cdsPos_z);
	      G4RotationMatrix* rotCounterCDCcell = new G4RotationMatrix;
	      rotCounterCDCcell->rotateZ(twopi/(double)N_CDC_CELL[i+5]*(j+0.5*CDC_CELL_XXP[i+5]));
	      G4Transform3D posCounterCDCcell(*rotCounterCDCcell, xyzCounter);

	      CDC_phys[i+5][j] = new G4PVPlacement(posCounterCDCcell,CDC_log[i+5], name_phy, CDC_body_log, false, n);
	      n++;
	    }
	}

      //X, X'
      for (G4int i=0; i<2; i++)
	{
	  rmin = cdc_radius[i+7]-cdc_off;
	  rmax = cdc_radius[i+7]+cdc_off;
      
	  sprintf(name_sol,"CDC_sol_%2.2d",i+7);
	  sprintf(name_log,"CDC_log_%2.2d",i+7);
      
	  CDC_tube[i+7] = new G4Tubs(name_sol, rmin, rmax, cdc_z, -0.5*twopi/(double)N_CDC_CELL[i+7], twopi/(double)N_CDC_CELL[i+7]);
	  CDC_log [i+7] = new G4LogicalVolume(CDC_tube[i+7],ChamberGas,name_log,0,0,0);
	  n=0;
	  for (G4int j=0; j<N_CDC_CELL[i+7]; j++)
	    {
	      sprintf(name_phy,"CDC_phys%2.2d%03d",i+8,j+1);

	      G4ThreeVector xyzCounter(cdsPos_x, cdsPos_y, cdsPos_z);
	      G4RotationMatrix* rotCounterCDCcell = new G4RotationMatrix;
	      rotCounterCDCcell->rotateZ(twopi/(double)N_CDC_CELL[i+7]*(j+0.5*CDC_CELL_XXP[i+7]));
	      G4Transform3D posCounterCDCcell(*rotCounterCDCcell, xyzCounter);

	      CDC_phys[i+7][j] = new G4PVPlacement(posCounterCDCcell, CDC_log[i+7], name_phy, CDC_body_log, false, n);
	      n++;
	    }
	}

      //U, U'
      for (G4int i=0; i<2; i++)
	{
	  rmin = cdc_radius[i+9]-cdc_off;
	  rmax = cdc_radius[i+9]+cdc_off;
      
	  sprintf(name_sol,"CDC_sol_%2.2d",i+9);
	  sprintf(name_log,"CDC_log_%2.2d",i+9);
      
	  CDC_twist_tube[i+9] = new G4TwistedTubs(name_sol, cdc_cell_twist[i+9], rmin, rmax, cdc_z, twopi/(double)N_CDC_CELL[i+9]);
	  CDC_log [i+9] = new G4LogicalVolume(CDC_twist_tube[i+9],ChamberGas,name_log,0,0,0);
	  n=0;
	  for (G4int j=0; j<N_CDC_CELL[i+9]; j++)
	    {
	      sprintf(name_phy,"CDC_phys%2.2d%03d",i+10,j+1);

	      G4ThreeVector xyzCounter(cdsPos_x, cdsPos_y, cdsPos_z);
	      G4RotationMatrix* rotCounterCDCcell = new G4RotationMatrix;
	      rotCounterCDCcell->rotateZ(twopi/(double)N_CDC_CELL[i+9]*(j+0.5*CDC_CELL_XXP[i+9]));
	      G4Transform3D posCounterCDCcell(*rotCounterCDCcell, xyzCounter);

	      CDC_phys[i+9][j] = new G4PVPlacement(posCounterCDCcell, CDC_log[i+9], name_phy, CDC_body_log, false, n);
	      n++;
	    }
	}

      //V, V'  
      for (G4int i=0; i<2; i++)
	{
	  rmin = cdc_radius[i+11]-cdc_off;
	  rmax = cdc_radius[i+11]+cdc_off;
      
	  sprintf(name_sol,"CDC_sol_%2.2d",i+11);
	  sprintf(name_log,"CDC_log_%2.2d",i+11);
      
	  CDC_twist_tube[i+11] = new G4TwistedTubs(name_sol, cdc_cell_twist[i+11], rmin, rmax,cdc_z, twopi/(double)N_CDC_CELL[i+11]);
	  CDC_log [i+11] = new G4LogicalVolume(CDC_twist_tube[i+11],ChamberGas,name_log,0,0,0);
	  n=0;
	  for (G4int j=0; j<N_CDC_CELL[i+11]; j++)
	    {
	      sprintf(name_phy,"CDC_phys%2.2d%03d",i+12,j+1);
	      G4ThreeVector xyzCounter(cdsPos_x, cdsPos_y, cdsPos_z);
	      G4RotationMatrix* rotCounterCDCcell = new G4RotationMatrix;
	      rotCounterCDCcell->rotateZ(twopi/(double)N_CDC_CELL[i+11]*(j+0.5*CDC_CELL_XXP[i+11]));
	      G4Transform3D posCounterCDCcell(*rotCounterCDCcell, xyzCounter);
	      CDC_phys[i+11][j] = new G4PVPlacement(posCounterCDCcell,CDC_log[i+11], name_phy, CDC_body_log, false, n);
	      n++;
	    }
	}
  
      if( CDCType !="B" )
	{
	  G4cout << "!!! CDC Type A " << G4endl; 
	  for (G4int i=0; i<2; i++)
	    {
	      rmin = cdc_radius[i+13]-cdc_off;
	      rmax = cdc_radius[i+13]+cdc_off;
	  
	      sprintf(name_sol,"CDC_sol_%2.2d",i+13);
	      sprintf(name_log,"CDC_log_%2.2d",i+13);
	  
	      CDC_tube[i+13] = new G4Tubs(name_sol, rmin, rmax, cdc_z,-0.5*twopi/(double)N_CDC_CELL[i+13], twopi/(double)N_CDC_CELL[i+13]);
	      CDC_log [i+13] = new G4LogicalVolume(CDC_tube[i+13],ChamberGas,name_log,0,0,0);
	      n=0;
	      for (G4int j=0; j<N_CDC_CELL[i+13]; j++)
		{
		  sprintf(name_phy,"CDC_phys%2.2d%03d",i+14,j+1);

		  G4ThreeVector xyzCounter(cdsPos_x, cdsPos_y, cdsPos_z);
		  G4RotationMatrix* rotCounterCDCcell = new G4RotationMatrix;
		  rotCounterCDCcell->rotateZ(twopi/(double)N_CDC_CELL[i+13]*(j+0.5*CDC_CELL_XXP[i+13]));
		  G4Transform3D posCounterCDCcell(*rotCounterCDCcell, xyzCounter);

		  CDC_phys[i+13][j] = new G4PVPlacement(posCounterCDCcell,CDC_log[i+13], name_phy, CDC_body_log, false, n);
		  n++;
		}
	    }
	}

      //--- Visualization ---//
      G4VisAttributes *CDC_body_att = new G4VisAttributes(LightBlue);
      CDC_body_att->SetForceWireframe(true);
      CDC_body_log->SetVisAttributes(CDC_body_att);
      G4VisAttributes *CDC_att = new G4VisAttributes(Yellow);
      CDC_att->SetForceWireframe(true);
      //CDC_att->SetForceSolid(true);

      for (G4int i=0; i<15; i++)
	{
	  //CDC_log[i]->SetSensitiveDetector(chamberSD);
	  //CDC_log[i]-> SetUserLimits( new G4UserLimits(1.0*mm) );
	  CDC_log[i]->SetVisAttributes(G4VisAttributes::Invisible);
	  //CDC_log[i]->SetVisAttributes(CDC_att);
	}


      //******************//
      //*** CDC window ***//
      //******************//
      G4double cdc_win_rmin = cdc_rmin-ThicknessOfCDCwCFRP*mm; 
      G4double cdc_win_rmax = cdc_rmin;
  
      G4VSolid* CDC_win_solid = new G4Tubs("CDC_sol",cdc_win_rmin,cdc_win_rmax,cdc_z,0.0, twopi);
      G4LogicalVolume* CDC_win_log   = new G4LogicalVolume(CDC_win_solid,Vacuum,"CDC_win_log",0,0,0); // Material CFRP
      G4PVPlacement* CDC_win_phys  = new G4PVPlacement(0, G4ThreeVector(cdsPos_x, cdsPos_y, cdsPos_z),CDC_win_log, "CDC_window",CDS_log, false, 0);

      G4double cdc_win_out_rmin = cdc_rmax; 
      G4double cdc_win_out_rmax = cdc_rmax+0.07*mm;

      G4VSolid* CDC_win_out_solid = new G4Tubs("CDC_out_sol",cdc_win_out_rmin,cdc_win_out_rmax,cdc_z,0.0, twopi);
      G4LogicalVolume* CDC_win_out_log   = new G4LogicalVolume(CDC_win_out_solid,Mylar,"CDC_win_out_log",0,0,0);
      G4PVPlacement* CDC_win_out_phys  = new G4PVPlacement(0, G4ThreeVector(cdsPos_x, cdsPos_y, cdsPos_z),CDC_win_out_log, "CDC_out_window",CDS_log, false, 0);

      //--- Visualization ---//
      G4VisAttributes *CDC_win_att = new G4VisAttributes(Purple);
      CDC_win_att->SetForceWireframe(true);
      CDC_win_log->SetVisAttributes(CDC_win_att);
      CDC_win_out_log->SetVisAttributes(CDC_win_att);
  

      int S_WIRE_VIS=0;
      int F_WIRE_VIS=0;

      //*****************//
      //*** CDC wires ***//
      //*****************//
      G4double r0 = 17.5*cm;

      G4double  radius[67];
      G4int     nwires[67];
      G4double  tilt[67];
      for(int i=0; i<67; i++)
	{
	  if(i==0)
	    radius[i] = r0; 
	  else
	    radius[i] = radius[i-1]+cdc_wire_dist;

	  if(i==1 || i==12 || i==13 || i==21 || i==22 || i==30 || i==31 || i==39 || i==40 || i==48 || i==49 || i==57 || i==58 || i==66)
	    radius[i] += 0.2*cm;
    
	  if(i<12)
	    {
	      nwires[i] =  72;
	      tilt[i] = 0;
	    }
	  else if (i>=12 && i<21)
	    {
	      nwires[i] =  90;
	      //tilt[i] = i==12 ? 0 : twopi/4-atan2(2*cdc_z, 2*radius[i]*sin(twopi/(double)nwires[i]));
	      tilt[i] = twopi/4-atan2(2*cdc_z, 2*radius[i]*sin(twopi/(double)nwires[i]*CDC_OFFSET/2.0));
	    }
	  else if (i>=21 && i<30)
	    { 
	      nwires[i] = 100;
	      //tilt[i] = i==21 ? 0 : -(twopi/4-atan2(2*cdc_z, 2*radius[i]*sin(twopi/(double)nwires[i])));
	      tilt[i] = twopi/4-atan2(2*cdc_z, 2*radius[i]*sin(twopi/(double)nwires[i]*CDC_OFFSET/2.0));
	    }
	  else if (i>=30 && i<39)
	    { 
	      nwires[i] = 120;
	      tilt[i] = 0;
	    }
	  else if (i>=39 && i<48)
	    { 
	      nwires[i] = 150;
	      //tilt[i] = i==39 ? 0 : twopi/4-atan2(2*cdc_z, 2*radius[i]*sin(twopi/(double)nwires[i]));
	      tilt[i] = twopi/4-atan2(2*cdc_z, 2*radius[i]*sin(twopi/(double)nwires[i]*CDC_OFFSET/2.0));
	    }
	  else if (i>=48 && i<57)
	    { 
	      nwires[i] = 160;
	      //tilt[i] = i==48 ? 0 : -(twopi/4-atan2(2*cdc_z, 2*radius[i]*sin(twopi/(double)nwires[i])));
	      tilt[i] = twopi/4-atan2(2*cdc_z, 2*radius[i]*sin(twopi/(double)nwires[i]*CDC_OFFSET/2.0));
	    }
	  else if (i>=57 && i<=66 )
	    {
	      nwires[i] = 180;
	      tilt[i] = 0;
	    }
	  else
	    {
	      nwires[i] =   0;
	    }
	}
      //for (int i=0; i<67; i++){
      //G4cout << radius[i] << " " << nwires[i] << " " << tilt[i]/deg << G4endl;
      //}   
 
      G4int total_sense_wires = 0;
      G4int total_field_wires = 0;

      G4int   idd = 0;
      G4int   IsSenseWire = 0;

      G4VSolid* CDC_sense_Wire[67];
      G4VSolid* CDC_field_Wire[67];
      G4LogicalVolume* CDCSenseWirelog[67];
      G4LogicalVolume* CDCFieldWirelog[67];

      std::vector<std::vector<G4PVPlacement*> > CDCSenseWirephy(67),CDCFieldWirephy(67);
      for(size_t i=0;i<CDCSenseWirephy.size();++i)
	{
	  CDCSenseWirephy[i].resize(nwires[i],nullptr);
	  CDCFieldWirephy[i].resize(nwires[i],nullptr);
	}
      for(G4int i=0; i<67; i++)
	{
	  if ( i== 0 || i== 1 || i== 3 || i== 5 || i== 7 || i== 9|| i==11 || i==12 || i==13 || i==15 || i==17 ||
	       i==19 || i==21 || i==22 || i==24 || i==26 || i==28|| i==30 || i==31 || i==33 || i==35 || i==37 ||
	       i==39 || i==40 || i==42 || i==44 || i==46 || i==48|| i==49 || i==51 || i==53 || i==55 || i==57 ||
	       i==58 || i==60 || i==62 || i==64 || i==66 )
	    {
	      idd = 0;
	    }
	  else
	    {
	      idd = 1;
	    }

	  if (i== 3 || i== 6 || i== 9 || i==15 || i==18 || i==24 || i==27 || i==33 || i==36 || 
	      i==42 || i==45 || i==51 || i==54 || i==60 || i==63             )
	    IsSenseWire = 1; 
	  else 
	    IsSenseWire = 0;

	  if(S_WIRE_VIS == 1)
	    CDC_sense_Wire[i] = new G4Tubs("CDC_sense_Wire", 0.0, 1*mm, cdc_z/cos(fabs(tilt[i])), 0.0, twopi);
	  else
	    CDC_sense_Wire[i] = new G4Tubs("CDC_sense_Wire", 0.0, 0.03*mm, cdc_z/cos(fabs(tilt[i])), 0.0, twopi);
      
	  if(F_WIRE_VIS == 1)
	    CDC_field_Wire[i] = new G4Tubs("CDC_field_Wire", 0.0, 1*mm, cdc_z/cos(fabs(tilt[i])), 0.0, twopi);
	  else
	    CDC_field_Wire[i] = new G4Tubs("CDC_field_Wire", 0.0, 0.10*mm, cdc_z/cos(fabs(tilt[i])), 0.0, twopi);

      
	  CDCSenseWirelog[i] =  new G4LogicalVolume(CDC_sense_Wire[i], W,  "CDCSWire_log", 0,0,0);
	  CDCFieldWirelog[i] =  new G4LogicalVolume(CDC_field_Wire[i], Al, "CDCFWire_log", 0,0,0);


	  G4double unit_ang = twopi/(G4double)nwires[i];
	  G4double Rwire    = radius[i];

	  for (int j=0; j<nwires[i]; j++)
	    {
	      G4double xxx = 0.0;
	      G4double yyy = 0.0;
	      G4double ang = 0.0;
	      if (idd==0)
		{
		  ang = unit_ang*(G4double)j;
		  xxx = Rwire*cos(ang);
		  yyy = Rwire*sin(ang); 
		}
	      else if (idd==1)
		{
		  ang = unit_ang*(G4double)j + unit_ang/2.0;
		  xxx = Rwire*cos(ang);
		  yyy = Rwire*sin(ang); 
		}
	
	      G4ThreeVector wire_pos(xxx, yyy, 0.0);
	      G4RotationMatrix* wire_rot = new G4RotationMatrix;
	      wire_rot->rotateX(-tilt[i]);
	      wire_rot->rotateZ(ang);
	      G4Transform3D wire_geom(*wire_rot, wire_pos);
	
	      //  /**
	      if (IsSenseWire == 1 )
		{
		  sprintf(name_phy,"CDCSWire_phy%06d", total_sense_wires);
		  //CDCSenseWirephy[total_sense_wires] =  new G4PVPlacement(0, wire_pos, CDCSenseWirelog, 
		  CDCSenseWirephy[i][j] =  new G4PVPlacement(wire_geom, CDCSenseWirelog[i], name_phy, CDC_body_log, false, 0 ); 	
		  total_sense_wires++;
		}
	      else if (IsSenseWire == 0 ) 
		{
		  sprintf(name_phy,"CDCFWire_phy%06d", total_field_wires);
		  //CDCFieldWirephy[total_field_wires] =  new G4PVPlacement(0, wire_pos, CDCFieldWirelog, 
		  CDCFieldWirephy[i][j] =  new G4PVPlacement(wire_geom, CDCFieldWirelog[i], name_phy, CDC_body_log, false, 0 );
		  total_field_wires++;
		}
	      else
		{
		  G4cout << "Wrong IsSenseWire assigned " << IsSenseWire << G4endl;
		}
	    }
	}  
  
      //G4cout << total_sense_wires << " " << total_field_wires << G4endl;
      //G4cout << "Number of wires generated = " << total_wires << G4endl;

      //--- Visualization ---//
      G4VisAttributes *CDCSenseWire_att = new G4VisAttributes(S_WIRE_VIS,Red);
      //CDCSenseWire_att->SetForceWireframe(true);
      CDCSenseWire_att->SetForceSolid(true);
      G4VisAttributes *CDCFieldWire_att = new G4VisAttributes(F_WIRE_VIS,Green);
      //CDCFieldWire_att->SetForceWireframe(true);
      CDCFieldWire_att->SetForceSolid(true);
      for (G4int i=0; i<67; i++)
	{  
	  CDCSenseWirelog[i]->SetVisAttributes(CDCSenseWire_att);
	  CDCFieldWirelog[i]->SetVisAttributes(CDCFieldWire_att);
	  //CDCSenseWirelog[i]->SetVisAttributes(G4VisAttributes::Invisible);
	  //CDCFieldWirelog[i]->SetVisAttributes(G4VisAttributes::Invisible);
	}
      //#endif


    }

  if(DoCDH==true)
    {
      //***********//
      //*** CDH ***//
      //***********//
      G4double cdh_rmin = CDH_RMAX*mm;
      G4double cdh_rmax = CDH_RMIN*mm;
      G4double cdh_z = CDH_Z/2.0*mm;

      char phys_name[100];

      G4Tubs* CDH_tube= new G4Tubs("CDH_tube", cdh_rmin, cdh_rmax, cdh_z, 0, twopi/(double)Ncdh);
      G4LogicalVolume* CDH_log = new G4LogicalVolume(CDH_tube, Scinti,"CDH_log", 0,0,0);
      G4int n=0;
      std::vector<G4PVPlacement*> CDH_phys(Ncdh,nullptr);
      for (G4int i=0; i<Ncdh; i++)
	{
	  sprintf(phys_name,"CDH_phys%02d", i+1);
	  G4ThreeVector xyzCounter(cdsPos_x, cdsPos_y, cdsPos_z);
	  G4RotationMatrix* rotCounterCDH = new G4RotationMatrix;
	  rotCounterCDH->rotateZ(twopi/(double)Ncdh*(i));
      
	  G4Transform3D posCounterCDH(*rotCounterCDH, xyzCounter);
      
	  CDH_phys[i] = new G4PVPlacement(posCounterCDH, CDH_log, phys_name, CDS_log, false, n);
	  n++;
	}
  
      //CDH_log->SetSensitiveDetector(counterSD);

      //--- Visualization ---//
      G4VisAttributes *CDH_att = new G4VisAttributes(Red);
      CDH_att->SetForceWireframe(true);
      //CDH_att->SetForceSolid(true);
      CDH_log->SetVisAttributes(CDH_att);

    }
  

  const G4double tarCham_rmax = TARGET_CHM_RMAX*mm;
  const G4double tarCham_z = 260.0/2.0*mm;
  const G4double tarChamPos_x = 0.0*m;
  const G4double tarChamPos_y = 0.0*m;
  const G4double tarChamPos_z = 0.0*m;

  if(DoTargetChamber == true)
    {
      //**********************//
      //*** Target-Chamber ***//
      //**********************//

      G4VSolid* TarCham_tube = 0;
      G4LogicalVolume* TarCham_log = 0; 
      G4PVPlacement* TarCham_phys = 0;

      if (TargetMaterialID!=2)
	{
	  TarCham_tube = new G4Tubs("TarCham_tube", 0.0, tarCham_rmax, tarCham_z, 0.0, twopi);
	  TarCham_log  = new G4LogicalVolume(TarCham_tube, Vacuum, "TarCham_log", 0,0,0);
	  TarCham_phys = new G4PVPlacement(0,G4ThreeVector(tarChamPos_x, tarChamPos_y, tarChamPos_z), TarCham_log, "Target_Chamber", CDS_log, false,0);
	}

      //**************//
      //*** Target ***//
      //**************//
      G4double target_r = TARGET_RMAX*mm;
      G4double target_z = TargetLength/2.0*mm;

      G4VSolid* Target_tube= new G4Tubs("Target_tube", 0.0, target_r, target_z, 0.0, twopi);
      G4LogicalVolume*  Target_log = nullptr;
      if (TargetMaterialID==0)
	Target_log = new G4LogicalVolume(Target_tube, Vacuum,   "Target_log", 0,0,0); // Material LHe3
      else if (TargetMaterialID==1)
	Target_log = new G4LogicalVolume(Target_tube, Vacuum,    "Target_log", 0,0,0); // Material CH2
      else if (TargetMaterialID==2)
	Target_log = new G4LogicalVolume(Target_tube, Vacuum,    "Target_log", 0,0,0); // Material D2Gas    
      else if (TargetMaterialID==999)
	Target_log = new G4LogicalVolume(Target_tube, Vacuum, "Target_log", 0,0,0);

      G4double targetPos_x = 0.0*m;
      G4double targetPos_y = 0.0*m;
      G4double targetPos_z = 0.0*m;

      G4PVPlacement* Target_phys;
      if (TargetMaterialID!=2)
	Target_phys = new G4PVPlacement(0, G4ThreeVector(targetPos_x, targetPos_y, targetPos_z), Target_log, "Target", TarCham_log, false,0);
      else
	Target_phys = new G4PVPlacement(0, G4ThreeVector(targetPos_x, targetPos_y, targetPos_z), Target_log, "Target", CDS_log, false,0);


      G4VSolid* TarPET_tube =0; 	
      G4LogicalVolume* TarPET_log  =0; 
      G4PVPlacement* TarPET_phys =0;

      G4VSolid* TarAl_tube =0; 
      G4LogicalVolume* TarAl_log =0; 
      G4PVPlacement* TarAl_phys=0; 

      G4VSolid* TarCFRP_tube = 0;
      G4LogicalVolume* TarCFRP_log =0; 
      G4PVPlacement* TarCFRP_phys  =0; 

      G4VSolid* TarChmEndCapU_solid = 0;
      G4LogicalVolume* TarChmEndCapU_log =0; 
      G4PVPlacement* TarChmEndCapU_phys  =0; 
  
      G4VSolid* TarChmEndCapD_solid =0;
      G4LogicalVolume* TarChmEndCapD_log=0;
      G4PVPlacement* TarChmEndCapD_phys  =0;

      G4VSolid* TarSUS_tube = 0;
      G4LogicalVolume* TarSUS_log =0;  
      G4PVPlacement* TarSUS_phys=0;

      G4VSolid* TarSUSEndCapU_solid =0;
      G4LogicalVolume* TarSUSEndCapU_log =0; 
      G4PVPlacement* TarSUSEndCapU_phys=0;
                   
      G4VSolid* TarSUSEndCapD_solid=0;
      G4LogicalVolume* TarSUSEndCapD_log=0;  
      G4PVPlacement* TarSUSEndCapD_phys =0;

  
      if (TargetMaterialID!=2)
	{
	  //******************//
	  //*** Target-PET ***//
	  //******************//
	  G4double tarPET_rmin = target_r;
	  G4double tarPET_rmax = target_r+ThicknessOfTgtCell*mm;
      
	  TarPET_tube = new G4Tubs("TarPET_tube", tarPET_rmin, tarPET_rmax, target_z, 0.0, twopi);
	  TarPET_log  = new G4LogicalVolume(TarPET_tube, Be, "TarPET_log", 0,0,0);
	  TarPET_phys = new G4PVPlacement(0, G4ThreeVector(targetPos_x, targetPos_y, targetPos_z), TarPET_log, "PET", TarCham_log, false,0);

	  //*****************//
	  //*** Target-Al ***//
	  //*****************//
	  G4double tarAl_rmin = 60.0*mm;
	  G4double tarAl_rmax = tarAl_rmin+ThicknessOfTgtAl*mm;
	  G4double tarAl_z    = 20./2.0*cm;

	  TarAl_tube= new G4Tubs("TarAl_tube", tarAl_rmin, tarAl_rmax, tarAl_z, 0.0, twopi);
	  TarAl_log = new G4LogicalVolume(TarAl_tube, Al, "TarAl_log", 0,0,0);
	  TarAl_phys = new G4PVPlacement(0,G4ThreeVector(targetPos_x, targetPos_y, targetPos_z),TarAl_log, "Al", TarCham_log, false,0);


	  //*******************//
	  //*** Target-CFRP ***//
	  //*******************//
	  G4double tarCFRP_rmax = tarCham_rmax;
	  G4double tarCFRP_rmin = tarCFRP_rmax-ThicknessOfTgtCFRP*mm;
	  G4double tarCFRP_z    = 26/2.0*cm;

	  TarCFRP_tube= new G4Tubs("TarCFRP_tube", tarCFRP_rmin, tarCFRP_rmax, tarCFRP_z, 0.0, twopi);
	  TarCFRP_log = new G4LogicalVolume(TarCFRP_tube, Vacuum, "TarCFRP_log", 0,0,0); // Material CFRP
	  TarCFRP_phys = new G4PVPlacement(0,G4ThreeVector(targetPos_x, targetPos_y, targetPos_z),TarCFRP_log, "CFRP", TarCham_log, false,0);

	  TarCFRP_log->SetUserLimits(new G4UserLimits(1.0*mm));

	  //*****************************//
	  //*** Target Chamber endcap ***//
	  //*****************************//
	  TarChmEndCapU_solid     = new G4Tubs("TarChmEndCapU_solid", 0.0, tarCFRP_rmax,0.25*mm, 0.0, twopi);
	  TarChmEndCapU_log       = new G4LogicalVolume(TarChmEndCapU_solid, Al, "TarChmEndCapU_log", 0, 0, 0);
	  TarChmEndCapU_phys      = new G4PVPlacement(0,G4ThreeVector(targetPos_x, targetPos_y, targetPos_z-tarCFRP_z-0.25*mm),
						      TarChmEndCapU_log, "TarChmEndCapU", CDS_log,false,0);
	  //TarChmEndCapU_log->SetSensitiveDetector(chamberSD);

	  TarChmEndCapD_solid     = new G4Tubs("TarChmEndCapD_solid", 0.0, tarCFRP_rmax,0.25*mm, 0.0, twopi);
	  TarChmEndCapD_log       = new G4LogicalVolume(TarChmEndCapD_solid, Al, "TarChmEndCapU_log", 0, 0, 0);
	  TarChmEndCapD_phys      = new G4PVPlacement(0,G4ThreeVector(targetPos_x, targetPos_y, targetPos_z+tarCFRP_z+0.25*mm),
						      TarChmEndCapD_log, "TarChmEndCapD", CDS_log,false,0);
	  //TarChmEndCapD_log->SetSensitiveDetector(chamberSD);

	  //--- Visualization ---//
	  TarCham_log->SetVisAttributes(G4VisAttributes::Invisible);
	  G4VisAttributes *Target_att = new G4VisAttributes(Green);
	  Target_att->SetForceWireframe(true);
	  Target_att->SetForceSolid(true);
	  Target_log->SetVisAttributes(Target_att);
	  G4VisAttributes *TarPET_att = new G4VisAttributes(Blue);
	  TarPET_att->SetForceWireframe(true);
	  TarPET_log->SetVisAttributes(TarPET_att);
	  G4VisAttributes *TarAl_att = new G4VisAttributes(Orange);
	  TarAl_att->SetForceWireframe(true);
	  TarAl_log->SetVisAttributes(TarAl_att);
	  G4VisAttributes *TarCFRP_att = new G4VisAttributes(Purple);
	  TarCFRP_att->SetForceWireframe(true);
	  TarCFRP_log->SetVisAttributes(TarCFRP_att);
	  G4VisAttributes *TarChmEndCap_att = new G4VisAttributes(LightBlue);
	  TarChmEndCap_att->SetForceWireframe(true);
	  TarChmEndCapU_log->SetVisAttributes(TarChmEndCap_att);
	  TarChmEndCapD_log->SetVisAttributes(TarChmEndCap_att);

	}
      else
	{
	  //******************//
	  //*** Target-SUS ***//
	  //******************//
	  G4double tarSUS_rmin = target_r;
	  G4double tarSUS_rmax = target_r+0.05*mm;

	  TarSUS_tube = new G4Tubs("TarSUS_tube", tarSUS_rmin, tarSUS_rmax, target_z, 0.0, twopi);
	  TarSUS_log  = new G4LogicalVolume(TarSUS_tube, Vacuum, "TarSUS_log", 0,0,0); // Material SUS 
	  TarSUS_phys = new G4PVPlacement(0,G4ThreeVector(targetPos_x, targetPos_y, targetPos_z),
					  TarSUS_log, "SUS", CDS_log, false,0);
    
	  //*****************************//
	  //*** Target Chamber endcap ***//
	  //*****************************//
	  TarSUSEndCapU_solid     = new G4Tubs("TarSUSEndCapU_solid", 0.0, tarSUS_rmax,0.5*mm, 0.0, twopi);
	  TarSUSEndCapU_log       = new G4LogicalVolume(TarSUSEndCapU_solid, Vacuum, "TarSUSEndCapU_log", 0, 0, 0); // Material SUS
	  TarSUSEndCapU_phys      = new G4PVPlacement(0,G4ThreeVector(targetPos_x, targetPos_y, targetPos_z-target_z-0.5*mm),
						      TarSUSEndCapU_log, "TarSUSEndCapU", CDS_log,false,0);
	  //TarChmEndCapU_log->SetSensitiveDetector(chamberSD);

	  TarSUSEndCapD_solid     = new G4Tubs("TarSUSEndCapD_solid", 0.0, tarSUS_rmax,0.5*mm, 0.0, twopi);
	  TarSUSEndCapD_log       = new G4LogicalVolume(TarSUSEndCapD_solid, Vacuum, "TarSUSEndCapU_log", 0, 0, 0); // Material SUS
	  TarSUSEndCapD_phys      = new G4PVPlacement(0,G4ThreeVector(targetPos_x, targetPos_y, targetPos_z+target_z+0.5*mm),
						      TarSUSEndCapD_log, "TarSUSEndCapD", CDS_log,false,0);
	  //TarChmEndCapD_log->SetSensitiveDetector(chamberSD);

	  //--- Visualization ---//
	  G4VisAttributes *Target_att = new G4VisAttributes(Green);
	  Target_att->SetForceWireframe(true);
	  Target_att->SetForceSolid(true);
	  Target_log->SetVisAttributes(Target_att);
	  G4VisAttributes *TarSUS_att = new G4VisAttributes(Blue);
	  TarSUS_att->SetForceWireframe(true);
	  TarSUS_log->SetVisAttributes(TarSUS_att);
	  G4VisAttributes *TarSUSEndCap_att = new G4VisAttributes(LightBlue);
	  TarSUSEndCap_att->SetForceWireframe(true);
	  TarSUSEndCapU_log->SetVisAttributes(TarSUSEndCap_att);
	  TarSUSEndCapD_log->SetVisAttributes(TarSUSEndCap_att);
	}
      //***************************//
      //*** Charge Veto Counter ***//
      //***************************//
      G4double CVC_Length = 10.0/2.0*mm;
      G4double CVC_Zpos   = target_z+CVC_Length;
      G4double cvc_r = 40.0*mm;
      //G4double cvc_r = 60.0*mm;

      G4VSolid* CVC_solid = new G4Tubs("CVC_solid",0.0, cvc_r, CVC_Length, 0.0, twopi);
      G4LogicalVolume* CVC_log = new G4LogicalVolume(CVC_solid, Scinti, "CVC_log", 0, 0, 0);
      G4PVPlacement* CVC_phys = new G4PVPlacement(0,G4ThreeVector(cdsPos_x, cdsPos_y, CVC_Zpos),
						  CVC_log, "CVC_phys", TarCham_log,false,0);  

      //CVC_log->SetSensitiveDetector(counterSD); 

      //--- Visualization ---//
      G4VisAttributes *CVC_att = new G4VisAttributes(Red);
      //CVC_att->SetForceWireframe(true);
      CVC_log->SetVisAttributes(CVC_att);

      //#if 1 // by fujioka
      //***********//
      //*** ZVC ***//
      //***********//
      G4double zvc_z = ZVCLength/2.0*mm;

      //#if 1
      int VertexType = 1;
      ZVertexChamber = 1;
      if(ZVertexChamber == 1)
	{
      

	  if(VertexType == 1)
	    {
	      G4double zvc_rmin    = 85.0*mm;
	      G4double zvc_rmax    = 140.0*mm;
      
	      G4Tubs* ZVC_body_tube[2];
	      G4LogicalVolume* ZVC_body_log[2] ;
	      G4PVPlacement* ZVC_body_phys[2];

	      ZVC_body_tube[0]= new G4Tubs("ZVC_body_tube1", zvc_rmin, zvc_rmax, zvc_z, 0.0, twopi);
	      ZVC_body_log[0] = new G4LogicalVolume(ZVC_body_tube[0], Vacuum, "ZVC_body_log1", 0,0,0); // Material ArCH4_90_10
	      ZVC_body_phys[0] = new G4PVPlacement(0,G4ThreeVector(cdsPos_x, cdsPos_y, cdsPos_z),ZVC_body_log[0], "ZVC(TPC)", CDS_log, false,0);
      
	      char name[100];

	      G4VSolid* ZVC_tube [1][8];
	      G4LogicalVolume* ZVC_log[1][8]; 
	      G4PVPlacement* ZVC_phys[1][8];      
	      G4double rmin = 0.0;
	      G4double rmax = 0.0;
	      // Kapton
	      rmin = zvc_rmin;
	      rmax = rmin+0.4*mm;
	      sprintf(name, "ZVC_inner_Kapton");
	      ZVC_tube[0][0] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
	      ZVC_log [0][0] = new G4LogicalVolume(ZVC_tube[0][0],Kapton,name,0,0,0);
	      ZVC_phys[0][0] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[0][0], name, ZVC_body_log[0], false, 0);  
	      // Kapton
	      rmin = 95.0*mm;
	      rmax = rmin+0.4*mm;
	      sprintf(name, "ZVC_field_in_Kapton");
	      ZVC_tube[0][1] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
	      ZVC_log [0][1] = new G4LogicalVolume(ZVC_tube[0][1],Kapton,name,0,0,0);
	      ZVC_phys[0][1] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[0][1], name, ZVC_body_log[0], false, 0);  
	      // Cu
	      rmin = rmax;
	      rmax = rmin+0.018*mm;
	      sprintf(name, "ZVC_field_in_Cu");
	      ZVC_tube[0][2] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
	      ZVC_log [0][2] = new G4LogicalVolume(ZVC_tube[0][2],Cu,name,0,0,0);
	      ZVC_phys[0][2] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[0][2], name, ZVC_body_log[0], false, 0);  
	      // Kapton
	      rmin = 130.0*mm;
	      rmax = rmin+0.4*mm;
	      sprintf(name, "ZVC_field_out_Kapton");
	      ZVC_tube[0][3] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
	      ZVC_log [0][3] = new G4LogicalVolume(ZVC_tube[0][3],Kapton,name,0,0,0);
	      ZVC_phys[0][3] = new G4PVPlacement(0, G4ThreeVector(0,0,0), ZVC_log[0][3], name, ZVC_body_log[0], false, 0);  
	      // Cu
	      rmin = rmax;
	      rmax = rmin+0.018*mm;
	      sprintf(name, "ZVC_field_out_Cu");
	      ZVC_tube[0][4] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
	      ZVC_log [0][4] = new G4LogicalVolume(ZVC_tube[0][4],Cu,name,0,0,0);
	      ZVC_phys[0][4] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[0][4], name, ZVC_body_log[0], false, 0);
      
	      // Kapton
	      rmin = zvc_rmax-1.0*mm-0.4*mm;
	      rmax = rmin+0.4*mm;
	      sprintf(name, "ZVC_outer_Kapton");
	      ZVC_tube[0][5] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
	      ZVC_log [0][5] = new G4LogicalVolume(ZVC_tube[0][5],Kapton,name,0,0,0);
	      ZVC_phys[0][5] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[0][5], name, ZVC_body_log[0], false, 0);  
	      // Al
	      rmin = rmax;
	      rmax = zvc_rmax;
	      sprintf(name, "ZVC_outer_Al");
	      ZVC_tube[0][6] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
	      ZVC_log [0][6] = new G4LogicalVolume(ZVC_tube[0][6],Al,name,0,0,0);
	      ZVC_phys[0][6] = new G4PVPlacement(0, G4ThreeVector(0,0,0),
						 ZVC_log[0][6], name, ZVC_body_log[0], false, 0);

	      // readout-cell
	      rmin = 95.0*mm+0.4*mm+0.018*mm;
	      rmax = rmin+0.001*mm;
	      const G4double ZVC_cell_size = ZVertexChamberCell; // full width
	      G4int ZVC_n_cell = (G4int)(2*zvc_z/ZVC_cell_size);
	      G4VSolid* ZVC_cell_tube[2];
	      G4LogicalVolume* ZVC_cell_log[2];
	      std::vector< std::vector<G4PVPlacement*> > ZVC_cell_phys(2,std::vector<G4PVPlacement*>(ZVC_n_cell,nullptr));

	      ZVC_cell_tube[0] = new G4Tubs(name, rmin, rmax, ZVC_cell_size/2, 0.0, twopi);
	      ZVC_cell_log[0]  = new G4LogicalVolume(ZVC_cell_tube[0],Vacuum,name,0,0,0);

	      char phys_name[100];
	      G4int n=0;
	      
	      for (G4int i=1; i<=ZVC_n_cell; i++)
		{
		  sprintf(phys_name,"ZVC_phys1%04d", i);
		  G4double zvcPos_z = -zvc_z+ZVC_cell_size*0.5+ZVC_cell_size*(i-1);
		  G4ThreeVector ZVC_Pos(cdsPos_x, cdsPos_y, zvcPos_z);
		  ZVC_cell_phys[0][i-1] = new G4PVPlacement(0, ZVC_Pos,ZVC_cell_log[0], phys_name,ZVC_body_log[0], false, n);
		  n++;
		}

	      // readout-cell
	      rmin = 130.0*mm-0.001*mm;
	      rmax = rmin+0.001*mm;
	      ZVC_cell_tube[1] = new G4Tubs(name, rmin, rmax, ZVC_cell_size/2, 0.0, twopi);
	      ZVC_cell_log[1]  = new G4LogicalVolume(ZVC_cell_tube[1],Vacuum,name,0,0,0);
	      n=0;
	      for (G4int i=1; i<=ZVC_n_cell; i++)
		{
		  sprintf(phys_name,"ZVC_phys2%04d", i);
		  G4double zvcPos_z = -zvc_z+ZVC_cell_size*0.5+ZVC_cell_size*(i-1);
		  G4ThreeVector ZVC_Pos(cdsPos_x, cdsPos_y, zvcPos_z);
		  ZVC_cell_phys[1][i-1] = new G4PVPlacement(0, ZVC_Pos,ZVC_cell_log[1], phys_name,ZVC_body_log[0], false, n);
		  n++;
		}

	      // drift-region
	      rmin = 95.0*mm+0.4*mm+0.018*mm+0.001*mm;
	      rmax = 130.0*mm-0.001*mm;
	      sprintf(name, "ZVC_drift_region");
	      ZVC_tube[0][7] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
	      ZVC_log [0][7] = new G4LogicalVolume(ZVC_tube[0][7],Vacuum,name,0,0,0); // Material ArCH4_90_10
	      ZVC_phys[0][7] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[0][7], name, ZVC_body_log[0], false, 0);

	
	      //ZVC_cell_log[0]->SetSensitiveDetector(chamberSD);
	      //ZVC_cell_log[1]->SetSensitiveDetector(chamberSD);
	
	      //ZVC_log[0][7]->SetSensitiveDetector(counterSD);

	      //--- Visualization ---//
	      G4VisAttributes *ZVC_body_att = new G4VisAttributes(LightBlue);
	      ZVC_body_att->SetForceWireframe(true);
	      //ZVC_body_att->SetForceSolid(true);
	      ZVC_body_log[0]->SetVisAttributes(ZVC_body_att);
	      for (G4int i=0; i<7; i++)
		ZVC_log [0][i]->SetVisAttributes(G4VisAttributes::Invisible);

	      G4VisAttributes *ZVC_cell_att = new G4VisAttributes(Blue);
	      ZVC_cell_att->SetForceWireframe(true);
	      //ZVC_cell_att->SetForceSolid(true);
	      //ZVC_cell_log[0]->SetVisAttributes(ZVC_cell_att);
	      //ZVC_cell_log[1]->SetVisAttributes(ZVC_cell_att);
	      ZVC_cell_log[0]->SetVisAttributes(G4VisAttributes::Invisible);
	      ZVC_cell_log[1]->SetVisAttributes(G4VisAttributes::Invisible);
	    }
    
	  else
	    {
	      //*******************//
	      //*** outer layer ***//
	      //*******************//
	      G4double zvc_rmin    = (ZVC_R-6.0)*mm;
	      G4double zvc_rmax    = (ZVC_R+1.0)*mm;
      	      
	      G4Tubs* ZVC_body_tube[2];
	      G4LogicalVolume* ZVC_body_log[2] ;
	      G4PVPlacement* ZVC_body_phys[2];
	      ZVC_body_tube[0]= new G4Tubs("ZVC_body_tube1", zvc_rmin, zvc_rmax, zvc_z, 0.0, twopi);
	      ZVC_body_log[0] = new G4LogicalVolume(ZVC_body_tube[0], ChamberGas, "ZVC_body_log1", 0,0,0);
	      ZVC_body_phys[0] = new G4PVPlacement(0, G4ThreeVector(cdsPos_x, cdsPos_y, cdsPos_z),ZVC_body_log[0], "ZVC1", CDS_log, false,0);
	      
	      char name[100];

	      G4VSolid* ZVC_tube [2][11];
	      G4LogicalVolume* ZVC_log[2][11]; 
	      G4PVPlacement* ZVC_phys[2][11];      

	      G4double rmin = 0.0;
	      G4double rmax = 0.0;
	      
	      // cathode-Kapton
	      rmin = (ZVC_R-5.196)*mm;
	      rmax = rmin+0.01*mm;
	      sprintf(name, "ZVC_cathode-Kapton1");
	      ZVC_tube[0][0] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
	      ZVC_log [0][0] = new G4LogicalVolume(ZVC_tube[0][0],Kapton,name,0,0,0);
	      ZVC_phys[0][0] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[0][0], name, ZVC_body_log[0], false, 0);  
	      // cathode-Al
	      rmin = rmax;
	      rmax = rmin+0.01*mm;
	      sprintf(name, "ZVC_cathode-Al1");
	      ZVC_tube[0][1] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
	      ZVC_log [0][1] = new G4LogicalVolume(ZVC_tube[0][1],Al,name,0,0,0);
	      ZVC_phys[0][1] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[0][1], name, ZVC_body_log[0], false, 0);  
	      // GEM1-Cu
	      rmin = rmax+2.0*mm;
	      rmax = rmin+0.008*mm;
	      sprintf(name, "ZVC_GEM1-Cu1");
	      ZVC_tube[0][2] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
	      ZVC_log [0][2] = new G4LogicalVolume(ZVC_tube[0][2],Cu,name,0,0,0);
	      ZVC_phys[0][2] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[0][2], name, ZVC_body_log[0], false, 0);  
	      // GEM1-Kapton
	      rmin = rmax;
	      rmax = rmin+0.1*mm;
	      sprintf(name, "ZVC_GEM1-Kapton1");
	      ZVC_tube[0][3] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
	      ZVC_log [0][3] = new G4LogicalVolume(ZVC_tube[0][3],Kapton,name,0,0,0);
	      ZVC_phys[0][3] = new G4PVPlacement(0, G4ThreeVector(0,0,0), ZVC_log[0][3], name, ZVC_body_log[0], false, 0);  
	      // GEM1-Cu
	      rmin = rmax;
	      rmax = rmin+0.008*mm;
	      sprintf(name, "ZVC_GEM1-Cu1");
	      ZVC_tube[0][4] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
	      ZVC_log [0][4] = new G4LogicalVolume(ZVC_tube[0][4],Cu,name,0,0,0);
	      ZVC_phys[0][4] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[0][4], name, ZVC_body_log[0], false, 0);
	      // GEM2-Cu
	      rmin = rmax+1.0*mm;
	      rmax = rmin+0.005*mm;
	      sprintf(name, "ZVC_GEM2-Cu1");
	      ZVC_tube[0][5] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
	      ZVC_log [0][5] = new G4LogicalVolume(ZVC_tube[0][5],Cu,name,0,0,0);
	      ZVC_phys[0][5] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[0][5], name, ZVC_body_log[0], false, 0);
	      // GEM2-Kapton
	      rmin = rmax;
	      rmax = rmin+0.05*mm;
	      sprintf(name, "ZVC_GEM1-Kapton1");
	      ZVC_tube[0][6] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
	      ZVC_log [0][6] = new G4LogicalVolume(ZVC_tube[0][6],Kapton,name,0,0,0);
	      ZVC_phys[0][6] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[0][6], name, ZVC_body_log[0], false, 0);  
	      // GEM2-Cu
	      rmin = rmax;
	      rmax = rmin+0.005*mm;
	      sprintf(name, "ZVC_GEM2-Cu1");
	      ZVC_tube[0][7] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
	      ZVC_log [0][7] = new G4LogicalVolume(ZVC_tube[0][7],Cu,name,0,0,0);
	      ZVC_phys[0][7] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[0][7], name, ZVC_body_log[0], false, 0);
	      // readout-Al
	      rmin = rmax+2.0*mm;
	      rmax = rmin+0.01*mm;
	      sprintf(name, "ZVC_readout-Al1");
	      ZVC_tube[0][8] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
	      ZVC_log [0][8] = new G4LogicalVolume(ZVC_tube[0][8],Al,name,0,0,0);
	      ZVC_phys[0][8] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[0][8], name, ZVC_body_log[0], false, 0);  
	      // readout-Al-cell
	      const G4double     ZVC_cell_size = ZVertexChamberCell; // full width
	      G4int ZVC_n_cell = (G4int)(2*zvc_z/ZVC_cell_size);

	      G4VSolid* ZVC_cell_tube[2];
	      G4LogicalVolume* ZVC_cell_log[2];
	      std::vector<std::vector<G4PVPlacement*> > ZVC_cell_phys(2,std::vector<G4PVPlacement*>(ZVC_n_cell,nullptr));
	      ZVC_cell_tube[0] = new G4Tubs(name, rmin, rmax, ZVC_cell_size/2, 0.0, twopi);
	      ZVC_cell_log[0]  = new G4LogicalVolume(ZVC_cell_tube[0],Al,name,0,0,0);
	      G4int n=0;
	      char phys_name[100];
	      for (G4int i=0; i<ZVC_n_cell; i++)
		{
		  sprintf(phys_name,"ZVC_phys1%04d", i+1);
		  G4double zvcPos_z = -zvc_z+ZVC_cell_size*0.5+ZVC_cell_size*(i);
		  G4ThreeVector ZVC_Pos(cdsPos_x, cdsPos_y, zvcPos_z);
		  ZVC_cell_phys[0][i] = new G4PVPlacement(0, ZVC_Pos,ZVC_cell_log[0], phys_name,ZVC_log[0][8], false, n);
		  n++;
		}
	      
	      // readout-Kapton
	      rmin = rmax;
	      rmax = rmin+0.01*mm;
	      sprintf(name, "ZVC_readout-Kapton1");
	      ZVC_tube[0][9] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
	      ZVC_log [0][9] = new G4LogicalVolume(ZVC_tube[0][9],Kapton,name,0,0,0);
	      ZVC_phys[0][9] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[0][9], name, ZVC_body_log[0], false, 0);    
	      // readout-Al
	      rmin = rmax;
	      rmax = rmin+0.01*mm;
	      sprintf(name, "ZVC_readout-Al1");
	      ZVC_tube[0][10] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
	      ZVC_log [0][10] = new G4LogicalVolume(ZVC_tube[0][10],Al,name,0,0,0);
	      ZVC_phys[0][10] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[0][10], name, ZVC_body_log[0], false, 0);  
	      
	      
	      //--- Visualization ---//
	      G4VisAttributes *ZVC_body_att = new G4VisAttributes(LightBlue);
	      ZVC_body_att->SetForceWireframe(true);
	      //ZVC_body_att->SetForceSolid(true);
	      ZVC_body_log[0]->SetVisAttributes(ZVC_body_att);
	      for (G4int i=0; i<11; i++)
		{
		  ZVC_log [0][i]->SetVisAttributes(G4VisAttributes::Invisible);
		}
	      G4VisAttributes *ZVC_cell_att = new G4VisAttributes(Blue);
	      ZVC_cell_att->SetForceWireframe(true);
	      //ZVC_cell_att->SetForceSolid(true);
	      //ZVC_cell_log[0]->SetVisAttributes(ZVC_cell_att);
	      ZVC_cell_log[0]->SetVisAttributes(G4VisAttributes::Invisible);
	      
       
	      //#define ZVC2
	      //#ifdef ZVC2
	      bool Def_ZVC2 = true;
	      if(Def_ZVC2 == true)
		{
	      
		  //*******************//
		  //*** inner layer ***//
		  //*******************//
		  zvc_rmin    = (ZVC_R-35.0-6.0)*mm;
		  zvc_rmax    = (ZVC_R-35.0+1.0)*mm;
	      
		  ZVC_body_tube[1]= new G4Tubs("ZVC_body_tube2", zvc_rmin, zvc_rmax, zvc_z, 0.0, twopi);
		  ZVC_body_log[1] = new G4LogicalVolume(ZVC_body_tube[1], ChamberGas, "ZVC_body_log2", 0,0,0);	      
		  ZVC_body_phys[1] = new G4PVPlacement(0,G4ThreeVector(cdsPos_x, cdsPos_y, cdsPos_z),ZVC_body_log[1], "ZVC2", CDS_log, false,0);
	      
		  // cathode-Kapton
		  rmin = (ZVC_R-35.0-5.196)*mm;
		  rmax = rmin+0.01*mm;
		  sprintf(name, "ZVC_cathode-Kapton2");
		  ZVC_tube[1][0] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
		  ZVC_log [1][0] = new G4LogicalVolume(ZVC_tube[1][0],Kapton,name,0,0,0);
		  ZVC_phys[1][0] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[1][0], name, ZVC_body_log[1], false, 0);  
		  // cathode-Al
		  rmin = rmax;
		  rmax = rmin+0.01*mm;
		  sprintf(name, "ZVC_cathode-Al2");
		  ZVC_tube[1][1] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
		  ZVC_log [1][1] = new G4LogicalVolume(ZVC_tube[1][1],Al,name,0,0,0);
		  ZVC_phys[1][1] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[1][1], name, ZVC_body_log[1], false, 0);  
		  // GEM1-Cu
		  rmin = rmax+2.0*mm;
		  rmax = rmin+0.008*mm;
		  sprintf(name, "ZVC_GEM1-Cu2");
		  ZVC_tube[1][2] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
		  ZVC_log [1][2] = new G4LogicalVolume(ZVC_tube[1][2],Cu,name,0,0,0);
		  ZVC_phys[1][2] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[1][2], name, ZVC_body_log[1], false, 0);  
		  // GEM1-Kapton
		  rmin = rmax;
		  rmax = rmin+0.1*mm;
		  sprintf(name, "ZVC_GEM1-Kapton2");
		  ZVC_tube[1][3] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
		  ZVC_log [1][3] = new G4LogicalVolume(ZVC_tube[1][3],Kapton,name,0,0,0);
		  ZVC_phys[1][3] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[1][3], name, ZVC_body_log[1], false, 0);  
		  // GEM1-Cu
		  rmin = rmax;
		  rmax = rmin+0.008*mm;
		  sprintf(name, "ZVC_GEM1-Cu2");
		  ZVC_tube[1][4] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
		  ZVC_log [1][4] = new G4LogicalVolume(ZVC_tube[1][4],Cu,name,0,0,0);
		  ZVC_phys[1][4] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[1][4], name, ZVC_body_log[1], false, 0);
		  // GEM2-Cu
		  rmin = rmax+1.0*mm;
		  rmax = rmin+0.005*mm;
		  sprintf(name, "ZVC_GEM2-Cu2");
		  ZVC_tube[1][5] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
		  ZVC_log [1][5] = new G4LogicalVolume(ZVC_tube[1][5],Cu,name,0,0,0);
		  ZVC_phys[1][5] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[1][5], name, ZVC_body_log[1], false, 0);
		  // GEM2-Kapton
		  rmin = rmax;
		  rmax = rmin+0.05*mm;
		  sprintf(name, "ZVC_GEM1-Kapton2");
		  ZVC_tube[1][6] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
		  ZVC_log [1][6] = new G4LogicalVolume(ZVC_tube[1][6],Kapton,name,0,0,0);
		  ZVC_phys[1][6] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[1][6], name, ZVC_body_log[1], false, 0);  
		  // GEM2-Cu
		  rmin = rmax;
		  rmax = rmin+0.005*mm;
		  sprintf(name, "ZVC_GEM2-Cu2");
		  ZVC_tube[1][7] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
		  ZVC_log [1][7] = new G4LogicalVolume(ZVC_tube[1][7],Cu,name,0,0,0);
		  ZVC_phys[1][7] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[1][7], name, ZVC_body_log[1], false, 0);
		  // readout-Al
		  rmin = rmax+2.0*mm;
		  rmax = rmin+0.01*mm;
		  sprintf(name, "ZVC_readout-Al2");
		  ZVC_tube[1][8] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
		  ZVC_log [1][8] = new G4LogicalVolume(ZVC_tube[1][8],Al,name,0,0,0);
		  ZVC_phys[1][8] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[1][8], name, ZVC_body_log[1], false, 0);  
		  // readout-Al-cell
		  ZVC_cell_tube[1] = new G4Tubs(name, rmin, rmax, ZVC_cell_size/2, 0.0, twopi);
		  ZVC_cell_log[1]  = new G4LogicalVolume(ZVC_cell_tube[1],Al,name,0,0,0);
		  n=0;
		  for (G4int i=1; i<=ZVC_n_cell; i++)
		    {
		      sprintf(phys_name,"ZVC_phys2%04d", i);
		      G4double zvcPos_z = -zvc_z+ZVC_cell_size*0.5+ZVC_cell_size*(i-1);
		      G4ThreeVector ZVC_Pos(cdsPos_x, cdsPos_y, zvcPos_z);
		      ZVC_cell_phys[1][i-1] = new G4PVPlacement(0, ZVC_Pos,ZVC_cell_log[1], phys_name,ZVC_log[1][8], false, n);
		      n++;
		    }
    
		  // readout-Kapton
		  rmin = rmax;
		  rmax = rmin+0.01*mm;
		  sprintf(name, "ZVC_readout-Kapton2");
		  ZVC_tube[1][9] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
		  ZVC_log [1][9] = new G4LogicalVolume(ZVC_tube[1][9],Kapton,name,0,0,0);
		  ZVC_phys[1][9] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[1][9], name, ZVC_body_log[1], false, 0);    
		  // readout-Al
		  rmin = rmax;
		  rmax = rmin+0.01*mm;
		  sprintf(name, "ZVC_readout-Al2");
		  ZVC_tube[1][10] = new G4Tubs(name, rmin, rmax, zvc_z, 0.0, twopi);
		  ZVC_log [1][10] = new G4LogicalVolume(ZVC_tube[1][10],Al,name,0,0,0);
		  ZVC_phys[1][10] = new G4PVPlacement(0, G4ThreeVector(0,0,0),ZVC_log[1][10], name, ZVC_body_log[1], false, 0);  
	      
	      
		  //--- Visualization ---//
		  ZVC_body_log[1]->SetVisAttributes(ZVC_body_att);
		  for (G4int i=0; i<11; i++)
		    {
		      ZVC_log [1][i]->SetVisAttributes(G4VisAttributes::Invisible);
		    }
		  //ZVC_cell_log[1]->SetVisAttributes(ZVC_cell_att);
		  ZVC_cell_log[1]->SetVisAttributes(G4VisAttributes::Invisible);
		  //#endif
		}
	    }
	}
    }
  //#endif
  //#endif
  bool Kaon_VetoCounter=false;
  if(Kaon_VetoCounter == true)
    {
      //#if 0 // by fujioka
      //*******************************//
      //*** Kaon Decay Veto Counter ***//
      //*******************************//
      G4VSolid* KDV_solid = 0;
#if 0 // phi300mm
      G4double KDV_Length = (47.5*cm-zvc_z)/2.0;//30.0cm
      G4double KDV_Zpos   = -(zvc_z+KDV_Length);
      KDV_solid              = new G4Tubs("KDV_solid",cdc_rmin-1.0*cm, cdc_rmin, KDV_Length, 0.0, twopi);
#else // phi150mm
      G4double KDV_Length = (47.5*cm-13.0*cm)/2.0;//34.5cm
      G4double KDV_Zpos   = -(13.0*cm+KDV_Length);
      KDV_solid              = new G4Tubs("KDV_solid",tarCham_rmax-1.0*cm, tarCham_rmax, KDV_Length, 0.0, twopi);
#endif
      
      
      G4LogicalVolume* KDV_log = new G4LogicalVolume(KDV_solid, Scinti, "KDV_log", 0, 0, 0);
      G4PVPlacement* KDV_phys  = new G4PVPlacement(0,G4ThreeVector(cdsPos_x, cdsPos_y, KDV_Zpos), KDV_log, "KDV_phys", CDS_log,false,0);  
      

      //--- Visualization ---//
      G4VisAttributes *KDV_att = new G4VisAttributes(Red);
      KDV_att->SetForceWireframe(true);
      //KDV_att->SetForceSolid(true);
      KDV_log->SetVisAttributes(KDV_att);
    }
  //#endif

  const G4double CDS_AC_space = 10.0*cm;
  const G4double AC_STC_space = 2.0*cm;
  const G4double STC_BLC_space = 2.0*cm;
  const G4double BLC_BLC_space = 2.0*cm;
  
  if(DoAC == true)
    {  
      // ==============================================================
      // AC
      // ==============================================================
	
      //**********//
      //*** AC ***//
      //**********//
      G4double AC_block_xy = 24./2.*cm;
      G4double AC_block_z = 5.0/2.*cm;
      G4VSolid* AC_box = new G4Box("AC_box",AC_block_xy, AC_block_xy,AC_block_z);
      G4LogicalVolume* AC_log = new G4LogicalVolume(AC_box,Vacuum,"AC_log",0,0,0);
      G4double AC_blockPos_x = 0.0*m;
      G4double AC_blockPos_y = 0.0*m;
      G4double AC_blockPos_z = -(cds_z+CDS_AC_space+AC_block_z);
      G4PVPlacement* AC_phys = new G4PVPlacement(0, G4ThreeVector(AC_blockPos_x,AC_blockPos_y,AC_blockPos_z),AC_log,"AC",experimentalHall_log,false,0);

      //********************************************//
      //*** sAC ( counter element in AC volume ) ***//
      //********************************************//
      G4double AC_seg_x = AC_block_xy;
      G4double AC_seg_y = AC_block_xy;
      G4double AC_seg_z = AC_block_z;
      G4VSolid* sAC_box = new G4Box("sAC_box",AC_seg_x,AC_seg_y,AC_seg_z);
      G4LogicalVolume* sAC_log = new G4LogicalVolume(sAC_box, Air, //<-material is temporal
						     "sAC_log", 0,0,0);
      G4PVPlacement* sAC_phys = new G4PVPlacement(0,G4ThreeVector(0,0,0), sAC_log, "sAC_phys", AC_log, false, 0);

      //--- Visualization ---//
      AC_log->SetVisAttributes(G4VisAttributes::Invisible);
      G4VisAttributes *sAC_att = new G4VisAttributes(Pink);
      sAC_att->SetForceWireframe(true);
      sAC_log->SetVisAttributes(sAC_att);


      // ==============================================================
      // STC
      // ==============================================================
	
      //***********//
      //*** STC ***//
      //***********//
      G4double STC_block_xy = 24./2.*cm;
      G4double STC_block_z = 1.0/2.*cm;
      G4VSolid* STC_box = new G4Box("STC_box",STC_block_xy, STC_block_xy,STC_block_z);
      G4LogicalVolume* STC_log = new G4LogicalVolume(STC_box,Vacuum,"STC_log",0,0,0);
      G4double STC_blockPos_x = 0.0*m;
      G4double STC_blockPos_y = 0.0*m;
      G4double STC_blockPos_z = -(cds_z+CDS_AC_space+AC_block_z*2+AC_STC_space+STC_block_z);
      G4PVPlacement* STC_phys = new G4PVPlacement(0, G4ThreeVector(STC_blockPos_x,STC_blockPos_y,STC_blockPos_z),
						  STC_log,"STC",experimentalHall_log,false,0);
	
      //**********************************************//
      //*** sSTC ( counter element in STC volume ) ***//
      //**********************************************//
      G4double STC_seg_x = STC_block_xy;
      G4double STC_seg_y = STC_block_xy;
      G4double STC_seg_z = STC_block_z;
      G4VSolid* sSTC_box = new G4Box("sSTC_box",STC_seg_x,STC_seg_y,STC_seg_z);
      G4LogicalVolume* sSTC_log = new G4LogicalVolume(sSTC_box, Scinti,"sSTC_log", 0,0,0);
      G4PVPlacement* sSTC_phys = new G4PVPlacement(0, G4ThreeVector(0,0,0),sSTC_log, "sSTC_phys", STC_log, false, 0);
	
      //--- Visualization ---//
      STC_log->SetVisAttributes(G4VisAttributes::Invisible);
      G4VisAttributes *sSTC_att = new G4VisAttributes(Red);
      sSTC_att->SetForceWireframe(true);
      sSTC_log->SetVisAttributes(sSTC_att);
	
	
      // ==============================================================
      // BLC (Beam-Line drift Chamber)
      // ==============================================================

      //***********//
      //*** BLC ***//
      //***********//
      G4double BLC_block_xy = 24./2.*cm;
      G4double BLC_block_z = 5.0/2.*cm;
      G4VSolid* BLC_box = new G4Box("BLC_box",BLC_block_xy, BLC_block_xy,BLC_block_z);
      G4LogicalVolume* BLC_log[4] = {0,0,0,0};
      char log_name[100];
      for (G4int i=0; i<4; i++)
	{
	  sprintf(log_name, "BLC_log%d", i+1);
	  BLC_log[i] = new G4LogicalVolume(BLC_box,Vacuum,log_name,0,0,0);
	}
      G4double BLC_pos_x = 0.0*m;
      G4double BLC_pos_y = 0.0*m;
      G4double BLC_pos_z_base = -(cds_z+CDS_AC_space+AC_block_z*2+AC_STC_space+STC_block_z*2+STC_BLC_space+BLC_block_z);
      G4double BLC_pos_z[4] = {BLC_pos_z_base-BLC_block_z*6-BLC_BLC_space*3,
			       BLC_pos_z_base-BLC_block_z*4-BLC_BLC_space*2,
			       BLC_pos_z_base-BLC_block_z*2-BLC_BLC_space,
			       BLC_pos_z_base};

      G4PVPlacement* BLC_phys[4] = {0,0,0,0} ;
      char phys_name[100];
      for (G4int i=0; i<4; i++)
	{
	  sprintf(phys_name, "BLC_phys%d", i+1);
	  BLC_phys[i] = new G4PVPlacement(0, G4ThreeVector(BLC_pos_x,BLC_pos_y,BLC_pos_z[i]),
					  BLC_log[i], phys_name, experimentalHall_log, false, 0);
	}
      const G4double BLC_cell_size = 8.0*mm; // full width
      G4double BLC_seg_x = BLC_cell_size/2;
      G4double BLC_seg_y = BLC_block_xy;
      G4double BLC_seg_z = 0.000001/2*mm;
      G4VSolid* sBLC_box = new G4Box("sBLC_box",BLC_seg_x,BLC_seg_y,BLC_seg_z);
      G4LogicalVolume* sBLC_log = new G4LogicalVolume(sBLC_box, ChamberGas, "sBLC_log", 0,0,0);
      G4double BLC_layer_pos[4]={-1.5*cm,-0.5*cm,0.5*cm,1.5*cm};
      G4int BLC_layer_n[4]={30,29,30,29};
      G4double BLC_layer_rot[4]={0*deg,0*deg,90*deg,90*deg};
      G4double BLC_layer_sta[4]={-BLC_block_xy+BLC_cell_size/2,-BLC_block_xy+BLC_cell_size,
				 -BLC_block_xy+BLC_cell_size/2,-BLC_block_xy+BLC_cell_size};

      std::vector< std::vector< std::vector<G4PVPlacement*> > > sBLC_phys(4, std::vector<std::vector<G4PVPlacement*> >(4) );
      for (G4int i=0; i<4; i++)
	{
	  for (G4int j=0; j<4; j++)
	    {
	      G4int n=0;
	      sBLC_phys[i][j].resize(BLC_layer_n[j],nullptr);
	      for (G4int k=0; k<BLC_layer_n[j]; k++)
		{
		  sprintf(phys_name, "sBLC_phys%d%d%02d", i+1, j+1, k+1);
		  G4ThreeVector xyzCounter;
		  if (j<2)
		    xyzCounter=G4ThreeVector(BLC_layer_sta[j]+BLC_cell_size*(k),0,BLC_layer_pos[j]);
		  else
		    xyzCounter=G4ThreeVector(0,BLC_layer_sta[j]+BLC_cell_size*(k),BLC_layer_pos[j]);
		  G4RotationMatrix* rotCounterBLC = new G4RotationMatrix;
		  rotCounterBLC->rotateZ(BLC_layer_rot[j]);
		  G4Transform3D posCounterBLC(*rotCounterBLC, xyzCounter);
		  sBLC_phys[i][j][k] = new G4PVPlacement(posCounterBLC,sBLC_log, phys_name, BLC_log[i], false, n);
		  n++;
		}
	    }
	}


      //--- Visualization ---//
      G4VisAttributes *BLC_att = new G4VisAttributes(Green);
      BLC_att->SetForceWireframe(true);
      for (G4int i=0; i<4; i++)
	{
	  //BLC_log[i-1]->SetVisAttributes(BLC_att);
	  BLC_log[i]->SetVisAttributes(G4VisAttributes::Invisible);
	}
      G4VisAttributes *sBLC_att = new G4VisAttributes(Yellow);
      sBLC_att->SetForceWireframe(true);
      sBLC_log->SetVisAttributes(sBLC_att);      



	
    }


    G4cout << "KnuclDetectorConstruction completed" << G4endl;
	
    return experimentalHall_phys;
}


void KnuclDetectorConstruction::ConstructSDandField()
{ 


  // ==============================================================
  // Definition of sensitive detectors
  // ==============================================================

  // G4SDManager* SDman = G4SDManager::GetSDMpointer();
  // counterSD   = new KnuclCounterSD  ("counterSD",   this);
  // chamberSD   = new KnuclChamberSD  ("chamberSD",   this, ChamberSmear, ZVertexChamberCell*mm, CDSLength*mm, CDCResolution*mm);
  // cherenkovSD = new KnuclCherenkovSD("cherenkovSD", this);
  // SDman->AddNewDetector(counterSD);
  // SDman->AddNewDetector(chamberSD);
  // SDman->AddNewDetector(cherenkovSD)

  //spTOF_log->SetSensitiveDetector(counterSD);
  // for (G4int i=0; i<15; i++)
  //   {
  //     CDC_log[i]->SetSensitiveDetector(chamberSD);
  //     //CDC_log[i]-> SetUserLimits( new G4UserLimits(1.0*mm) );
  //   }
  //CDH_log->SetSensitiveDetector(counterSD);
  //CVC_log->SetSensitiveDetector(counterSD); 
  //ZVC_cell_log[0]->SetSensitiveDetector(chamberSD);
  //ZVC_cell_log[1]->SetSensitiveDetector(chamberSD);
	
  //ZVC_log[0][7]->SetSensitiveDetector(counterSD);

  //ZVC_cell_log[0] -> SetSensitiveDetector(chamberSD);
  //ZVC_cell_log[1] -> SetSensitiveDetector(chamberSD);

  //KDV_log->SetSensitiveDetector(counterSD); 
  //sAC_log->SetSensitiveDetector(counterSD);

  //sSTC_log->SetSensitiveDetector(counterSD);
  //sBLC_log->SetSensitiveDetector(chamberSD);

  //sSDC1_log->SetSensitiveDetector(chamberSD);
  //sSDC2_log->SetSensitiveDetector(chamberSD);

  // ==============================================================
  // magnetic field
  // ==============================================================

   G4double fCDSField    = -FieldInCDC   *tesla;  
   G4double fKuramaField =  FieldInKurama*tesla;  

   G4ThreeVector fKURAMA(0.0, fKuramaField, 0.0      );
   G4ThreeVector fCDC   (0.0, 0.0,          fCDSField);  
// #ifdef NOYET
//   fEmFieldSetupKURAMA = new KnuclFieldSetup(fKURAMA) ;
//   fEmFieldSetupCDC    = new KnuclFieldSetup(fCDC) ;
// #endif


   fMagneticField = new G4SolSimpleMagneticField();
   fFieldMgr = new G4FieldManager();
   fFieldMgr->SetDetectorField(fMagneticField);
   fFieldMgr->CreateChordFinder(fMagneticField);

   G4bool forceToAllDaughters = true;
   //KuramaAperture_log->SetFieldManager(fFieldMgr,forceToAllDaughters);
   CDS_log->SetFieldManager(fFieldMgr,forceToAllDaughters);
   CDS_endcap_log->SetFieldManager(fFieldMgr,forceToAllDaughters);

   G4AutoDelete::Register(fMagneticField);
   G4AutoDelete::Register(fFieldMgr);
    
   
}
