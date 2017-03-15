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
// ------------------------------------------------------- 
// Implementation of the G4SolRunData class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------------

#include "G4SolRunData.hh"
//#include "G4SolAnalysis.hh"

#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include "G4SystemOfUnits.hh"

#include "G4SolHit.hh"
#include "G4SolStackingAction.hh"

#include "TG4Sol_Hit.hh"

#include <map>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4SolRunData::G4SolRunData(const G4String& name) : G4Run(),namefile(name),fileOut(nullptr),Tree(nullptr),LookCheckFile(false),LookCheckTree(false),CloseDone(false)
{
  std::cout<<"G4SolRunData Ctr :"<<name<<" "<<namefile<<std::endl;
  if(LookCheckFile==false)
    {
      fileOut = new TFile(namefile,"RECREATE");
      LookCheckFile = true;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4SolRunData::~G4SolRunData()
{
  Close();
}

void G4SolRunData::Close()
{
  if(CloseDone==false)
    {
      fileOut->cd();
      Tree->Write();
      fileOut->Close();
      fileOut->Delete();
      CloseDone=true;
    }
}

void G4SolRunData::InitTree(const std::vector<G4String>& nameDet, const G4SolConfig& config)
{
  if(LookCheckTree==false)
    {
      if(LookCheckFile==false)
	{
	  fileOut = new TFile(namefile,"RECREATE");
	  LookCheckFile = true;
	}

      std::map<std::string, double> parameterToFile;
      parameterToFile.insert(std::make_pair("Target_Size",config.Get<double>("Target_Size")/cm));
      parameterToFile.insert(std::make_pair("Target_PosX",config.Get<double>("Target_PosX")/cm));
      parameterToFile.insert(std::make_pair("Target_PosY",config.Get<double>("Target_PosY")/cm));
      parameterToFile.insert(std::make_pair("Target_PosZ",config.Get<double>("Target_PosZ")/cm));
      parameterToFile.insert(std::make_pair("Field_CDS_Bz",config.Get<double>("Field_CDS_Bz")/tesla));

      fileOut->cd();
      fileOut->WriteObjectAny(&nameDet, "std::vector<std::string>", "nameDet");
      fileOut->WriteObjectAny(&parameterToFile, "std::map<std::string,double>", "simParameters");
      Tree = new TTree("G4Tree","Geant4 Tree");
      
      fEvent = new TG4Sol_Event; 
      Tree->Branch("TG4Sol_Event",&fEvent,3200,2);
      for(auto& nameBranch : nameDet)
	{
	  //std::cout<<" Branch Construction : "<<nameBranch<<std::endl;
	  TClonesArray* SubEvent_Layer = new TClonesArray("TG4Sol_Hit",20);
	  Tree->Branch(nameBranch,&SubEvent_Layer,3200,2);
	  addrCloneArray.push_back(SubEvent_Layer);
	}

      LookCheckTree = true;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SolRunData::FillPerEvent(const G4Event* event)
{
  // // get analysis manager
  // G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // //accumulate statistic
  // //

  // for ( G4int i=0; i<kDim; i++) {
  //   // fill histograms
  //   analysisManager->FillH1(i+1, fEdep[i]);
  //   analysisManager->FillH1(kDim+i+1, fTrackLength[i]);

  //   // fill ntuple
  //   analysisManager->FillNtupleDColumn(i, fEdep[i]);
  //   analysisManager->FillNtupleDColumn(kDim+i, fTrackLength[i]);
  // }  

  // analysisManager->AddNtupleRow();  
  std::vector<int> prim_id;
  
  unsigned int Nprimary = event->GetNumberOfPrimaryVertex();
  if(Nprimary > 0)
    {
      //G4cout<<" Write Event str :"<<Nprimary<<G4endl;
      for(unsigned int iPrimary = 0; iPrimary < Nprimary; ++iPrimary)
	{
	  G4PrimaryVertex* PrimVertex = event->GetPrimaryVertex(iPrimary);
	  
	  if(PrimVertex!=nullptr)
	    {
	      unsigned int Nparticle = PrimVertex->GetNumberOfParticle();
	      fEvent->InteractionPoint_X=PrimVertex->GetX0()/cm;
	      fEvent->InteractionPoint_Y=PrimVertex->GetY0()/cm;
	      fEvent->InteractionPoint_Z=PrimVertex->GetZ0()/cm;
	      
	      //G4cout<< "Vtx#"<<iPrimary<<" Nparticle :"<<Nparticle<<G4endl;

	      for(unsigned int iParticle = 0; iParticle<Nparticle;++iParticle)
		{
		  G4PrimaryParticle* PrimParticle = PrimVertex->GetPrimary(iParticle);
		  
		  if(PrimParticle!=nullptr)
		    {
		      //G4cout<<" Particle#"<<iParticle<<" "<<PrimParticle->GetParticleDefinition()->GetParticleName()
		      //<<" "<<PrimParticle->GetTotalEnergy()<<G4endl;
		      //G4cout<<"Daugther"<<PrimParticle->GetDaughter()<<G4endl;//->Print();
		      //PrimParticle->GetUserInformation();
		      fEvent->BeamNames.push_back(PrimParticle->GetParticleDefinition()->GetParticleName());
		      fEvent->BeamMasses.push_back(PrimParticle->GetMass()/GeV);
		      fEvent->BeamCharges.push_back(PrimParticle->GetCharge());
		      fEvent->BeamMomentums_X.push_back(PrimParticle->GetPx()/GeV);
		      fEvent->BeamMomentums_Y.push_back(PrimParticle->GetPy()/GeV);
		      fEvent->BeamMomentums_Z.push_back(PrimParticle->GetPz()/GeV);
		      fEvent->BeamTrackID.push_back(PrimParticle->GetTrackID());
		    }	  
		}
	    }
	}
    }

  const G4SolStackingAction* stackingUser = dynamic_cast<const G4SolStackingAction*> (G4RunManager::GetRunManager()->GetUserStackingAction());
  for(unsigned int idPr = 0 ; idPr < fEvent->BeamTrackID.size() ;++idPr)
    {
      if(stackingUser->Get_MotherInfo(fEvent->BeamTrackID[idPr]))
	prim_id.push_back(idPr);
    }

  if(prim_id.size()!=0)
    {
      if(prim_id.size()!=1)
	{
	  G4cout<<" G4RunData : too much mother decaying ! "<<G4endl;
	}
      else
	{
	  int MotherID = prim_id[0];
	  const G4SolStacking::Daugthers_Info tempDaugther = stackingUser->Get_DaugthersInfo(fEvent->BeamTrackID[MotherID]);
	  if(tempDaugther.name_daughters.size()!=0)
	    {
	      fEvent->WasDecay = 1;
	      fEvent->MotherTrackID = MotherID;
	      fEvent->MotherName = fEvent->BeamNames[MotherID];
	      fEvent->MotherMass = fEvent->BeamMasses[MotherID];
	      fEvent->DecayTime = tempDaugther.decaytime/ns;
	      fEvent->DecayVertex_X = tempDaugther.secondary_vertex.x()/cm;
	      fEvent->DecayVertex_Y = tempDaugther.secondary_vertex.y()/cm;
	      fEvent->DecayVertex_Z = tempDaugther.secondary_vertex.z()/cm;

	      for(unsigned int idDaug = 0 ; idDaug < tempDaugther.name_daughters.size(); ++idDaug)
		{
		  int tempTrackID = tempDaugther.trackID_daughters[idDaug];
		  bool isAChain = stackingUser->Get_MotherInfo(tempTrackID);
		  if(isAChain==false)
		    {
		      fEvent->DaughterNames.push_back(tempDaugther.name_daughters[idDaug]);
		      fEvent->DaughterMasses.push_back(tempDaugther.mass_daughters[idDaug]/GeV);
		      fEvent->DaughterCharges.push_back(tempDaugther.charge_daughters[idDaug]);
		      fEvent->DaughterMomentums_X.push_back(tempDaugther.mom_daughters[idDaug].x()/GeV);
		      fEvent->DaughterMomentums_Y.push_back(tempDaugther.mom_daughters[idDaug].y()/GeV);
		      fEvent->DaughterMomentums_Z.push_back(tempDaugther.mom_daughters[idDaug].z()/GeV);
		      fEvent->DaughterTrackID.push_back(tempDaugther.trackID_daughters[idDaug]);
		  
		      fEvent->MotherMomentumAtDecay_X += tempDaugther.mom_daughters[idDaug].x()/GeV;
		      fEvent->MotherMomentumAtDecay_Y += tempDaugther.mom_daughters[idDaug].y()/GeV;
		      fEvent->MotherMomentumAtDecay_Z += tempDaugther.mom_daughters[idDaug].z()/GeV;
		    }
		  else
		    {
		      const G4SolStacking::Daugthers_Info tempDaugther2 = stackingUser->Get_DaugthersInfo(tempTrackID);
		      for(unsigned int idDaug2 = 0 ; idDaug2 < tempDaugther2.name_daughters.size(); ++idDaug2)
			{
			  fEvent->DaughterNames.push_back(tempDaugther2.name_daughters[idDaug2]);
			  fEvent->DaughterMasses.push_back(tempDaugther2.mass_daughters[idDaug2]/GeV);
			  fEvent->DaughterCharges.push_back(tempDaugther2.charge_daughters[idDaug2]);
			  fEvent->DaughterMomentums_X.push_back(tempDaugther2.mom_daughters[idDaug2].x()/GeV);
			  fEvent->DaughterMomentums_Y.push_back(tempDaugther2.mom_daughters[idDaug2].y()/GeV);
			  fEvent->DaughterMomentums_Z.push_back(tempDaugther2.mom_daughters[idDaug2].z()/GeV);
			  fEvent->DaughterTrackID.push_back(tempDaugther2.trackID_daughters[idDaug2]);
			  
			  fEvent->MotherMomentumAtDecay_X += tempDaugther2.mom_daughters[idDaug2].x()/GeV;
			  fEvent->MotherMomentumAtDecay_Y += tempDaugther2.mom_daughters[idDaug2].y()/GeV;
			  fEvent->MotherMomentumAtDecay_Z += tempDaugther2.mom_daughters[idDaug2].z()/GeV;
			}
		    }
		}
	    }
	}
    }

  G4HCofThisEvent* hce = event->GetHCofThisEvent();
  if (!hce) 
    {
      G4ExceptionDescription msg;
      msg << "No hits collection of this event found." << G4endl; 
      G4Exception("B5EventAction::EndOfEventAction()","B5Code001", JustWarning, msg);
      return;
    }   

  for( unsigned int idCol = 0;idCol<addrCloneArray.size();++idCol)
    {
      G4SolHitsCollection* TempCol = dynamic_cast<G4SolHitsCollection*>(hce->GetHC(idCol));
      TClonesArray* TempCArray = addrCloneArray[idCol];
      if(TempCol != nullptr)
  	{
  	  G4int nhits = TempCol->entries();
	  //std::cout<<" Col "<<TempCol->GetName()<<" "<<nhits<<std::endl;
  	  for(G4int ihit = 0 ; ihit<nhits;++ihit)
  	    {
  	      G4SolHit* TempHit = (*TempCol)[ihit];
	      TG4Sol_Hit *RootHit = dynamic_cast<TG4Sol_Hit*>(TempCArray->ConstructedAt(TempCArray->GetEntries()));
  	    
	      RootHit->TrackID = TempHit->TrackID;
	      RootHit->LayerID = TempHit->LayerID;
	      RootHit->HitPosX = TempHit->HitPosX/cm; 
	      RootHit->HitPosY = TempHit->HitPosY/cm; 
	      RootHit->HitPosZ = TempHit->HitPosZ/cm; 
	      RootHit->ExitPosX = TempHit->ExitPosX/cm; 
	      RootHit->ExitPosY = TempHit->ExitPosY/cm; 
	      RootHit->ExitPosZ = TempHit->ExitPosZ/cm;
	      RootHit->MomX = TempHit->MomX/GeV;
	      RootHit->MomY = TempHit->MomY/GeV; 
	      RootHit->MomZ = TempHit->MomZ/GeV;
	      RootHit->Energy = TempHit->Energy/MeV;
	      RootHit->Time = TempHit->Time/ns; 
	      RootHit->Pname = TempHit->Pname;
	      RootHit->Mass = TempHit->Mass/GeV;
	    }
  	}
      else
  	{
  	  G4ExceptionDescription msg;
  	  msg << "Some of hits collections of this event not found." << G4endl; 
  	  G4Exception("B5EventAction::EndOfEventAction()","B5Code001", JustWarning, msg);
  	  return;
  	}
    }

  Tree->Fill();
  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SolRunData::Reset()
{ 
  fEvent->Zero();
  for( unsigned int i = 0; i<addrCloneArray.size(); ++i )
    addrCloneArray[i]->Clear("C");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
