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
// Implementation of the G4SolRunAction class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------------

#include "G4SolRunAction.hh"
#include "G4SolRunData.hh"
//#include "G4SolAnalysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
//#include "G4UnitsTable.hh"
//#include "G4SystemOfUnits.hh"

//#include "G4SolSensetiveD.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4SolRunAction::G4SolRunAction(const G4String& name, const std::vector<G4String>& nameSD_Det) : G4UserRunAction(),OutputFileName(name),NameDetectorsSD(nameSD_Det)
{ 
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     
  std::cout<<"!> G4SolRunAction Ctr :"<<OutputFileName<<" "<<NameDetectorsSD.size()<<" "<<nameSD_Det.size()<<std::endl;
  for(auto& nameD : NameDetectorsSD)
    std::cout<<"-> "<<nameD<<std::endl;
  
  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in G4SolAnalysis.hh
  
  // G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // G4cout << "Using " << analysisManager->GetType() << G4endl;

  // // Create directories 
  // //analysisManager->SetHistoDirectoryName("histograms");
  // //analysisManager->SetNtupleDirectoryName("ntuple");
  // analysisManager->SetVerboseLevel(1);
  // analysisManager->SetFirstHistoId(1);

  // // Book histograms, ntuple
  // //
  
  // // Creating histograms
  // analysisManager->CreateH1("1","Edep in absorber", 100, 0., 800*MeV);
  // analysisManager->CreateH1("2","Edep in gap", 100, 0., 100*MeV);
  // analysisManager->CreateH1("3","trackL in absorber", 100, 0., 1*m);
  // analysisManager->CreateH1("4","trackL in gap", 100, 0., 50*cm);

  // // Creating ntuple
  // //
  // analysisManager->CreateNtuple("G4Sol", "Edep and TrackL");

  // analysisManager->CreateNtupleIColumn("EventID");
  // analysisManager->CreateNtupleIColumn("LayerID");
  // analysisManager->CreateNtupleIColumn("HitID");
  // analysisManager->CreateNtupleIColumn("TrackID");
  // analysisManager->CreateNtupleDColumn("posX");
  // analysisManager->CreateNtupleDColumn("posY");
  // analysisManager->CreateNtupleDColumn("posZ");
  // analysisManager->CreateNtupleDColumn("momX");
  // analysisManager->CreateNtupleDColumn("momY");
  // analysisManager->CreateNtupleDColumn("momZ");
  // analysisManager->CreateNtupleDColumn("mass");

  // analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// void G4SolRunAction::Close()
// {
//   OutputFile->cd();
//   OutTree->Write();
//   OutputFile->Close();
//   OutputFile->Delete();
//   OutputFile=0;
// }

// void G4SolRunAction::Fill() const
// {
//   int ret = OutTree->Fill();
//   G4cout<<" Fill Tree "<<ret<<G4endl;

// }

G4SolRunAction::~G4SolRunAction()
{
  // if(OutputFile!=nullptr)
  //   Close();

  //  delete G4AnalysisManager::Instance();  
}

G4Run* G4SolRunAction::GenerateRun()
{
  return (new G4SolRunData(OutputFileName));
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SolRunAction::BeginOfRunAction(const G4Run*)
{ 
  G4SolRunData* hyprun = dynamic_cast<G4SolRunData*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  
  //const G4SolRunData* hyprun = dynamic_cast<const G4SolRunData*>(run);
  hyprun->InitTree(NameDetectorsSD);

  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Get analysis manager
  // G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // // Open an output file
  // //
  // G4String fileName = "G4SolOut";
  // analysisManager->OpenFile(fileName);

  // G4SDManager *SDman = G4SDManager::GetSDMpointer();

  // OutputFile = new TFile(OutputFileName,"RECREATE");
  // OutputFile->cd();

  // OutTree = new TTree("G4Tree","Geant4 Simulation Output Tree");
  
  // h1 = new TH1F("h_field","h_field",1000,0,1000);
  // h1->SetDirectory(OutputFile);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SolRunAction::EndOfRunAction(const G4Run*)
{

  G4cout<<" End Run : Closing root file ";
  // print histogram statistics
  //
  G4SolRunData* hyprun = dynamic_cast<G4SolRunData*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  //const G4SolRunData* hyprun = dynamic_cast<const G4SolRunData*>(run);
  hyprun->Close();

  //G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  // if ( analysisManager->GetH1(1) ) {
  //   G4cout << "\n ----> print histograms statistic ";
  //   if(isMaster) {
  //     G4cout << "for the entire run \n" << G4endl; 
  //   }
  //   else {
  //     G4cout << "for the local thread \n" << G4endl; 
  //   }
    
  //   G4cout << " EAbs : mean = " 
  //      << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy") 
  //      << " rms = " 
  //      << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") << G4endl;
    
  //   G4cout << " EGap : mean = " 
  //      << G4BestUnit(analysisManager->GetH1(2)->mean(), "Energy") 
  //      << " rms = " 
  //      << G4BestUnit(analysisManager->GetH1(2)->rms(),  "Energy") << G4endl;
    
  //   G4cout << " LAbs : mean = " 
  //     << G4BestUnit(analysisManager->GetH1(3)->mean(), "Length") 
  //     << " rms = " 
  //     << G4BestUnit(analysisManager->GetH1(3)->rms(),  "Length") << G4endl;

  //   G4cout << " LGap : mean = " 
  //     << G4BestUnit(analysisManager->GetH1(4)->mean(), "Length") 
  //     << " rms = " 
  //     << G4BestUnit(analysisManager->GetH1(4)->rms(),  "Length") << G4endl;
  // }

  // // save histograms & ntuple
  // //
  // analysisManager->Write();
  // analysisManager->CloseFile();

  //Close();
  
  G4cout<<" done !"<<G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
