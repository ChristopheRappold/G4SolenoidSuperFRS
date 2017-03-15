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
// Implementation of the G4SolActionInitialization class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------------

#include "G4SolActionInitialization.hh"
#include "HypHIPrimaryGeneratorAction.hh"
#include "G4SolSimplePrimaryGeneratorAction.hh"
#include "G4SolRunAction.hh"
#include "G4SolEventAction.hh"
#include "G4SolStackingAction.hh"

#include <memory>
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4SolActionInitialization::G4SolActionInitialization (G4SolGeometryController* geoControl, const G4SolConfig& ConfFile) : G4VUserActionInitialization(),fGeoController(geoControl),Conf(ConfFile)
{
  OutputFile = Conf.Get<std::string>("Output_Namefile");
  std::cout<<"ActionInit : done ! "<<OutputFile<<std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4SolActionInitialization::~G4SolActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SolActionInitialization::BuildForMaster() const
{
  std::vector<G4String> nameD = fGeoController->GetNameDetectors();
  std::cout<<"!> G4SolActionInitialization BuildForMaster:"<<nameD.size()<<" "<<fGeoController->GetNameDetectors().size()<<std::endl;
  
  SetUserAction(new G4SolRunAction(OutputFile, fGeoController->GetNameDetectors(),Conf));

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SolActionInitialization::Build() const
{
  std::vector<G4String> nameD = fGeoController->GetNameDetectors();
  std::cout<<"!> G4SolActionInitialization Build:"<<nameD.size()<<" "<<fGeoController->GetNameDetectors().size()<<std::endl;

  std::unique_ptr<HypHIPrimaryGeneratorAction> newPrimGen (new HypHIPrimaryGeneratorAction(Conf));//= std::make_unique<HypHIPrimaryGeneratorAction>(Conf);
  if(newPrimGen->GetStatus()!=0)
    SetUserAction(new G4SolSimplePrimaryGeneratorAction(Conf));
  else
    SetUserAction(newPrimGen.release());
    
  SetUserAction(new G4SolRunAction(OutputFile, fGeoController->GetNameDetectors(),Conf));
  
  G4SolEventAction* eventAction = new G4SolEventAction( fGeoController->GetNameDetectors());
  SetUserAction(eventAction);

  SetUserAction(new G4SolStackingAction());
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
