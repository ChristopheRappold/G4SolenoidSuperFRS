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
// -------------------------------------------------------------
// Implementation of the G4SolSimpleActionInitialization class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------------------

#include "G4SolSimpleActionInitialization.hh"

#include "G4SolSimplePrimaryGeneratorAction.hh"
//#include "G4SolSimpleRunAction.hh"
//#include "G4SolSimpleEventAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4SolSimpleActionInitialization::G4SolSimpleActionInitialization(const G4SolConfig& conf)
    : G4VUserActionInitialization(), Conf(conf)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4SolSimpleActionInitialization::~G4SolSimpleActionInitialization() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SolSimpleActionInitialization::BuildForMaster() const
{
  // G4SolSimpleEventAction* eventAction = 0;
  // SetUserAction(new G4SolSimpleRunAction(eventAction));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SolSimpleActionInitialization::Build() const
{
  SetUserAction(new G4SolSimplePrimaryGeneratorAction(Conf));

  // G4SolSimpleEventAction* eventAction = new G4SolSimpleEventAction;
  // SetUserAction(eventAction);

  // SetUserAction(new G4SolSimpleRunAction(eventAction));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
