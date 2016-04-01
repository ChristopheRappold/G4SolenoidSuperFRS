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
// Implementation of the G4SolSimpleMagneticField class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------------------

#include "G4SolSimpleMagneticField.hh"

#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4SolSimpleMagneticField::G4SolSimpleMagneticField() : G4MagneticField(), fMessenger(0), fBx(0.),fBy(0.),fBz(1.0*tesla)
{
    // define commands for this class
    DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4SolSimpleMagneticField::~G4SolSimpleMagneticField()
{ 
    delete fMessenger; 
}

void G4SolSimpleMagneticField::GetFieldValue(const G4double [4],double *bField) const
{
    bField[0] = fBx;
    bField[1] = fBy;
    bField[2] = fBz;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4SolSimpleMagneticField::DefineCommands()
{
    // Define /G4SolSimple/field command directory using generic messenger class
    fMessenger = new G4GenericMessenger(this, "/G4SolSimple/field/", "Field control");

    // fieldValue command 
    G4GenericMessenger::Command& valueCmd = fMessenger->DeclareMethodWithUnit("valueBx","tesla", &G4SolSimpleMagneticField::SetFieldBx, "Set fieldBx strength.");
    valueCmd.SetParameterName("fieldBx", true);
    valueCmd.SetDefaultValue("0.");

    // fieldValue command 
    G4GenericMessenger::Command& valueCmd2 = fMessenger->DeclareMethodWithUnit("valueBy","tesla", &G4SolSimpleMagneticField::SetFieldBy, "Set fieldBy strength.");
    valueCmd2.SetParameterName("fieldBy", true);
    valueCmd2.SetDefaultValue("0.");

    // fieldValue command 
    G4GenericMessenger::Command& valueCmd3 = fMessenger->DeclareMethodWithUnit("valueBz","tesla", &G4SolSimpleMagneticField::SetFieldBz, "Set fieldBz strength.");
    valueCmd3.SetParameterName("fieldBz", true);
    valueCmd3.SetDefaultValue(".5");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
