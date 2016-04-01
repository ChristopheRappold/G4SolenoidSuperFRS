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
// ------------------------------------------------- 
// Definition of the G4SolSimpleMagneticField class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------

#ifndef G4SolSimpleMagneticField_H
#define G4SolSimpleMagneticField_H 1

#include "globals.hh"
#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"

class G4GenericMessenger;

/// Magnetic field

class G4SolSimpleMagneticField : public G4MagneticField
{
public:
  G4SolSimpleMagneticField();
  virtual ~G4SolSimpleMagneticField();
    
  virtual void GetFieldValue(const G4double point[4],double* bField ) const;
  
  void SetField(const G4ThreeVector& val) { fBx = val.x(); fBx = val.y(); fBx = val.z(); fField = val;}
  void SetFieldBx(double val) { fBx = val;}
  void SetFieldBy(double val) { fBy = val;}
  void SetFieldBz(double val) { fBz = val;}
  G4ThreeVector GetField() const { return fField; }
    
private:
  void DefineCommands();

  G4GenericMessenger* fMessenger;
  G4double fBx;
  G4double fBy;
  G4double fBz;
  G4ThreeVector fField;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
