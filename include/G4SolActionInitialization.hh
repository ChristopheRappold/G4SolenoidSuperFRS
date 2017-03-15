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
// 
//
// ------------------------------------------------- 
// Definition of the G4SolActionInitialization class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------

#ifndef G4SolActionInitialization_h
#define G4SolActionInitialization_h 1

#include "globals.hh"
#include "G4VUserActionInitialization.hh"
//#include "KnuclDetectorConstruction.hh"
//#include "G4SolVDetectorConstruction.hh"
#include "G4SolGeometryController.hh"
#include "G4SolConfig.hh"


//class KnuclDetectorConstruction;
//class G4SolVDetectorConstruction;
class G4SolGeometryController;
class THypHi_Par;

class G4SolActionInitialization : public G4VUserActionInitialization
{
public:
  G4SolActionInitialization(G4SolGeometryController*, const G4SolConfig&);
  virtual ~G4SolActionInitialization();
  
  virtual void BuildForMaster() const;
  virtual void Build() const;
 
private:
  //KnuclDetectorConstruction* fGeoController;
  //G4SolVDetectorConstruction* fGeoController;
  G4SolGeometryController* fGeoController;
  const G4SolConfig& Conf;
  std::string OutputFile;

};

#endif

    
