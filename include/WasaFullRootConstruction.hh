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
// $Id: WasaFullRootConstruction.hh 75215 2013-10-29 16:07:06Z gcosmo $
//
/// \file WasaFullRootConstruction.hh
/// \brief Definition of the WasaFullRootConstruction class

#ifndef WasaFullRootConstruction_h
#define WasaFullRootConstruction_h 1

#include "G4SolSimpleMagneticField.hh"
#include "G4SolVDetectorConstruction.hh"
//#include "THypHi_Par.hh"
#include "G4ChordFinder.hh"
#include "G4ClassicalRK4.hh"
#include "G4Color.hh"
#include "G4FieldManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4NystromRK4.hh"
#include "G4SolConfig.hh"
#include "G4VPhysicalVolume.hh"

class G4VSolid;
class G4PVPlacement;
class G4VPhysicalVolume;
class G4Material;
class G4VSensitiveDetector;
class G4VisAttributes;
class G4GenericMessenger;

class G4SolSimpleMagneticField;
class G4SolWASAMapMagneticField;
class G4FieldManager;

/// Detector construction class to define materials and geometry.
/// The calorimeter is a box made of a given number of layers. A layer consists
/// of an absorber plate and of a detection gap. The layer is replicated.
///
/// Four parameters define the geometry of the calorimeter :
///
/// - the thickness of an absorber plate,
/// - the thickness of a gap,
/// - the number of layers,
/// - the transverse size of the calorimeter (the input face is a square).
///
/// In addition a transverse uniform magnetic field is defined
/// via G4GlobalMagFieldMessenger class.

class WasaFullRootConstruction : public G4SolVDetectorConstruction
{
public:
  explicit WasaFullRootConstruction(G4SolConfig& conf);
  virtual ~WasaFullRootConstruction();

  virtual G4VPhysicalVolume* Construct() override;
  virtual void ConstructSDandField() override;

private:
  static G4ThreadLocal G4MagneticField* fMagneticField;
  static G4ThreadLocal G4FieldManager* fFieldMgr;

  G4ChordFinder* fChordFinder;
  G4Mag_UsualEqRhs* fEquation;
  G4MagIntegratorStepper* fStepper;

  G4LogicalVolume* experimentalHall_log;
  G4VPhysicalVolume* experimentalHall_phys;

  G4LogicalVolume* experimentalHall_logOutRoot;
  G4VPhysicalVolume* experimentalHall_physOutRoot;

  G4VPhysicalVolume* FindVolPhys(const G4String& name);

  // methods
  //
  void DefineMaterials();
  G4VPhysicalVolume* DefineVolumes();

};

#endif
