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
// $Id: WasaDetectorConstruction.hh 75215 2013-10-29 16:07:06Z gcosmo $
//
/// \file WasaDetectorConstruction.hh
/// \brief Definition of the WasaDetectorConstruction class

#ifndef WasaDetectorConstruction_h
#define WasaDetectorConstruction_h 1

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

class WasaDetectorConstruction : public G4SolVDetectorConstruction
{
public:
  explicit WasaDetectorConstruction(G4SolConfig& conf);
  virtual ~WasaDetectorConstruction();

  virtual G4VPhysicalVolume* Construct() override;
  virtual void ConstructSDandField() override;

private:
  static G4ThreadLocal G4MagneticField* fMagneticField;
  static G4ThreadLocal G4FieldManager* fFieldMgr;

  G4ChordFinder* fChordFinder;
  G4Mag_UsualEqRhs* fEquation;
  G4MagIntegratorStepper* fStepper;

  const G4Colour Blue        = {0.0, 0.0, 1.0};
  const G4Colour Gray        = {0.5, 0.5, 0.5};
  const G4Colour Red         = {1.0, 0.0, 0.0};
  const G4Colour LightBlue   = {0.0, 0.0, 1.0, 0.7};
  const G4Colour Yellow      = {1.0, 1.0, 0.0};
  const G4Colour Purple      = {.5, 0.0, .5};
  const G4Colour LightPurple = {106. / 255., 27. / 255., 189. / 255.};
  const G4Colour Green       = {0.0, 1.0, 0.0};
  const G4Colour Orange      = {1.0, 0.647, 0};
  const G4Colour Pink        = {1.0, 0.753, 0.796};

  const G4Colour  color_x     = {1., 1., 0.8};
  const G4Colour  color_u     = {1., 0.8, 1.};
  const G4Colour  color_v     = {0.8, 1., 1.};
  
  // G4Colour  white   (1.0, 1.0, 1.0) ;
  // G4Colour  grey    (0.5, 0.5, 0.5) ;
  // G4Colour  lgrey   (.85, .85, .85) ;
  // G4Colour  red     (1.0, 0.0, 0.0) ;
  // G4Colour  blue    (0.0, 0.0, 1.0) ;
  // G4Colour  cyan    (0.0, 1.0, 1.0) ;
  // G4Colour  magenta (1.0, 0.0, 1.0) ;
  // G4Colour  yellow  (1.0, 1.0, 0.0) ;
  // G4Colour  orange  (.75, .55, 0.0) ;
  // G4Colour  lblue   (0.0, 0.0, .75) ;
  // G4Colour  lgreen  (0.0, .75, 0.0) ;
  // G4Colour  green   (0.0, 1.0, 0.0) ;
  // G4Colour  brown   (0.7, 0.4, 0.1) ;

  const G4Colour PaletteOrange[5] = {{254 / 255., 237 / 255., 222 / 255.},
                                     {253 / 255., 190 / 255., 133 / 255.},
                                     {253 / 255., 141 / 255., 060 / 255.},
                                     {230 / 255., 85 / 255., 13 / 255.},
                                     {166 / 255., 54 / 255., 3 / 255.}};
  const G4Colour PaletteGreen[5]  = {{237 / 255., 248 / 255., 233 / 255.},
                                    {186 / 255., 228 / 255., 179 / 255.},
                                    {116 / 255., 196 / 255., 118 / 255.},
                                    {49 / 255., 163 / 255., 84 / 255.},
                                    {0, 109 / 255., 44 / 255.}};
  const G4Colour PaletteRed[5]    = {{254 / 255., 229 / 255., 217 / 255.},
				     {252 / 255., 174 / 255., 145 / 255.},
				     {251 / 255., 106 / 255., 74 / 255.},
				     {222 / 255., 45 / 255., 38 / 255.},
				     {165 / 255., 15 / 255., 21 / 255.}};

  const G4Colour ColorCDC[17] = {
      PaletteOrange[0], PaletteOrange[0], PaletteOrange[1], PaletteOrange[1], PaletteOrange[1], PaletteOrange[2],
      PaletteOrange[2], PaletteOrange[2], PaletteOrange[3], PaletteOrange[3], PaletteOrange[3], PaletteOrange[3],
      PaletteOrange[4], PaletteOrange[4], PaletteOrange[4], PaletteOrange[4], PaletteOrange[4]};

  G4LogicalVolume* experimentalHall_log;
  G4VPhysicalVolume* experimentalHall_phys;

  G4LogicalVolume* MFLD_log;
  G4VPhysicalVolume* MFLD_phys;

  G4LogicalVolume* INNER_log;
  G4VPhysicalVolume* INNER_phys;

  G4ThreeVector transMFLD;
  G4RotationMatrix rotMFLD;

  std::vector<G4PVPlacement*> AllPlacements;

  G4VPhysicalVolume* FindVolPhys(const G4String& name);

  // methods
  //
  void DefineMaterials();
  G4VPhysicalVolume* DefineVolumes();

  G4bool fCheckOverlaps; // option to activate checking of volumes overlaps
};

#endif
