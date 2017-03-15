// --------------------------------------------------------- 
// Definition of the G4SolVDetectorConstruction class
// Created by C.Rappold (c.rappold@gsi.de)
//----------------------------------------------------------


#ifndef G4SolVDetectorConstruction_h
#define G4SolVDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "globals.hh"

#include "G4SolConfig.hh"


class G4VPhysicalVolume;


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

class G4SolVDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  explicit G4SolVDetectorConstruction(const G4SolConfig& _par);
  virtual ~G4SolVDetectorConstruction();

public:
  virtual G4VPhysicalVolume* Construct() = 0;
  virtual void ConstructSDandField() = 0;

  const std::vector<G4String>& GetNameDetectors() const {  return NameDetectorsSD; };

  G4LogicalVolume* experimentalHall_logOutRoot;
  G4VPhysicalVolume*   experimentalHall_physOutRoot;

protected:

  G4LogicalVolume* FindVolume(const G4String& name);
  const G4SolConfig& Par;
  std::vector<G4String> NameDetectorsSD;


};

     

#endif

