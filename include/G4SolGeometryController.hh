#ifndef G4SolGeometryController_hh
#define G4SolGeometryController_hh 1

#include "globals.hh"
#include "G4String.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4SolVDetectorConstruction.hh"

#include "G4SolConfig.hh"

/**
 * Controller for geometry selection
 *
 * This controller is called by the geometry messenger and used to
 * select the geometry. Each available geometry must have unique name
 * and it must be known by the geometry controller.
 */
class G4SolGeometryController
{
public:
  explicit G4SolGeometryController(const G4SolConfig& _par);
  ~G4SolGeometryController();

  /**
   * Select a geometry by name.
   */
  int SetGeometry(G4String);
  const std::vector<G4String>& GetNameDetectors();
  void ConvertG4toRoot(const std::string& nameConvertRoot);

private:
  void registerGeometry(G4VUserDetectorConstruction *det);

  const G4SolConfig& Par;
  G4String nameGeometry;
  G4SolVDetectorConstruction *detectorBuilder;
};


#endif
