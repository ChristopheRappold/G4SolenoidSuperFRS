#include "G4SolGeometryController.hh"

#include "G4RunManager.hh"
#include "G4SolConvertGeo.hh"
#include "G4VUserParallelWorld.hh"
#include "KnuclDetectorConstruction.hh"
#include "WasaDetectorConstruction.hh"
#include "WasaFullRootConstruction.hh"
#include "WasaSimpleDetectorConstruction.hh"
#include "HIAFSimpleDetectorConstruction.hh"

G4SolGeometryController::G4SolGeometryController(G4SolConfig& _par) : Par(_par) {}

G4SolGeometryController::~G4SolGeometryController() {}

int G4SolGeometryController::SetGeometry(G4String nameGeo)
{
  G4cout << "Activating geometry " << nameGeo << G4endl;
  // if(name == "G4Sol")
  //   {
  //     detectorBuilder = new G4SolDetectorConstruction(Par);
  //     registerGeometry(detectorBuilder);
  //   }
  // else if(name == "Phase0")
  //   {
  //     detectorBuilder = new HypHIPhase0DetectorConstruction(Par);
  //     registerGeometry(detectorBuilder);
  //   }
  // else if(name == "Phase0.5")
  //   {
  //     detectorBuilder = new HypHIPhase0DetectorConstruction(Par);
  //     registerGeometry(detectorBuilder);
  //   }
  // else
  //   {
  //     G4cout <<"Unknown geometry: " << name << ". Geometry not changed." << G4endl;
  //   }
  if(nameGeo == "CDS")
    {
      detectorBuilder = new KnuclDetectorConstruction(Par);
      // Mandatory user initialization classes
      registerGeometry(detectorBuilder);
    }
  else if(nameGeo == "WasaSimple")
    {
      detectorBuilder = new WasaSimpleDetectorConstruction(Par);
      // Mandatory user initialization classes
      registerGeometry(detectorBuilder);
    }
  else if(nameGeo == "Wasa")
    {
      detectorBuilder = new WasaDetectorConstruction(Par);
      registerGeometry(detectorBuilder);
    }
  else if(nameGeo == "WasaFullRoot")
    {
      detectorBuilder = new WasaFullRootConstruction(Par);
      registerGeometry(detectorBuilder);
    }
  else if(nameGeo == "HIAFSimple")
    {
      detectorBuilder = new HIAFSimpleDetectorConstruction(Par);
      registerGeometry(detectorBuilder);
    }
  else
    {
      std::cout << "E> No Geometry selected !" << nameGeo << "\n";
      return -2;
    }

  nameGeometry = nameGeo;
  return 0;
}

void G4SolGeometryController::registerGeometry(G4VUserDetectorConstruction* detector)
{
  G4RunManager* runManager = G4RunManager::GetRunManager();
  runManager->SetUserInitialization(detector);
  runManager->GeometryHasBeenModified();
}

const std::vector<G4String>& G4SolGeometryController::GetNameDetectors() { return detectorBuilder->GetNameDetectors(); }

void G4SolGeometryController::ConvertG4toRoot(const std::string& nameConvertRoot)
{

  auto* logi = detectorBuilder->experimentalHall_physOutRoot->GetLogicalVolume();
  if(logi != nullptr)
    {
      std::cout << logi->GetName() << std::endl;
      G4SolConvertGeo convertor(detectorBuilder->experimentalHall_physOutRoot);
      convertor.Convert(nameConvertRoot, nameGeometry);
    }
  else
    {
      std::cout << "E> Convertion G4toRoot failed!\n";
      return;
    }
}
