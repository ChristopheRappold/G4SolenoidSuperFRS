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


#include "KnuclDetectorConstruction.hh"
#include "G4SolSimpleActionInitialization.hh"
#include "G4SolActionInitialization.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "G4StepLimiterPhysics.hh"

#include "G4SolConfig.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{

  std::cout<<" Config :"<<std::endl;
  G4SolConfig config(argc,argv);
  if(config.ProperConf()!=0)
    return -1;
  config.CheckConfig();

  int guimode = config.Get<int>("Gui"); 
  
#ifdef G4UI_USE
  // Detect interactive mode (if no arguments) and define UI session
  //
  G4UIExecutive* ui = 0;
  if ( guimode )
    {
      ui = new G4UIExecutive(argc, argv);
    }
#endif

  // Construct the default run manager
  //
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
#else
  G4RunManager* runManager = new G4RunManager;
#endif

  // std::string nameP = config.Get<std::string>("Particle");
  // double sizeTarget = config.Get<double>("Target_Size");
  // std::string unit = config.Get<std::string>("Target_Size.unit");
  // std::string nameUnit("Target_Size.unit."+unit);

  // double unitVal = config.Get<double>(nameUnit);
  // std::cout<<" get:"<<nameP<<" "<<sizeTarget/unitVal<<" "<<unit<<" ("<<unitVal<<" "<<sizeTarget<<")"<<std::endl;
  // std::cout<<" done !"<<std::endl;

  KnuclDetectorConstruction* Geometry = new KnuclDetectorConstruction(config); 
  // Mandatory user initialization classes
  runManager->SetUserInitialization(Geometry);

  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  //physicsList->SetVerboseLevel(0);
  runManager->SetUserInitialization(physicsList);

  // User action initialization
  if(config.Get<bool>("SimpleGeo")==true)
    runManager->SetUserInitialization(new G4SolSimpleActionInitialization(config));
  else
    runManager->SetUserInitialization(new G4SolActionInitialization(Geometry,config));
 
#ifdef G4VIS_USE
  // Visualization manager construction
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif
    
  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  boost::optional<std::string> macro = config.Get<boost::optional<std::string> >("MacroFile"); 
  
  if ( macro )
    {
      // execute an argument macro file if exist
      G4String command = "/control/execute ";
      G4String fileName(*macro);
      UImanager->ApplyCommand(command+fileName);
    }
  else
    {
#ifdef G4UI_USE
#ifdef G4VIS_USE
      UImanager->ApplyCommand("/control/execute init_vis.mac"); 
#else
      UImanager->ApplyCommand("/control/execute init.mac"); 
#endif
      if (ui->IsGUI())
	{
	  UImanager->ApplyCommand("/control/execute gui.mac");
	}     
      // start interactive session
      ui->SessionStart();
      delete ui;
#endif
    }
  
  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
