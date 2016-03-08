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
// $Id: HypHIPrimaryGeneratorAction.hh 68058 2013-03-13 14:47:43Z gcosmo $
// 
/// \file HypHIPrimaryGeneratorAction.hh
/// \brief Definition of the HypHIPrimaryGeneratorAction class

#ifndef HypHIPrimaryGeneratorAction_h
#define HypHIPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4SolConfig.hh"
#include "G4GenericMessenger.hh"

#include <fstream>

class G4ParticleGun;
class G4GenericMessenger;
class G4Event;
class G4ParticleDefinition;

/// The primary generator action class with particle gum.
///
/// It defines a single particle which hits the calorimeter 
/// perpendicular to the input face. The type of the particle
/// can be changed via the G4 build-in commands of G4ParticleGun class 
/// (see the macros provided with this example).

class HypHIPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  explicit HypHIPrimaryGeneratorAction(const G4SolConfig& conf);
  virtual ~HypHIPrimaryGeneratorAction();

  virtual void GeneratePrimaries(G4Event* event);
  int GetStatus() const {return status;}

private:
  void DefineCommands();

  int status;

  G4ParticleGun*  fParticleGun; // G4 particle gun
  G4GenericMessenger* fMessenger;
  G4String nameInputFile;
  std::ifstream InStream;
  G4ParticleDefinition* ConstParticle;

  G4double fPosX;
  G4double fPosY;
  G4double fPosZ;
  G4double fSpotSizeSigma;
  G4double fTargetSize;
  
  G4long nEvents;
  
  G4bool fRandomizePrimary[3] = {true,true,true};

  const G4SolConfig& Par;
  
  G4ParticleDefinition* GetParticle(const G4String& particleName);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif


