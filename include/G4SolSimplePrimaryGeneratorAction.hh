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
// ----------------------------------------------------------- 
// Definition of the G4SolSimplePrimaryGeneratorAction class
// Created by C.Rappold (c.rappold@gsi.de)
//------------------------------------------------------------

#ifndef G4SolSimplePrimaryGeneratorAction_h
#define G4SolSimplePrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "Randomize.hh"

#include "G4SolConfig.hh"
#include "G4ThreeVector.hh"

#include <functional>
#include <unordered_map>
#include <tuple>

class G4ParticleGun;
class G4GenericMessenger;
class G4Event;
class G4ParticleDefinition;

/// Primary generator
///
/// A single particle is generated.
/// User can select 
/// - the initial momentum and angle
/// - the momentum and angle spreads
/// - random selection of a particle type from proton, kaon+, pi+, muon+, e+ 


class G4SolSimplePrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  explicit G4SolSimplePrimaryGeneratorAction(const G4SolConfig& conf);
  virtual ~G4SolSimplePrimaryGeneratorAction();
    
  virtual void GeneratePrimaries(G4Event*);

  void SetParticle(G4String& name) { nameParticle =  name;} 

  void SetMomentum(G4double val) { fMomentum = val; }
  G4double GetMomentum() const { return fMomentum; }

  void SetSigmaMomentum(G4double val) { fSigmaMomentum = val; fRandomizePrimary[1] = true; }
  G4double GetSigmaMomentum() const { return fSigmaMomentum; }

  void SetSigmaAngle(G4double val) { fSigmaAngle = val; fRandomizePrimary[2] = true; }
  G4double GetSigmaAngle() const { return fSigmaAngle; }

  void SetRandomize(G4bool val, int i) { fRandomizePrimary[i] = val; }
  G4bool GetRandomize(int i = 0) const { return fRandomizePrimary[i]; }

  G4ParticleDefinition* GetParticle(const G4String& particleName);

private:
  void DefineCommands();

  G4ParticleGun* fParticleGun;
  G4GenericMessenger* fMessenger;
  G4String nameParticle;
  G4ParticleDefinition* ConstParticle;
  G4double fMomentum;
  G4double fKineticE;
  G4double fSigmaMomentum;
  G4double fSigmaAngle;
  G4double fDirX;
  G4double fDirY;
  G4double fDirZ;
  G4ThreeVector fDir;
  G4double fDirXSigma;
  G4double fDirYSigma;
  G4double fDirZSigma;
  G4double fPosX;
  G4double fPosY;
  G4double fPosZ;
  G4double fSpotSizeSigma;
  G4double fTargetSize;

  G4double fDirTheta_1;
  G4double fDirTheta_2;

  G4double fDirPhi_1;
  G4double fDirPhi_2;

  G4int fRandomizePrimary[5] = {0,0,0,0,0}; // [0] : Position / [1] : Total Mom / [2] : Beam Dir / [3] : Beam Theta / [4] : Beam Phi 

  std::unordered_map<int, std::function<double(double, double)> > RandTable {
									     {1, [](double mean, double sigma) { return G4RandGauss::shoot(mean, sigma);} },
									     {2, [](double min, double max) { return G4RandFlat::shoot(min, max); } },
									     {3, [](double a, double) { return a;} },};

  G4bool setBinRand;
  G4long eventPerBin;
  G4long currentEvent;
  std::tuple<unsigned int,unsigned int> currentBin;
  G4int NbBin;
  std::vector<std::tuple<double,double> > ThetaBins;
  std::vector<std::tuple<double,double> > PhiBins;
  

  const G4SolConfig& Par;

  
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
