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

#ifndef G4SolSimplePrimaryGeneratorAction_h
#define G4SolSimplePrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4SolConfig.hh"

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
  G4double fSigmaMomentum;
  G4double fSigmaAngle;
  G4bool fRandomizePrimary[3] = {true,true,true};

  const G4SolConfig& Par;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
