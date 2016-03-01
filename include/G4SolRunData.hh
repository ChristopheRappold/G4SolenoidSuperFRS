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
// $Id: G4SolRunData.hh 69223 2013-04-23 12:36:10Z gcosmo $
// 
/// \file G4SolRunData.hh
/// \brief Definition of the G4SolRunData class

#ifndef G4SolRunData_h
#define G4SolRunData_h 1

#include "G4Run.hh"
#include "globals.hh"
#include "G4Event.hh"

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"

#include "TG4Sol_Event.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

enum {
  kAbs = 0,
  kGap = 1,
  kDim = 2
};  

///  Run data class
///
/// It defines data members to hold the energy deposit and track lengths
/// of charged particles in Absober and Gap layers.
/// 
/// In order to reduce the number of data members a 2-dimensions array 
/// is introduced for each quantity:
/// - fEdep[], fTrackLength[].
///
/// The data are collected step by step in G4SolSteppingAction, and
/// the accumulated values are filled in histograms and entuple
/// event by event in B4EventAction.

class G4SolRunData : public G4Run
{
public:
  explicit G4SolRunData(const G4String& namefile);
  virtual ~G4SolRunData();

  void FillPerEvent(const G4Event* event);

  void Reset();
  void Close();
  void InitTree(const std::vector<G4String>& nameDet);

private:

  const G4String& namefile;
  TFile* fileOut;
  TTree* Tree;

  TG4Sol_Event* fEvent;
  std::vector<TClonesArray*> addrCloneArray;

  bool LookCheckFile;
  bool LookCheckTree;
  bool CloseDone;

};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

