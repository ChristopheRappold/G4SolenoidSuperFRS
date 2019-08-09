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
//
// ------------------------------------------------- 
// Definition of the G4SolStackingAction class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------

#ifndef G4SolStackingAction_H
#define G4SolStackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"
#include "G4ThreeVector.hh"
#include "G4String.hh"

#include <vector>
#include <unordered_map>
#include <tuple>
#include <set>
#include <functional>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
namespace G4SolStacking
{
struct Daugthers_Info
{

  G4ThreeVector secondary_vertex;
  G4double decaytime = 0.0;
  std::vector<G4ThreeVector> mom_daughters;
  std::vector<G4double> mass_daughters;
  std::vector<G4double> charge_daughters;
  std::vector<G4String> name_daughters;
  std::vector<G4int> trackID_daughters;

  void Print() const
  {
    G4cout<<" Daugther :"<<G4endl;
    G4cout<<" Secondary vertex :"<<secondary_vertex<<G4endl;

    for(unsigned int i=0; i< name_daughters.size();++i)
      {
	G4cout<<" name :"<<name_daughters[i]<<G4endl;
	G4cout<<"  '--> mom"<<mom_daughters[i]<<G4endl;
      }
  }
};

// class HashVec
// {
//   std::hash<std::string> hashFunc;
// public :
//   size_t operator() (const std::vector<int>& a) const
//   {
//     std::string nameHash("");
//     for( auto i : a )
//       {
// 	std::string nameHashTemp = std::to_string(i);
// 	nameHash += nameHashTemp;
//       }
//     return hashFunc(nameHash);
//   }

// };

class HashTuple
{
  std::hash<int> hashFunc;
public :
  size_t operator() (const std::tuple<int,int>& a) const
  {
    return hashFunc(std::get<0,int>(a));
  }

};

}


class G4SolStackingAction : public G4UserStackingAction
{
public:
  G4SolStackingAction();
  virtual ~G4SolStackingAction();

public:
  virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
  virtual void NewStage();
  virtual void PrepareNewEvent();

  bool Get_MotherInfo(G4int ) const;
  const G4SolStacking::Daugthers_Info Get_DaugthersInfo(G4int ) const;
private:
  
  std::unordered_map<G4int,G4SolStacking::Daugthers_Info> mother_daugthersInfo;

  std::set<G4int> list_PDGSelected;
  std::unordered_map<std::tuple<int,int>,std::vector<std::tuple<double,double,double,double,double> >, G4SolStacking::HashTuple > optic_lines; // PID -> graph : pathLength, X, Y, angleXZ, angleYZ 

  
  const G4ThreeVector beamAxis;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
