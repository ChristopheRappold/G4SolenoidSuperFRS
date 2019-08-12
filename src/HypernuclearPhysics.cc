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
// -------------------------------------------------------------
// Implementation of the HypernuclearPhysics class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------------------

#include "HypernuclearPhysics.hh"

#include "TH3L.hh"
#include "TH4L.hh"
#include "THe4L.hh"
#include "THe5L.hh"
#include "TTritonStar.hh"
#include "TnnL.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HypernuclearPhysics::HypernuclearPhysics(const G4String& name, const G4SolConfig& _par)
    : G4VPhysicsConstructor(name), Par(_par)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HypernuclearPhysics::~HypernuclearPhysics() {}

void HypernuclearPhysics::ConstructParticle()
{
  TTritonStar::TritonStarDefinition();

  TH3L* h3l   = TH3L::H3LDefinition();
  TH4L* h4l   = TH4L::H4LDefinition();
  THe4L* he4l = THe4L::He4LDefinition();
  THe5L* he5l = THe5L::He5LDefinition();
  TnnL* nnL   = TnnL::nnLDefinition();

  std::vector<double> vecBR;
  vecBR.resize(5);
  vecBR[0] = Par.Get<double>("HyperNuclei_H4L_br_mode1");
  vecBR[1] = Par.Get<double>("HyperNuclei_H4L_br_mode2");
  vecBR[2] = Par.Get<double>("HyperNuclei_H4L_br_mode3");
  vecBR[3] = Par.Get<double>("HyperNuclei_H4L_br_mode4");
  vecBR[4] = Par.Get<double>("HyperNuclei_H4L_br_mode5");
  h4l->setBR(vecBR);
  h4l->setLT(Par.Get<double>("HyperNuclei_H4LT12"));

  vecBR.resize(5);
  vecBR[0] = Par.Get<double>("HyperNuclei_H3L_br_mode1");
  vecBR[1] = Par.Get<double>("HyperNuclei_H3L_br_mode2");
  vecBR[2] = Par.Get<double>("HyperNuclei_H3L_br_mode3");
  vecBR[3] = Par.Get<double>("HyperNuclei_H3L_br_mode4");
  vecBR[4] = Par.Get<double>("HyperNuclei_H3L_br_mode5");
  h3l->setBR(vecBR);
  h3l->setLT(Par.Get<double>("HyperNuclei_H3L_T12"));

  vecBR.resize(7);
  vecBR[0] = Par.Get<double>("HyperNuclei_He4L_br_mode1");
  vecBR[1] = Par.Get<double>("HyperNuclei_He4L_br_mode2");
  vecBR[2] = Par.Get<double>("HyperNuclei_He4L_br_mode3");
  vecBR[3] = Par.Get<double>("HyperNuclei_He4L_br_mode4");
  vecBR[4] = Par.Get<double>("HyperNuclei_He4L_br_mode5");
  vecBR[5] = Par.Get<double>("HyperNuclei_He4L_br_mode6");
  vecBR[6] = Par.Get<double>("HyperNuclei_He4L_br_mode7");
  he4l->setBR(vecBR);
  he4l->setLT(Par.Get<double>("HyperNuclei_He4L_T12"));

  vecBR.resize(2);
  vecBR[0] = Par.Get<double>("HyperNuclei_He5L_br_mode1");
  vecBR[1] = Par.Get<double>("HyperNuclei_He5L_br_mode2");
  he5l->setBR(vecBR);
  he5l->setLT(Par.Get<double>("HyperNuclei_He5L_T12"));

  vecBR.resize(2);
  vecBR[0] = Par.Get<double>("HyperNuclei_nnL_br_mode1");
  vecBR[1] = Par.Get<double>("HyperNuclei_nnL_br_mode2");
  nnL->setBR(vecBR);
  nnL->setLT(Par.Get<double>("HyperNuclei_nnL_T12"));
}

void HypernuclearPhysics::ConstructProcess() {}
