// -------------------------------------------------------------
// Implementation of the THe4L class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------------------

#include "THe4L.hh"

#include "G4DecayTable.hh"
#include "G4ParticleTable.hh"
#include "G4PhaseSpaceDecayChannel.hh"
// #include "G4Lambda.hh"
// #include "G4Triton.hh"
// #include "G4Alpha.hh"
// #include "G4Proton.hh"
// #include "G4PionMinus.hh"
// #include "G4Triton.hh"
// #include "G4Deuteron.hh"

#include "G4SystemOfUnits.hh"

//#include "THypHi_ParCreator.hh"

/*******************************************************************/
THe4L* THe4L::theInstance = 0;

/*******************************************************************/
THe4L* THe4L::Definition()
{
  if(theInstance != 0)
    return theInstance;
  //
  G4int lepton_number = 0;
  G4int baryon_number = 4;
  G4double mass       = 3921.19 * MeV;
  G4double width      = 0.0 * GeV;
  G4double charge     = 2.0 * eplus;
  G4int spin          = 0;
  G4int parity        = 1;
  G4int C_conjugation = 1;
  G4int Isospin       = 0;
  G4int IsospinZ      = 0;
  G4int G_parity      = 0;
  G4String pType      = "nucleus";
  G4int encoding      = 0;
  G4bool stable       = false;
  G4double lifetime   = 0.256 * ns;
  G4bool shortlived   = false;
  //
  // THypHi_Par *Par = THypHi_ParCreator::GetParameter();
  // lifetime = Par->Get_HyperNuclei_He4L_T12();
  //
  const G4String name = "He4L";
  // search in particle table
  G4ParticleTable* pTable          = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if(anInstance == 0)
    {
      anInstance = new G4ParticleDefinition(name, mass, width, charge, spin, parity, C_conjugation, Isospin, IsospinZ,
                                            G_parity, pType, lepton_number, baryon_number, encoding, stable, lifetime,
                                            NULL, shortlived, "hypernucleus");
      // anInstance->SetAtomicNumber(1);
      // anInstance->SetAtomicMass(4);
      // Decay brunches
      // G4Alpha::AlphaDefinition();
      // G4Proton::ProtonDefinition();
      // G4PionMinus::PionMinusDefinition();
      // G4Triton::TritonDefinition();
      // G4Deuteron::DeuteronDefinition();
      G4DecayTable* decayTable = new G4DecayTable();

      G4VDecayChannel* mode_1 = new G4PhaseSpaceDecayChannel("He4L", 0.56 * 0.690, 2, "pi0", "alpha");
      G4VDecayChannel* mode_2 = new G4PhaseSpaceDecayChannel("He4L", 0.27 * 0.840, 3, "pi-", "Z2A3", "proton");
      G4VDecayChannel* mode_3 = new G4PhaseSpaceDecayChannel("He4L", 0.56 * 0.260, 3, "pi0", "proton", "triton");
      G4VDecayChannel* mode_4 = new G4PhaseSpaceDecayChannel("He4L", 0.190 * 0.340, 2, "deuteron", "deuteron");
      G4VDecayChannel* mode_5 = new G4PhaseSpaceDecayChannel("He4L", 0.27 * 0.160, 4, "pi-", "proton", "proton", "deuteron");
      G4VDecayChannel* mode_6 = new G4PhaseSpaceDecayChannel("He4L", 0.56 * 0.049, 3, "pi0", "deuteron", "deuteron");
      G4VDecayChannel* mode_7 = new G4PhaseSpaceDecayChannel("He4L", 0.56 * 0.001, 3, "pi0", "Z2A3", "neutron");
      //---------
      //-------
      decayTable->Insert(mode_1);
      decayTable->Insert(mode_2);
      decayTable->Insert(mode_3);
      decayTable->Insert(mode_4);
      decayTable->Insert(mode_5);
      decayTable->Insert(mode_6);
      decayTable->Insert(mode_7);

      anInstance->SetDecayTable(decayTable);
    }
  theInstance = reinterpret_cast<THe4L*>(anInstance);
  return theInstance;
}

/*******************************************************************/
THe4L* THe4L::He4LDefinition() { return Definition(); }

/*******************************************************************/
THe4L* THe4L::He4L() { return Definition(); }

void THe4L::setBR(const std::vector<double>& vecBR)
{
  G4DecayTable* decayTable = theInstance->GetDecayTable();
  if(vecBR.size() == (unsigned int)(decayTable->entries()))
    {
      for(unsigned int idBR = 0; idBR < vecBR.size(); ++idBR)
        decayTable->GetDecayChannel(idBR)->SetBR(vecBR[idBR]);
    }
}

void THe4L::setLT(double LT) { theInstance->SetPDGLifeTime(LT); }
