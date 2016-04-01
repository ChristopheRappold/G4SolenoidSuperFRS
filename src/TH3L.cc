// ------------------------------------------------------------- 
// Implementation of the TH3L class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------------------

#include "TH3L.hh"
#include "G4ParticleTable.hh"
#include "G4DecayTable.hh"
//#include "G4PhaseSpaceDecayChannel.hh"
// # include "G4Lambda.hh"
// #include "G4Triton.hh"
// #include "G4Alpha.hh"
// #include "G4He3.hh"
// #include "G4Proton.hh"
// #include "G4PionMinus.hh"
// #include "G4Triton.hh"
// #include "G4Deuteron.hh"

#include "G4PhaseSpaceDecayChannel.hh"

#include "G4SystemOfUnits.hh"

//#include "THypHi_ParCreator.hh"

/*******************************************************************/
TH3L* TH3L::theInstance = 0;

/*******************************************************************/
TH3L* TH3L::Definition()
{
  if (theInstance !=nullptr) return theInstance;
  //
  G4int    lepton_number = 0;
  G4int    baryon_number = 3;
  G4double mass          = 2991.14*MeV;
  G4double width         = 0.0*GeV;
  G4double charge        = 1.0*eplus;
  G4int    spin          = 1;
  G4int    parity        = 1;
  G4int    C_conjugation = 1;
  G4int    Isospin       = 0;
  G4int    IsospinZ      = 0;
  G4int    G_parity      = 0;
  G4String pType         = "nucleus";
  G4int    encoding      = 0;
  G4bool   stable        = false;
  G4double lifetime      = 0.246*ns;
  G4bool   shortlived    = false;
  //

  const G4String name = "H3L";
  // search in particle table
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance==nullptr)
    {
      anInstance = new G4ParticleDefinition(name,mass,width,charge,
					    spin,parity,C_conjugation,
					    Isospin,IsospinZ,G_parity,
					    pType,lepton_number,baryon_number,encoding,
					    stable,lifetime,NULL,
					    shortlived,"hypernucleus");
    //anInstance->SetAtomicNumber(1);
    //anInstance->SetAtomicMass(4);
    // Decay brunches
    // G4Alpha::AlphaDefinition();
    // G4He3::He3Definition(); //27Oct2012
    // G4Proton::ProtonDefinition();
    // G4PionMinus::PionMinusDefinition();
    // G4Triton::TritonDefinition();
    // G4Deuteron::DeuteronDefinition();
    G4DecayTable * decayTable = new G4DecayTable();
    G4VDecayChannel* mode_1 = new G4PhaseSpaceDecayChannel("H3L",0.40891,3,"pi-","proton","deuteron");
    G4VDecayChannel* mode_2 = new G4PhaseSpaceDecayChannel("H3L",0.25349,2,"pi-","He3"); //"pi-","Z2A3"); //27Oct2012
    G4VDecayChannel* mode_3 = new G4PhaseSpaceDecayChannel("H3L",0.20446,3,"pi0","deuteron","neutron");
    G4VDecayChannel* mode_4 = new G4PhaseSpaceDecayChannel("H3L",0.12674,2,"pi0","triton");
    G4VDecayChannel* mode_5 = new G4PhaseSpaceDecayChannel("H3L",0.00640,4,"pi-","neutron","proton","proton");
    decayTable->Insert(mode_1);
    decayTable->Insert(mode_2);
    decayTable->Insert(mode_3);
    decayTable->Insert(mode_4);
    decayTable->Insert(mode_5);
    anInstance->SetDecayTable(decayTable);
    }
  theInstance = reinterpret_cast<TH3L*>(anInstance);
  return theInstance;
}

/*******************************************************************/
TH3L*  TH3L::H3LDefinition()
{
  return Definition();
}

/*******************************************************************/
TH3L*  TH3L::H3L()
{
  return Definition();
}

void TH3L::setBR(const std::vector<double>& vecBR)
{
  G4DecayTable* decayTable = theInstance->GetDecayTable();
  if(vecBR.size()==(unsigned int)(decayTable->entries()))
    {
      for(unsigned int idBR = 0 ; idBR < vecBR.size(); ++idBR)
	decayTable->GetDecayChannel(idBR)->SetBR(vecBR[idBR]);
    }
}

void TH3L::setLT(double LT)
{
  theInstance->SetPDGLifeTime(LT);
}
