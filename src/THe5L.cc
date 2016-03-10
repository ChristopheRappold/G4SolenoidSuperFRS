#include "THe5L.hh"
#include "G4ParticleTable.hh"
#include "G4DecayTable.hh"
#include "G4PhaseSpaceDecayChannel.hh"

#include "G4SystemOfUnits.hh"
// // #include "G4Lambda.hh"
// #include "G4Triton.hh"
// #include "G4Alpha.hh"
// #include "G4Proton.hh"
// #include "G4PionMinus.hh"
// #include "G4Triton.hh"

//#include "THypHi_ParCreator.hh"

/*******************************************************************/
THe5L* THe5L::theInstance = 0;

/*******************************************************************/
THe5L* THe5L::Definition()
{
  if (theInstance !=0) return theInstance;
  //
  G4int    lepton_number = 0;
  G4int    baryon_number = 5;
  G4double mass          = 4839.9*MeV;
  G4double width         = 0.0*GeV;
  G4double charge        = 2.0*eplus;
  G4int    spin          = 1;
  G4int    parity        = -1;
  G4int    C_conjugation = 1;
  G4int    Isospin       = 0;
  G4int    IsospinZ      = 0;
  G4int    G_parity      = 0;
  G4String pType         = "nucleus";
  G4int    encoding      = 0;
  G4bool   stable        = false;
  G4double lifetime      = 0.256*ns;
  G4bool   shortlived    = false;
  //
  //THypHi_Par *Par = THypHi_ParCreator::GetParameter();
  //lifetime = Par->Get_HyperNuclei_He5L_T12();
  //
  const G4String name = "He5L";
  // search in particle table
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance ==0)
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
      // G4Proton::ProtonDefinition();
      // G4PionMinus::PionMinusDefinition();
      // G4Triton::TritonDefinition();
      G4DecayTable * decayTable = new G4DecayTable();
      G4VDecayChannel* mode_1 = new G4PhaseSpaceDecayChannel("He5L",0.67,3,"pi-","proton","alpha");
      G4VDecayChannel* mode_2 = new G4PhaseSpaceDecayChannel("He5L",0.33,3,"pi0","neutron","alpha");
      decayTable->Insert(mode_1);
      decayTable->Insert(mode_2);
      anInstance->SetDecayTable(decayTable);
    }
  theInstance = reinterpret_cast<THe5L*>(anInstance);
  return theInstance;
}

/*******************************************************************/
THe5L*  THe5L::He5LDefinition()
{
  return Definition();
}

/*******************************************************************/
THe5L*  THe5L::He5L()
{
  return Definition();
}


void THe5L::setBR(const std::vector<double>& vecBR)
{
  G4DecayTable* decayTable = theInstance->GetDecayTable();
  if(vecBR.size()==(unsigned int)(decayTable->entries()))
    {
      for(unsigned int idBR = 0 ; idBR < vecBR.size(); ++idBR)
	decayTable->GetDecayChannel(idBR)->SetBR(vecBR[idBR]);
    }
}

void THe5L::setLT(double LT)
{
  theInstance->SetPDGLifeTime(LT);
}
