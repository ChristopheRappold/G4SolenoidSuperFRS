#include "TTritonStar.hh"
#include "G4ParticleTable.hh"
#include "G4DecayTable.hh"
#include "G4PhaseSpaceDecayChannel.hh"
// // #include "G4Lambda.hh"
// #include "G4Triton.hh"
// #include "G4Alpha.hh"
// #include "G4He3.hh"
// #include "G4Proton.hh"
// #include "G4PionMinus.hh"
// #include "G4Triton.hh"
// #include "G4Deuteron.hh"
#include "G4SystemOfUnits.hh"

//#include "THypHi_ParCreator.hh"

/*******************************************************************/
TTritonStar* TTritonStar::theInstance = 0;

/*******************************************************************/
TTritonStar* TTritonStar::Definition()
{
  if (theInstance !=0) return theInstance;
  //
  G4int    lepton_number = 0;
  G4int    baryon_number = 3;
  G4double mass          = (1875.613+939.56563+0.1)*MeV;
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
  G4double lifetime      = 0.0*ns;
  G4bool   shortlived    = true;
  //
  //THypHi_Par *Par = THypHi_ParCreator::GetParameter();
  //lifetime = Par->Get_HyperNuclei_H3L_T12();
  //
  const G4String name = "triton*";
  // search in particle table
  G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* anInstance = pTable->FindParticle(name);
  if (anInstance ==0)
    {
      anInstance = new G4ParticleDefinition(name,mass,width,charge,
					    spin,parity,C_conjugation,
					    Isospin,IsospinZ,G_parity,
					    pType,lepton_number,baryon_number,encoding,
					    stable,lifetime,NULL,shortlived);
      //anInstance->SetAtomicNumber(1);
      //anInstance->SetAtomicMass(4);
      // Decay brunches

      G4DecayTable * decayTable = new G4DecayTable();
      G4VDecayChannel* mode_1 = new G4PhaseSpaceDecayChannel("triton*",1,2,"neutron","deuteron");
      decayTable->Insert(mode_1);
      anInstance->SetDecayTable(decayTable);
    }
  theInstance = reinterpret_cast<TTritonStar*>(anInstance);
  return theInstance;
}

/*******************************************************************/
TTritonStar*  TTritonStar::TritonStarDefinition()
{
  return Definition();
}

/*******************************************************************/
TTritonStar*  TTritonStar::TritonStar()
{
  return Definition();
}
