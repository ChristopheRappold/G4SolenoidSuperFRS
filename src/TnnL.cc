// ------------------------------------------------------------- 
// Implementation of the TnnL class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------------------

#include "TnnL.hh"
#include "G4ParticleTable.hh"
#include "G4DecayTable.hh"
//#include "THyperNuclDecChan.hh"
// // #include "G4Lambda.hh"
// #include "G4Triton.hh"
// #include "G4Alpha.hh"
// #include "G4He3.hh"
// #include "G4Proton.hh"
// #include "G4PionMinus.hh"
// #include "G4Triton.hh"
// #include "G4Deuteron.hh"
#include "G4PhaseSpaceDecayChannel.hh"

#include "TTritonStar.hh"

#include "G4SystemOfUnits.hh"

//#include "THypHi_ParCreator.hh"

/*******************************************************************/
TnnL* TnnL::theInstance = 0;

/*******************************************************************/
TnnL* TnnL::Definition()
{
  if (theInstance !=0) return theInstance;
  //
  G4int    lepton_number = 0;
  G4int    baryon_number = 3;
  G4double mass          = 2993.7*MeV;
  G4double width         = 0.0*GeV;
  G4double charge        = 0.0*eplus;
  G4int    spin          = 1;
  G4int    parity        = 1;
  G4int    C_conjugation = 1;
  G4int    Isospin       = 1;
  G4int    IsospinZ      = 0;
  G4int    G_parity      = 0;
  G4String pType         = "nucleus";
  G4int    encoding      = 0;
  G4bool   stable        = false;
  G4double lifetime      = 0.180*ns;
  G4bool   shortlived    = false;
  //
  //THypHi_Par *Par = THypHi_ParCreator::GetParameter();
  //lifetime = Par->Get_HyperNuclei_H3L_T12();
  //
  const G4String name = "nnL";
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
      // G4He3::He3Definition(); //27Oct2012
      // G4Proton::ProtonDefinition();
      // G4PionMinus::PionMinusDefinition();
      // G4Triton::TritonDefinition();
      // G4Deuteron::DeuteronDefinition();
      TTritonStar::TritonStarDefinition();
      G4DecayTable * decayTable = new G4DecayTable();

      G4VDecayChannel* mode_1 = new G4PhaseSpaceDecayChannel("nnL",.5,2,"pi-","triton*");
      G4VDecayChannel* mode_2 = new G4PhaseSpaceDecayChannel("nnL",.5,2,"pi-","triton");
      //mass_dpi = 2059.3*MeV;
      decayTable->Insert(mode_1);
      decayTable->Insert(mode_2);
      anInstance->SetDecayTable(decayTable);
    }
  theInstance = reinterpret_cast<TnnL*>(anInstance);
  return theInstance;
}

/*******************************************************************/
TnnL*  TnnL::nnLDefinition()
{
  return Definition();
}

/*******************************************************************/
TnnL*  TnnL::nnL()
{
  return Definition();
}


void TnnL::setBR(const std::vector<double>& vecBR)
{
  G4DecayTable* decayTable = theInstance->GetDecayTable();
  if(vecBR.size()==(unsigned int)(decayTable->entries()))
    {
      for(unsigned int idBR = 0 ; idBR < vecBR.size(); ++idBR)
	decayTable->GetDecayChannel(idBR)->SetBR(vecBR[idBR]);
    }
}

void TnnL::setLT(double LT)
{
  theInstance->SetPDGLifeTime(LT);
}
