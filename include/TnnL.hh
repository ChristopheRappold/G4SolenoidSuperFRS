// -----------------------------------------------------
// Definition of the TnnL class
// Created by C.Rappold (c.rappold@gsi.de)
//------------------------------------------------------

#ifndef TnnL_h
#define TnnL_h 1

#include "G4ParticleDefinition.hh"
#include "G4ios.hh"
#include "globals.hh"
//
class G4VDecayChannel;

class TnnL : public G4ParticleDefinition
{
public:
  static TnnL* Definition();
  static TnnL* nnLDefinition();
  static TnnL* nnL();
  void setBR(const std::vector<double>&);
  void setLT(double);

private:
  // G4double mass_dpi;// = 2059.3*MeV;
  TnnL() {}
  ~TnnL() {}
  //
  static TnnL* theInstance;
};

#endif
