// -----------------------------------------------------
// Definition of the TTritonStar class
// Created by C.Rappold (c.rappold@gsi.de)
//------------------------------------------------------

#ifndef TTritonStar_h
#define TTritonStar_h 1

#include "G4ParticleDefinition.hh"
#include "G4ios.hh"
#include "globals.hh"

// ######################################################################
// ###                 TritonStar with decay vertex tracking              ###
// ######################################################################

class TTritonStar : public G4ParticleDefinition
{
private:
  static TTritonStar* theInstance;
  TTritonStar() {}
  ~TTritonStar() {}

public:
  static TTritonStar* Definition();
  static TTritonStar* TritonStarDefinition();
  static TTritonStar* TritonStar();
};

#endif
