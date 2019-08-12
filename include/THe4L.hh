// -----------------------------------------------------
// Definition of the THe4L class
// Created by C.Rappold (c.rappold@gsi.de)
//------------------------------------------------------

#ifndef THe4L_h
#define THe4L_h 1

#include "G4ParticleDefinition.hh"
#include "G4ios.hh"
#include "globals.hh"
//
class G4VDecayChannel;

class THe4L : public G4ParticleDefinition
{
public:
  static THe4L* Definition();
  static THe4L* He4LDefinition();
  static THe4L* He4L();
  void setBR(const std::vector<double>&);
  void setLT(double);

private:
  THe4L() {}
  ~THe4L() {}
  //
  static THe4L* theInstance;
};

#endif
