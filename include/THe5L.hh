#ifndef THe5L_h
#define THe5L_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
//
class G4VDecayChannel;

class THe5L : public G4ParticleDefinition
{
public:
  static THe5L* Definition();
  static THe5L* He5LDefinition();
  static THe5L* He5L();
  void setBR(const std::vector<double>&);
  void setLT(double);
private:
  THe5L(){}
  ~THe5L(){}
  //
  static THe5L* theInstance;
};

#endif

