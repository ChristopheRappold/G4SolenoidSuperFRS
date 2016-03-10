#ifndef TH3L_h
#define TH3L_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
//
class G4VDecayChannel;

class TH3L : public G4ParticleDefinition{
public:
  static TH3L* Definition();
  static TH3L* H3LDefinition();
  static TH3L* H3L();

  void setBR(const std::vector<double>&);
  void setLT(double);
private:
  TH3L(){}
  ~TH3L(){}
  //
  static TH3L* theInstance;
};

#endif

