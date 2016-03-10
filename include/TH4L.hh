#ifndef TH4L_h
#define TH4L_h 1

#include "globals.hh"
#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
//
class G4VDecayChannel;

class TH4L : public G4ParticleDefinition{
public:
  static TH4L* Definition();
  static TH4L* H4LDefinition();
  static TH4L* H4L();

  void setBR(const std::vector<double>&);
  void setLT(double);
private:
  TH4L(){}
  ~TH4L(){}
  //
  static TH4L* theInstance;
};

#endif

