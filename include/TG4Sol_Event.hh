// ----------------------------------------------------- 
// Definition of the TG4Sol_Event class
// Created by C.Rappold (c.rappold@gsi.de)
//------------------------------------------------------

#ifndef TG4Sol_Event_header
#define TG4Sol_Event_header 1

#include <iostream>
#include <fstream>
#include <string>

#include "TObject.h"
#include <vector>

/********************************************************************/
class TG4Sol_Event:public TObject
{
public: // Primary decay vertex
  Int_t WasDecay;
  Int_t MotherTrackID;
  std::string MotherName;
  Double_t MotherMass;
  Double_t DecayTime;
  
  Double_t InteractionPoint_X;
  Double_t InteractionPoint_Y;
  Double_t InteractionPoint_Z;
  
  Double_t DecayVertex_X;
  Double_t DecayVertex_Y;
  Double_t DecayVertex_Z;
  
  Double_t MotherMomentumAtDecay_X;
  Double_t MotherMomentumAtDecay_Y;
  Double_t MotherMomentumAtDecay_Z;
  
  std::vector<std::string> DaughterNames;
  std::vector<Double_t> DaughterMasses;
  std::vector<Double_t> DaughterCharges;
  std::vector<Double_t> DaughterMomentums_X;
  std::vector<Double_t> DaughterMomentums_Y;
  std::vector<Double_t> DaughterMomentums_Z;
  std::vector<Int_t> DaughterTrackID;
  //
  // beam composition
  std::vector<std::string> BeamNames;
  std::vector<Double_t> BeamMasses;
  std::vector<Double_t> BeamCharges;
  std::vector<Double_t> BeamMomentums_X;
  std::vector<Double_t> BeamMomentums_Y;
  std::vector<Double_t> BeamMomentums_Z;
  std::vector<Int_t> BeamTrackID;
  //***********************************************************
          
public:
  TG4Sol_Event();
  ~TG4Sol_Event();
  //
  void Display();
  void Zero();
  /*------------------------------------------------------------*/
  ClassDef(TG4Sol_Event,1); // TG4Sol_Event
};

#endif
