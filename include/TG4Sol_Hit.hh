// ----------------------------------------------------- 
// Definition of the TG4Sol_Hit class
// Created by C.Rappold (c.rappold@gsi.de)
//------------------------------------------------------

#ifndef TG4SOL_HIT_H
#define TG4SOL_HIT_H

#include "TObject.h"
#include <string>

class TG4Sol_Hit : public TObject
{
  public:
  TG4Sol_Hit();
  virtual ~TG4Sol_Hit();

  virtual void Print(const Option_t* = "") const;
  
  Int_t TrackID;  
  Int_t LayerID;
  Double32_t HitPosX;
  Double32_t HitPosY;
  Double32_t HitPosZ;

  Double32_t ExitPosX; 
  Double32_t ExitPosY; 
  Double32_t ExitPosZ; 

  Double32_t MomX;   
  Double32_t MomY;   
  Double32_t MomZ;   
  Double32_t Mass;   

  Double32_t Energy; 
  Double32_t Time;   

  std::string Pname;
    
  ClassDef(TG4Sol_Hit,1)  
};

#endif // TG4SOL_HIT_H


