// --------------------------------------------------------- 
// Definition of the G4SolConvertGeo class
// Created by C.Rappold (c.rappold@gsi.de)
//----------------------------------------------------------

#ifndef G4SolConvertGeo_h
#define G4SolConvertGeo_h 1

#include "G4VPhysicalVolume.hh"

class G4SolConvertGeo
{
  G4VPhysicalVolume* physiWorld;
  
public :
  G4SolConvertGeo(G4VPhysicalVolume* w);
  ~G4SolConvertGeo();
  int Convert(const std::string& nameOut);
   
};


#endif
