// ------------------------------------------------------------- 
// Implementation of the TG4Sol_Hit class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------------------

#include "TG4Sol_Hit.hh"
#include <iostream>

TG4Sol_Hit::TG4Sol_Hit():TrackID(-1), LayerID(-1),
			 HitPosX(-9999.) , HitPosY(-9999.) , HitPosZ(-9999.) , 
			 ExitPosX(-9999.), ExitPosY(-9999.), ExitPosZ(-9999.),
			 MomX(-9999.)    , MomY(-9999.)    , MomZ(-9999.)    ,Mass(-9999.),
			 Energy(-9999.), Time(-9999.), TrackLength(-9999.),
			 Pname("")
			   
{}

TG4Sol_Hit::~TG4Sol_Hit() {}

void TG4Sol_Hit::Print(const Option_t*) const
{
  std::cout<<"The TTG4Sol_Hit #" <<TrackID<<" Layer:"<<LayerID<<std::endl;
  std::cout<<"Particle name:"<<Pname<<" Edep:"<<Energy<<" Time:"<<Time<<std::endl;
  std::cout<<"Position :"<<HitPosX<<" "<<HitPosY<<" "<<HitPosZ<<std::endl;
  std::cout<<"Exit Position :"<<ExitPosX<<" "<<ExitPosY<<" "<<ExitPosZ<<std::endl;
  std::cout<<"Momentum :"<<MomX<<" "<<MomY<<" "<<MomZ<<std::endl;
} 
