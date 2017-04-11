// ------------------------------------------------------- 
// Implementation of the G4SolHit class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------------

#include "G4SolHit.hh"

#include "G4VVisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

G4ThreadLocal G4Allocator<G4SolHit>* G4SolHitAllocator;


G4SolHit::G4SolHit():TrackID(-1), 
		     HitPosX(-9999.) , HitPosY(-9999.) , HitPosZ(-9999.) , 
		     ExitPosX(-9999.), ExitPosY(-9999.), ExitPosZ(-9999.),
		     MomX(-9999.)    , MomY(-9999.)    , MomZ(-9999.)    ,Mass(-9999),
		     Energy(-9999.), Time(-9999.), TrackLength(-9999.), 
		     Pname(""),LayerID(-1)
			   
{}

G4SolHit::G4SolHit(G4int z):TrackID(-1), 
			    HitPosX(-9999.) , HitPosY(-9999.) , HitPosZ(-9999.) , 
			    ExitPosX(-9999.), ExitPosY(-9999.), ExitPosZ(-9999.),
			    MomX(-9999.)    , MomY(-9999.)    , MomZ(-9999.)    ,Mass(-9999),
			    Energy(-9999.), Time(-9999.), TrackLength(-9999.),
			    Pname(""),LayerID(z)
			   
{}

G4SolHit::G4SolHit(const G4SolHit& hit):G4VHit()
{
  TrackID = hit.TrackID; 
  HitPosX = hit.HitPosX; 
  HitPosY = hit.HitPosY; 
  HitPosZ = hit.HitPosZ; 
  ExitPosX = hit.ExitPosX; 
  ExitPosY = hit.ExitPosY; 
  ExitPosZ = hit.ExitPosZ;
  MomX = hit.MomX;
  MomY = hit.MomY; 
  MomZ = hit.MomZ;
  Mass = hit.Mass;
  Energy = hit.Energy;
  Time = hit.Time;
  TrackLength = hit.TrackLength;
  Pname = hit.Pname;
  LayerID = hit.LayerID;
  
}

const G4SolHit& G4SolHit::operator=(const G4SolHit& hit)
{
  TrackID = hit.TrackID; 
  HitPosX = hit.HitPosX; 
  HitPosY = hit.HitPosY; 
  HitPosZ = hit.HitPosZ; 
  ExitPosX = hit.ExitPosX; 
  ExitPosY = hit.ExitPosY; 
  ExitPosZ = hit.ExitPosZ;
  MomX = hit.MomX;
  MomY = hit.MomY; 
  MomZ = hit.MomZ;
  Mass = hit.Mass;
  Energy = hit.Energy;
  Time = hit.Time; 
  TrackLength = hit.TrackLength;
  Pname = hit.Pname;
  LayerID = hit.LayerID;
  
  return *this;
}

int G4SolHit::operator==(const G4SolHit& hit) const
{
  return LayerID==hit.LayerID && TrackID==hit.TrackID;
}

void G4SolHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
    {
      G4Circle circle(G4Point3D(HitPosX,HitPosY,HitPosZ));
      circle.SetScreenSize(2);
      circle.SetFillStyle(G4Circle::filled);
      G4Colour colour(1.,1.,0.);
      G4VisAttributes attribs(colour);
      circle.SetVisAttributes(attribs);
      pVVisManager->Draw(circle);
    }
}



G4SolHit::~G4SolHit() {}

void G4SolHit::Print()
{
  std::cout<<"The TG4SolHit #" <<TrackID<<std::endl;
  std::cout<<"Particle name:"<<Pname<<" Edep:"<<Energy<<" Time:"<<Time<<std::endl;
  std::cout<<"Position :"<<HitPosX<<" "<<HitPosY<<" "<<HitPosZ<<std::endl;
  std::cout<<"Exit Position :"<<ExitPosX<<" "<<ExitPosY<<" "<<ExitPosZ<<std::endl;
  std::cout<<"Momentum :"<<MomX<<" "<<MomY<<" "<<MomZ<<std::endl;
} 
