// ------------------------------------------------- 
// Definition of the G4SolHit class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------


#ifndef G4SOL_HIT_H
#define G4SOL_HIT_H

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"

//#include <string>

class G4SolHit : public G4VHit
{
  public:
  G4SolHit();
  explicit G4SolHit(G4int z);
  G4SolHit(const G4SolHit& hit);
  virtual ~G4SolHit();

  const G4SolHit& operator=(const G4SolHit& right);
  int operator==(const G4SolHit &right) const;

  inline void *operator new(size_t);
  inline void operator delete(void *aHit);

  virtual void Draw();
  // virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
  // virtual std::vector<G4AttValue>* CreateAttValues() const;
  virtual void Print();

  
  G4int TrackID;   ///< Track Id

  G4double HitPosX;
  G4double HitPosY;
  G4double HitPosZ;

  G4double ExitPosX; 
  G4double ExitPosY; 
  G4double ExitPosZ; 

  G4double MomX;   
  G4double MomY;   
  G4double MomZ;   
  G4double Mass;   

  G4double Energy; 
  G4double Time;   
  G4double TrackLength;
  
  G4String Pname;

  G4int LayerID; 

};
typedef G4THitsCollection<G4SolHit> G4SolHitsCollection;

extern G4ThreadLocal G4Allocator<G4SolHit>* G4SolHitAllocator;

inline void* G4SolHit::operator new(size_t)
{
    if (!G4SolHitAllocator)
        G4SolHitAllocator = new G4Allocator<G4SolHit>;
    return (void*)G4SolHitAllocator->MallocSingle();
}

inline void G4SolHit::operator delete(void* aHit)
{
    G4SolHitAllocator->FreeSingle((G4SolHit*) aHit);
}

#endif // CALHIT_H


