// ------------------------------------------------- 
// Definition of the G4Sol_SD_Det class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------

#ifndef G4Sol_SensetiveD_H
#define G4Sol_SensetiveD_H 1

#include "G4VSensitiveDetector.hh"
#include "G4SolHit.hh"
#include <unordered_map>

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class G4Sol_SD_Det: public G4VSensitiveDetector
{
public:
  //
  G4Sol_SD_Det(const G4String& nameVol);
  virtual ~G4Sol_SD_Det();
  
  G4bool ProcessHits(G4Step*aStep,G4TouchableHistory* rohist);
  virtual void Initialize(G4HCofThisEvent*HCE);
  virtual void EndOfEvent(G4HCofThisEvent*);
  virtual void clear();

  G4SolHitsCollection* fHitsCollection;
  G4int fHCID;

  std::unordered_map< int, std::unordered_map<int,int> > mapTrackID_Hits;

};

// class G4Sol_SD_CDC: public G4VSensitiveDetector
// {
// public:
//   //
//   G4Sol_SD_CDC(const G4String& name, const G4String& nameVolPhysm);
//   virtual ~G4Sol_SD_CDC();
  
//   G4bool ProcessHits(G4Step*aStep,G4TouchableHistory* rohist);
//   virtual void Initialize(G4HCofThisEvent*HCE);
//   virtual void EndOfEvent(G4HCofThisEvent*);
//   virtual void clear();

//   G4SolHitsCollection* fHitsCollection;
//   G4int fHCID;

//   std::map<int,int> mapTrackID_Hits;

// };

// class G4Sol_SD_Cal: public G4VSensitiveDetector
// {
// public:
//   //
//   G4Sol_SD_Cal(const G4String& name, const G4String& nameVolPhysm);
//   virtual ~G4Sol_SD_Cal();
  
//   G4bool ProcessHits(G4Step*aStep,G4TouchableHistory* rohist);
//   virtual void Initialize(G4HCofThisEvent*HCE);
//   virtual void EndOfEvent(G4HCofThisEvent*);
//   virtual void clear();

//   G4SolHitsCollection* fHitsCollection;
//   G4int fHCID;

//   std::map<int,int> mapTrackID_Hits;

// };



#endif
