// -------------------------------------------------------
// Implementation of the G4Sol_SD_Det class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------------

#include "G4SolSensitiveD.hh"

#include "G4Material.hh"
#include "G4ParticleTypes.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4VSolid.hh"
#include "G4VTouchable.hh"

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

G4Sol_SD_Det::G4Sol_SD_Det(const G4String& Vol_name)
    : G4VSensitiveDetector(Vol_name), fHitsCollection(nullptr), fHCID(-1)
{
  collectionName.insert("G4SolColl");
}

G4Sol_SD_Det::~G4Sol_SD_Det() {}

void G4Sol_SD_Det::Initialize(G4HCofThisEvent* hce)
{
  fHitsCollection = new G4SolHitsCollection(SensitiveDetectorName, collectionName[0]);
  if(fHCID < 0)
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);

  hce->AddHitsCollection(fHCID, fHitsCollection);
  if(mapTrackID_Hits.size() != 0)
    mapTrackID_Hits.clear();
}

void G4Sol_SD_Det::clear() { mapTrackID_Hits.clear(); }

void G4Sol_SD_Det::EndOfEvent(G4HCofThisEvent*)
{
  //   Hits->Clear("C");
  mapTrackID_Hits.clear();
}

G4bool G4Sol_SD_Det::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{

  G4double energ_depos = aStep->GetTotalEnergyDeposit();
  //
  // std::cout<<" Current SD :"<<SensitiveDetectorName<<" "<<energ_depos<<std::endl;
  if(energ_depos < 1e-4)
    return true;

  G4String PhysName(aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName());
  // std::cout<<" ->"<<PhysName<<std::endl;

  // if(PhysName!=SensitiveDetectorName)
  //   return true;

  const G4TouchableHistory* touchable =
      dynamic_cast<const G4TouchableHistory*>(aStep->GetPreStepPoint()->GetTouchable());
  G4int idVolume = touchable->GetHistoryDepth() > 3
                       ? 1
                       : 0; // more than 3 meaning it is Central drifttube | less or equal to 3 meaning the rest

  // G4VPhysicalVolume* motherPhysical = touchable->GetVolume(1); // mother
  G4VPhysicalVolume* currentPhysical = touchable->GetVolume(idVolume);
  // G4int copyNo1 = motherPhysical->GetCopyNo();
  G4int copyNo = currentPhysical->GetCopyNo();

  // G4int replicat1 = touchable->GetReplicaNumber(1);
  // G4int replicat = touchable->GetReplicaNumber(0);

  int CurrentTrack = aStep->GetTrack()->GetTrackID();

  // std::cout<<" -> mother :"<<motherPhysical->GetName()<<" "<<copyNo1<<" / "<<replicat1<<" | current
  // :"<<currentPhysical->GetName()<<" "<<copyNo<<" / "<<replicat<<" # "<<touchable->GetHistoryDepth()<<std::endl;
  // std::cout<<" current :"<<currentPhysical->GetName()<<" "<<copyNo<<" # "<<touchable->GetHistoryDepth()<<std::endl;

  auto it_layer = mapTrackID_Hits.find(copyNo);
  if(it_layer == mapTrackID_Hits.end())
    {
      int IdHit        = fHitsCollection->GetSize();
      G4SolHit* newHit = new G4SolHit(fHCID);

      newHit->Pname       = aStep->GetTrack()->GetDefinition()->GetParticleName();
      newHit->TrackID     = aStep->GetTrack()->GetTrackID();
      newHit->Energy      = energ_depos;
      newHit->Time        = aStep->GetPreStepPoint()->GetGlobalTime();
      newHit->TrackLength = aStep->GetTrack()->GetTrackLength();
      newHit->HitPosX     = aStep->GetTrack()->GetPosition().x();
      newHit->HitPosY     = aStep->GetTrack()->GetPosition().y();
      newHit->HitPosZ     = aStep->GetTrack()->GetPosition().z();
      newHit->ExitPosX    = aStep->GetTrack()->GetPosition().x();
      newHit->ExitPosY    = aStep->GetTrack()->GetPosition().y();
      newHit->ExitPosZ    = aStep->GetTrack()->GetPosition().z();
      newHit->MomX        = aStep->GetTrack()->GetMomentum().x();
      newHit->MomY        = aStep->GetTrack()->GetMomentum().y();
      newHit->MomZ        = aStep->GetTrack()->GetMomentum().z();
      newHit->Mass        = aStep->GetTrack()->GetDefinition()->GetPDGMass();
      newHit->LayerID     = copyNo;

      mapTrackID_Hits.insert(std::pair<int, std::unordered_map<int, int> >(copyNo, {{CurrentTrack, IdHit}}));

      fHitsCollection->insert(newHit);
    }
  else
    {
      int IdHit = -1;
      // for(int i=(int)Event->id.size()-1;i>-1;i--)

      // G4cout<<" --> Track#"<<CurrentTrack<<G4endl;

      auto it_track = it_layer->second.find(CurrentTrack);

      if(it_track != it_layer->second.end())
        IdHit = it_track->second;

      // G4cout<<" --> IdHit:"<<IdHit;

      if(IdHit == -1)
        {
          // G4cout<<" NewHit !"<<G4endl;
          IdHit            = fHitsCollection->GetSize();
          G4SolHit* newHit = new G4SolHit(fHCID);

          newHit->Pname       = aStep->GetTrack()->GetDefinition()->GetParticleName();
          newHit->TrackID     = aStep->GetTrack()->GetTrackID();
          newHit->Energy      = energ_depos;
          newHit->Time        = aStep->GetPreStepPoint()->GetGlobalTime();
          newHit->TrackLength = aStep->GetTrack()->GetTrackLength();
          newHit->HitPosX     = aStep->GetTrack()->GetPosition().x();
          newHit->HitPosY     = aStep->GetTrack()->GetPosition().y();
          newHit->HitPosZ     = aStep->GetTrack()->GetPosition().z();
          newHit->ExitPosX    = aStep->GetTrack()->GetPosition().x();
          newHit->ExitPosY    = aStep->GetTrack()->GetPosition().y();
          newHit->ExitPosZ    = aStep->GetTrack()->GetPosition().z();
          newHit->MomX        = aStep->GetTrack()->GetMomentum().x();
          newHit->MomY        = aStep->GetTrack()->GetMomentum().y();
          newHit->MomZ        = aStep->GetTrack()->GetMomentum().z();
          newHit->Mass        = aStep->GetTrack()->GetDefinition()->GetPDGMass();
          newHit->LayerID     = copyNo;
          it_layer->second.insert(std::pair<int, int>(CurrentTrack, IdHit));

          fHitsCollection->insert(newHit);
        }
      else
        {
          // G4cout<<" Hit#"<<IdHit<<G4endl;
          G4SolHit* CurrentHit = dynamic_cast<G4SolHit*>((*fHitsCollection)[IdHit]);
          if(CurrentHit != nullptr)
            {
              CurrentHit->Energy += energ_depos;
              CurrentHit->ExitPosX = aStep->GetTrack()->GetPosition().x();
              CurrentHit->ExitPosY = aStep->GetTrack()->GetPosition().y();
              CurrentHit->ExitPosZ = aStep->GetTrack()->GetPosition().z();
            }
        }
    }

  return true;
}
