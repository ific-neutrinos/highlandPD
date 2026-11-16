#include "EventBoxKaon.hxx"
//#include "CutUtils.hxx"
#include "AnalysisUtils.hxx"
#include "pdDataClasses.hxx"

//********************************************************************
EventBoxKaon::EventBoxKaon():EventBoxB(){
//******************************************************************** 

}

//********************************************************************
EventBoxKaon::~EventBoxKaon(){
//********************************************************************

} 

//********************************************************************
void boxUtils::FillKaonCandidatesAndDaughters(AnaEventB& event, SubDetId::SubDetEnum det){
//********************************************************************

  EventBoxB* EventBox = event.EventBoxes[EventBoxId::kEventBoxKaon];

  EventBox->nRecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters]=0;
  anaUtils::CreateArray(EventBox->RecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters], NMAXPARTICLES);

  AnaParticleB** parts = static_cast<AnaEventB*>(&event)->Particles;
  int nParts           = static_cast<AnaEventB*>(&event)->nParticles;
    
  //loop over particles
  int nCandidatesAndDaughters = 0;
  for(Int_t i = 0; i < nParts; i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    //skip beam particle
    if(part->isPandora)continue;
    //look for candidates
    if(part->DaughtersIDs.size() == 1){
      //check this particle hasn't been added before (e.g. as the daughter of another candidate)
      bool added = false;
      for(int j = 0; j < EventBox->nRecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters]; j++){
	if(EventBox->RecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters][j] == part){
	  added = true;
	  break;
	}
      }
      //if it hasn't, add it
      if(!added){
	EventBox->RecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters][EventBox->nRecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters]++] = part;
	nCandidatesAndDaughters++;
      }
      //now go for the daughter
      if(part->Daughters.empty())continue; //daughter not in the event
      AnaParticleB* dau = static_cast<AnaParticleB*>(part->Daughters[0]);
      if(!dau)continue;
      added = false;
      for(int j = 0; j < EventBox->nRecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters]; j++){
	if(EventBox->RecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters][j] == dau){
	  added = true;
	  break;
	}
      }
      if(!added){
	EventBox->RecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters][EventBox->nRecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters]++] = dau;
	nCandidatesAndDaughters++;
      }
    }
  }
  anaUtils::ResizeArray(EventBox->RecObjectsInGroup [EventBoxKaon::kCandidatesAndDaughters], 
                        EventBox->nRecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters],
                        nCandidatesAndDaughters);
}

//********************************************************************
void boxUtils::FillKaonXS(AnaEventB& event, SubDetId::SubDetEnum det){
//********************************************************************

  EventBoxB* EventBox = event.EventBoxes[EventBoxId::kEventBoxKaon];

  EventBox->nRecObjectsInGroup[EventBoxKaon::kKaonXS]=0;
  anaUtils::CreateArray(EventBox->RecObjectsInGroup[EventBoxKaon::kKaonXS], 10*NMAXPARTICLES);

  AnaParticleB** parts = static_cast<AnaEventB*>(&event)->Particles;
  int nParts           = static_cast<AnaEventB*>(&event)->nParticles;
    
  //loop over particles
  int nCandidatesAndDaughters = 0;
  for(Int_t i = 0; i < nParts; i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    EventBox->RecObjectsInGroup[EventBoxKaon::kKaonXS][EventBox->nRecObjectsInGroup[EventBoxKaon::kKaonXS]++] = part;
    nCandidatesAndDaughters++;
  }
  anaUtils::ResizeArray(EventBox->RecObjectsInGroup [EventBoxKaon::kKaonXS], 
                        EventBox->nRecObjectsInGroup[EventBoxKaon::kKaonXS],
                        nCandidatesAndDaughters);
}


//********************************************************************
void boxUtils::FillProtonCandidates(AnaEventB& event, SubDetId::SubDetEnum det){
//********************************************************************

  EventBoxB* EventBox = event.EventBoxes[EventBoxId::kEventBoxKaon];

  EventBox->nRecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters]=0;
  anaUtils::CreateArray(EventBox->RecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters], NMAXPARTICLES);

  AnaParticleB** parts = static_cast<AnaEventB*>(&event)->Particles;
  int nParts           = static_cast<AnaEventB*>(&event)->nParticles;
    
  //loop over particles
  int nCandidates = 0;
  for(Int_t i = 0; i < nParts; i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    //skip beam particle
    if(part->isPandora)continue;
    //look for candidates
    if(part->DaughtersIDs.size() == 0){
      //check this particle hasn't been added before (e.g. as the daughter of another candidate)
      bool added = false;
      for(int j = 0; j < EventBox->nRecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters]; j++){
	if(EventBox->RecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters][j] == part){
	  added = true;
	  break;
	}
      }
      //if it hasn't, add it
      if(!added){
	EventBox->RecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters][EventBox->nRecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters]++] = part;
	nCandidates++;
      }
    }
  }
  anaUtils::ResizeArray(EventBox->RecObjectsInGroup [EventBoxKaon::kCandidatesAndDaughters], 
                        EventBox->nRecObjectsInGroup[EventBoxKaon::kCandidatesAndDaughters],
                        nCandidates);
}
