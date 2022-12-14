#include "pandoraPreselection.hxx"
#include "EventBoxPD.hxx"
#include "pdAnalysisUtils.hxx"
#include "pionAnalysisUtils.hxx"

//********************************************************************
pandoraPreselection::pandoraPreselection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxPD) {
//********************************************************************

}

//********************************************************************
void pandoraPreselection::DefineSteps(){
//********************************************************************

  // Steps must be added in the right order
  // if "true" is added to the constructor of the step,
  // the step sequence is broken if cut is not passed (default is "false")
  AddStep(StepBase::kAction, "find Pandora track",         new FindPandoraTrackAction());  
  AddStep(StepBase::kCut,    "candidate exists",           new CandidateExistsCut());
  AddStep(StepBase::kCut,    "beam quality",               new BeamQualityCut());
  //AddStep(StepBase::kCut,    "pandora reco worked",        new CandidateIsBeamCut());

  SetBranchAlias(0,"pandora");
  
  // No preselection for the moment
  SetPreSelectionAccumLevel(-1);
}

//**************************************************
bool FindPandoraTrackAction::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  // Cast the ToyBox to the appropriate type
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB); 

  // This action fills box.MainTrack 
  box.MainTrack = NULL;
  
  // Get the array of parts from the event
  AnaParticleB** parts = static_cast<AnaEventB*>(&event)->Particles;
  int nParts           = static_cast<AnaEventB*>(&event)->nParticles;

  //look over the particles in the event
  for(int i = 0; i < nParts; i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    if(part->isPandora){
      box.MainTrack = part;
      break;
    }
  }
  return true;
}

//**************************************************
bool CandidateExistsCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   
  if (!box.MainTrack) return false;
  else return true;
}

//**************************************************
bool BeamQualityCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  // Cast the candidate
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);
  AnaParticlePD* part = static_cast<AnaParticlePD*>(box.MainTrack);
  //Get the beam from the event
  AnaBeamPD* beam = static_cast<AnaBeamPD*>(static_cast<AnaEventB*>(&event)->Beam);

  bool candidateIsBeam=false;
  if (useIsBeamLike)      
    candidateIsBeam = pdAnaUtils::isBeamLike(part,beam);
  else
    candidateIsBeam = (bool)part->isBeamPart;
  
  return candidateIsBeam;
}

//**************************************************
bool CandidateIsBeamCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  // Cast the candidate
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);
  AnaParticlePD* part = static_cast<AnaParticlePD*>(box.MainTrack);
  //Get the beam from the event
  AnaBeamPD* beam = static_cast<AnaBeamPD*>(static_cast<AnaEventB*>(&event)->Beam);

  bool candidateIsBeam=false;
  if (useIsBeamLike)      
    candidateIsBeam = pdAnaUtils::isBeamLike(part,beam);
  else
    candidateIsBeam = (bool)part->isBeamPart;
  
  return candidateIsBeam;
}

//**************************************************
bool BeamPDGCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  //Get the beam from the event and the beam particle
  AnaBeamPD* beam = static_cast<AnaBeamPD*>(static_cast<AnaEventB*>(&event)->Beam);
  AnaParticlePD* beamPart = static_cast<AnaParticlePD*>(beam->BeamParticle);
  if(!beamPart)return false;

  if(event.GetIsMC()){
    AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(beamPart->TrueObject);
    if(!truePart)return false;
    if(abs(_PDG)==abs(truePart->PDG))return true;
    else return false;
  }
  else{
    if(std::find(beam->PDGs.begin(),beam->PDGs.end(),_PDG) != beam->PDGs.end())return true;
    else return false;
  }

  return false;
}

//**************************************************
void pandoraPreselection::InitializeEvent(AnaEventC& eventBB){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventBB); 

  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxPD])
    event.EventBoxes[EventBoxId::kEventBoxPD] = new EventBoxPD();

  boxUtils::FillCandidateAndDaughters(event);
}
