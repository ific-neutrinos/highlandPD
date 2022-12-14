#include "pionSelection.hxx"
#include "EventBoxPD.hxx"
#include "pdAnalysisUtils.hxx"
#include "pionAnalysisUtils.hxx"
#include "pdBaseSelection.hxx"
#include "CNNUtils.hxx"


/*
  This file contains:
  - the method DefineSteps in which actions and cuts are added to the selection in the right order. 
  - the implementation of each of the steps
  - the method InitializeEvent, in which the EventBox is created. This box is used to make a subsample 
    of event objects (Particles, hits, etc), which are frequently needed, such that we don't have to create 
    the subsample for every toy experiment dusing systematic propagation
 */


//********************************************************************
pionSelection::pionSelection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxPD) {
//********************************************************************

}

//********************************************************************
void pionSelection::DefineSteps(){
//********************************************************************

  // Steps must be added in the right order
  // if "true" is added to the constructor of the step,
  // the step sequence is broken if cut is not passed (default is "false")

  //copy steps from pandoraPreselecion
  AddStep(StepBase::kAction, "find Pandora track",         new FindBeamTrackAction()); // in pdBaseAnalysis/src/pandoraPreselection  
  AddStep(StepBase::kCut,    "candidate exists",           new BeamTrackExistsCut());  // in pdBaseAnalysis/src/pandoraPreselection  
  AddStep(StepBase::kCut,    "beam is pion",               new BeamPionCut());
  //AddStep(StepBase::kCut,    "pandora reco worked",        new CandidateIsBeamCut());              // in pdBaseAnalysis/src/pandoraPreselection  
  AddStep(StepBase::kCut,    "candidate is track",         new CandidateIsTrackCut());
  AddStep(StepBase::kCut,    "pion track end",             new PionEndsAPA3Cut());
  //AddStep(StepBase::kCut,    "seltrk chi2 cut",            new PionPassChi2Cut());
  AddStep(StepBase::kCut,    "no pion daughter",           new NoPionDaughterCut());

  //Add a split to the trunk with 2 branches.
  AddSplit(2);

  //First branch is for pi0 showers
  AddStep(0, StepBase::kCut, "pi0 showers",        new PionCexCut());

  //Second branch is for pi absortion
  AddStep(1, StepBase::kCut, "pi absorption",      new PionAbsorptionCut());

  // Set the branch aliases to the three branches 
  SetBranchAlias(0,"pi0 shower",     0);
  SetBranchAlias(1,"pi absorption",  1);

  // This number tells the selection to stop when a given accum level (accumulated cut level) is not passed. For example, 2 as argument would mean
  // that selection is stopped for the current toy experiment when the third cut (numbering starts at 0) is not passed. This is a way of saving time.
  // -1 means no preselection
  SetPreSelectionAccumLevel(5);   // TODO: No need to Apply dEdx systematics
}

//**************************************************
bool BeamPionCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  BeamPDGCut* beamCut1 = new BeamPDGCut(211);
  BeamPDGCut* beamCut2 = new BeamPDGCut(13);

  if(beamCut1->Apply(event,boxB) || beamCut2->Apply(event,boxB))return true;
  else return false;
}

//**************************************************
bool CandidateIsTrackCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   

  // Main track shouldmust exist
  if (!box.MainTrack) return false;
  
  //if the seltrk is a track, accept it
  if (box.MainTrack->Type == AnaParticlePD::kTrack) return true;
  else return false;
}

//**************************************************
bool PionEndsAPA3Cut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  Float_t cutAPA3_Z = 226.;

  // Cast the ToyBox to the appropiate type
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   
  if (!box.MainTrack) return false;

  if(box.MainTrack->PositionEnd[2] < cutAPA3_Z)return true;
  else return false;  
}

//**************************************************
bool PionPassChi2Cut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  Float_t cut_chi2 = 140;

  //Cast the ToyBox to the appropiate type
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);
  if (!box.MainTrack) return false;
  
  if(box.MainTrack->Chi2Proton/box.MainTrack->Chi2ndf > cut_chi2) return true;
  else return false;
}

//**************************************************
bool NoPionDaughterCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  Float_t cut_CNNTrackScore = 0.3;
  double cut_dEdx = 3.8;
  
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   
  if (!box.MainTrack) return false;

  // Look for pion daughters (track hypothesis)
  bool noPion = true;
  for(UInt_t i = 0; i < box.MainTrack->Daughters.size(); i++){
    AnaParticlePD* daughter = static_cast<AnaParticlePD*>(box.MainTrack->Daughters[i]);

    //    daughter->Print();
    
    if(daughter->Type           == AnaParticlePD::kTrack &&
       daughter->UniqueID       != -999 &&
       daughter->truncLibo_dEdx <= cut_dEdx){

      // To minimize the number of times CNN is called, we only call it when needed, that is
      // when the selection gets to this point and the CNN is around the cut value
      if (_cnnSystematicEnabled){
        if (fabs(daughter->CNNscore[0]-cut_CNNTrackScore)<_cnnRecomputeCut){
          event.Print();
          pdAnaUtils::ComputeParticleCNN((static_cast<AnaEventPD*>(&event))->CNNwires,*daughter);
          //          pdAnaUtils::ComputeParticleCNN(*daughter);
        }
      }            
      if (daughter->CNNscore[0]    > cut_CNNTrackScore){        
        noPion = false;
        break;
      }
    }
  }

  return noPion;  
}

//**************************************************
//bool NoPionDaughterCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
/*
  (void)event;
  Float_t cut_daughter_track_distance = 10;
  Float_t cut_CNNTrackScore = 0.35;
  Float_t cut_chi2 = 50;
   
  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   
  if (!box.MainTrack) return false;

  // Look for pion daughters (track hypothesis)
  bool noPion = true;
  for(UInt_t i = 0; i < box.MainTrack->Daughters.size()/2; i++){
    AnaParticlePD* daughter = static_cast<AnaParticlePD*>(box.MainTrack->Daughters[i]);
    if(daughter->UniqueID != -999 &&
       daughter->CNNscore[0]                  > cut_CNNTrackScore &&
       daughter->Chi2Proton/daughter->Chi2ndf > cut_chi2 &&
       box.DaughterDistanceToVertex[i]        < cut_daughter_track_distance){
      noPion = false;
      break;
    }
  }
  return noPion;  
}*/

//**************************************************
bool PionCexCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  int cut_nHits_shower_low  = 40;
  Float_t cut_CNNTrackScore = 0.3;

  ToyBoxPD& box = *static_cast<ToyBoxPD*>(&boxB);   
  if (!box.MainTrack) return false;

  // loop over daughters assuming shower hypothesis
  bool cex = false;
  for(UInt_t i = 0; i < box.MainTrack->Daughters.size(); i++){
    AnaParticlePD* daughter = static_cast<AnaParticlePD*>(box.MainTrack->Daughters[i]);    
    if(daughter->Type       == AnaParticlePD::kShower &&
       daughter->NHits       > cut_nHits_shower_low &&
       daughter->CNNscore[0] != -999){

      // To minimize the number of times CNN is called, we only call it when needed, that is
      // when the selection gets to this point and the CNN is around the cut value
      if (_cnnSystematicEnabled){
        if (fabs(daughter->CNNscore[0]-cut_CNNTrackScore)<_cnnRecomputeCut){
          pdAnaUtils::ComputeParticleCNN(*daughter);
        }
      }            

      if (daughter->CNNscore[0] < cut_CNNTrackScore){
        cex = true;
        break;
      }
    }
  }
  return cex;
}

//**************************************************
bool PionAbsorptionCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  // TODO: make sure CNN is computed only once
  return !_cut.Apply(event, boxB);
}

//**************************************************
void pionSelection::InitializeEvent(AnaEventC& eventC){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventC); 

  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxPD])
    event.EventBoxes[EventBoxId::kEventBoxPD] = new EventBoxPD();

  boxUtils::FillCandidateAndDaughters(event);
  boxUtils::FillTrueCandidateAndDaughters(event);
}

