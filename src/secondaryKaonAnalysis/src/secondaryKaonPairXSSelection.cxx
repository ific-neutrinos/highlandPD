#include "TPad.h"
#include "TMultiGraph.h"
#include "secondaryKaonPairXSSelection.hxx"
#include "secondaryKaonSelection.hxx"
#include "EventBoxKaon.hxx"
#include "Parameters.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdBaseSelection.hxx"
#include "muCandidateUtils.hxx"

//********************************************************************
secondaryKaonPairXSSelection::secondaryKaonPairXSSelection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxKaon) {
//********************************************************************
  
}

//********************************************************************
void secondaryKaonPairXSSelection::DefineSteps(){
//********************************************************************

  //copy steps from pdBaseSelection
  AddStep(StepBase::kAction,   "find beam track",            new FindBeamTrackAction(), true   );// in pdBaseAnalysis/src/pdBaseSelection  
  AddStep(StepBase::kCut,      "Pandora track exists",       new BeamTrackExistsCut(), true    );// in pdBaseAnalysis/src/pdBaseSelection  
  // Hasta el 25-nov exigía un haz tipo pion (BeamPDGCut(211)); para compararla con secondaryKaonProductionXS, que no
  // corta en la especie del haz (B se mide por especie después), es configurable y por defecto no corta
  const int beamPDG = ND::params().GetParameterI("secondaryKaonAnalysis.pairXSAnalysis.BeamPDG");
  AddStep(StepBase::kCut,    beamPDG ? Form("beam particle PDG %d", beamPDG) : "beam particle (any species)", new PairXSBeamPDGCut(beamPDG));
  //select events using beam instrumentation

  AddStep(StepBase::kAction,   "get a vector of kaon candidates",      new GetKaonsForPairXSAction(),true   );
  AddStep(StepBase::kCut,      "we have a kaon",             new EventHasKaonCut(),  true);
  
  AddStep(StepBase::kAction,   "get forced daughter",        new GetForcedDaughterPairXSAction()   );

  SetPreSelectionAccumLevel(-1);

  SetBranchAlias(0,"NO BRANCH :)");
}

//**************************************************
bool GetKaonsForPairXSAction::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  
  // Cast the ToyBox to the appropriate type
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);

  // Get the beam track
  AnaParticlePD* beamTrack = box.MainTrack;
    
  //look over daughters and look for candidates
  for(int i = 0; i < (int)beamTrack->Daughters.size(); i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(beamTrack->Daughters[i]);
    if((int)part->Daughters.size()>1)continue;
    std::pair<double,int> pid = pdAnaUtils::Chi2PID(*part,321);
    if(pid.first<0)continue;
    double chi_kaon = pid.first/pid.second;
    pid = pdAnaUtils::Chi2PID(*part,2212);
    if(pid.first<0)continue;
    double chi_prot = pid.first/pid.second;
    if(chi_prot-chi_kaon<50)
      box.Candidates.push_back(part);
  }
  
  return true;  
}

  

//   //cast the box

//   //get the current branch  

//   //if there is kaon, cut in its chi2
//   else return false;

//**************************************************
bool GetForcedDaughterPairXSAction::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  
  // Cast the ToyBox to the appropriate type
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);

  // Get the array of parts from the event
  AnaParticleB** parts = static_cast<AnaEventB*>(&event)->Particles;
  int nParts           = static_cast<AnaEventB*>(&event)->nParticles;

  // vista normal: fuera los PFParticles que el módulo no guardaría en un evento normal
  const AnaEventInfoPD* info = static_cast<AnaEventInfoPD*>(static_cast<AnaEventB*>(&event)->EventInfo);
  const bool fullDump = info && info->FullDump;
  bool hasSlice = false;
  for(int ipart = 0; ipart < nParts; ipart++)if(static_cast<AnaParticlePD*>(parts[ipart])->BeamOrigin){ hasSlice = true; break; }
  auto inNormalView = [&](const AnaParticlePD* p){ return !fullDump || muCandUtils::NormalSaveMode(p, hasSlice) != muCandUtils::kNotSaved; };

  // loop over candidates
  for(int ican = 0; ican < (int)box.Candidates.size(); ican++){
    AnaParticlePD* candidate = box.Candidates[ican];
    // skip candidates with daughters
    if(!candidate->Daughters.empty())continue;
    double dis_min = 10000;
    int i_closest = -1;
    //loop over particles and look for the one starting closer to the end point of the candidate
    for(int ipart = 0; ipart < nParts; ipart++){
      AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[ipart]);
      if(candidate->UniqueID == part->UniqueID)continue; //skip itself
      if(part->isPandora)continue; //skip beam particle
      if(!inNormalView(part))continue;
      if(part->Type!=2)continue; //skip showers maybe not needed?
      double mom = pdAnaUtils::ComputeRangeMomentum(part->Length,13);
      if(mom < 0.2 || mom > 0.3)continue;
      double dis = 0;
      for(int i = 0; i < 3; i++)dis += pow(part->PositionStart[i]-candidate->PositionEnd[i],2);
      dis = sqrt(dis);
      if(dis < dis_min){
	dis_min = dis;
	i_closest = ipart;
      }
    }
    if(i_closest == -1 || dis_min>3){ //above three almost everything are protons
      //if here and no forced daughter, remove candidate
      if(candidate->Daughters.empty())
	box.Candidates.erase(box.Candidates.begin()+ican);
      ican--;
      continue;
    }
    //associate the closest one as forced daughter
    AnaParticlePD* forced_dau = static_cast<AnaParticlePD*>(parts[i_closest]);
    candidate->forced_daughter = true;
    candidate->Daughters.push_back(forced_dau);
    
    //check true matching if possible
    AnaTrueParticlePD* truecandidate = static_cast<AnaTrueParticlePD*>(candidate->TrueObject);
    AnaTrueParticlePD* trueforced_dau = static_cast<AnaTrueParticlePD*>(forced_dau->TrueObject);
    if(!truecandidate || !trueforced_dau)continue;
    for(int i = 0; i < (int)truecandidate->Daughters.size(); i++){
      if(truecandidate->Daughters[i] == trueforced_dau->ID){
	candidate->forced_daughter_matched = true;
	break;
      }
    }
  }
    
  return true;  
}

//**************************************************
bool PairXSBeamPDGCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  if(_pdg == 0)return true;
  BeamPDGCut cut(_pdg);
  return cut.Apply(event, boxB);
}

//**************************************************
void secondaryKaonPairXSSelection::InitializeEvent(AnaEventC& eventC){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventC); 

  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxKaon])
    event.EventBoxes[EventBoxId::kEventBoxKaon] = new EventBoxKaon();

  boxUtils::FillKaonXS(event);
}

