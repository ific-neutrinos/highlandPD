#include "secondaryKaonSelection.hxx"
#include "EventBoxPD.hxx"
#include "pdAnalysisUtils.hxx"
#include "pandoraPreselection.hxx"

//********************************************************************
secondaryKaonSelection::secondaryKaonSelection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxPD) {
//********************************************************************

}

//********************************************************************
void secondaryKaonSelection::DefineSteps(){
//********************************************************************

  // Steps must be added in the right order
  // if "true" is added to the constructor of the step,
  // the step sequence is broken if cut is not passed (default is "false")

  //secondary kaon selection
  AddStep(StepBase::kAction,   "get a vector of kaons",      new GetKaonsAction()        );
  AddStep(StepBase::kCut,      "we have a kaon",             new EventHasKaonCut(),  true);
  //split the selection in branches, one for each possible candidate
  AddSplit(secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES);
  //next cuts have to be applied to each branch
  for(int i = 0; i < (int)secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES; i++){
    AddStep(i, StepBase::kCut, "kaon daughter is a track",   new MuonIsTrackCut(),             true);
    AddStep(i, StepBase::kCut, "kaon daughter chi2 cut",     new MuonChi2Cut(),                true);
    AddStep(i, StepBase::kCut, "kaon daughter CNN cut",      new MuonCNNCut(),                 true);
    AddStep(i, StepBase::kCut, "kaon daughter mom cut",      new MuonRangeMomCut(),            true);
    AddStep(i, StepBase::kCut, "kaon CNN cut",               new KaonCNNCut(),                 true);
    AddStep(i, StepBase::kCut, "kaon-muon angle cut",        new KaonMuonAngleCut(),           true);
    AddStep(i, StepBase::kCut, "kaon-muon distance cut",     new KaonMuonDistanceCut(),        true);
  }

  //Set the branch aliases to the different branches 
  for(int i = 0; i < (int)secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES; i++){
    std::stringstream ssi;
    ssi << i;
    SetBranchAlias(i,("possible candidate "+ssi.str()+"").c_str(),i);
  }  

  // This number tells the selection to stop when a given accum level (accumulated cut level) is not passed. 
  // For example, 2 as argument would mean that selection is stopped for the current toy experiment when the 
  // third cut (numbering starts at 0) is not passed. This is a way of saving time.
  // -1 means no preselection
  SetPreSelectionAccumLevel(-1);
}

//**************************************************
bool GetKaonsAction::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  
  // Cast the ToyBox to the appropriate type
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB); 
    
  // Get the array of parts from the event
  AnaParticleB** parts = static_cast<AnaEventB*>(&event)->Particles;
  int nParts           = static_cast<AnaEventB*>(&event)->nParticles;

  //look over the particles in the event
  for(int i = 0; i < nParts; i++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[i]);
    if(part->Daughters.size() == 1 && part->ParentID != -1)box.Candidates.push_back(part);
  }
  
  return true;  
}

//**************************************************
bool EventHasKaonCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);   

  //if the main track has at least a single daughter with a single daughter, pass
  if (box.Candidates.size()>0) return true;
  else return false;
}

//**************************************************
bool MuonIsTrackCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;

  //cast the box
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);   

  //get the current branch  
  std::vector<UInt_t> branchesIDs = GetBranchUniqueIDs();
  
  //make sure there is a candidate for this branch
  if(branchesIDs[0] > box.Candidates.size()-1)return false;

  AnaParticlePD* dau = static_cast<AnaParticlePD*>(box.Candidates[branchesIDs[0]]->Daughters[0]);
  if(dau->Type == 2)return true;
  else return false; 
}

//**************************************************
bool MuonChi2Cut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);   

  //get the current branch  
  std::vector<UInt_t> branchesIDs = GetBranchUniqueIDs();

  //if there is kaon, cut in its daughter chi2
  AnaParticlePD* dau = static_cast<AnaParticlePD*>(box.Candidates[branchesIDs[0]]->Daughters[0]);
  /*std::pair<double,int> prot_result = pdAnaUtils::Chi2PID(*dau,2212);
  std::pair<double,int> muon_result = pdAnaUtils::Chi2PID(*dau,13);
  dau->Chi2Proton = prot_result.first;
  dau->Chi2Muon   = muon_result.first;
  dau->Chi2ndf    = muon_result.second;*/
  if(dau->Chi2Muon / dau->Chi2ndf < 10)return true;
  else return false; 
}

//**************************************************
bool MuonCNNCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);   

  //get the current branch  
  std::vector<UInt_t> branchesIDs = GetBranchUniqueIDs();

  //if there is kaon, cut in its daughter CNN score
  AnaParticlePD* dau = static_cast<AnaParticlePD*>(box.Candidates[branchesIDs[0]]->Daughters[0]);
  if(dau->CNNscore[0] > 0.6)return true;
  else return false; 
}

//**************************************************
bool MuonRangeMomCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);   

  //get the current branch  
  std::vector<UInt_t> branchesIDs = GetBranchUniqueIDs();

  //if there is kaon, cut in its daughter range mom
  AnaParticlePD* dau = static_cast<AnaParticlePD*>(box.Candidates[branchesIDs[0]]->Daughters[0]);
  double mom = pdAnaUtils::ComputeRangeMomentum(dau->Length,13);
  if(fabs(mom-0.22) < 0.02)return true;
  else return false; 
}

//**************************************************
bool KaonCNNCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);   

  //get the current branch  
  std::vector<UInt_t> branchesIDs = GetBranchUniqueIDs();

  //if there is kaon, cut on its CNN
  if(box.Candidates[branchesIDs[0]]->CNNscore[0] > 0.8)return true;
  else return false; 
}

//**************************************************
bool KaonMuonAngleCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);   

  //get the current branch  
  std::vector<UInt_t> branchesIDs = GetBranchUniqueIDs();

  //if there is kaon, cut on the angle between it and its daughter muon
  AnaParticlePD* kaon = box.Candidates[branchesIDs[0]];
  AnaParticlePD* muon = static_cast<AnaParticlePD*>(box.Candidates[branchesIDs[0]]->Daughters[0]);
  double cos = 0;
  for(int i = 0; i < 3; i++)cos = cos + kaon->DirectionEnd[i] * muon->DirectionStart[i];
  if(cos < 0.1 )return true;
  else return false; 
}

//**************************************************
bool KaonMuonDistanceCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************
  
  (void)event;

  //cast the box
  ToyBoxKaon& box = *static_cast<ToyBoxKaon*>(&boxB);   

  //get the current branch  
  std::vector<UInt_t> branchesIDs = GetBranchUniqueIDs();

  //if there is kaon, cut on the distance between it and its daughter muon
  AnaParticlePD* kaon = box.Candidates[branchesIDs[0]];
  AnaParticlePD* muon = static_cast<AnaParticlePD*>(box.Candidates[branchesIDs[0]]->Daughters[0]);
  double dis = 0;
  for(int i = 0; i < 3; i++)dis = dis + pow(kaon->PositionEnd[i] - muon->PositionStart[i],2);
  dis = sqrt(dis);
  if(dis < 10)return true;
  else return false; 
}

//**************************************************
void secondaryKaonSelection::InitializeEvent(AnaEventC& eventC){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventC); 

  // Create the appropriate EventBox if it does not exist yet
  if (!event.EventBoxes[EventBoxId::kEventBoxPD])
    event.EventBoxes[EventBoxId::kEventBoxPD] = new EventBoxPD();

  boxUtils::FillCandidateAndDaughters(event);
  boxUtils::FillTrueCandidateAndDaughters(event);
}
