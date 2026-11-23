#include "secondaryKaonProductionXSSelection.hxx"
#include "EventBoxKaon.hxx"
#include "pdBaseSelection.hxx"
#include "kaonAnalysisUtils.hxx"

//********************************************************************
std::string muCandCuts::Name(int id){
//********************************************************************

  switch(id){
    case kPopulation:    return "poblaci#acute{o}n";
    case kFaces35:       return Form("extremos a > %.0f cm de las caras", kFacesMin);
    case kNoClearCosmic: return "no IsClearCosmic";
    case kAxis60:        return Form("eje < %.0f cm", kAxisMax);
    case kVertex15:      return Form("v#acute{e}rtice > %.0f cm", kVertexMin);
    case kChi2Mu:        return Form("#chi^{2}_{#mu} < %.0f", kChi2MuMax);
    case kdEdxStop:      return Form("dE/dx parada > %.0f", kdEdxStopMin);
    case kChi2p:         return Form("#chi^{2}_{p} > %.0f", kChi2pMin);
    case kdEdxStart:     return Form("dE/dx inicio > %.1f", kdEdxStartMin);
    case kGaps:          return "sin extremos en los huecos";
    case kWindow:        return Form("p rango en [%.0f, %.0f]", kWindowMin, kWindowMax);
    case kMichelHits:    return Form("hits Michel #geq %d", kMichelHitsMin);
  }
  return "?";
}

//********************************************************************
bool muCandCuts::Pass(int id, const muCandUtils::Candidate& c, int pop){
//********************************************************************

  switch(id){
    case kPopulation:    return c.inSlice == (pop == 0);
    case kFaces35:       return std::min(FacesStop(c), FacesStart(c)) > kFacesMin;
    case kNoClearCosmic: return !c.clearCosmic;
    case kAxis60:        return AxisMin(c) < kAxisMax;
    case kVertex15:      return VertexMin(c) < 0 || VertexMin(c) > kVertexMin;
    case kChi2Mu:        return Chi2MuStop(c) >= 0 && Chi2MuStop(c) < kChi2MuMax;
    case kdEdxStop:      return dEdxStop(c) > kdEdxStopMin;
    case kChi2p:         return Chi2pStop(c) > kChi2pMin;
    case kdEdxStart:     return dEdxStart(c) > kdEdxStartMin;
    case kGaps:          return GapDist(c.startZ) >= kGapHalfWidth && GapDist(c.stopZ) >= kGapHalfWidth;
    case kWindow:        return c.pRange >= kWindowMin && c.pRange <= kWindowMax;
    case kMichelHits:    return c.arNMichel[1] < -900 ? c.cnnMichel > 0.02 : c.arNMichel[1] >= kMichelHitsMin;
  }
  return false;
}

//********************************************************************
std::vector<int> muCandCuts::Chain(int pop){
//********************************************************************

  if(pop == 0)return {kPopulation, kFaces35, kVertex15, kChi2Mu, kdEdxStop, kGaps, kWindow, kMichelHits};
  return {kPopulation, kNoClearCosmic, kAxis60, kChi2Mu, kdEdxStop, kChi2p, kdEdxStart, kGaps, kWindow, kMichelHits};
}

//********************************************************************
secondaryKaonProductionXSSelection::secondaryKaonProductionXSSelection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxKaon) {
//********************************************************************

}

//********************************************************************
void secondaryKaonProductionXSSelection::DefineSteps(){
//********************************************************************

  // El primario de Pandora, para las variables del haz (no corta: siempre devuelve true)
  AddStep(StepBase::kAction, "find beam track",     new FindBeamTrackAction());// in pdBaseAnalysis/src/pdBaseSelection
  AddStep(StepBase::kAction, "get mu candidates",   new GetMuCandidatesAction());
  AddStep(StepBase::kCut,    "mu candidate exists", new MuCandidateExistsCut(), true);

  // una rama por población
  AddSplit(2);
  for(int pop = 0; pop < 2; pop++){
    const std::vector<int> chain = muCandCuts::Chain(pop);
    for(int i = 0; i < (int)chain.size(); i++){
      const std::string title = chain[i] == muCandCuts::kPopulation ? (pop == 0 ? "en la beam slice" : "fuera de la beam slice") : muCandCuts::Name(chain[i]);
      AddStep(pop, StepBase::kCut, title, new MuCandidateCut(pop, i + 1, chain[i]), true);
    }
  }
  SetBranchAlias(0, "en la beam slice", 0);
  SetBranchAlias(1, "fuera de la beam slice", 1);

  SetPreSelectionAccumLevel(-1);
}

//**************************************************
bool GetMuCandidatesAction::Apply(AnaEventC& eventC, ToyBoxB& boxB) const{
//**************************************************

  ToyBoxMuCand& box = *static_cast<ToyBoxMuCand*>(&boxB);
  AnaEventPD& event = *static_cast<AnaEventPD*>(&eventC);

  std::vector<AnaParticlePD*> parts;
  for(int i = 0; i < event.nParticles; i++)parts.push_back(static_cast<AnaParticlePD*>(event.Particles[i]));
  std::vector<AnaTrueParticleB*> trueParts;
  for(int i = 0; i < event.nTrueParticles; i++)trueParts.push_back(event.TrueParticles[i]);
  const AnaBeamPD* beam = static_cast<AnaBeamPD*>(event.Beam);
  const AnaTrueParticleB* beamTrue = (beam && beam->BeamParticle) ? static_cast<AnaTrueParticleB*>(beam->BeamParticle->TrueObject) : nullptr;
  const int isMC = event.GetIsMC() ? 1 : 0;

  muCandUtils::EventContext ctx;
  muCandUtils::MakeContext(ctx, parts, event.UnassociatedHits, trueParts, beamTrue, isMC);
  box.HasVertex = ctx.hasVtx;
  if(ctx.hasVtx){ box.VertexPos[0] = ctx.vtx.X(); box.VertexPos[1] = ctx.vtx.Y(); box.VertexPos[2] = ctx.vtx.Z(); }

  for(const AnaParticlePD* p : parts){
    muCandUtils::Candidate c;
    if(!muCandUtils::Compute(p, ctx, c))continue;
    box.NPreselected++;
    // preselección geométrica del módulo: un extremo a >= 10 cm de las caras y uno a <= 150 cm del eje
    if(std::max(c.facesA, c.facesB) < 10 || muCandCuts::AxisMin(c) > 150)continue;
    box.Candidates.push_back(c);
    box.Population.push_back(c.inSlice ? 0 : 1);
    box.Level.push_back(0);
  }
  return true;
}

//**************************************************
bool MuCandidateExistsCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  return !static_cast<ToyBoxMuCand*>(&boxB)->Candidates.empty();
}

//**************************************************
bool MuCandidateCut::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  (void)event;
  ToyBoxMuCand& box = *static_cast<ToyBoxMuCand*>(&boxB);
  bool any = false;
  for(size_t i = 0; i < box.Candidates.size(); i++){
    if(box.Population[i] != _pop || box.Level[i] != _index - 1)continue;
    if(muCandCuts::Pass(_id, box.Candidates[i], _pop))box.Level[i] = _index;
    if(box.Level[i] == _index)any = true;
  }
  return any;
}

//**************************************************
void secondaryKaonProductionXSSelection::InitializeEvent(AnaEventC& eventC){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventC);

  // EventBox de los kaones (grupos de objetos para los sistemáticos; ahora no se usan)
  if (!event.EventBoxes[EventBoxId::kEventBoxKaon])
    event.EventBoxes[EventBoxId::kEventBoxKaon] = new EventBoxKaon();

  boxUtils::FillKaonXS(event);
}
