#include "secondaryKaonProductionXSSelection.hxx"
#include "EventBoxKaon.hxx"
#include "pdBaseSelection.hxx"

//********************************************************************
secondaryKaonProductionXSSelection::secondaryKaonProductionXSSelection(bool forceBreak): SelectionBase(forceBreak,EventBoxId::kEventBoxKaon) {
//********************************************************************

}

//********************************************************************
void secondaryKaonProductionXSSelection::DefineSteps(){
//********************************************************************

  // Pasos comunes a las dos ramas
  // - el primario de Pandora (box.MainTrack), solo para las variables del haz del árbol; no corta (siempre true)
  AddStep(StepBase::kAction, "find beam track",     new FindBeamTrackAction());// in pdBaseAnalysis/src/pdBaseSelection
  // - los candidatos a mu, con sus variables, en la caja (nivel 0)
  AddStep(StepBase::kAction, "get mu candidates",   new GetMuCandidatesAction());
  // - corte de evento: hay algún candidato (accum_level 1)
  AddStep(StepBase::kCut,    "mu candidate exists", new MuCandidateExistsCut(), true);

  // Una rama por población: 0 en la beam slice, 1 fuera
  AddSplit(2);

  // Cortes de cada rama. Cada corte sobre candidatos necesita saber su rama (pop) y su posición en la cadena
  // (index: 1 el primero, 2 el segundo...), para la contabilidad de los niveles de MuCandidateCutBase; esta función
  // los pone en orden, así que la cadena se escribe como cualquier selección de highland.
  int index[2] = {0, 0};
  auto add = [&](int pop, const std::string& title, MuCandidateCutBase* cut){
    cut->SetPopulationAndIndex(pop, ++index[pop]);
    AddStep(pop, StepBase::kCut, title, cut, true);   // true: si ningún candidato pasa, la rama se para aquí
  };

  // rama 0: en la beam slice (fondo: trozos de piones y protones del haz)
  add(0, "en la beam slice",                      new MuCandPopulationCut());
  add(0, "extremos a > 35 cm de las caras",       new MuCandFacesCut(muCandCuts::kFacesMin));
  add(0, "v#acute{e}rtice > 15 cm",               new MuCandVertexCut(muCandCuts::kVertexMin));
  add(0, "#chi^{2}_{#mu} < 6",                    new MuCandChi2MuCut(muCandCuts::kChi2MuMax));
  add(0, "dE/dx parada > 3",                      new MuCandStopdEdxCut(muCandCuts::kdEdxStopMin));
  add(0, "sin extremos en los huecos",            new MuCandGapsCut(muCandCuts::kGapHalfWidth));
  add(0, "p rango en [205, 240]",                 new MuCandRangeMomCut(muCandCuts::kWindowMin, muCandCuts::kWindowMax));
  add(0, "hits Michel #geq 2",                    new MuCandMichelHitsCut(muCandCuts::kMichelHitsMin));

  // rama 1: fuera de la beam slice (fondo: mu cósmicos, trozos EM, piones y protones de interacciones secundarias)
  add(1, "fuera de la beam slice",                new MuCandPopulationCut());
  add(1, "no IsClearCosmic",                      new MuCandNoClearCosmicCut());
  add(1, "eje < 60 cm",                           new MuCandAxisCut(muCandCuts::kAxisMax));
  add(1, "#chi^{2}_{#mu} < 6",                    new MuCandChi2MuCut(muCandCuts::kChi2MuMax));
  add(1, "dE/dx parada > 3",                      new MuCandStopdEdxCut(muCandCuts::kdEdxStopMin));
  add(1, "#chi^{2}_{p} > 80",                     new MuCandChi2pCut(muCandCuts::kChi2pMin));
  add(1, "dE/dx inicio > 1.8",                    new MuCandStartdEdxCut(muCandCuts::kdEdxStartMin));
  add(1, "sin extremos en los huecos",            new MuCandGapsCut(muCandCuts::kGapHalfWidth));
  add(1, "p rango en [205, 240]",                 new MuCandRangeMomCut(muCandCuts::kWindowMin, muCandCuts::kWindowMax));
  add(1, "hits Michel #geq 2",                    new MuCandMichelHitsCut(muCandCuts::kMichelHitsMin));

  SetBranchAlias(0, "en la beam slice", 0);
  SetBranchAlias(1, "fuera de la beam slice", 1);

  // -1: sin preselección (no se para la selección antes de tiempo)
  SetPreSelectionAccumLevel(-1);
}

//**************************************************
bool GetMuCandidatesAction::Apply(AnaEventC& eventC, ToyBoxB& boxB) const{
//**************************************************

  ToyBoxMuCand& box = *static_cast<ToyBoxMuCand*>(&boxB);
  AnaEventPD& event = *static_cast<AnaEventPD*>(&eventC);

  // Lo que necesita muCandUtils del evento: todas las partículas, la truth y la partícula true del haz
  std::vector<AnaParticlePD*> parts;
  for(int i = 0; i < event.nParticles; i++)parts.push_back(static_cast<AnaParticlePD*>(event.Particles[i]));
  std::vector<AnaTrueParticleB*> trueParts;
  for(int i = 0; i < event.nTrueParticles; i++)trueParts.push_back(event.TrueParticles[i]);
  const AnaBeamPD* beam = static_cast<AnaBeamPD*>(event.Beam);
  const AnaTrueParticleB* beamTrue = (beam && beam->BeamParticle) ? static_cast<AnaTrueParticleB*>(beam->BeamParticle->TrueObject) : nullptr;
  const int isMC = event.GetIsMC() ? 1 : 0;

  // Lo que se calcula una vez por evento (beam slice, vértice, nube de puntos 3D, mapas de partículas)
  muCandUtils::EventContext ctx;
  muCandUtils::MakeContext(ctx, parts, event.UnassociatedHits, trueParts, beamTrue, isMC);
  box.HasVertex = ctx.hasVtx;
  if(ctx.hasVtx){ box.VertexPos[0] = ctx.vtx.X(); box.VertexPos[1] = ctx.vtx.Y(); box.VertexPos[2] = ctx.vtx.Z(); }

  // Cada partícula: Compute dice si es candidato (longitud y extremos válidos, p de rango de mu en [100, 600] MeV/c)
  // y calcula sus variables; después, la preselección geométrica del módulo
  for(const AnaParticlePD* p : parts){
    muCandUtils::Candidate c;
    if(!muCandUtils::Compute(p, ctx, c))continue;
    box.NPreselected++;
    // un extremo a >= 10 cm de las caras y uno a <= 150 cm del eje (lo que el módulo guarda completo)
    if(std::max(c.facesA, c.facesB) < 10 || muCandCuts::AxisMin(c) > 150)continue;
    box.Candidates.push_back(c);
    box.Population.push_back(c.inSlice ? 0 : 1);   // a qué rama pertenece
    box.Level.push_back(0);                          // todavía no ha pasado ningún corte de su rama
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
bool MuCandidateCutBase::Apply(AnaEventC& event, ToyBoxB& boxB) const{
//**************************************************

  // Contabilidad común a todos los cortes sobre candidatos. Este corte ocupa la posición _index en la cadena de la
  // rama _pop. Un candidato "llega" a este corte si es de esa rama y ha pasado todos los anteriores, es decir, si su
  // nivel es _index - 1. Si además pasa este corte, su nivel sube a _index. El evento pasa el corte si al menos un
  // candidato llega a _index. (Los candidatos que se quedaron atrás conservan su nivel: es lo que se guarda en
  // cand_level.)
  (void)event;
  ToyBoxMuCand& box = *static_cast<ToyBoxMuCand*>(&boxB);
  bool anyPasses = false;
  for(size_t i = 0; i < box.Candidates.size(); i++){
    if(box.Population[i] != _pop)continue;           // candidato de la otra rama
    if(box.Level[i] != _index - 1)continue;          // no pasó algún corte anterior de esta rama
    if(PassCandidate(box.Candidates[i])){            // el corte concreto (clase derivada)
      box.Level[i] = _index;
      anyPasses = true;
    }
  }
  return anyPasses;
}

//**************************************************
void secondaryKaonProductionXSSelection::InitializeEvent(AnaEventC& eventC){
//**************************************************

  AnaEventB& event = *static_cast<AnaEventB*>(&eventC);

  // EventBox de los kaones: grupos de objetos para los sistemáticos (hoy apagados)
  if (!event.EventBoxes[EventBoxId::kEventBoxKaon])
    event.EventBoxes[EventBoxId::kEventBoxKaon] = new EventBoxKaon();

  boxUtils::FillKaonXS(event);
}
