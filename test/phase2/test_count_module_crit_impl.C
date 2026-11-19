// Implementación de test_count_module_crit.C. Cuatro eventos de MC:
//   1. normal con beam slice: una jerarquía completa por cercanía (bit 1), otra completa solo por el
//      candidato a mu (bit 4, dos PFParticles) y otra en resumen (0)
//   2. normal sin beam slice: una jerarquía con el bit 4, completa, con un mu de K+
//   3. normal sin beam slice, sin cósmicos guardados (NBeamRegionSpacePoints rellenado)
//   4. con volcado y sin beam slice: dos mu de K+ de la cadena del haz; el primero con su K en
//      una jerarquía con el bit 4, el segundo en una sin él; una tercera jerarquía con el bit 4
//      sin mu
#include "count_module_minitrees_impl.C"
#include <iostream>

namespace {
  int nErr = 0;
  void Check(const char* what, double got, double expected){
    if(std::fabs(got - expected) < 1e-9)return;
    std::cout << "FALLA " << what << ": " << got << " (esperado " << expected << ")" << std::endl;
    nErr++;
  }
  AnaTrueParticlePD* True(AnaSpillPD* s, int id, int parent, int pdg, AnaTrueParticleB::ProcessEnum end = AnaTrueParticleB::Unknown){
    AnaTrueParticlePD* t = new AnaTrueParticlePD();
    t->ID = id; t->ParentID = parent; t->PDG = pdg; t->ProcessEnd = end;
    s->TrueParticles.push_back(t);
    return t;
  }
  // PFParticle; complete: con un punto de trayectoria
  AnaParticlePD* Reco(AnaSpillPD* s, int uid, int parent, bool slice, AnaTrueParticleB* match, int crit, bool complete){
    AnaParticlePD* p = new AnaParticlePD();
    p->UniqueID = uid; p->ParentID = parent; p->BeamOrigin = slice; p->isPandora = slice && parent == -1;
    p->TrueObject = match; p->NHits = 10;
    p->CosmicCriteria = slice ? -999 : crit;
    if(complete)p->TrjPoints.push_back(AnaTrajectoryPointPD());
    static_cast<AnaBunchPD*>(s->Bunches[0])->Particles.push_back(p);
    return p;
  }
  AnaSpillPD* Spill(bool fullDump, AnaTrueParticlePD*& beam){
    AnaSpillPD* s = new AnaSpillPD();
    AnaEventInfoPD* info = new AnaEventInfoPD(); info->IsMC = true; info->FullDump = fullDump;
    info->NBeamRegionSpacePoints = 50; s->EventInfo = info;
    s->Bunches.push_back(new AnaBunchPD());
    beam = True(s, 1, 0, 211);
    return s;
  }
}

int test_count_module_crit_impl(){
  const int kNear = AnaParticlePD::kCosmicNearBeamSlice, kMu = AnaParticlePD::kCosmicMuCandidate;
  CritCounter cc;
  AnaTrueParticlePD* beam = nullptr;

  AnaSpillPD* s1 = Spill(false, beam);
  Reco(s1, 0, -1, true, beam, -999, true);
  Reco(s1, 10, 1000, false, nullptr, kNear, true);
  Reco(s1, 20, 1000, false, nullptr, kMu, true);
  Reco(s1, 21, 20, false, nullptr, kMu, true);
  Reco(s1, 30, 1000, false, nullptr, 0, false);
  cc.Process(s1, beam);

  AnaSpillPD* s2 = Spill(false, beam);
  True(s2, 2, 1, 321, AnaTrueParticleB::Decay);
  AnaTrueParticlePD* mu2 = True(s2, 3, 2, -13);
  Reco(s2, 10, 1000, false, mu2, kMu, true);
  cc.Process(s2, beam);

  AnaSpillPD* s3 = Spill(false, beam);
  cc.Process(s3, beam);

  AnaSpillPD* s4 = Spill(true, beam);
  AnaTrueParticlePD* k1  = True(s4, 2, 1, 321, AnaTrueParticleB::Decay);
  AnaTrueParticlePD* mu1 = True(s4, 3, 2, -13);
  True(s4, 4, 1, 321, AnaTrueParticleB::Decay);
  AnaTrueParticlePD* muB = True(s4, 5, 4, -13);
  Reco(s4, 10, 1000, false, k1, kMu, true);
  Reco(s4, 11, 10, false, mu1, kMu, true);
  Reco(s4, 20, 1000, false, muB, 0, true);
  Reco(s4, 30, 1000, false, nullptr, kMu | kNear, true);
  cc.Process(s4, beam);

  cc.Print();
  Check("eventos normales", cc.events, 3);
  Check("con beam slice: jerarquías completas", cc.compH.sum, 2);
  Check("con beam slice: PFParticles completos", cc.compP.sum, 3);
  Check("con beam slice: solo por el candidato a mu", cc.muOnlyH.sum, 1);
  Check("con beam slice: PFParticles solo por el candidato a mu", cc.muOnlyP.sum, 2);
  Check("sin beam slice: eventos", cc.nsNormal, 2);
  Check("sin beam slice: con alguna jerarquía", cc.nsWith, 1);
  Check("sin beam slice: incoherencias", cc.nsBadCrit + cc.nsIncomplete, 0);
  Check("sin beam slice: jerarquías con mu de K+", cc.nsHierMuK, 1);
  Check("volcado: eventos", cc.fdNS, 1);
  Check("volcado: mu de K+", cc.fdMu, 2);
  Check("volcado: mu con el bit 4", cc.fdMuBit4, 1);
  Check("volcado: con su K en la misma jerarquía", cc.fdMuBit4K, 1);
  Check("volcado: jerarquías con el bit 4", cc.fdBit4H, 2);
  Check("volcado: de ellas con mu de K+", cc.fdBit4HMuK, 1);

  std::cout << (nErr ? "TEST CON FALLOS" : "TEST OK") << std::endl;
  return nErr;
}
