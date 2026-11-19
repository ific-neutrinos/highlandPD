// Implementación de test_kmu_pairing.C. Un evento por caso, con un pi+ del haz (traza primaria
// de z = 0 a 100 cm), un K+ secundario reconstruido como hija del pi (de z = 100 a 130) que se
// desintegra en un mu+ de 236 MeV/c, y el mu en distintas situaciones:
//   1. hija única de la K, empieza en su final                 -> "hija" y "vecino" lo encuentran
//   2. hija del pi, empieza a 3 cm del final de la K           -> solo "vecino"
//   3. fuera de la beam slice (resumen), al revés, fin a 2 cm  -> solo "vecino" (no "solo slice")
//   4. hija única de la K pero a 20 cm                         -> nadie
//   5. la K es el primario de la beam slice                    -> no cuenta como señal
//   6. sus hits dentro de la traza de la K                     -> nadie
// En el evento 1 hay además un protón hija del pi con una traza de 236 MeV/c (como mu) que empieza
// en su final: fondo para "vecino".
#include "kmu_pairing_minitrees_impl.C"
#include <iostream>

namespace {
  int nErr = 0;
  void Check(const char* what, long got, long expected){
    if(got == expected)return;
    std::cout << "FALLA " << what << ": " << got << " (esperado " << expected << ")" << std::endl;
    nErr++;
  }
  AnaTrueParticlePD* True(AnaSpillPD* s, int id, int parent, int pdg, AnaTrueParticleB::ProcessEnum end = AnaTrueParticleB::Unknown){
    AnaTrueParticlePD* t = new AnaTrueParticlePD();
    t->ID = id; t->ParentID = parent; t->PDG = pdg; t->ProcessEnd = end;
    t->NHitsTrue[0] = t->NHitsTrue[1] = t->NHitsTrue[2] = 10;
    s->TrueParticles.push_back(t);
    for(AnaTrueParticleB* p : s->TrueParticles)if(p->ID == parent)static_cast<AnaTrueParticle*>(p)->Daughters.push_back(id);
    return t;
  }
  AnaParticlePD* Reco(AnaSpillPD* s, int uid, int parent, bool slice, AnaTrueParticleB* match, TVector3 a, TVector3 b){
    AnaParticlePD* p = new AnaParticlePD();
    p->UniqueID = uid; p->ParentID = parent; p->BeamOrigin = slice; p->isPandora = slice && parent == -1;
    p->Type = AnaParticlePD::kTrack;
    p->TrueObject = match;
    for(int i = 0; i < 3; i++){ p->PositionStart[i] = a[i]; p->PositionEnd[i] = b[i]; }
    p->Length = (b - a).Mag();
    p->NHits = 10;
    for(AnaParticleB* q : static_cast<AnaBunchPD*>(s->Bunches[0])->Particles)if(q->UniqueID == parent)
      static_cast<AnaParticlePD*>(q)->DaughtersIDs.push_back(uid);
    static_cast<AnaBunchPD*>(s->Bunches[0])->Particles.push_back(p);
    return p;
  }
}

int test_kmu_pairing_impl(){
  using namespace kmu;
  ND::params().SetReadParamOverrideFilePointPassed();

  double L = 0;
  for(double l = 1; l < 200; l += 0.1)if(pdAnaUtils::ComputeRangeMomentum(l, 13) >= 0.233){ L = l; break; }
  std::cout << "Longitud de un mu de 233 MeV/c: " << L << " cm" << std::endl;

  KMuPairing km;
  for(int c = 1; c <= 6; c++){
    AnaSpillPD* s = new AnaSpillPD();
    AnaEventInfoPD* info = new AnaEventInfoPD(); info->IsMC = true; s->EventInfo = info;
    AnaBeamPD* beam = new AnaBeamPD(); beam->BeamParticle = new AnaParticlePD(); s->Beam = beam;
    s->Bunches.push_back(new AnaBunchPD());
    AnaTrueParticlePD* pi = True(s, 1, 0, 211, AnaTrueParticleB::Unknown);
    beam->BeamParticle->TrueObject = pi;
    AnaTrueParticlePD* k  = True(s, 2, 1, 321, AnaTrueParticleB::Decay);
    AnaTrueParticlePD* mu = True(s, 3, 2, -13);
    AnaTrueParticlePD* pr = True(s, 4, 1, 2212);
    const TVector3 kEnd(0,0,130);
    const TVector3 dir(1,0,0);
    if(c == 5){
      Reco(s, 0, -1, true, k, TVector3(0,0,100), kEnd);
      Reco(s, 1, 0, true, mu, kEnd, kEnd + L*dir);
    }
    else{
      Reco(s, 0, -1, true, pi, TVector3(0,0,0), TVector3(0,0,100));
      AnaParticlePD* kr = Reco(s, 1, 0, true, k, TVector3(0,0,100), kEnd);
      if(c == 1){
        Reco(s, 2, 1, true, mu, kEnd, kEnd + L*dir);
        Reco(s, 3, 0, true, pr, TVector3(0,0,100), TVector3(0,L,100));
        Reco(s, 4, 0, true, nullptr, TVector3(0,L,100), TVector3(0,L,100) + L*dir);
      }
      if(c == 2)Reco(s, 2, 0, true, mu, kEnd + TVector3(0,3,0), kEnd + TVector3(0,3,0) + L*dir);
      if(c == 3)Reco(s, 2, -1, false, mu, kEnd + TVector3(0,2,0) + L*dir, kEnd + TVector3(0,2,0));
      if(c == 4)Reco(s, 2, 1, true, mu, kEnd + TVector3(0,20,0), kEnd + TVector3(0,20,0) + L*dir);
      if(c == 6)for(int i = 0; i < 10; i++){ AnaHitPD h; h.TrueID = mu->ID; kr->Hits[2].push_back(h); }
    }
    km.Process(s);
  }
  km.Print();

  Check("pares K-mu con la K en la beam slice", km.nPairs, 6);
  Check("K primaria", km.nKPrimary, 1);
  Check("mu hija única", km.rel.n[kOnlyDaughter], 2);   // casos 1 y 4
  Check("mu en la slice, otra relación", km.rel.n[kSliceOther], 1);
  Check("mu fuera en resumen", km.rel.n[kRelOutSummary], 1);
  Check("mu fuera, al revés", km.rel.nReversed[kRelOutSummary], 1);
  Check("mu dentro de la K", km.rel.n[kInK], 1);
  Check("K de señal posibles", km.nSignalK, 5);
  Check("hija: señal", km.sigFound[kHija], 1);
  Check("vecino: señal", km.sigFound[kVecino], 3);
  Check("vecino solo slice: señal", km.sigFound[kVecinoSlice], 2);
  Check("hija: fondo", km.bkgPairs[kHija], 0);
  Check("vecino: fondo (protón con la traza de su final)", km.bkgPairs[kVecino], 1);

  std::cout << (nErr ? "TEST CON FALLOS" : "TEST OK") << std::endl;
  return nErr;
}
