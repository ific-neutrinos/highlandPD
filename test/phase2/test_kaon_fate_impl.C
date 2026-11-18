// Implementación de test_kaon_fate.C: eventos de MC construidos a mano con un pi+ del haz (o un
// K+), un K+ secundario que se desintegra y su mu+, en distintas situaciones de reconstrucción.
#include "pdDataClasses.hxx"
#include <iostream>
#include <string>

#include "kaon_fate_minitrees_impl.C"

namespace {
  AnaTrueParticlePD* True(AnaSpillPD* s, int id, int parent, int pdg, int nhits2, AnaTrueParticleB::ProcessEnum end = AnaTrueParticleB::Unknown){
    AnaTrueParticlePD* t = new AnaTrueParticlePD();
    t->ID = id; t->ParentID = parent; t->PDG = pdg; t->ProcessEnd = end;
    t->NHitsTrue[0] = t->NHitsTrue[1] = t->NHitsTrue[2] = nhits2;
    s->TrueParticles.push_back(t);
    for(AnaTrueParticleB* p : s->TrueParticles)if(p->ID == parent)static_cast<AnaTrueParticle*>(p)->Daughters.push_back(id);
    return t;
  }
  // Partícula reco con nhits hits en el plano 2 de la partícula true "hitsOf"
  AnaParticlePD* Reco(AnaSpillPD* s, int uid, int parent, int gen, bool slice, bool primary, AnaTrueParticleB* match,
                      int hitsOf, int nhits, bool clearCosmic = false){
    AnaParticlePD* p = new AnaParticlePD();
    p->UniqueID = uid; p->ParentID = parent; p->Generation = gen;
    p->BeamOrigin = slice; p->isPandora = primary; p->IsClearCosmic = clearCosmic;
    p->TrueObject = match;
    p->NHits = nhits;
    for(int i = 0; i < nhits; i++){ AnaHitPD h; h.TrueID = hitsOf; p->Hits[2].push_back(h); }
    static_cast<AnaBunchPD*>(s->Bunches[0])->Particles.push_back(p);
    return p;
  }
  AnaSpillPD* Spill(int beamPDG){
    AnaSpillPD* s = new AnaSpillPD();
    AnaEventInfoPD* info = new AnaEventInfoPD();
    info->IsMC = true;
    s->EventInfo = info;
    AnaBeamPD* beam = new AnaBeamPD();
    beam->BeamParticle = new AnaParticlePD();
    s->Beam = beam;
    s->Bunches.push_back(new AnaBunchPD());
    beam->BeamParticle->TrueObject = True(s, 1, 0, beamPDG, 50);
    return s;
  }
  int nErr = 0;
  void Check(const char* what, long got, long expected){
    if(got == expected)return;
    std::cout << "FALLA " << what << ": " << got << " (esperado " << expected << ")" << std::endl;
    nErr++;
  }
}

int test_kaon_fate_impl(){
  using namespace kfate;
  KaonFate kf;
  const AnaTrueParticleB::ProcessEnum dec = AnaTrueParticleB::Decay;

  // A: K hija del primario y mu hija de la K, los dos en la beam slice (40 % de los hits del mu)
  { AnaSpillPD* s = Spill(211);
    AnaTrueParticlePD* k = True(s, 2, 1, 321, 20, dec); AnaTrueParticlePD* mu = True(s, 3, 2, -13, 10);
    Reco(s, 0, -1, 0, true, true,  s->TrueParticles[0], 1, 30);
    Reco(s, 1,  0, 1, true, false, k, 2, 20);
    Reco(s, 2,  1, 2, true, false, mu, 3, 4);
    kf.Process(s); }
  // B: K en la beam slice; mu en un cósmico en resumen (sin hits) marcado IsClearCosmic
  { AnaSpillPD* s = Spill(211);
    AnaTrueParticlePD* k = True(s, 2, 1, 321, 20, dec); AnaTrueParticlePD* mu = True(s, 3, 2, -13, 10);
    Reco(s, 0, -1, 0, true, true, s->TrueParticles[0], 1, 30);
    Reco(s, 1,  0, 1, true, false, k, 2, 20);
    Reco(s, 7, -1, -1, false, false, mu, 3, 0, true);
    kf.Process(s); }
  // C: K y mu unidos en una traza de la beam slice emparejada con el mu: la K sale "dentro de otra"
  { AnaSpillPD* s = Spill(211);
    True(s, 2, 1, 321, 20, dec); AnaTrueParticlePD* mu = True(s, 3, 2, -13, 10);
    AnaParticlePD* p = Reco(s, 0, -1, 0, true, true, s->TrueParticles[0], 1, 30);
    AnaParticlePD* m = Reco(s, 1,  0, 1, true, false, mu, 3, 10);
    for(int i = 0; i < 5; i++){ AnaHitPD h; h.TrueID = 2; m->Hits[2].push_back(h); }
    (void)p;
    kf.Process(s); }
  // D: evento sin beam slice (solo un cósmico)
  { AnaSpillPD* s = Spill(211);
    AnaTrueParticlePD* k = True(s, 2, 1, 321, 20, dec); True(s, 3, 2, -13, 10);
    Reco(s, 7, -1, -1, false, false, k, 2, 20);
    kf.Process(s); }
  // D2: evento sin beam slice con volcado completo: la K en un cósmico completo
  { AnaSpillPD* s = Spill(211);
    static_cast<AnaEventInfoPD*>(s->EventInfo)->FullDump = true;
    AnaTrueParticlePD* k = True(s, 2, 1, 321, 20, AnaTrueParticleB::kaonplusInelastic);
    Reco(s, 7, -1, -1, false, false, k, 2, 20);
    kf.Process(s); }
  // E: K sin reconstruir (deja hits) y mu sin hits; y un cósmico completo con un K de fuera de la cadena
  { AnaSpillPD* s = Spill(211);
    True(s, 2, 1, 321, 20, dec); True(s, 3, 2, -13, 0);
    AnaTrueParticlePD* kc = True(s, 99, 0, 321, 20, dec);   // no desciende del haz: no cuenta
    Reco(s, 0, -1, 0, true, true, s->TrueParticles[0], 1, 30);
    Reco(s, 8, -1, -1, false, false, kc, 99, 20);
    kf.Process(s); }
  // F: el haz es un K+ y es el primario; inelástico, sin mu
  { AnaSpillPD* s = Spill(321);
    Reco(s, 0, -1, 0, true, true, s->TrueParticles[0], 1, 50);
    kf.Process(s); }

  kf.Print();

  Check("K del haz", kf.beamKaon.all, 1);
  Check("K del haz primario", kf.beamKaon.fate[kSlicePrimary], 1);
  Check("K secundarios", kf.secKaon.all, 6);
  Check("K secundarios que se desintegran", kf.secKaonDecay.all, 5);
  Check("K hija del primario", kf.secKaon.fate[kSliceDaughter], 2);
  Check("K dentro de otra (slice)", kf.secKaon.fate[kSliceInOther], 1);
  Check("K evento sin slice", kf.secKaon.fate[kNoSlice], 2);
  Check("K sin slice con volcado completo", kf.secKaonNoSlice.all, 1);
  Check("  en un cósmico completo", kf.secKaonNoSlice.fate[kOutComplete], 1);
  Check("K no reconstruida", kf.secKaon.fate[kNotReco], 1);
  Check("mu", kf.muon.all, 5);
  Check("mu generación >= 2", kf.muon.fate[kSliceDeeper], 1);
  Check("mu hija del primario", kf.muon.fate[kSliceDaughter], 1);
  Check("mu cósmico en resumen", kf.muon.fate[kOutSummary], 1);
  Check("mu IsClearCosmic", kf.muon.fate[kOutClearCosmic], 1);
  Check("mu evento sin slice", kf.muon.fate[kNoSlice], 1);
  Check("mu sin hits", kf.muon.fate[kNoHits], 1);
  Check("K y mu en la slice", kf.bothInSlice, 2);
  Check("mu hija de la K", kf.muDaughterOfKaonReco, 1);
  Check("tabla slice x fuera", kf.joint[kCSlice][kCOut], 1);
  Check("tabla no reco x sin hits", kf.joint[kCNotReco][kCNoHits], 1);
  Check("fracción de hits del mu en A", std::lround(100*kf.muon.hitFracSum[kSliceDeeper]), 40);

  std::cout << (nErr ? "TEST CON FALLOS" : "TEST OK") << std::endl;
  return nErr;
}
