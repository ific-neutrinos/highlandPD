// Implementación de test_cosmic_criteria.C. Un evento con una traza del haz a lo largo de z
// (0-200 cm) y, fuera de la beam slice:
//   slice 5: el mu de una K+ (empieza a 2 cm de la mitad de la traza, longitud de 236 MeV/c) y
//            otro PFParticle lejos, en otra jerarquía del mismo slice
//   slice 6: un cósmico vertical que pasa a 15 cm de la traza (puntuación 0,05)
//   slice 7: un pi+ de la cadena del haz lejos de la traza, con puntuación 0,9
//   slices 10-19: cósmicos lejos, puntuaciones 0,00-0,09
#include "cosmic_criteria_minitrees_impl.C"
#include <iostream>

namespace {
  int nErr = 0;
  void Check(const char* what, double got, double expected, double tol = 1e-6){
    if(std::fabs(got - expected) <= tol)return;
    std::cout << "FALLA " << what << ": " << got << " (esperado " << expected << ")" << std::endl;
    nErr++;
  }
  AnaTrueParticlePD* True(AnaSpillPD* s, int id, int parent, int pdg, AnaTrueParticleB::ProcessEnum end = AnaTrueParticleB::Unknown){
    AnaTrueParticlePD* t = new AnaTrueParticlePD();
    t->ID = id; t->ParentID = parent; t->PDG = pdg; t->ProcessEnd = end;
    s->TrueParticles.push_back(t);
    return t;
  }
  AnaParticlePD* Reco(AnaSpillPD* s, int uid, bool slice, AnaTrueParticleB* match, TVector3 a, TVector3 b,
                      int sliceID = -1, float score = -999){
    AnaParticlePD* p = new AnaParticlePD();
    p->UniqueID = uid; p->ParentID = -1; p->BeamOrigin = slice; p->isPandora = slice;
    p->TrueObject = match;
    p->SliceID = sliceID;
    p->BeamCosmicScore = score;
    for(int i = 0; i < 3; i++){ p->PositionStart[i] = a[i]; p->PositionEnd[i] = b[i]; }
    p->Length = (b - a).Mag();
    p->NHits = 10;
    p->Type = AnaParticlePD::kTrack;
    static_cast<AnaBunchPD*>(s->Bunches[0])->Particles.push_back(p);
    return p;
  }
}

int test_cosmic_criteria_impl(){
  using namespace ccrit;
  ND::params().SetReadParamOverrideFilePointPassed();   // MuRangeP lee parámetros

  // geometría
  Check("punto-segmento", PointSeg(TVector3(3,4,5), {TVector3(0,0,0), TVector3(0,0,10)}), 5);
  Check("punto-segmento fuera del extremo", PointSeg(TVector3(0,0,13), {TVector3(0,0,0), TVector3(0,0,10)}), 3);
  Check("segmentos que se cruzan a 2 cm", SegSeg({TVector3(-5,2,5), TVector3(5,2,5)}, {TVector3(0,0,0), TVector3(0,0,10)}), 2);
  Check("segmentos paralelos", SegSeg({TVector3(1,0,2), TVector3(1,0,8)}, {TVector3(0,0,0), TVector3(0,0,10)}), 1);
  Check("segmentos alejados por los extremos", SegSeg({TVector3(0,0,13), TVector3(0,0,20)}, {TVector3(0,0,0), TVector3(0,0,10)}), 3);

  // longitud del mu de 236 MeV/c
  double L = 0;
  for(double l = 1; l < 200; l += 0.1)if(MuRangeP(l) >= kMuP){ L = l; break; }
  std::cout << "Longitud de un mu de 236 MeV/c (momento por rango): " << L << " cm" << std::endl;
  if(L < 30 || L > 80){ std::cout << "FALLA longitud del mu fuera de lo esperado" << std::endl; nErr++; }

  AnaSpillPD* s = new AnaSpillPD();
  AnaEventInfoPD* info = new AnaEventInfoPD(); info->IsMC = true; s->EventInfo = info;
  AnaBeamPD* beam = new AnaBeamPD(); beam->BeamParticle = new AnaParticlePD(); s->Beam = beam;
  s->Bunches.push_back(new AnaBunchPD());
  AnaTrueParticlePD* pi = True(s, 1, 0, 211);
  beam->BeamParticle->TrueObject = pi;
  True(s, 2, 1, 321, AnaTrueParticleB::Decay);
  AnaTrueParticlePD* mu  = True(s, 3, 2, -13);
  AnaTrueParticlePD* pi2 = True(s, 4, 1, 211);
  AnaTrueParticlePD* cos = True(s, 50, 0, 13);

  AnaParticlePD* tr = Reco(s, 0, true, pi, TVector3(0,0,0), TVector3(0,0,200), 1);
  for(int z = 0; z <= 200; z++){ AnaTrajectoryPointPD tp; tp.Position_NoSCE.SetXYZ(0,0,z); tr->TrjPoints.push_back(tp); }
  Reco(s, 1, false, mu,  TVector3(2,0,100), TVector3(2+L,0,100), 5);
  Reco(s, 2, false, cos, TVector3(300,300,60), TVector3(300,-300,60), 5);
  Reco(s, 3, false, cos, TVector3(15,300,60), TVector3(15,-300,60), 6, 0.05);
  Reco(s, 4, false, pi2, TVector3(250,0,300), TVector3(250,0,320), 7, 0.9);
  for(int i = 0; i < 10; i++)Reco(s, 10+i, false, cos, TVector3(-300,300,30*i), TVector3(-300,-300,30*i), 10+i, 0.01*i);
  // (g): un trozo corto en otro slice que acaba a 5 cm del extremo lejano del mu (x = 2+L)
  Reco(s, 30, false, cos, TVector3(2+L,8,100), TVector3(2+L,5,100), 30);

  CosmicCriteria cc;
  cc.Process(s);
  cc.Print();

  // eficiencia y coste de cada criterio
  std::map<std::string, std::pair<int,int>> res;   // nombre -> (objetivos que pasan, jerarquías extra)
  for(const Crit& c : cc.Criteria()){
    int tp = 0, xh = 0;
    for(const Event& e : cc.events)for(const Slice& sl : e.slices)for(const Hier& h : sl.hiers){
      const bool pass = c.pass(sl, h);
      if(pass && !h.complete)xh++;
      for(const PartV& p : h.parts)if(p.target && pass)tp++;
    }
    res[c.name] = {tp, xh};
  }
  auto get = [&res](const std::string& prefix){
    for(auto& r : res)if(r.first.compare(0, prefix.size(), prefix) == 0)return r.second;
    std::cout << "FALLA no está el criterio " << prefix << std::endl; nErr++;
    return std::make_pair(-1,-1);
  };
  Check("slices de fuera", cc.events[0].slices.size(), 14);
  Check("slices con objetivos y puntuación (el pi)", cc.scoreTarget.size(), 1);
  Check("slices con objetivos sin puntuación (el mu)", cc.nScoreNone[1], 1);
  Check("(a) R = 20: objetivos", get("(a) cilindro R = 20").first, 1);
  Check("(a) R = 20: mu y cósmico cercano", get("(a) cilindro R = 20").second, 2);
  Check("(b) R = 20: solo el mu", get("(b) extremos R = 20").second, 1);
  Check("(c): el mu", get("(c) p rango").first, 1);
  Check("(d) p90: solo el pi", get("(d) puntuación > p90").first, 1);
  Check("(d) p90: el pi y un slice de fondo (el 10 % de 11)", get("(d) puntuación > p90").second, 2);
  Check("(e) R = 20: el slice del mu entero", get("(e) slice cercano R = 20").second, 2);
  Check("(d) p95 o (e): mu y pi", get("(d) p95 o (e)").first, 2);
  Check("(d) p95 o (e): 3 jerarquías", get("(d) p95 o (e)").second, 3);
  Check("(f) R = 20: el mu", get("(f) candidato a mu R = 20").second, 1);
  Check("(f) R = 20: objetivo (el mu)", get("(f) candidato a mu R = 20").first, 1);
  Check("(g) R = 10: el mu y el trozo a 5 cm", get("(g) entorno del mu R = 10").second, 2);

  // Emulación de CosmicsMode 2 (puntos a < 50 cm de un extremo de la beam slice, z = 0 o 200):
  // un cósmico con un space point a 40 cm del final y otro con trayectoria a 60 cm, los dos
  // completos en el minitree. En un evento normal el segundo sale como desacuerdo (completo pero
  // emulado en resumen); con FullDump la emulación deja completo solo el primero.
  for(bool fd : {false, true}){
    AnaSpillPD* s2 = new AnaSpillPD();
    AnaEventInfoPD* i2 = new AnaEventInfoPD(); i2->IsMC = true; i2->FullDump = fd; s2->EventInfo = i2;
    AnaBeamPD* b2 = new AnaBeamPD(); b2->BeamParticle = new AnaParticlePD(); s2->Beam = b2;
    s2->Bunches.push_back(new AnaBunchPD());
    AnaTrueParticlePD* p2 = True(s2, 1, 0, 211);
    b2->BeamParticle->TrueObject = p2;
    Reco(s2, 0, true, p2, TVector3(0,0,0), TVector3(0,0,200), 1);
    AnaParticlePD* c1 = Reco(s2, 1, false, nullptr, TVector3(40,300,200), TVector3(40,-300,200), 5);
    AnaParticlePD* c2 = Reco(s2, 2, false, nullptr, TVector3(60,300,200), TVector3(60,-300,200), 6);
    AnaSpacePointPD sp; sp.Position[0] = 40; sp.Position[1] = 0; sp.Position[2] = 200;
    c1->SpacePoints.push_back(sp);
    AnaTrajectoryPointPD tp; tp.Position_NoSCE.SetXYZ(60,0,200);
    c2->TrjPoints.push_back(tp);
    CosmicCriteria c;
    c.Process(s2);
    int nEmul = 0, nComplete = 0;
    for(const Slice& sl : c.events[0].slices)for(const Hier& h : sl.hiers){ nEmul += h.emulated; nComplete += h.complete; }
    Check("emulación: completas en el minitree", nComplete, 2);
    Check("emulación: completas emuladas", nEmul, 1);
    Check("emulación: marca de volcado", c.events[0].fullDump, fd);
    Check("validación: completas y emuladas", c.emulAgree[1][1], fd ? 0 : 1);
    Check("validación: completas emuladas en resumen", c.emulAgree[1][0], fd ? 0 : 1);
    if(fd)c.Print();
  }

  std::cout << (nErr ? "TEST CON FALLOS" : "TEST OK") << std::endl;
  return nErr;
}
