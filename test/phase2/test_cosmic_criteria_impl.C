// Implementación de test_cosmic_criteria.C. Un evento: traza del haz a lo largo de z (0-200 cm), mu
// de K+ fuera de la beam slice que empieza a 2 cm de la mitad de la traza con la longitud de
// 236 MeV/c, un cósmico vertical que pasa a 15 cm de la traza y otro lejos.
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
  AnaParticlePD* Reco(AnaSpillPD* s, int uid, bool slice, AnaTrueParticleB* match, TVector3 a, TVector3 b){
    AnaParticlePD* p = new AnaParticlePD();
    p->UniqueID = uid; p->ParentID = -1; p->BeamOrigin = slice; p->isPandora = slice;
    p->TrueObject = match;
    for(int i = 0; i < 3; i++){ p->PositionStart[i] = a[i]; p->PositionEnd[i] = b[i]; }
    p->Length = (b - a).Mag();
    p->NHits = 10;
    static_cast<AnaBunchPD*>(s->Bunches[0])->Particles.push_back(p);
    return p;
  }
}

int test_cosmic_criteria_impl(){
  using namespace ccrit;

  // geometría
  Check("punto-segmento", PointSeg(TVector3(3,4,5), {TVector3(0,0,0), TVector3(0,0,10)}), 5);
  Check("punto-segmento fuera del extremo", PointSeg(TVector3(0,0,13), {TVector3(0,0,0), TVector3(0,0,10)}), 3);
  Check("segmentos que se cruzan a 2 cm", SegSeg({TVector3(-5,2,5), TVector3(5,2,5)}, {TVector3(0,0,0), TVector3(0,0,10)}), 2);
  Check("segmentos paralelos", SegSeg({TVector3(1,0,2), TVector3(1,0,8)}, {TVector3(0,0,0), TVector3(0,0,10)}), 1);
  Check("segmentos alejados por los extremos", SegSeg({TVector3(0,0,13), TVector3(0,0,20)}, {TVector3(0,0,0), TVector3(0,0,10)}), 3);

  ND::params().SetReadParamOverrideFilePointPassed();   // ComputeRangeMomentum lee parámetros

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
  AnaTrueParticlePD* mu = True(s, 3, 2, -13);
  AnaTrueParticlePD* cos = True(s, 50, 0, 13);

  AnaParticlePD* tr = Reco(s, 0, true, pi, TVector3(0,0,0), TVector3(0,0,200));
  for(int z = 0; z <= 200; z++){ AnaTrajectoryPointPD tp; tp.Position_NoSCE.SetXYZ(0,0,z); tr->TrjPoints.push_back(tp); }
  Reco(s, 1, false, mu,  TVector3(2,0,100), TVector3(2+L,0,100));
  Reco(s, 2, false, cos, TVector3(15,300,60), TVector3(15,-300,60));
  Reco(s, 3, false, cos, TVector3(250,300,60), TVector3(250,-300,60));

  CosmicCriteria cc;
  cc.Process(s);
  cc.Print();

  Check("mu de K fuera", cc.muons.size(), 1);
  Check("mu: distancia de extremos", cc.muons[0].dEnd, 2);
  Check("mu: momento por rango", cc.muons[0].p, kMuP, 1.);
  Check("jerarquías en resumen", cc.summaryHier, 3);
  auto idx = [&cc](const std::string& n){ for(size_t c = 0; c < cc.crit.size(); c++)if(cc.crit[c].name == n)return (int)c; return -1; };
  Check("(a) R = 5: solo el mu", cc.extraHier[idx("(a) cilindro R = 5 cm")], 1);
  Check("(a) R = 20: mu y cósmico cercano", cc.extraHier[idx("(a) cilindro R = 20 cm")], 2);
  Check("(b) R = 50: solo el mu", cc.extraHier[idx("(b) extremos R = 50 cm")], 1);
  Check("(c) +-10 %: solo el mu", cc.extraHier[idx("(c) p rango = 236 +- 10 %")], 1);

  std::cout << (nErr ? "TEST CON FALLOS" : "TEST OK") << std::endl;
  return nErr;
}
