// Implementación de test_beam_axis.C: eventos con la instrumentación, el primario de la beam
// slice y la trayectoria true del haz en posiciones y direcciones conocidas.
#include "pdDataClasses.hxx"
#include <cmath>
#include <iostream>

#include "beam_axis_minitrees_impl.C"

namespace {
  int nErr = 0;
  void Check(const char* what, double got, double expected, double tol = 1e-3){
    if(std::fabs(got - expected) <= tol)return;
    std::cout << "FALLA " << what << ": " << got << " (esperado " << expected << ")" << std::endl;
    nErr++;
  }
  void Set(Float_t* v, double x, double y, double z){ v[0] = x; v[1] = y; v[2] = z; }

  // Dirección d (sin normalizar); la instrumentación acaba en z = -50 y el primario empieza en
  // zStart, los dos sobre la recta que pasa por (x0, y0, 0). primTrue: el primario es el haz true
  // primDir: dirección del primario si no es d (para probar el corte del cos)
  AnaSpillPD* Spill(bool mc, float mom, double x0, double y0, TVector3 d, double zStart, bool primTrue, bool withTraj,
                    TVector3 primDir = TVector3()){
    d = d.Unit();
    if(primDir.Mag() == 0)primDir = d;
    primDir = primDir.Unit();
    auto at = [&](double z){ return TVector3(x0, y0, 0) + d*(z/d.Z()); };
    AnaSpillPD* s = new AnaSpillPD();
    AnaEventInfoPD* info = new AnaEventInfoPD(); info->IsMC = mc; info->NominalBeamMom = mom; s->EventInfo = info;
    AnaBeamPD* beam = new AnaBeamPD(); beam->BeamParticle = new AnaParticlePD(); s->Beam = beam;
    const TVector3 e = at(-50);
    Set(beam->BeamParticle->PositionEnd, e.X(), e.Y(), e.Z());
    Set(beam->BeamParticle->DirectionEnd, d.X(), d.Y(), d.Z());
    s->Bunches.push_back(new AnaBunchPD());
    AnaTrueParticlePD* bt = nullptr;
    if(mc){
      bt = new AnaTrueParticlePD(); bt->ID = 1; s->TrueParticles.push_back(bt);
      beam->BeamParticle->TrueObject = bt;
      if(withTraj)for(double z : {-100., -10., 30.}){
        AnaTrueTrajPointPD tp; const TVector3 p = at(z);
        Set(tp.Position, p.X(), p.Y(), p.Z()); Set(tp.Momentum, 6*d.X(), 6*d.Y(), 6*d.Z());
        bt->TrjPoints.push_back(tp);
      }
    }
    AnaParticlePD* prim = new AnaParticlePD();
    prim->BeamOrigin = true;
    const TVector3 ps = at(zStart);
    Set(prim->PositionStart, ps.X(), ps.Y(), ps.Z());
    Set(prim->DirectionStart, primDir.X(), primDir.Y(), primDir.Z());
    if(mc && primTrue)prim->TrueObject = bt;
    static_cast<AnaBunchPD*>(s->Bunches[0])->Particles.push_back(prim);
    return s;
  }
}

int test_beam_axis_impl(){
  using namespace beamaxis;
  BeamAxis ba;
  const TVector3 d(-0.18, -0.20, 0.96);
  // MC 6 GeV: tres eventos sobre el eje del fcl desplazados en x; el tercero con el primario
  // que no es el haz y empieza a z = 150 (fuera del corte)
  ba.Process(Spill(true, 6, -31, 423, d, 5, true, true));
  ba.Process(Spill(true, 6, -30, 423, d, 5, true, true));
  ba.Process(Spill(true, 6, -29, 423, d, 150, false, false));
  // MC 6 GeV: primario en otra dirección (falla el cos con la instrumentación)
  ba.Process(Spill(true, 6, -30, 423, d, 5, true, true, TVector3(0.5, 0, 0.866)));   // cos 0,74
  // Datos 7 GeV: un evento
  ba.Process(Spill(false, 7, -27, 425, d, 2, false, false));

  Check("grupos", ba.groups.size(), 2);
  const Group& mc = ba.groups[{6, 1}];
  Check("MC: eventos", mc.nEvents, 4);
  Check("MC: instrumentación", mc.instr.x0.size(), 4);
  Check("MC: instrumentación x0 mediana", Quantile(mc.instr.x0, .5), -30);
  Check("MC: instrumentación y0", Quantile(mc.instr.y0, .5), 423);
  Check("MC: primario = haz", mc.primBeam.x0.size(), 3);
  Check("MC: primario con el corte", mc.primCut.x0.size(), 2);
  Check("MC: primario con el corte, ángulo con el fcl", Quantile(mc.primCut.angle, .5), 0, 1e-2);
  Check("MC: distancia al eje del fcl (x0 = -30)", Quantile(mc.primCut.dist, .84), 0.885, 0.01);   // -31: 0,9 cm en x del eje, 0,9·sin(ángulo con x) = 0,885
  Check("MC: truth", mc.truth.x0.size(), 3);
  Check("MC: truth z = 0 interpolado", Quantile(mc.truth.zref, .5), 0, 1e-4);
  Check("MC: truth dirección z", Quantile(mc.truth.uz, .16), d.Unit().Z(), 1e-4);
  const Group& data = ba.groups[{7, 0}];
  Check("Datos: eventos", data.nEvents, 1);
  Check("Datos: primario con el corte", data.primCut.x0.size(), 1);
  Check("Datos: x0", Quantile(data.primCut.x0, .5), -27);
  Check("Datos: sin truth", data.truth.x0.size(), 0);
  ba.Print();
  std::cout << std::endl << (nErr ? "TEST FALLA" : "TEST OK") << std::endl;
  return nErr;
}
