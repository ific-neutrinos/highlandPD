#ifndef secondaryKaonProductionXSSelection_h
#define secondaryKaonProductionXSSelection_h

#include "SelectionBase.hxx"
#include "Parameters.hxx"
#include "ToyBoxMuCand.hxx"
#include "EventBoxKaon.hxx"
#include "EventBoxId.hxx"
#include "SystId.hxx"
#include "SubDetId.hxx"
#include <string>
#include <vector>

/*
  Selección de la medida de producción de K+ (secondaryKaonProductionXS) por el mu+ de K+ -> mu+ nu en reposo
  (23-nov-2026).

  QUÉ SE SELECCIONA
  Candidatos a mu: PFParticles con traza (o sin traza, con su traza forzada) con p de rango de mu en [100, 600] MeV/c,
  un extremo a >= 10 cm de las caras y uno a <= 150 cm del eje del haz (los cortes del módulo que guarda los
  minitrees). Sus variables se calculan en muCandidateUtils.hxx (las mismas que en las macros de test/phase2).

  CÓMO ESTÁ ORGANIZADA (diferencia con una selección habitual de highland)
  En una selección habitual cada corte mira UN objeto (la traza del haz, o el candidato i de la rama i) y dice si el
  evento pasa. Aquí hay muchos candidatos por evento (decenas), así que cada corte se aplica a TODOS los candidatos a
  la vez:
    - cada candidato lleva un contador, su "nivel" (ToyBoxMuCand::Level): cuántos cortes de su cadena ha pasado
      seguidos;
    - cada corte solo mira los candidatos que han llegado hasta él (nivel = cortes anteriores) y sube el nivel de los
      que lo pasan;
    - el evento pasa el corte si al menos un candidato lo pasa.
  Así el accum_level del evento es el del mejor candidato, y en el árbol se guarda el nivel de cada candidato
  (cand_level) para poder dibujar cualquier variable en cualquier etapa de la cadena.
  Esa contabilidad está en una sola clase base (MuCandidateCutBase). Cada corte es una clase pequeña derivada, como en
  highland, que solo dice si un candidato pasa (PassCandidate), con el umbral en el constructor.

  RAMAS: una por población, porque las reconstruyen pasadas distintas de Pandora y tienen fondos distintos.
    rama 0, en la beam slice: población, caras > 35, vértice > 15, chi2 mu < 6, dE/dx parada > 3, huecos,
                              Michel por hits, ventana [205, 240]
    rama 1, fuera:            población, no IsClearCosmic, eje < 80, chi2 mu < 6, dE/dx parada > 3, chi2 p > 80,
                              dE/dx inicio > 1,8, huecos, Michel por hits, ventana [205, 240]
  La ventana va la última: en la medida B es el pico que se ajusta (la muestra del ajuste es la del penúltimo nivel,
  7 dentro y 9 fuera; el rango del ajuste lo pone la macro del ajuste), y en la de conteo, el último corte.
  (la cadena corta de test/phase2/mu_cand_common.h con el veto de los huecos y el Michel por hits; los umbrales,
  abajo en muCandCuts).
*/

namespace muCandCuts {

  // ---- Umbrales por defecto de la cadena (revisión del 23-nov) ----
  const double kFacesMin      = 35.;     // cm, los dos extremos a las caras del volumen activo (en la beam slice)
  const double kAxisMax       = 80.;     // cm (60 hasta el 26-nov-2026: el barrido con v1_00-v1_05 da el óptimo en 80-100), distancia mínima de un extremo al eje del haz (fuera)
  const double kVertexMin     = 15.;     // cm, los dos extremos al vértice del haz (en la beam slice; sin vértice pasa)
  const double kChi2MuMax     = 6.;      // chi2/ndf de mu con la parada en el extremo de mayor dE/dx
  const double kdEdxStopMin   = 3.;      // MeV/cm, media en los 5 cm de la parada (pico de Bragg)
  const double kChi2pMin      = 80.;     // chi2/ndf de protón, misma orientación (fuera)
  const double kdEdxStartMin  = 1.8;     // MeV/cm, media en los 5 cm del inicio (fuera)
  const double kGapZ[2]       = {230., 462.};   // cm, huecos entre APAs (electron diverters)
  const double kGapHalfWidth  = 5.;      // cm
  const double kWindowMin     = 205.;    // MeV/c, p de rango (hipótesis mu): pico del mu de K en reposo
  const double kWindowMax     = 240.;
  const int    kMichelHitsMin = 2;       // hits con CNN de Michel > 0,5 a < 10 cm de la parada

  // Niveles de la cadena de cada rama (0 dentro, 1 fuera; ver DefineSteps): la ventana de p de rango es el último
  // corte; la muestra del ajuste de la medida B es la de todos los cortes menos ella. Se guardan en el microtree
  // (sel_fitlevel, sel_countlevel) para que las macros no dependan del orden. Si se cambia la cadena, cambiarlos.
  const int kFitLevel[2]   = {7, 9};
  const int kCountLevel[2] = {8, 10};

  // ---- Variables derivadas de un candidato, referidas al extremo de parada (el de mayor dE/dx: stopAtB) ----
  inline double dEdxStop(const muCandUtils::Candidate& c){ return c.stopAtB ? c.dEdxB : c.dEdxA; }
  inline double dEdxStart(const muCandUtils::Candidate& c){ return c.stopAtB ? c.dEdxA : c.dEdxB; }
  inline double Chi2MuStop(const muCandUtils::Candidate& c){ return c.stopAtB ? c.chi2muF : c.chi2muR; }
  inline double Chi2MuOther(const muCandUtils::Candidate& c){ return c.stopAtB ? c.chi2muR : c.chi2muF; }
  inline double Chi2pStop(const muCandUtils::Candidate& c){ return c.stopAtB ? c.chi2pF : c.chi2pR; }
  inline double FacesStop(const muCandUtils::Candidate& c){ return c.stopAtB ? c.facesB : c.facesA; }
  inline double FacesStart(const muCandUtils::Candidate& c){ return c.stopAtB ? c.facesA : c.facesB; }
  inline double AxisMin(const muCandUtils::Candidate& c){ return std::min(c.axisA, c.axisB); }
  inline double VertexMin(const muCandUtils::Candidate& c){ return std::min(c.vtxA, c.vtxB); }
  /// Distancia de z al hueco entre APAs más cercano (1e9 si no hay posición)
  inline double GapDist(double z){ return z < -900 ? 1e9 : std::min(std::abs(z - kGapZ[0]), std::abs(z - kGapZ[1])); }
}

class secondaryKaonProductionXSSelection: public SelectionBase{
 public:
  secondaryKaonProductionXSSelection(bool forceBreak=true);
  virtual ~secondaryKaonProductionXSSelection(){}

  //---- These are mandatory functions
  void DefineSteps();
  void DefineDetectorFV(){SetDetectorFV(SubDetId::kSubdet1_1);} //dummy (not needed for this particular selection)
  ToyBoxB* MakeToyBox(){return new ToyBoxMuCand();}
  void InitializeEvent(AnaEventC&);

  // These ones are also mandatory, although only used in some cases. A dummy implementation is enough if many cases
  bool FillEventSummary(AnaEventC&, Int_t*){return false;}
  SampleId_h GetSampleId(){return UNASSIGNEDID;}
  Int_t GetRelevantRecObjectGroupsForSystematic(SystId_h, Int_t* IDs, Int_t) const{ IDs[0] = EventBoxKaon::kCandidatesAndDaughters;return 1;}
  Int_t GetRelevantTrueObjectGroupsForSystematic(SystId_h, Int_t* IDs, Int_t) const{IDs[0] = EventBoxKaon::kCandidatesAndDaughters;return 1;}
};

//------------------------------------------------------------------------------------------------------------------
// Acciones y cortes de evento
//------------------------------------------------------------------------------------------------------------------

/// Acción: calcula los candidatos a mu del evento (muCandUtils::Compute) y aplica la preselección geométrica del
/// módulo. Los guarda en la caja (ToyBoxMuCand) con nivel 0.
class GetMuCandidatesAction: public StepBase{
 public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new GetMuCandidatesAction();}
};

/// Corte de evento: hay algún candidato
class MuCandidateExistsCut: public StepBase{
 public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new MuCandidateExistsCut();}
};

//------------------------------------------------------------------------------------------------------------------
// Cortes sobre los candidatos (uno por clase, como en highland)
//------------------------------------------------------------------------------------------------------------------

/// Base de los cortes sobre candidatos. Hace la contabilidad de los niveles (ver la explicación arriba); las clases
/// derivadas solo implementan PassCandidate. pop: rama/población (0 en la beam slice, 1 fuera); index: posición del
/// corte en la cadena de su rama (1 el primero). Los pone DefineSteps.
class MuCandidateCutBase: public StepBase{
 public:
  MuCandidateCutBase(int pop = 0, int index = 0): _pop(pop), _index(index){}
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  /// ¿Pasa este candidato el corte?
  virtual bool PassCandidate(const muCandUtils::Candidate& c) const = 0;
  void SetPopulationAndIndex(int pop, int index){ _pop = pop; _index = index; }
 protected:
  int _pop, _index;
};

/// El candidato es de la población de esta rama (primer corte de cada rama)
class MuCandPopulationCut: public MuCandidateCutBase{
 public:
  using MuCandidateCutBase::MuCandidateCutBase;
  bool PassCandidate(const muCandUtils::Candidate& c) const { return c.inSlice == (_pop == 0); }
  StepBase* MakeClone(){ return new MuCandPopulationCut(_pop, _index); }
};

/// Los dos extremos a más de "min" cm de las caras del volumen activo
class MuCandFacesCut: public MuCandidateCutBase{
 public:
  MuCandFacesCut(double min = muCandCuts::kFacesMin, int pop = 0, int index = 0): MuCandidateCutBase(pop, index), _min(min){}
  bool PassCandidate(const muCandUtils::Candidate& c) const { return std::min(muCandCuts::FacesStop(c), muCandCuts::FacesStart(c)) > _min; }
  StepBase* MakeClone(){ return new MuCandFacesCut(_min, _pop, _index); }
 private:
  double _min;
};

/// Pandora no lo marca como cósmico claro (IsClearCosmic)
class MuCandNoClearCosmicCut: public MuCandidateCutBase{
 public:
  using MuCandidateCutBase::MuCandidateCutBase;
  bool PassCandidate(const muCandUtils::Candidate& c) const { return !c.clearCosmic; }
  StepBase* MakeClone(){ return new MuCandNoClearCosmicCut(_pop, _index); }
};

/// Algún extremo a menos de "max" cm del eje del haz
class MuCandAxisCut: public MuCandidateCutBase{
 public:
  MuCandAxisCut(double max = muCandCuts::kAxisMax, int pop = 0, int index = 0): MuCandidateCutBase(pop, index), _max(max){}
  bool PassCandidate(const muCandUtils::Candidate& c) const { return muCandCuts::AxisMin(c) < _max; }
  StepBase* MakeClone(){ return new MuCandAxisCut(_max, _pop, _index); }
 private:
  double _max;
};

/// Los dos extremos a más de "min" cm del vértice del haz (sin vértice, pasa). Quita trozos de piones y protones que
/// salen del vértice; ojo, también las K cortas de la primera interacción (medida A)
class MuCandVertexCut: public MuCandidateCutBase{
 public:
  MuCandVertexCut(double min = muCandCuts::kVertexMin, int pop = 0, int index = 0): MuCandidateCutBase(pop, index), _min(min){}
  bool PassCandidate(const muCandUtils::Candidate& c) const { return muCandCuts::VertexMin(c) < 0 || muCandCuts::VertexMin(c) > _min; }
  StepBase* MakeClone(){ return new MuCandVertexCut(_min, _pop, _index); }
 private:
  double _min;
};

/// chi2/ndf de la hipótesis mu (parada en el extremo de mayor dE/dx) menor que "max"
class MuCandChi2MuCut: public MuCandidateCutBase{
 public:
  MuCandChi2MuCut(double max = muCandCuts::kChi2MuMax, int pop = 0, int index = 0): MuCandidateCutBase(pop, index), _max(max){}
  bool PassCandidate(const muCandUtils::Candidate& c) const { return muCandCuts::Chi2MuStop(c) >= 0 && muCandCuts::Chi2MuStop(c) < _max; }
  StepBase* MakeClone(){ return new MuCandChi2MuCut(_max, _pop, _index); }
 private:
  double _max;
};

/// dE/dx medio en los 5 cm de la parada mayor que "min" (hay pico de Bragg)
class MuCandStopdEdxCut: public MuCandidateCutBase{
 public:
  MuCandStopdEdxCut(double min = muCandCuts::kdEdxStopMin, int pop = 0, int index = 0): MuCandidateCutBase(pop, index), _min(min){}
  bool PassCandidate(const muCandUtils::Candidate& c) const { return muCandCuts::dEdxStop(c) > _min; }
  StepBase* MakeClone(){ return new MuCandStopdEdxCut(_min, _pop, _index); }
 private:
  double _min;
};

/// chi2/ndf de la hipótesis protón (misma orientación) mayor que "min" (no es un protón)
class MuCandChi2pCut: public MuCandidateCutBase{
 public:
  MuCandChi2pCut(double min = muCandCuts::kChi2pMin, int pop = 0, int index = 0): MuCandidateCutBase(pop, index), _min(min){}
  bool PassCandidate(const muCandUtils::Candidate& c) const { return muCandCuts::Chi2pStop(c) > _min; }
  StepBase* MakeClone(){ return new MuCandChi2pCut(_min, _pop, _index); }
 private:
  double _min;
};

/// dE/dx medio en los 5 cm del inicio mayor que "min" (quita trozos EM)
class MuCandStartdEdxCut: public MuCandidateCutBase{
 public:
  MuCandStartdEdxCut(double min = muCandCuts::kdEdxStartMin, int pop = 0, int index = 0): MuCandidateCutBase(pop, index), _min(min){}
  bool PassCandidate(const muCandUtils::Candidate& c) const { return muCandCuts::dEdxStart(c) > _min; }
  StepBase* MakeClone(){ return new MuCandStartdEdxCut(_min, _pop, _index); }
 private:
  double _min;
};

/// Ningún extremo a menos de "halfWidth" cm de un hueco entre APAs (trazas cortadas en los electron diverters)
class MuCandGapsCut: public MuCandidateCutBase{
 public:
  MuCandGapsCut(double halfWidth = muCandCuts::kGapHalfWidth, int pop = 0, int index = 0): MuCandidateCutBase(pop, index), _halfWidth(halfWidth){}
  bool PassCandidate(const muCandUtils::Candidate& c) const { return muCandCuts::GapDist(c.startZ) >= _halfWidth && muCandCuts::GapDist(c.stopZ) >= _halfWidth; }
  StepBase* MakeClone(){ return new MuCandGapsCut(_halfWidth, _pop, _index); }
 private:
  double _halfWidth;
};

/// p de rango (hipótesis mu) en [min, max] MeV/c: la ventana del pico del mu de K en reposo
class MuCandRangeMomCut: public MuCandidateCutBase{
 public:
  MuCandRangeMomCut(double min = muCandCuts::kWindowMin, double max = muCandCuts::kWindowMax, int pop = 0, int index = 0): MuCandidateCutBase(pop, index), _min(min), _max(max){}
  bool PassCandidate(const muCandUtils::Candidate& c) const { return c.pRange >= _min && c.pRange <= _max; }
  StepBase* MakeClone(){ return new MuCandRangeMomCut(_min, _max, _pop, _index); }
 private:
  double _min, _max;
};

/// Al menos "min" hits con CNN de Michel > 0,5 a menos de 10 cm de la parada (de cualquier PFParticle o no asociados).
/// Ficheros sin esa variable (-999): se usa la CNN de Michel media de la propia partícula > 0,02 (corte antiguo)
class MuCandMichelHitsCut: public MuCandidateCutBase{
 public:
  MuCandMichelHitsCut(int min = muCandCuts::kMichelHitsMin, int pop = 0, int index = 0): MuCandidateCutBase(pop, index), _min(min){}
  // Sin hit de calorimetría en la parada (arNMichel = -1: showers, partículas en resumen) no se pueden contar los hits
  // de alrededor: se usa la CNN de Michel media de la partícula (hasta el 24-nov solo con -999, que no se daba)
  bool PassCandidate(const muCandUtils::Candidate& c) const { return c.arNMichel[1] < 0 ? c.cnnMichel > 0.02 : c.arNMichel[1] >= _min; }
  StepBase* MakeClone(){ return new MuCandMichelHitsCut(_min, _pop, _index); }
 private:
  int _min;
};

#endif
