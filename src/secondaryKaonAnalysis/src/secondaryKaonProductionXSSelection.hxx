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

// Selección de la medida de producción de K+ (secondaryKaonProductionXS) por el mu+ de K+ -> mu+ nu en reposo
// (23-nov-2026). Candidatos a mu: PFParticles con traza (o sin traza, con su traza forzada) con p de rango de mu en
// [100, 600] MeV/c, un extremo a >= 10 cm de las caras y uno a <= 150 cm del eje del haz (los cortes del módulo).
// Variables en muCandidateUtils.hxx. Dos ramas, una por población (las reconstruyen pasadas distintas de Pandora):
//   rama 0, en la beam slice:  caras35, vtx>15, chi2mu<6, dedxStop>3, huecos, ventana [205, 240], Michel por hits
//   rama 1, fuera:             noclear, eje<60, chi2mu<6, dedxStop>3, chi2p>80, dedxStart>1.8, huecos, ventana, Michel
// (la cadena corta de test/phase2/mu_cand_common.h con gGapVeto y gMichelHits). Cada corte quita candidatos: el evento
// pasa si queda alguno. Cada candidato guarda hasta qué corte de su rama llega (ToyBoxMuCand::Level).
namespace muCandCuts {

  enum CutId { kPopulation = 0, kFaces35, kNoClearCosmic, kAxis60, kVertex15, kChi2Mu, kdEdxStop, kChi2p, kdEdxStart, kGaps, kWindow, kMichelHits, kNCuts };

  // Umbrales (los de la revisión del 23-nov)
  const double kFacesMin      = 35.;     // cm, los dos extremos a las caras (en la beam slice)
  const double kAxisMax       = 60.;     // cm, distancia mínima al eje del haz (fuera)
  const double kVertexMin     = 15.;     // cm, los dos extremos al vértice del haz (en la beam slice; sin vértice pasa)
  const double kChi2MuMax     = 6.;      // chi2/ndf de mu con la parada en el extremo de mayor dE/dx
  const double kdEdxStopMin   = 3.;      // MeV/cm, media en los 5 cm de la parada
  const double kChi2pMin      = 80.;     // chi2/ndf de protón, misma orientación (fuera)
  const double kdEdxStartMin  = 1.8;     // MeV/cm, media en los 5 cm del inicio (fuera)
  const double kGapZ[2]       = {230., 462.};   // cm, huecos entre APAs
  const double kGapHalfWidth  = 5.;      // cm
  const double kWindowMin     = 205.;    // MeV/c, p de rango (hipótesis mu)
  const double kWindowMax     = 240.;
  const int    kMichelHitsMin = 2;       // hits con CNN de Michel > 0,5 a < 10 cm de la parada

  /// Nombre del corte (para los títulos de los pasos)
  std::string Name(int id);
  /// ¿Pasa el candidato el corte? (pop: 0 en la beam slice, 1 fuera; solo lo usa kPopulation)
  bool Pass(int id, const muCandUtils::Candidate& c, int pop);
  /// Cadena de cada población
  std::vector<int> Chain(int pop);

  // Variables derivadas de un candidato, referidas al extremo de parada
  inline double dEdxStop(const muCandUtils::Candidate& c){ return c.stopAtB ? c.dEdxB : c.dEdxA; }
  inline double dEdxStart(const muCandUtils::Candidate& c){ return c.stopAtB ? c.dEdxA : c.dEdxB; }
  inline double Chi2MuStop(const muCandUtils::Candidate& c){ return c.stopAtB ? c.chi2muF : c.chi2muR; }
  inline double Chi2MuOther(const muCandUtils::Candidate& c){ return c.stopAtB ? c.chi2muR : c.chi2muF; }
  inline double Chi2pStop(const muCandUtils::Candidate& c){ return c.stopAtB ? c.chi2pF : c.chi2pR; }
  inline double FacesStop(const muCandUtils::Candidate& c){ return c.stopAtB ? c.facesB : c.facesA; }
  inline double FacesStart(const muCandUtils::Candidate& c){ return c.stopAtB ? c.facesA : c.facesB; }
  inline double AxisMin(const muCandUtils::Candidate& c){ return std::min(c.axisA, c.axisB); }
  inline double VertexMin(const muCandUtils::Candidate& c){ return std::min(c.vtxA, c.vtxB); }
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

/// Calcula los candidatos a mu del evento (muCandUtils::Compute) y aplica la preselección geométrica del módulo
class GetMuCandidatesAction: public StepBase{
 public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new GetMuCandidatesAction();}
};

/// El evento tiene algún candidato
class MuCandidateExistsCut: public StepBase{
 public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new MuCandidateExistsCut();}
};

/// Corte sobre los candidatos de una población: el candidato sube a nivel "index" si estaba en index-1 y pasa el
/// corte; el evento pasa si algún candidato de la población llega a index
class MuCandidateCut: public StepBase{
 public:
  MuCandidateCut(int pop = 0, int index = 0, int id = 0): _pop(pop), _index(index), _id(id){}
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new MuCandidateCut(_pop, _index, _id);}
 private:
  int _pop, _index, _id;
};

#endif
