#ifndef secondaryKaonPairXSSelection_h
#define secondaryKaonPairXSSelection_h

#include "secondaryKaonAnalysis.hxx"
#include "secondaryKaonSelection.hxx"
#include "SelectionBase.hxx"
#include "Parameters.hxx"
#include "ToyBoxKaon.hxx"
#include "EventBoxKaon.hxx"
#include "EventBoxId.hxx"
#include "SystId.hxx"
#include "SubDetId.hxx"

class secondaryKaonPairXSSelection: public SelectionBase{
 public:
  secondaryKaonPairXSSelection(bool forceBreak=true);
  virtual ~secondaryKaonPairXSSelection(){}

  //---- These are mandatory functions
  void DefineSteps();
  void DefineDetectorFV(){SetDetectorFV(SubDetId::kSubdet1_1);} //dummy (not needed for this particular selection)
  ToyBoxB* MakeToyBox(){return new ToyBoxKaon();}
  void InitializeEvent(AnaEventC&);

  // These ones are also mandatory, although only used in some cases. A dummy implementation is enough if many cases  
  bool FillEventSummary(AnaEventC&, Int_t*){return false;}
  SampleId_h GetSampleId(){return UNASSIGNEDID;}
  Int_t GetRelevantRecObjectGroupsForSystematic(SystId_h, Int_t* IDs, Int_t) const{ IDs[0] = EventBoxKaon::kCandidatesAndDaughters;return 1;}
  Int_t GetRelevantTrueObjectGroupsForSystematic(SystId_h, Int_t* IDs, Int_t) const{IDs[0] = EventBoxKaon::kCandidatesAndDaughters;return 1;}

  //------------------

protected:

private:


};



class GetKaonsForPairXSAction: public StepBase{
 public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new GetKaonsForPairXSAction();}   // devolvía GetKaonsAction (la de la tesis): corregido el 25-nov
};

// class KaonChi2Cut: public StepBase{
//  public:
//   StepBase* MakeClone(){return new KaonChi2Cut();}
// };

/// Corte de especie del haz (BeamPDGCut de pdBaseSelection) configurable: pdg = 0 deja pasar todo. El paso está siempre,
/// para que los niveles de accum_level no cambien con el parámetro.
class PairXSBeamPDGCut: public StepBase{
 public:
  PairXSBeamPDGCut(int pdg = 0): _pdg(pdg){}
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new PairXSBeamPDGCut(_pdg);}
 private:
  int _pdg;
};

/// Busca la forced daughter entre lo que el módulo guardaría en un evento normal (muCandUtils::NormalSaveMode), como
/// secondaryKaonProductionXS: en los eventos con volcado completo (MC con un K+ en la cadena del haz) no se usan los
/// PFParticles que un evento normal no tendría (25-nov). La selección solo usa tipo, longitud e inicio, que las copias
/// en resumen conservan, así que se usa la partícula original.
class GetForcedDaughterPairXSAction: public StepBase{
 public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new GetForcedDaughterPairXSAction();}
};




#endif
