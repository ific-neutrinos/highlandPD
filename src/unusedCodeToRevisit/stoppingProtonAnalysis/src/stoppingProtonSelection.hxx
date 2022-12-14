#ifndef stoppingProtonSelection_h
#define stoppingProtonSelection_h

#include "SelectionBase.hxx"
#include "Parameters.hxx"
#include "ToyBoxPD.hxx"
#include "EventBoxId.hxx"
#include "EventBoxPD.hxx"
#include "SystId.hxx"
#include "SubDetId.hxx"

class stoppingProtonSelection: public SelectionBase{
 public:
  stoppingProtonSelection(bool forceBreak=true);
  virtual ~stoppingProtonSelection(){}

  //---- These are mandatory functions
  void DefineSteps();
  void DefineDetectorFV(); //dummy (not needed for this particular selection)
  ToyBoxB* MakeToyBox(){return new ToyBoxPD();}
  void InitializeEvent(AnaEventC&);

  // These ones are also mandatory, although only used in some cases. A dummy implementation is enough if many cases  
  bool FillEventSummary(AnaEventC&, Int_t*){return false;}
  SampleId_h GetSampleId(){return UNASSIGNEDID;}
  Int_t GetRelevantRecObjectGroupsForSystematic(SystId_h, Int_t* IDs, Int_t) const{ IDs[0] = EventBoxPD::kLongTracks;return 1;}
  Int_t GetRelevantTrueObjectGroupsForSystematic(SystId_h, Int_t*, Int_t) const{return 0;}

  //------------------

  bool IsRelevantTrueObjectForSystematicInToy(const AnaEventC&, const ToyBoxB&, AnaTrueObjectC*, SystId_h systId, Int_t branch=0) const;
  
protected:

  Int_t _ProtonRangeCutIndex;
  Int_t _ProtonRangeStepIndex;
  Int_t _FindMainTrackStepIndex;
  Int_t _TotalMultiplicityCutIndex;
};


class BeamProtonAngleCut: public StepBase{
 public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new BeamProtonAngleCut();}
};

class ProtonCSDARangeCut: public StepBase{
 public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new ProtonCSDARangeCut();}
};


class ProtonPIDACut: public StepBase{
 public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new ProtonPIDACut();}
};


#endif
