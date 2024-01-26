#ifndef QPIXSelection_h
#define QPIXSelection_h

#include "SelectionBase.hxx"
#include "Parameters.hxx"
#include "ToyBoxPD.hxx"
#include "EventBoxId.hxx"
#include "EventBoxPD.hxx"
#include "SystId.hxx"
#include "SubDetId.hxx"

class QPIXSelection: public SelectionBase{
 public:
  QPIXSelection(bool forceBreak=true);
  virtual ~QPIXSelection(){}

  //---- These are mandatory functions
  void DefineSteps();
  void DefineDetectorFV(); //dummy (not needed for this particular selection)
  ToyBoxB* MakeToyBox(){return new ToyBoxPD();}
  void InitializeEvent(AnaEventC&);

  // These ones are also mandatory, although only used in some cases. A dummy implementation is enough in many cases  
  bool FillEventSummary(AnaEventC&, Int_t*){return false;}
  SampleId_h GetSampleId(){return UNASSIGNEDID;}
  Int_t GetRelevantRecObjectGroupsForSystematic(SystId_h, Int_t* IDs, Int_t) const{ IDs[0] = EventBoxPD::kLongTracks;return 1;}
  Int_t GetRelevantTrueObjectGroupsForSystematic(SystId_h, Int_t*, Int_t) const{return 0;}

  //------------------

protected:
};


class DummyCut: public StepBase{
 public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new DummyCut();}
};

#endif
