#ifndef pdBaseSelection_h
#define pdBaseSelection_h

#include "SelectionBase.hxx"
#include "Parameters.hxx"
#include "ToyBoxPD.hxx"
#include "EventBoxPD.hxx"
#include "EventBoxId.hxx"
#include "SystId.hxx"
#include "SubDetId.hxx"

/// This file defines steps that are common to all analyses.


class FindBeamTrackAction: public StepBase{
 public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new FindBeamTrackAction();}
};

class BeamTrackExistsCut: public StepBase{
public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new BeamTrackExistsCut();}
};

class BeamPDGCut: public StepBase{
public:
  using StepBase::StepBase;
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new BeamPDGCut();}
};

class BeamQualityCut: public StepBase{
 public:
  using StepBase::Apply;
  bool Apply(AnaEventC& event, ToyBoxB& box) const;
  StepBase* MakeClone(){return new BeamQualityCut();}
};

#endif
