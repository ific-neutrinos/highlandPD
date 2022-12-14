#ifndef BrokenTrackWeight_h
#define BrokenTrackWeight_h

#include "EventWeightBase.hxx"
#include "BinnedParams.hxx"


class BrokenTrackWeight: public EventWeightBase, public BinnedParams {
public:
  
  BrokenTrackWeight();
  virtual ~BrokenTrackWeight(){}
  
  using EventWeightBase::ComputeWeight;
  Weight_h ComputeWeight(const ToyExperiment& toy, const AnaEventC& event, const ToyBoxB& box);


  /// Check whether a AnaRecObject is relevant for this systematic or not
  virtual bool IsRelevantRecObject(const AnaEventC&, const AnaRecObjectC&) const;
  
};

#endif
