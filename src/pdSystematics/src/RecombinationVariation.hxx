#ifndef RecombinationVariation_h
#define RecombinationVariation_h

#include "EventVariationBase.hxx"
#include "BinnedParams.hxx"
#include "Calorimetry.hxx"

class RecombinationVariation: public EventVariationBase, public BinnedParams{
public:

  RecombinationVariation();
  virtual ~RecombinationVariation(){} 
  
  /// Apply the systematic
  virtual void Apply(const ToyExperiment& toy, AnaEventC& event);
  
  /// Undo  the systematic variations done by ApplyVariation. This is faster tha reseting the full Spill
  bool UndoSystematic(AnaEventC& event);
  
protected:
  
  /// Is this particle relevant for this systematic ?
  bool IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const;
  
  BinnedParams* _BP_A;
  BinnedParams* _BP_B;  
  
  Calorimetry* _cal;
};

#endif
