#ifndef SCEVariation_h
#define SCEVariation_h

#include "EventVariationBase.hxx"
#include "BinnedParams.hxx"
#include "Calorimetry.hxx"

class SCEVariation: public EventVariationBase, public BinnedParams{
public:

  SCEVariation();
  virtual ~SCEVariation(); 

  /// Apply the systematic
  virtual void Apply(const ToyExperiment& toy, AnaEventC& event);
  
  /// Undo  the systematic variations done by ApplyVariation. This is faster tha reseting the full Spill
  bool UndoSystematic(AnaEventC& event);

  /// Scaling of the nominal SCE map for this toy (global variation).
  /// Use it with the nominal map, pd::services().SCE().
  SCEScaling GetToyScaling(const ToyExperiment& toy);
  
protected:
  
  /// Is this particle relevant for this systematic ?
  bool IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const;
  
  const Calorimetry* _cal; ///< owned by pd::services()
  const SpaceCharge* _sce; ///< nominal map, owned by pd::services()
};

#endif
