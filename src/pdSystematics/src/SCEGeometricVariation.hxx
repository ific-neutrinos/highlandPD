#ifndef SCEGeometricVariation_h
#define SCEGeometricVariation_h

#include "EventVariationBase.hxx"
#include "BinnedParams.hxx"
#include "Calorimetry.hxx"

class SCEGeometricVariation: public EventVariationBase, public BinnedParams{
public:

  SCEGeometricVariation();
  virtual ~SCEGeometricVariation(); 

  
  /// Apply the systematic
  virtual void Apply(const ToyExperiment& toy, AnaEventC& event);
  
  /// Undo  the systematic variations done by ApplyVariation. This is faster tha reseting the full Spill
  bool UndoSystematic(AnaEventC& event);

  /// Scaling of the nominal SCE map for this toy (global variation)
  SCEScaling GetToyScaling(const ToyExperiment& toy);

protected:
  
  /// Is this particle relevant for this systematic ?
  bool IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const;
  
  const SpaceCharge* _sce; ///< nominal map, owned by pd::services()
};

#endif
