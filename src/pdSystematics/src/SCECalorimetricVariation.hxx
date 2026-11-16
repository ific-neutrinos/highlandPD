#ifndef SCECalorimetricVariation_h
#define SCECalorimetricVariation_h

#include "SCEVariationBase.hxx"
#include "Calorimetry.hxx"

/// SCE systematic on hit positions and on the pitch correction of dQ/dx
/// (dQdx_SCE), used by the dE/dx calibration analyses. It does not change
/// the track geometry (see SCEGeometricVariation).
/// See SCEVariationBase for the variation modes.
class SCECalorimetricVariation: public SCEVariationBase{
public:

  SCECalorimetricVariation();
  virtual ~SCECalorimetricVariation(); 

  /// Apply the systematic
  virtual void Apply(const ToyExperiment& toy, AnaEventC& event);
  
  /// Undo  the systematic variations done by ApplyVariation. This is faster tha reseting the full Spill
  bool UndoSystematic(AnaEventC& event);

  
protected:
  
  /// Is this particle relevant for this systematic ?
  bool IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const;
  
  const Calorimetry* _cal; ///< owned by pd::services()
};

#endif
