#ifndef SCEGeometricVariation_h
#define SCEGeometricVariation_h

#include "SCEVariationBase.hxx"
#include "Calorimetry.hxx"

/// SCE systematic on the track geometry: trajectory points, start/end
/// positions, length, hit positions and pitch, and the quantities derived
/// from them (residual range, chi2, range momentum). It does not change
/// dQ/dx or dE/dx: the calorimetric part is covered by dEdxVariation.
/// See SCEVariationBase for the variation modes.
class SCEGeometricVariation: public SCEVariationBase{
public:

  SCEGeometricVariation();
  virtual ~SCEGeometricVariation(); 

  
  /// Apply the systematic
  virtual void Apply(const ToyExperiment& toy, AnaEventC& event);
  
  /// Undo  the systematic variations done by ApplyVariation. This is faster tha reseting the full Spill
  bool UndoSystematic(AnaEventC& event);


protected:
  
  /// Is this particle relevant for this systematic ?
  bool IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const;
  
};

#endif
