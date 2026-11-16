#ifndef SCECorrection_h
#define SCECorrection_h

#include "CorrectionBase.hxx"
#include "BaseDataClasses.hxx"
#include "SpaceCharge.hxx"
#include "Calorimetry.hxx"

/// This class corrects hits and track end/start position by SCE

class SCECorrection: public CorrectionBase {

public:
  
  SCECorrection();

  virtual ~SCECorrection() {}

  /// Apply the sce correction
  void Apply(AnaSpillC& spill);
  
protected:

  const SpaceCharge* _sce; ///< owned by pd::services()
  const Calorimetry* _cal; ///< owned by pd::services()
  
};

#endif
