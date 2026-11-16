//Class with calorimetry calibration

#ifndef Calorimetry_h
#define Calorimetry_h

#include "TFile.h"
#include "TF1.h"
#include "TH2F.h"
#include "TSpline.h"
#include "TVector3.h"

#include "pdDataClasses.hxx"
#include "SpaceCharge.hxx"

/// What can change per event or per toy when applying calorimetry:
/// electron lifetime, ModBox recombination parameters and the space
/// charge map. Start from Calorimetry::NominalContext() and modify what
/// the correction or systematic needs.
struct CaloContext{
  double Lifetime;          ///< electron lifetime (us)
  double ModBoxA;           ///< ModBox alpha
  double ModBoxB;           ///< ModBox beta
  const SpaceCharge* SCE;   ///< space charge map (not owned)
  SCEScaling SCEScale;      ///< scaling of the space charge map (nominal by default)
};

/// Calorimetry calibration. It holds only fixed configuration (constants
/// from pdUtils.parameters.dat, calibration histograms and the nominal
/// space charge map), so a single instance is shared through
/// pd::services().Calo(). Methods without a CaloContext argument use
/// NominalContext().
class Calorimetry{
public:

  explicit Calorimetry(const SpaceCharge& sce);
  virtual ~Calorimetry();

  /// MC lifetime, nominal ModBox parameters and nominal SCE
  CaloContext NominalContext() const;

  /// Electron lifetime for the run of this event (MC lifetime for MC)
  double GetLifetime(const AnaEventPD& event) const;

  /// Nominal values
  double GetLifetime() const {return _Lifetime;}
  double GetModBoxA()  const {return _ModBoxA;}
  double GetModBoxB()  const {return _ModBoxB;}
  const SpaceCharge& GetSCE() const {return _sce;}

  void ApplySCECorrection(AnaParticlePD* part) const;
  void ApplySCECorrection(AnaHitPD& hit) const;
  void ApplySCECorrection(AnaParticlePD* part, const CaloContext& ctx) const;
  void ApplySCECorrection(AnaHitPD& hit, const CaloContext& ctx) const;

  void ApplyLifetimeCorrection(AnaParticlePD* part) const;
  void ApplyLifetimeCorrection(AnaHitPD& hit) const;
  void ApplyLifetimeCorrection(AnaParticlePD* part, const CaloContext& ctx) const;
  void ApplyLifetimeCorrection(AnaHitPD& hit, const CaloContext& ctx) const;
  void UndoLifetimeCorrection(AnaParticlePD* part) const;
  void UndoLifetimeCorrection(AnaHitPD& hit) const;
  void UndoLifetimeCorrection(AnaParticlePD* part, const CaloContext& ctx) const;
  void UndoLifetimeCorrection(AnaHitPD& hit, const CaloContext& ctx) const;

  void ApplyNormalization(AnaHitPD &hit) const;
  void ApplyXCalibration(AnaHitPD &hit) const;
  void ApplyYZCalibration(AnaHitPD &hit) const;

  void ApplyRecombination(AnaParticlePD* part) const;
  void ApplyRecombination(AnaHitPD &hit) const;
  void ApplyRecombination(AnaParticlePD* part, const CaloContext& ctx) const;
  void ApplyRecombination(AnaHitPD &hit, const CaloContext& ctx) const;

  void CalibratedQdx(AnaParticlePD* part) const; //think about naming
  void CalibratedQdx(AnaHitPD &hit) const; //think about naming

  double GetNormalization(AnaHitPD &hit) const;
  double GetXCalibration(AnaHitPD &hit) const;
  double GetYZCalibration(AnaHitPD &hit) const;

protected:

  Calorimetry(const Calorimetry&) = delete;
  Calorimetry& operator=(const Calorimetry&) = delete;

  void CreateNormHistogram();
  void CreateXCalHistogram();
  void CreateYZCalHistogram();

  double GetElectronsFromADCArea(double dQdx, int PlaneID) const {return dQdx/_CalAreaConstants[PlaneID];}

  TH1F* _h_norm;
  TH1F* _h_XCal[3];
  TH2F* _h_YZCal[3][2];

  double _Efield;
  double _ModBoxA;
  double _ModBoxB;

  double _Lifetime;
  double _Vdrift;
  double _APAx;

  double _CalAreaConstants[3];

  const SpaceCharge& _sce;

};

#endif
