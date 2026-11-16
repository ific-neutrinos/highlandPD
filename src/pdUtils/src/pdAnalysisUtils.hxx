#ifndef pdAnalysisUtils_h
#define pdAnalysisUtils_h

#include "ParticleId.hxx"
#include "AnalysisUtils.hxx"
#include "ToyBoxPD.hxx"
#include "EventBoxId.hxx"
#include "pdDataClasses.hxx"
#include "TProfile.h"
#include <limits>


namespace pdAnaUtils{

  /// Particle masses (MeV), PDG 2022. Use these instead of local values.
  constexpr double kMuonMass   = 105.6583755;
  constexpr double kPionMass   = 139.57039;
  constexpr double kKaonMass   = 493.677;
  constexpr double kProtonMass = 938.27208816;

  /// Name ("muon", "pion", "kaon", "proton") and mass (MeV) of a particle
  /// hypothesis (PDG 13, 211, 321, 2212). Returns false for any other PDG;
  /// callers must then report it and skip the calculation.
  bool GetParticleHypothesis(Int_t pdg, std::string& name, double& mass);


  /// Computes the range momentum
  Float_t ComputeRangeMomentum(double trkrange, int pdg);

  /// Computes the CSDARange
  Float_t ComputeCSDARange(double beammom, int pdg);

  /// Find the beam true particle 
  AnaTrueParticle* GetBeamTrueParticle(const AnaSpillB& spill);  
  
  /// The AnaTrueParticle with a given ID (anaUtils::GetTrueParticleByID with the PD type)
  AnaTrueParticlePD* GetTrueParticle(AnaEventB* event, Int_t ID);
  AnaTrueParticlePD* GetTrueParticle(const std::vector<AnaTrueParticleB*>& trueParticles, Int_t ID);

  AnaParticlePD* GetRecoParticleWithAssociatedTrueID(const std::vector<AnaParticleB*> particles, Int_t true_ID);

  // retreieve the BI particle
  AnaParticlePD* GetBeamParticle(const AnaEventC& event);
  AnaParticlePD* GetBeamParticleInTPC(const AnaEventPD* event);

  // retreieve the true BI particle
  AnaTrueParticlePD* GetTrueBeamParticle(const AnaEventC& event);

  /// PID chi2 and ndf for a pdg hypothesis (13, 321, 2212), using the hits up to a
  /// maximum residual range (all of them by default)
  std::pair< double, int > Chi2PID(const AnaParticlePD& part, const int pdg,
                                   const double maxRR = std::numeric_limits<double>::infinity());

  // Compute the track length using the trajectory points
  Float_t ComputeTrackLengthFromTrajectoryPoints(AnaParticlePD* part); 

  void ComputeParticlePositionAndDirection(AnaParticlePD* part); 

  // Compute the truncated mean of an std vector
  Float_t ComputeTruncatedMean(float truncate_low, float truncate_high, const std::vector<double>& dEdx);
  Float_t ComputeTruncatedMean(float truncate_low, float truncate_high, const std::vector<AnaHitPD>& hits);

  // Distance mother daughter
  Float_t ComputeDistanceMotherDaughter(AnaParticlePD* mother, AnaParticlePD* daughter);
  
  // Cos mother daughter
  Float_t ComputeCosMotherDaughter(AnaParticlePD* mother, AnaParticlePD* daughter);

  // Average dEdx
  Float_t ComputeAveragedEdxOverResRange(AnaParticlePD* part, double maxresrange = 9999);

  bool IsStoppingInFV(AnaParticlePD *part);

  int GetHitTPCid(AnaHitPD& hit);
  int GetPosTPCid(TVector3 pos);

  void EstimateHitsDirection(AnaParticlePD* part);

  void ComputeResidualRange(AnaParticlePD* part);

  Double_t ComputeDepositedEnergy(AnaParticlePD* part);
 
  Double_t EstimateTrueMomAtAPABorder(AnaParticlePD* part);

  Double_t ComputeDistanceToClosestParticle(AnaParticlePD* part, AnaParticleB** parts, const int nparts);

  bool GetBeamQuality(AnaEventPD* event);
  void GetBeamQualityCuts(AnaEventPD* event, 
			  double &mean_x, double &mean_y, double &mean_z,
			  double &sigma_x, double &sigma_y, double &sigma_z,
			  double &cos);

  Float_t GetdEdxLikelihood(AnaParticlePD* part, Int_t PDG);
  Float_t GetdEdxLikelihood_UpToRR(AnaParticlePD* part, Int_t PDG, const double maxRR);
  Float_t dEdxLikelihood(TGraph* tg, const TGraph* tg_ke, 
			 Float_t mass);
  std::pair<Float_t,Float_t> GetdEdxLikelihoodFreeRange(AnaParticlePD* part, Int_t PDG);
  std::pair<Float_t,Float_t> GetdEdxLikelihoodFreeRange_UpToRR(AnaParticlePD* part, Int_t PDG, const double maxRR);
  std::pair<Float_t,Float_t> dEdxLikelihoodFreeRange(TGraph* tg, const TGraph* tg_ke, 
				  Float_t mass);
  double GetDensityCorrection(double beta, double gamma);
  double GetdEdxBetheBloch(double KE, double mass);
  double GetWmax(double KE, double mass);
  double GetLandauXi(double KE, double dx, double mass);
  double dEdxPDF(double *x, double *par);
}

#endif
