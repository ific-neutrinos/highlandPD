//*********************************************************************
//minimal implementetation of the space charge effects in ProtoDUNE-SP
//only considering spline_th3 methods 
//and backwards (correction) implementation
//Get_Offsets give the effect forward of the SCE effect (it does not work) 
//GetCaL_Offsets give the effect backward (it works)
//
//M. García, migarpez@ific.uv.es
//*********************************************************************

#ifndef SpaceCharge_h
#define SpaceCharge_h

#include "TFile.h"
#include "TF1.h"
#include "TH3F.h"
#include "TSpline.h"
#include "TVector3.h"

#include "pdDataClasses.hxx"

/// Relative uncertainty of each node of the calibration maps, for the
/// voxel variation of the space charge map. Index 0 is the x<0 drift
/// volume and 1 the x>0 one; each vector is indexed by the TH3 global bin
/// of the maps of that volume (the six maps of a volume share binning).
/// Index = -1 means the node is not varied.
struct SCEVoxelSigmas{
  std::vector<Int_t>   Index[2];  ///< systematic parameter of each node
  std::vector<Float_t> Sigma[2];  ///< relative uncertainty of each node
};

/// How the nominal space charge map is scaled, e.g. for a systematic toy.
/// The map interpolation is linear in the map values, so scaling the map
/// is done on the fly, without copying it:
/// - Global: a single factor for the whole map (1 = nominal).
/// - Voxel: if Voxels and ToyVariations are set, each node of the map is
///   scaled by 1 + Sigma*ToyVariations[Index] (Global is then ignored).
struct SCEScaling{
  double Global = 1.;
  const SCEVoxelSigmas* Voxels = nullptr;
  const Float_t* ToyVariations = nullptr;

  bool IsVoxel() const {return Voxels && ToyVariations;}

  /// Scaling factor of a node (TH3 global bin) of drift volume dv (0: x<0, 1: x>0)
  double NodeFactor(int dv, int bin) const{
    Int_t i = Voxels->Index[dv][bin];
    if(i < 0) return 1.;
    return 1 + Voxels->Sigma[dv][bin]*ToyVariations[i];
  }
};

class SpaceCharge{
public:
  
  SpaceCharge();
  virtual ~SpaceCharge();
    
  void Initialize();

  TVector3 GetPosOffsets(TVector3 const& point) const;
  TVector3 GetEfieldOffsets(TVector3 const& point) const;
  /// The optional SCEScaling argument scales the map (default: nominal)
  TVector3 GetCalPosOffsets(TVector3 const& point, int const& TPCid, const SCEScaling& scale = SCEScaling()) const;
  TVector3 GetCalEfieldOffsets(TVector3 const& point, int const& TPCid, const SCEScaling& scale = SCEScaling()) const;

  void ApplyParticlePositionCorrection(AnaParticlePD* part, const SCEScaling& scale = SCEScaling()) const;

  void ApplyPositionCorrection(AnaParticlePD* part, const SCEScaling& scale = SCEScaling()) const;
  void ApplyPositionCorrection(AnaHitPD& hit, const SCEScaling& scale = SCEScaling()) const;
  void ApplyPitchCorrection(AnaParticlePD* part, const SCEScaling& scale = SCEScaling()) const;
  void ApplyPitchCorrection(AnaHitPD& hit, const SCEScaling& scale = SCEScaling()) const;

  void ApplyTrjPointPositionCorrection(AnaParticlePD* part, const SCEScaling& scale = SCEScaling()) const;
  void ApplyTrjPointPositionCorrection(AnaTrajectoryPointPD& tp, const SCEScaling& scale = SCEScaling()) const;
  void ApplyTrjPointDirectionCorrection(AnaParticlePD* part) const;

  UInt_t GetNbinsX() const;
  UInt_t GetNbinsY() const;
  UInt_t GetNbinsZ() const;

  double GetBinCenterX(UInt_t bin) const;
  double GetBinCenterY(UInt_t bin) const;
  double GetBinCenterZ(UInt_t bin) const;

  double GetEFieldXAtPoint(double x, double y, double z, const SCEScaling& scale = SCEScaling()) const;
  double GetEFieldYAtPoint(double x, double y, double z, const SCEScaling& scale = SCEScaling()) const;
  double GetEFieldZAtPoint(double x, double y, double z, const SCEScaling& scale = SCEScaling()) const;

  /// Calibration (backward) maps: 0-5 x>0 volume (Dx,Dy,Dz,Ex,Ey,Ez), 6-11 x<0 volume
  const TH3F* GetCalHistogram(int i) const {return CalSCEhistograms.at(i);}
  
protected:
  
  TVector3 GetOffsets(TVector3 const& point, TH3F* hX, TH3F* hY, TH3F* hZ, int maptype, int driftvol) const;

  void InitializeHistograms();
  void InitializeHistogram(TH3F** h, const char* histo, const char* newname);
  void InitializeSplines();

  void ClearSplines();
  void ClearVectorOfSplines(std::vector<std::vector<TSpline3*>> &splines);

  std::vector<TH3F*> SCEhistograms = std::vector<TH3F*>(12); //Histograms are Dx, Dy, Dz, dEx/E0, dEy/E0, dEz/E0 (positive; repeat for negative)
  std::vector<TH3F*> CalSCEhistograms = std::vector<TH3F*>(12); 
  
  bool IsInsideBoundaries(TVector3 const& point) const;
  bool IsTooFarFromBoundaries(TVector3 const& point) const;
  TVector3 PretendAtBoundary(TVector3 const& point) const;

  TVector3 ElectronDiverterPosOffsets(TVector3 const& point) const;
  bool   fEnableElectronDiverterDistortions1;
  bool   fEnableElectronDiverterDistortions2;
  double fEDZCenterMin;
  double fEDZCenterMax;
  double fEDAXPosOffs;
  double fEDBZPosOffs;
  double fEDs;
  double fEDChargeLossZLowMin;
  double fEDChargeLossZLowMax;
  double fEDChargeLossZHighMin;
  double fEDChargeLossZHighMax;

  TSpline3* MakeSpline(TH3F* spline_hist, int dim1, int dim2_bin, int dim3_bin, int maptype, int driftvol) const;
  double InterpolateSplines(TH3F* interp_hist, double xVal, double yVal, double zVal, int dim, int maptype, int driftvol) const;

  /// Voxel variation: same interpolation as InterpolateSplines, with each
  /// node of the map scaled by its factor. ihist is the index in
  /// CalSCEhistograms. Uses the cardinal basis of the spline (_basis).
  double InterpolateVoxel(int ihist, double xVal, double yVal, double zVal, int dim, const SCEScaling& scale) const;
  TVector3 GetOffsetsVoxel(TVector3 const& point, int ihist0, const SCEScaling& scale) const;

  /// Cardinal basis of the natural cubic spline along the spline axis of each
  /// calibration map: _basis[ihist][i] is the spline through the unit vector e_i.
  /// The map spline through values y_i is sum_i y_i*_basis[ihist][i].
  void InitializeBasis();
  std::vector<TSpline3*> _basis[12];
  /// Same basis as plain arrays, for speed: knot positions, and coefficients
  /// (y,b,c,d) of basis i in interval k at _basisCoef[ih][4*(i*n+k)]
  std::vector<double> _basisKnots[12];
  std::vector<double> _basisCoef[12];
  /// Map values by TH3 global bin, for speed
  std::vector<double> _mapContent[12];
  
  std::string fInputFilename;
  TFile* fInputFile;
    
  TH3F* hDx_cal_pos;
  TH3F* hDy_cal_pos;
  TH3F* hDz_cal_pos;
  TH3F* hEx_cal_pos;
  TH3F* hEy_cal_pos;
  TH3F* hEz_cal_pos;
  		   
  TH3F* hDx_cal_neg;
  TH3F* hDy_cal_neg;
  TH3F* hDz_cal_neg;
  TH3F* hEx_cal_neg;
  TH3F* hEy_cal_neg;
  TH3F* hEz_cal_neg;

  std::vector<std::vector<TSpline3*>> spline_dx_fwd_neg;
  std::vector<std::vector<TSpline3*>> spline_dy_fwd_neg;
  std::vector<std::vector<TSpline3*>> spline_dz_fwd_neg;
  
  std::vector<std::vector<TSpline3*>> spline_dx_bkwd_neg;
  std::vector<std::vector<TSpline3*>> spline_dy_bkwd_neg;
  std::vector<std::vector<TSpline3*>> spline_dz_bkwd_neg;
  
  std::vector<std::vector<TSpline3*>> spline_dEx_neg;
  std::vector<std::vector<TSpline3*>> spline_dEy_neg;
  std::vector<std::vector<TSpline3*>> spline_dEz_neg;
  
  std::vector<std::vector<TSpline3*>> spline_dx_fwd_pos;
  std::vector<std::vector<TSpline3*>> spline_dy_fwd_pos;
  std::vector<std::vector<TSpline3*>> spline_dz_fwd_pos;
  
  std::vector<std::vector<TSpline3*>> spline_dx_bkwd_pos;
  std::vector<std::vector<TSpline3*>> spline_dy_bkwd_pos;
  std::vector<std::vector<TSpline3*>> spline_dz_bkwd_pos;
  
  std::vector<std::vector<TSpline3*>> spline_dEx_pos;
  std::vector<std::vector<TSpline3*>> spline_dEy_pos;
  std::vector<std::vector<TSpline3*>> spline_dEz_pos;
};
  
#endif
