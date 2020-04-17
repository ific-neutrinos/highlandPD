#ifndef pdDataClasses_hxx
#define pdDataClasses_hxx

#include "DataClasses.hxx"
#include "ParticleId.hxx"


/// AnaParticle
class AnaParticlePD: public AnaParticle{
public :

  AnaParticlePD();
  virtual ~AnaParticlePD();

  /// Clone this object.
  virtual AnaParticlePD* Clone() {
    return new AnaParticlePD(*this);
  }

  /// Dump the object to screen.
  virtual void Print() const;

protected:

  /// Copy constructor is protected, as Clone() should be used to copy this object.
  AnaParticlePD(const AnaParticlePD& part);

public:

    enum PartTypeEnum {
      kUnknown=0, 
      kShower,
      kTrack
  };

  /// Track or Shower
  PartTypeEnum Type;

  /// Is this the beam particle according to Pandora
  bool isBeamPart;

  /// 
  bool isPandoraPart;


  /// Number of hits in each wire plane
  Int_t NHitsPerPlane[3];

  /// Residual range for each wire in each plane
  Float_t ResidualRange[3][NMAXHITSPERPLANE];

  /// dEdx for each wire in each plane
  Float_t dEdx[3][NMAXHITSPERPLANE];
  Float_t dEdx_corr[3][NMAXHITSPERPLANE];

  /// dQdx for each wire in each plane
  Float_t dQdx[3][NMAXHITSPERPLANE];  
  Float_t dQdx_corr[3][NMAXHITSPERPLANE];  

  /// dQdx for each wire in each plane
  Float_t HitX[3][NMAXHITSPERPLANE]; 
  Float_t HitY[3][NMAXHITSPERPLANE]; 
  Float_t HitZ[3][NMAXHITSPERPLANE]; 

    /// Particle ID hypothesis used in the fit (if any)
  Int_t FitPDG;
  
  /// PDG of the most probable particle hypothesis used at reconstruction level
  Int_t ReconPDG[3]; 

  /// PID variables
  Float_t PID[3][10];

  Float_t PIDA[3];

  Float_t Chi2Proton;
  Float_t Chi2ndf;

  Float_t CNNscore[3];
  
  /// CALO variables
  Float_t CALO[3][10];
  
  /// Momentum by range for muon and proton hypotheses
  Float_t RangeMomentum[2];
};

/// AnaTrueParticle
class AnaTrueParticlePD: public AnaTrueParticle{
public :

  AnaTrueParticlePD();
  virtual ~AnaTrueParticlePD();

  /// Clone this object.
  virtual AnaTrueParticlePD* Clone() {
    return new AnaTrueParticlePD(*this);
  }

  /// Dump the object to screen.
  virtual void Print() const;

protected:

  /// Copy constructor is protected, as Clone() should be used to copy this object.
  AnaTrueParticlePD(const AnaTrueParticlePD& truePart);

public:

  /// Vector Pi0 decays IDs
  std::vector<Int_t> Pi0_decay_ID; 

  Bool_t Matched;

  /// Origin
  Int_t  Origin;

  /// The particle length inside the TPC
  Float_t LengthInTPC;
  
  /// The true momentum at the TPC entrance
  Float_t MomentumInTPC;  
};


/// Representation of the beam information, including POT and quality.
class AnaBeamPD: public AnaBeam{
public :

  AnaBeamPD();
  virtual ~AnaBeamPD();

  /// Clone this object.
  virtual AnaBeamPD* Clone() {
    return new AnaBeamPD(*this);
  }

  /// Dump the object to screen.
  virtual void Print() const;
  
protected:

  /// Copy constructor is protected, as Clone() should be used to copy this object.
  AnaBeamPD(const AnaBeamPD& beam);

public:

  /// Other relevant beam info
  int BeamTrigger;
  double TOF;
  int    CerenkovStatus[2];
  double CerenkovTime[2];
  double CerenkovPressure[2];
  double BeamTrackTime;
  double BeamMomentum;
  double BeamMomentumInTPC;

  int nFibers[3];
  size_t nMomenta;  
  size_t nTracks;  
  
  std::vector< int > PDGs;

};



class PDCounters{

public:
  
  PDCounters(){
    ntrue_beamdaughter_piplus=0;
    ntrue_beamdaughter_piminus=0;
    ntrue_beamdaughter_pi0=0;
    ntrue_beamdaughter_proton=0;
    ntrue_beamdaughter_neutron=0;
    ntrue_beamdaughter_nucleus=0;    
  }
  virtual ~PDCounters(){}
  
  Int_t ntrue_beamdaughter_pi0;
  Int_t ntrue_beamdaughter_piplus;
  Int_t ntrue_beamdaughter_piminus;
  Int_t ntrue_beamdaughter_proton;
  Int_t ntrue_beamdaughter_neutron;
  Int_t ntrue_beamdaughter_nucleus;
};


#endif
