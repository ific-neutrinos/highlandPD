#ifndef pdDataClasses_hxx
#define pdDataClasses_hxx

#include "DataClasses.hxx"
#include "ParticleId.hxx"

class AnaWireID{
public:

  AnaWireID(){Wire=Plane=TPC=Cryostat=0;}
  AnaWireID(const AnaWireID& wireID);
  virtual ~AnaWireID(){}
  
  Int_t Wire;
  Int_t Plane;
  Int_t TPC;
  Int_t Cryostat;
};

class AnaHitPD{
public:
  AnaHitPD();
  virtual ~AnaHitPD(){}
  
  AnaHitPD(Int_t wire, Float_t integ, Float_t peakT, Float_t peakAmp, TVector3 pos){
    WireID.Plane =wire;
    Integral = integ;
    PeakTime = peakT;
    PeakAmplitude = peakAmp;
    Position = pos;
  }
    
  AnaHitPD(double integral){
    Integral = integral;
  }
  
  /// Dump the object to screen.
  virtual void Print() const;
  /*    
  Float_t PeakAmplitude() const {return fPeakAmplitude;}
  Float_t PeakTime() const {return fPeakTime;}
  AnaWireID  WireID()const {return fWireID;}
  Float_t Integral()const {return fIntegral;}

  UInt_t Channel()const {return fChannel;}
  Int_t View()const {return fView;}

  std::vector<Float_t> Signal() const {return fSignal;}
  */
  //protected:

  /// Copy constructor is protected, as Clone() should be used to copy this object.
  AnaHitPD(const AnaHitPD& part);
  
public:

  AnaWireID WireID;    
  Float_t Integral;
  Float_t PeakTime;
  Float_t PeakAmplitude;
  TVector3 Position;

  UInt_t StartTick;
  UInt_t EndTick;
  
  UInt_t Channel;
  Int_t  View;

  /// Residual range for each wire in each plane
  Float_t ResidualRange;

  /// dEdx for each wire in each plane
  Float_t dEdx;
  Float_t dEdx_corr;

  /// dQdx for each wire in each plane
  Float_t dQdx;  
  Float_t dQdx_corr;  

  /// wave form associated to this hit
  std::vector<Float_t> Signal;  
};


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

  /// Is this the beam particle according to Pandora ?
  bool isPandora;

  ///  Pandora beam particle that passes geometric cuts
  bool isBeamPart;

  /// Vector of hits for each plane
  std::vector<AnaHitPD> Hits[3];

  /// Total Number of hits in each wire plane (The vector of Hits above might be a subsample)
  Int_t NHitsPerPlane[3];
  
  Float_t truncLibo_dEdx;
  
  // PID variables
  Float_t Chi2Proton;
  Float_t Chi2ndf;

  Float_t CNNscore[3];
  
  /// CALO variables
  Float_t CALO[3][10];
  
  /// Momentum by range for muon and proton hypotheses
  Float_t RangeMomentum[2];
  Float_t RangeMomentum_alt[2];

  /// Alternate length
  Float_t Length_alt;

  // ---- OBSOLETE PID VARIABLES ----------
  
  /// Particle ID hypothesis used in the fit (if any)
  Int_t FitPDG;
  
  /// PDG of the most probable particle hypothesis used at reconstruction level
  Int_t ReconPDG[3]; 

  /// PID variables
  Float_t PID[3][10];

  /// PIDA
  Float_t PIDA[3];  
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

  /// True-reco matching flag, need for 
  Bool_t Matched;

  /// Origin
  Int_t  Origin;

  /// The particle length inside the TPC
  Float_t LengthInTPC;
  
  /// The true momentum at the TPC entrance
  Float_t MomentumInTPC;  
};

/// Extension of AnaEvent to include specific information of the ProtoDUNE beam line instrumentation
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


/// ProtoDUNE counters

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

// Extension of AnaEvent to include the APA wire wafeforms, needed to recompute the CNN 
class AnaSpillPD: public AnaSpill{
public :

  AnaSpillPD();
  virtual ~AnaSpillPD();

  /// Clone this object.
  virtual AnaSpillPD* Clone() {
    AnaSpillPD* spill = new AnaSpillPD(*this);
    spill->RedoLinks();
    spill->isClone=true;
    return spill;
  }

  /// Dump the object to screen.
  virtual void Print() const;

protected:

  /// Copy constructor is protected, as Clone() should be used to copy this object.
  AnaSpillPD(const AnaSpillPD& spill);

public:

  std::vector<std::vector<Float_t > > ADC;
  
};


// Extension of AnaEvent to include the APA wire wafeforms, needed to recompute the CNN 
class AnaEventPD: public AnaEvent{
public :

  AnaEventPD();
  AnaEventPD(const AnaSpillPD& spill, const AnaBunch& bunch);
  virtual ~AnaEventPD();

  /// Clone this object.
  virtual AnaEventPD* Clone() {
    AnaEventPD* spill = new AnaEventPD(*this);
    spill->isClone=true;
    return spill;
  }

  /// Dump the object to screen.
  virtual void Print() const;

protected:

  /// Copy constructor is protected, as Clone() should be used to copy this object.
  AnaEventPD(const AnaEventPD& event);

public:

  std::vector<std::vector<Float_t > > ADC;
};



#endif
