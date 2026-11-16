#ifndef pdDataClasses_hxx
#define pdDataClasses_hxx

#include "DataClasses.hxx"
#include "ParticleId.hxx"

const UInt_t NMAXHITSPERPLANE        = 300;
const UInt_t NMAXHITSPERPLANE_SELTRK = 2500;

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


  /// Dump the object to screen.
  virtual void Print() const;

  ///   Copy constructor
  AnaHitPD(const AnaHitPD& part);
  
public:

  // ---- Raw hit information
  AnaWireID WireID;

  Float_t Integral;
  Float_t PeakTime;
  Float_t PeakAmplitude;
  int TPCid;
  int PlaneID;

  UInt_t StartTick;
  UInt_t EndTick;
  
  UInt_t Channel;
  Int_t  View;

  /// Key of the recob::Hit in LArSoft
  Int_t UniqueID;

  /// ID of the true particle that deposits most energy in this hit (MC only)
  Int_t TrueID;

  /// Hit fit: width (ticks), goodness of fit and number of hits fitted together
  Float_t RMS;
  Float_t GoodnessOfFit;
  Int_t   Multiplicity;

  Float_t CNN[3];
  
  /// wave form associated to this hit
  std::vector<Float_t> Signal;  

  //------------------------------------------------------------

  /// Calorimetric information
  /// No SCE correction
  Float_t dQdx_NoSCE;
  Float_t dEdx_NoSCE;
  Float_t ResidualRange_NoSCE;
  Float_t Pitch_NoSCE;
  TVector3 Position_NoSCE;
  TVector3 Direction_NoSCE;


  /// SCE correction (intermediate step of the calibrations)
  Float_t dQdx_SCE;

  /// SCE correction + e⁻ lifetime correction (intermediate step of the calibrations)
  Float_t dQdx_elife;

  /// SCE correction + e⁻ lifetime correction + XYZT correction + dEdx calibration
  Float_t dQdx;
  Float_t dEdx;
  Float_t ResidualRange;
  Float_t Pitch;
  TVector3 Position;
  TVector3 Direction;
};

//trajectory point representation. Has position and direction information.
//this information is somehow redundant since each hit has a tjp associated 
//from which position and direction are obtained. This should be solved at
//some point. However, all tjp are needed to compute track lenght and sce
//systematics
class AnaTrajectoryPointPD{
public:
  AnaTrajectoryPointPD();
  virtual ~AnaTrajectoryPointPD(){}

  /// Dump the object to screen.
  virtual void Print() const;

  /// Copy constructor
  AnaTrajectoryPointPD(const AnaTrajectoryPointPD& tjp);
  
public:

  TVector3 Position;
  TVector3 Direction;

  TVector3 Position_NoSCE;
  TVector3 Direction_NoSCE;

  bool IsValid() {return (Position_NoSCE.X() != -999 && Position_NoSCE.Y() != -999 && Position_NoSCE.Z() != -999);}
};

/// 3D space point of a PFParticle (Michels, showers and PFParticles without track).
/// Value type without inheritance: there may be millions of them.
class AnaSpacePointPD{
public:
  AnaSpacePointPD();

  /// Dump the object to screen.
  void Print() const;

  Float_t Position[3];

  /// UniqueID of its hit: Pandora builds one space point per 2D hit
  Int_t HitID;
};

/// Track fit forced on a PFParticle without track (shower-like): pandora2Track and its calorimetry.
/// Value type kept apart from the AnaParticlePD members so that the analyses do not change
class AnaForcedTrackPD{
public:
  AnaForcedTrackPD();

  /// Dump the object to screen.
  void Print() const;

  Int_t   ID;
  Float_t PositionStart[3];
  Float_t PositionEnd[3];
  Float_t DirectionStart[3];
  Float_t DirectionEnd[3];
  Float_t Length;

  /// PID with the collection plane calorimetry
  Float_t Chi2Proton;
  Float_t Chi2Muon;
  Float_t Chi2ndf;

  /// Hits with calorimetry, for each plane
  std::vector<AnaHitPD> Hits[3];

  std::vector<AnaTrajectoryPointPD> TrjPoints;
};

/// True trajectory point, to follow the chain of the K+- (only for the configured PDGs)
class AnaTrueTrajPointPD{
public:
  AnaTrueTrajPointPD();

  /// Dump the object to screen.
  void Print() const;

  Float_t Position[3];
  Float_t Momentum[3];

  /// Geant4 process at this point (AnaTrueParticleB::ProcessEnum), -999 if none
  Int_t Process;
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

  /// MotherID
  Int_t ParentID;

  /// for kaon analysis
  bool IsCandidate;

  /// Track or Shower
  PartTypeEnum Type;

  /// ID of the Pandora track (pandoraTrack; the forced one is ForcedTrack[0].ID)
  Int_t TrackID;

  /// ID of the Pandora shower
  Int_t ShowerID;

  /// Is this the beam particle according to Pandora ?
  bool isPandora;

  /// Is this the beam particle or a descendant according to Pandora ?
  bool BeamOrigin;

  ///  Pandora beam particle that passes geometric cuts
  bool isBeamPart;

  double ThetaXZ;
  double ThetaYZ;

  /// Vector of hits for each plane, with calorimetry
  std::vector<AnaHitPD> Hits[3];

  /// Hits of the PFParticle without calorimetry, for each plane: those of tracks that are not
  /// in the calorimetry (mostly 2D hits without 3D point) and all those of showers and
  /// PFParticles without track
  std::vector<AnaHitPD> HitsNoCalo[3];

  /// Total Number of hits in each wire plane (The vector of Hits above might be a subsample)
  Int_t NHitsPerPlane[3];

  // Vector of trajectory points
  std::vector<AnaTrajectoryPointPD> TrjPoints;

  // PID variables
  Float_t Chi2Proton;
  Float_t Chi2Muon;
  Float_t Chi2ndf;

  Float_t CNNscore[3];

  /// Momentum by range for muon and proton hypotheses
  Float_t RangeMomentum[2];

  /// Generation of the particle (primary particle, daughter, granddaughter, etc.)
  Int_t Generation;

  /// Pandora slice: distinguishes the beam slice from the cosmic slices
  Int_t SliceID;

  /// Pandora metadata
  Bool_t  IsClearCosmic;
  Float_t PandoraTrackScore;

  /// T0 of the PFParticle (-999 if none)
  Float_t T0;

  /// Pandora beam/cosmic score and vertex
  Float_t BeamCosmicScore;
  Float_t PandoraVertex[3];

  /// Shower information (shower-like PFParticles), for each plane
  Float_t ShowerEnergy[3];
  Float_t ShowerdEdx[3];
  Float_t ShowerOpeningAngle;

  /// Track fit forced on PFParticles without track (0 or 1 elements)
  std::vector<AnaForcedTrackPD> ForcedTrack;

  /// 3D space points
  std::vector<AnaSpacePointPD> SpacePoints;

  /// Distance to the closest particle (computed by secondaryProtonSelection)
  Double_t Distance_to_closest_particle;

  /// Daughter forced by secondaryKaonXSSelection, and whether it matches the true one
  bool forced_daughter;
  bool forced_daughter_matched;
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

  /// Is this a truth signal object? Quite useful, potentially could be added in CoreDataClasses in AnaTrueObject
  Bool_t Signal;

  /// True-reco matching flag, need for 
  Bool_t Matched;

  /// Origin
  Int_t  Origin;

  /// Generation
  Int_t Generation;

  /// True trajectory (only for the configured PDGs)
  std::vector<AnaTrueTrajPointPD> TrjPoints;

  /// Hits in each plane where this particle deposits most energy, and energy (MeV) it deposits in
  /// all the hits of the event: to know whether it was reconstructed and where (MC only)
  Int_t   NHitsTrue[3];
  Float_t EDepTrue;
};

/// Data quality of ProtoDUNE-SP: what the FEMB and unstable HV filters check, saved instead of
/// filtering (data only)
class AnaDataQualityPD: public AnaDataQuality{
public :

  AnaDataQualityPD();
  virtual ~AnaDataQualityPD(){}

  /// Clone this object.
  virtual AnaDataQualityPD* Clone() {
    return new AnaDataQualityPD(*this);
  }

  /// Dump the object to screen.
  virtual void Print() const;

protected:

  /// Copy constructor is protected, as Clone() should be used to copy this object.
  AnaDataQualityPD(const AnaDataQualityPD& dq);

public:

  /// Number of active FEMBs in each APA (20 = all). -999 if the raw digits are not in the input
  Int_t NActiveFembs[6];

  /// Timestamps consistent among all the APAs and among the beam-side ones (0, 2, 4):
  /// 1 yes, 0 no, -999 if the raw timestamps are not in the input
  Int_t TimeStampConsistent;
  Int_t TimeStampConsistentBeamSide;

  /// Time from the timing system (s) and whether it is in a period of unstable HV
  Double_t TimingTimeStamp;
  Bool_t   UnstableHV;
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

  int nFibers[3];
  size_t nMomenta;  
  size_t nTracks;  
  
  std::vector< int > PDGs;
};


// Extension of AnaBunch with the hits without PFParticle
class AnaBunchPD: public AnaBunch{
public :
  
  AnaBunchPD();
  virtual ~AnaBunchPD();
  
  /// Clone this object.
  virtual AnaBunchPD* Clone() {
    return new AnaBunchPD(*this);
  }
  
  /// Dump the object to screen.
  virtual void Print() const;
  
protected:
  
  /// Copy constructor is protected, as Clone() should be used to copy this object.
  AnaBunchPD(const AnaBunchPD& bunch);
  
public:

  /// Hits without PFParticle near the ends of the beam slice
  std::vector<AnaHitPD> UnassociatedHits;
};

class AnaSpillPD: public AnaSpill{
public :

  AnaSpillPD(){}
  virtual ~AnaSpillPD(){}

  /// Clone this object.
  virtual AnaSpillPD* Clone() {
    AnaSpillPD* spill = new AnaSpillPD(*this);
    spill->RedoLinks();
    spill->isClone=true;
    return spill;
  }

protected:

public:

};

// Extension of AnaEvent with the hits without PFParticle
class AnaEventPD: public AnaEvent{
public :

  AnaEventPD();
  AnaEventPD(const AnaSpillPD& spill, const AnaBunchPD& bunch);
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

  /// Hits without PFParticle near the ends of the beam slice (copied from the bunch)
  std::vector<AnaHitPD> UnassociatedHits;
};

// Extension of AnaEventInfo
class AnaEventInfoPD: public AnaEventInfo{
public :

  AnaEventInfoPD();
  virtual ~AnaEventInfoPD();

  /// Clone this object.
  virtual AnaEventInfoPD* Clone() {
    AnaEventInfoPD* eventinfo = new AnaEventInfoPD(*this);
    return eventinfo;
  }

  /// Dump the object to screen.
  virtual void Print() const;

protected:

  /// Copy constructor is protected, as Clone() should be used to copy this object.
  AnaEventInfoPD(const AnaEventInfoPD& eventinfo);

public:

  Float_t NominalBeamMom;

  Bool_t EmptyEvent; //this should probably go in DataQuality
  Bool_t HasPandora; //this should probably go in DataQuality
};



#endif
