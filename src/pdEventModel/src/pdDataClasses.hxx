#ifndef pdDataClasses_hxx
#define pdDataClasses_hxx

// Event model of ProtoDUNE-SP: extensions of the highland classes (DataClasses.hxx) and
// ProtoDUNE-specific classes. The minitrees written by CreateHighlandMiniTree (highlandAnalysis)
// and read by the analyses are made of these objects.
//
// Conventions:
//  - Lengths and positions in cm, energies in MeV and momenta in GeV, unless stated otherwise.
//  - -999 means "not filled" (kIntUnassigned, kFloatUnassigned, kDoubleUnassigned in the .cxx).
//  - "SCE" means corrected for the space charge effect; "NoSCE", as reconstructed.
//  - Every new member must go to the copy constructor (or it is lost in Clone(), which is used
//    to write and read the minitrees) and to Print().

#include "DataClasses.hxx"
#include "ParticleId.hxx"

/// Maximum number of hits per plane saved in the output trees of the analyses (kaonTree) and
/// read from the pduneana trees (PDSPAnalyzerTreeConverter)
const UInt_t NMAXHITSPERPLANE        = 300;

/// Maximum number of hits per plane of the selected track saved in standardPDTree
const UInt_t NMAXHITSPERPLANE_SELTRK = 2500;

/// Wire of a hit (geo::WireID in LArSoft)
class AnaWireID{
public:

  AnaWireID(){Wire=Plane=TPC=Cryostat=0;}
  AnaWireID(const AnaWireID& wireID);
  virtual ~AnaWireID(){}

  Int_t Wire;     ///< Wire number within the plane
  Int_t Plane;    ///< Plane: 0 (U, induction), 1 (V, induction), 2 (collection)
  Int_t TPC;      ///< TPC number
  Int_t Cryostat; ///< Cryostat number (always 0 in ProtoDUNE-SP)
};

/// Hit on a wire (recob::Hit), with its calorimetry when it belongs to a track.
/// Value type: it is stored in vectors inside the particles.
class AnaHitPD{
public:
  AnaHitPD();
  virtual ~AnaHitPD(){}

  /// Dump the object to screen.
  virtual void Print() const;

  /// Copy constructor
  AnaHitPD(const AnaHitPD& part);

public:

  // ---- Raw hit information (from recob::Hit) ----

  AnaWireID WireID;        ///< Wire of the hit
  Float_t Integral;        ///< Integral of the fitted pulse (ADC x tick): the charge of the hit
  Float_t PeakTime;        ///< Time of the peak (ticks)
  Float_t PeakAmplitude;   ///< Amplitude of the peak (ADC)
  int TPCid;               ///< TPC of the hit (the same as WireID.TPC)
  int PlaneID;             ///< Plane of the hit (the same as WireID.Plane)
  UInt_t StartTick;        ///< First tick of the hit region
  UInt_t EndTick;          ///< Last tick of the hit region
  UInt_t Channel;          ///< Readout channel
  Int_t  View;             ///< Wire view (geo::View_t)
  Int_t UniqueID;          ///< Key of the recob::Hit in its LArSoft collection: identifies the hit in the event
  Float_t RMS;             ///< Width of the fitted pulse (ticks)
  Float_t GoodnessOfFit;   ///< Goodness of the pulse fit (chi2/ndf)
  Int_t   Multiplicity;    ///< Number of pulses fitted together in the same region

  /// ID of the true particle that deposits most energy in this hit (MC only, from the backtracker).
  /// Negative for EM shower particles not kept by Geant4
  Int_t TrueID;

  /// CNN scores of the hit (emtrkmichelid): 0 track-like, 1 EM-like, 2 Michel-like
  Float_t CNN[3];

  /// Waveform of the hit between StartTick and EndTick (ADC). Only with SaveHitSignal
  std::vector<Float_t> Signal;

  // ---- Calorimetry: only for the hits of the tracks (AnaParticlePD::Hits) ----
  // dQdx in ADC/cm, dEdx in MeV/cm, Pitch and ResidualRange in cm

  /// Without SCE corrections (pandoracalonosce)
  Float_t dQdx_NoSCE;
  Float_t dEdx_NoSCE;
  Float_t ResidualRange_NoSCE;   ///< Distance to the end of the track along it
  Float_t Pitch_NoSCE;           ///< Track length seen by the wire
  TVector3 Position_NoSCE;       ///< 3D position of the hit on the track
  TVector3 Direction_NoSCE;      ///< Track direction at the hit (from its trajectory point)

  /// Intermediate steps of the calibrations (dEdxCalibrationAnalysis, Calorimetry): dQdx with the
  /// SCE corrections, and with the SCE and electron lifetime corrections
  Float_t dQdx_SCE;
  Float_t dQdx_elife;

  /// Fully corrected and calibrated (pandoracali: SCE, electron lifetime, XYZT corrections and
  /// dE/dx calibration). These are the ones the analyses use; the corrections and systematics
  /// of highlandPD modify them
  Float_t dQdx;
  Float_t dEdx;
  Float_t ResidualRange;
  Float_t Pitch;
  TVector3 Position;    ///< 3D position of the hit, SCE corrected
  TVector3 Direction;   ///< Track direction at the hit, SCE corrected
};

/// Trajectory point of a track (recob::Track), with position and direction.
/// It is partly redundant with the hits (each hit has a trajectory point), but all the
/// trajectory points are needed to compute the track length and the SCE systematics.
class AnaTrajectoryPointPD{
public:
  AnaTrajectoryPointPD();
  virtual ~AnaTrajectoryPointPD(){}

  /// Dump the object to screen.
  virtual void Print() const;

  /// Copy constructor
  AnaTrajectoryPointPD(const AnaTrajectoryPointPD& tjp);

public:

  TVector3 Position;         ///< Position, SCE corrected (filled by ParticlePositionSCECorrection)
  TVector3 Direction;        ///< Direction, SCE corrected (filled by ParticlePositionSCECorrection)
  TVector3 Position_NoSCE;   ///< Position as reconstructed
  TVector3 Direction_NoSCE;  ///< Direction as reconstructed

  /// False for the points without a valid position (Pandora adds them for the hits without 3D point)
  bool IsValid() {return (Position_NoSCE.X() != -999 && Position_NoSCE.Y() != -999 && Position_NoSCE.Z() != -999);}
};

/// 3D space point of a PFParticle (recob::SpacePoint): the 3D cloud of the PFParticles without
/// track (Michels, showers), and of the tracks if SaveTrackSpacePoints.
/// Value type without inheritance: there may be millions of them.
class AnaSpacePointPD{
public:
  AnaSpacePointPD();

  /// Dump the object to screen.
  void Print() const;

  Float_t Position[3];  ///< Position as reconstructed (no SCE correction)

  /// UniqueID of its hit: Pandora builds one space point per 2D hit
  Int_t HitID;
};

/// Track fit forced on a PFParticle without track (shower-like), from pandora2Track and its
/// calorimetry (pandora2cali). Useful when a muon or a short track is reconstructed as a shower.
/// Kept apart from the members of AnaParticlePD so that the analyses do not change.
class AnaForcedTrackPD{
public:
  AnaForcedTrackPD();

  /// Dump the object to screen.
  void Print() const;

  Int_t   ID;                 ///< ID of the track in pandora2Track
  Float_t PositionStart[3];   ///< As reconstructed (no SCE correction)
  Float_t PositionEnd[3];
  Float_t DirectionStart[3];
  Float_t DirectionEnd[3];
  Float_t Length;

  /// PID with the collection plane calorimetry: chi2 with the proton and muon dE/dx templates,
  /// and number of degrees of freedom
  Float_t Chi2Proton;
  Float_t Chi2Muon;
  Float_t Chi2ndf;

  /// Hits with calorimetry, for each plane (only the collection plane unless SaveHitsAllPlanes)
  std::vector<AnaHitPD> Hits[3];

  /// Trajectory points of the forced track
  std::vector<AnaTrajectoryPointPD> TrjPoints;
};

/// True trajectory point (simb::MCTrajectory), to follow the chain of the K+- (scatterings,
/// inelastic interactions, decay). Only for the beam particle and the PDGs in TrueTrajPDGs.
class AnaTrueTrajPointPD{
public:
  AnaTrueTrajPointPD();

  /// Dump the object to screen.
  void Print() const;

  Float_t Position[3];  ///< True position
  Float_t Momentum[3];  ///< True momentum (GeV)

  /// Geant4 process at this point (AnaTrueParticleB::ProcessEnum), -999 if none.
  /// Only filled if the MC production saved the trajectory processes
  Int_t Process;
};

/// Reconstructed particle: a Pandora PFParticle, with its track or shower.
/// The base classes (AnaParticle, AnaParticleB, AnaParticleMomB, AnaParticleE) provide, among
/// others, UniqueID (the PFParticle index), PositionStart/End, DirectionStart/End, Length,
/// Momentum, NHits, DaughtersIDs, TrueObject (matched true particle), TrueEff and TruePur.
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

  // ---- Identity and hierarchy ----

  /// UniqueID of the parent PFParticle (-1 for the primary ones)
  Int_t ParentID;

  /// Generation in the Pandora hierarchy: 0 the beam particle, 1 its daughters, etc.
  /// -1 for the primary cosmics (their daughters have 0, 1, ...)
  Int_t Generation;

  /// Track or shower (if the PFParticle has both, the track)
  PartTypeEnum Type;

  /// ID of the Pandora track (pandoraTrack; the forced one is ForcedTrack[0].ID)
  Int_t TrackID;

  /// ID of the Pandora shower (pandoraShower)
  Int_t ShowerID;

  /// Is this the beam particle according to Pandora (primary PFParticle of the beam slice)?
  bool isPandora;

  /// Is this the beam particle or a descendant according to Pandora (beam slice)? False for cosmics
  bool BeamOrigin;

  /// Pandora beam particle that passes the beam quality cuts (only from the pduneana trees;
  /// used by BrokenTracksAnalysis)
  bool isBeamPart;

  // ---- Pandora information ----

  /// Pandora slice: distinguishes the beam slice from the cosmic slices
  Int_t SliceID;

  /// Pandora tags it as a clear cosmic
  Bool_t IsClearCosmic;

  /// Pandora track/shower score: close to 1 track-like, close to 0 shower-like
  Float_t PandoraTrackScore;

  /// T0 of the PFParticle (ns), for the PFParticles crossing the cathode or an anode. -999 if none
  Float_t T0;

  /// Pandora score of the slice being beam or cosmic
  Float_t BeamCosmicScore;

  /// Pandora vertex of the PFParticle (no SCE correction)
  Float_t PandoraVertex[3];

  // ---- Hits, trajectory and space points ----

  /// Hits with calorimetry (those of the tracks), for each plane. The analyses use these.
  /// Only the collection plane unless SaveHitsAllPlanes
  std::vector<AnaHitPD> Hits[3];

  /// Hits of the PFParticle without calorimetry, for each plane: those of tracks that are not
  /// in the calorimetry (mostly 2D hits without 3D point) and all those of showers and
  /// PFParticles without track
  std::vector<AnaHitPD> HitsNoCalo[3];

  /// Number of hits with calorimetry in each plane
  Int_t NHitsPerPlane[3];

  /// Trajectory points of the track
  std::vector<AnaTrajectoryPointPD> TrjPoints;

  /// 3D space points (PFParticles without track, and tracks if SaveTrackSpacePoints)
  std::vector<AnaSpacePointPD> SpacePoints;

  /// Track fit forced on PFParticles without track (0 or 1 elements)
  std::vector<AnaForcedTrackPD> ForcedTrack;

  // ---- PID ----

  /// Chi2 of the collection plane dE/dx with the proton and muon templates, and number of
  /// degrees of freedom (number of hits used). Recomputed by the corrections and systematics
  /// that change dE/dx (pdAnaUtils::RecomputeChi2PID)
  Float_t Chi2Proton;
  Float_t Chi2Muon;
  Float_t Chi2ndf;

  /// Average CNN score of the collection plane hits: 0 track-like, 1 EM-like, 2 Michel-like
  Float_t CNNscore[3];

  /// Momentum by range: 0 proton hypothesis, 1 muon hypothesis
  Float_t RangeMomentum[2];

  // ---- Shower information (shower-like PFParticles) ----

  Float_t ShowerEnergy[3];     ///< Energy in each plane (units of the shower reconstruction)
  Float_t ShowerdEdx[3];       ///< dE/dx at the start of the shower in each plane (MeV/cm)
  Float_t ShowerOpeningAngle;  ///< Opening angle (rad)

  // ---- Filled by the analyses ----

  /// Candidate of the kaon analysis (CreateMiniTreeCandidates, CreateMiniTreeEventCandidates)
  bool IsCandidate;

  /// Track angles used by the dE/dx calibration (only from the Michel removing trees)
  double ThetaXZ;
  double ThetaYZ;

  /// Distance to the closest particle (computed by secondaryProtonSelection)
  Double_t Distance_to_closest_particle;

  /// Daughter forced by secondaryKaonXSSelection, and whether it matches the true one
  bool forced_daughter;
  bool forced_daughter_matched;
};

/// True particle (simb::MCParticle).
/// The base classes (AnaTrueParticle, AnaTrueParticleB) provide, among others, ID (Geant4
/// TrackId), ParentID, PDG, ParentPDG, GParentPDG, Position/PositionEnd, Direction/DirectionEnd,
/// Momentum/MomentumEnd, ProcessStart/ProcessEnd and Daughters.
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

  /// Signal of the analysis (set by the analysis, e.g. secondaryKaonAnalysis::SetIsSignal)
  Bool_t Signal;

  /// For the true beam particle: is it the true particle matched to the Pandora beam particle?
  Bool_t Matched;

  /// Origin of the particle (simb::Origin_t): 0 unknown, 2 cosmic ray (CORSIKA), 4 single
  /// particle (the beam). Daughters with unknown origin take the one of their parent
  Int_t  Origin;

  /// Generation: 0 the beam particle (or a primary cosmic), 1 its daughters, etc.
  /// -999 for the particles added only because they are matched to a reconstructed one
  Int_t Generation;

  /// True trajectory (only for the beam particle and the PDGs in TrueTrajPDGs)
  std::vector<AnaTrueTrajPointPD> TrjPoints;

  /// Hits in each plane where this particle deposits most energy, and energy it deposits in the
  /// collection plane hits of the whole event (MC only, SaveTrueHitCounts). 0 and 0 if it is not
  /// seen in the TPC: to know whether it was reconstructed and where
  Int_t   NHitsTrue[3];
  Float_t EDepTrue;
};

/// Data quality of ProtoDUNE-SP (data only). Besides GoodDaq (beamline trigger, from the base
/// class), what the FEMB and unstable HV filters check, saved instead of filtering.
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

  /// Time from the timing system (s)
  Double_t TimingTimeStamp;

  /// Is the event in a period of unstable HV (same list as ProtoDUNEUnstableHVFilter)?
  Bool_t   UnstableHV;
};

/// Beam line instrumentation of ProtoDUNE-SP (beam::ProtoDUNEBeamEvent).
/// The base class provides BeamParticle (the beam track reconstructed by the instrumentation:
/// position, direction and momentum before the TPC; in MC, TrueObject is the true beam
/// particle), POTSincePreviousSavedSpill and SpillsSincePreviousSavedSpill (spills not saved
/// before this one, for the normalization; 1 and 0 in each event written by the module).
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

  int    BeamTrigger;          ///< Timing trigger word
  double TOF;                  ///< Time of flight between the TOF counters (ns)
  int    CerenkovStatus[2];    ///< Trigger status of the two Cherenkov counters
  double CerenkovTime[2];      ///< Time stamp of the two Cherenkov counters
  double CerenkovPressure[2];  ///< Pressure of the two Cherenkov counters
  double BeamTrackTime;        ///< Time of the beam trigger (s; ns only resolved to ~0.2 us)

  /// Momentum measured by the spectrometer (GeV), the same as BeamParticle->Momentum
  /// (BeamMomSmearingCorrection keeps them equal)
  double BeamMomentum;

  int    nFibers[3];           ///< Active fibers in the three profile monitors of the spectrometer
  size_t nMomenta;             ///< Number of reconstructed beam momenta
  size_t nTracks;              ///< Number of reconstructed beam tracks

  /// PDGs compatible with the beam line PID (TOF and Cherenkov) for the nominal momentum
  std::vector< int > PDGs;
};

/// Bunch of ProtoDUNE-SP (there is only one per event). The base class holds the reconstructed
/// particles: the beam slice (BeamOrigin) and the saved cosmics.
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

  /// Hits without PFParticle near the ends of the beam slice (2D distance in each plane,
  /// UnassociatedHitsRadius). No calorimetry
  std::vector<AnaHitPD> UnassociatedHits;
};

/// Spill of ProtoDUNE-SP: one per event. The base class holds EventInfo, DataQuality, Beam,
/// Bunches and the true particles. Clone() does not clone the truth (see CLAUDE.md).
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
};

/// Event seen by the analyses, built from the spill and its bunch
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

/// Event information of ProtoDUNE-SP. The base classes provide Run, SubRun, Event, IsMC and
/// EventTime (s, from the art event).
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

  /// Nominal beam momentum of the run (GeV)
  Float_t NominalBeamMom;

  /// No beam activity at the TPC entrance, as protoana::ProtoDUNEEmptyEventFinder (the
  /// reco_reconstructable_beam_event of PDSPAnalyzer, negated): fewer than NHitsThreshold space
  /// points of SpacePointSolver (reco3d, independent of Pandora) in a box at the beam entrance.
  /// False when there are no reco3d space points (NBeamRegionSpacePoints = -1)
  Bool_t EmptyEvent;

  /// Number of reco3d space points in the box of the empty event finder (-1 without reco3d,
  /// -999 not filled: minitrees from before 2026-11-19 or from the PDSPAnalyzer converter)
  Int_t NBeamRegionSpacePoints;

  /// The Pandora products are missing or there are no PFParticles (CreateHighlandMiniTree)
  Bool_t MissingReco;

  /// Pandora found a beam slice (it should probably go in DataQuality)
  Bool_t HasPandora;

  /// CreateHighlandMiniTree saved all the PFParticles of the event complete, also in events
  /// without beam slice (FullDumpPrescale, or FullDumpTruePDGs in MC). Selection biased by the
  /// truth when it comes from FullDumpTruePDGs
  Bool_t FullDump;
};

#endif
