#ifndef secondaryKaonProductionXSAnalysis_h
#define secondaryKaonProductionXSAnalysis_h

#include "baseAnalysis.hxx"
#include "ToyBoxKaon.hxx"
#include "standardPDTree.hxx"
#include "pdDataClasses.hxx"
#include "kaonTree.hxx"
#include "kaonAnalysisUtils.hxx"

#include "TH1F.h"

namespace secondaryKaonProductionXSAnalysisConstants{

  const UInt_t NMAXSAVEDCANDIDATES     = 10;
}

class secondaryKaonProductionXSAnalysis: public baseAnalysis {
 public:
  secondaryKaonProductionXSAnalysis(AnalysisAlgorithm* ana=NULL);
  virtual ~secondaryKaonProductionXSAnalysis(){}

  //---- These are mandatory functions
  void DefineSelections();
  void DefineCorrections();
  void DefineSystematics();
  void DefineConfigurations();
  void DefineMicroTrees(bool addBase=true);
  void DefineTruthTree();

  void FillMicroTrees(bool addBase=true);
  void FillToyVarsInMicroTrees(bool addBase=true);

  bool CheckFillTruthTree(const AnaTrueVertex& vtx){return true;}
  bool CheckFillTruthTree(const AnaTrueParticlePD* part);

  bool FinalizeConfiguration();
  void FillTruthTree();
  void FillTruthTree(const AnaTrueVertex& vtx){return;}
  void FillTruthTree(const AnaTrueParticlePD& part);
  //--------------------

  bool Initialize();
  void FillCategories();

  /// Signal of the production measurement (measure A): the first K+ of the chain (going up through the
  /// K+ -> K+ of kaonplusInelastic) is a direct daughter of the triggering beam particle. K+ from
  /// secondary interactions are background, and the beam K+ after K+ -> K+ is not production.
  bool IsProductionSignal(const AnaTrueParticlePD* kaon);
  /// Sets the Signal flag of the true particle of a candidate (used by the reduced categories, code 3210)
  void SetIsSignal(AnaParticlePD* part);
  void DefineInputConverters();


  /// Get a casted AnaSpillC to AnaSpillPD from the InputManager
  AnaSpillPD& GetSpill(){return *static_cast<AnaSpillPD*>(&input().GetSpill());}
  
  /// Get a casted AnaBunchBB to AnaBunchPD from the InputManager
  AnaBunchPD& GetBunch(){return *static_cast<AnaBunchPD*>(&input().GetBunch());}

  /// Get a casted AnaEventC to AnaEventPD 
  AnaEventPD& GetEvent(){return *static_cast<AnaEventPD*>(_event);}
  
  /// Create the appropriate event time from an Spill and a Bunch in that spill
  virtual AnaEventC* MakeEvent(){
    return new AnaEventPD(GetSpill(),GetBunch());
  }
  
  /// Returns the ToyBoxKaon
  virtual const ToyBoxKaon& box(Int_t isel=-1) const {return *static_cast<const ToyBoxKaon*>(&boxB(isel));}

  /// Returns the vertex for the ToyBoxKaon
  virtual AnaVertexB* GetVertex() const{return box().Vertex;}

  /// Returns the true vertex for the ToyBoxKaon
  virtual AnaTrueVertexB* GetTrueVertex() const {return box().TrueVertex;}


private:

  bool _FillBeamInstrumentationInfo;
  bool _FillBeamParticleInfo;
  bool _FillBeamParticleHitsInfo;
  bool _FillBeamParticleDaughtersInfo;
  bool _FillKaonCandidatesInfo;
  bool _FillKaonCandidatesHitInfo;
  bool _FillKaonCandidatesPIDInfo;
  bool _FillBestKaonCandidateInfo;
  bool _FillToyVariablesInfo;

  bool _ApplydQdxSystematic;
  bool _ApplyRecombinationSystematic;
  bool _ApplySCESystematic;
  bool _ApplyBrokenTracksSystematic;
  bool _ApplyBeamPIDEfficiencySystematic;
  bool _ApplyBeamPartWeightSystematic;
  bool _ApplyBeamMomWeightSystematic;
  bool _ApplyProtonBackgroundWeightSystematic;

  bool _ApplySCECorrection;
  bool _ApplydEdxCorrection;
  bool _ApplyRecombinationCorrection;

public:

  enum enumConf_secondaryKaonProductionXSAnalysis{
    detmass_syst=baseAnalysis::enumConfLast_baseAnalysis+1,    
    enumConfLast_secondaryKaonProductionXSAnalysis
  };

  enum enumSyst_secondaryKaonProductionXSAnalysis{
    kdEdxVariation=0,
    kRecombination,
    kSCEGeometric,
    kBrokenTracks,
    kBeamPIDEfficiency,
    kBeamPartWeight,
    kBeamMomWeight,
    kProtonBackgroundWeight,
    enumSystLast_secondaryKaonProductionXSAnalysis
  };  
};

#endif
