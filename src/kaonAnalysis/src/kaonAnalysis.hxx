#ifndef kaonAnalysis_h
#define kaonAnalysis_h

#include "baseAnalysis.hxx"
#include "ToyBoxPD.hxx"
#include "standardPDTree.hxx"


/* This is an example of analysis in ProtoDUNE-SP detector 
   This example contains several elements:
     - A simple event selection of kaons decaying at rest from a beam of single kaons at 1GeV. Than means 
       that this example should be run over that particular sample
          /pnfs/dune/scratch/dunepro/v06_05_00/mergeana/gen_protoDune_kaon_1GeV_mono/13405798_0/anahist.root
     - A example of propagation of event variation systematic: dE/dx scaling
     - Filling of a root tree containing summary information about the selection and systematic mentioned above
     - A macro to produce few plots using the root file produced by the executable and the DrawingTools
 */


namespace kaonAnalysisConstants{

  const UInt_t NMAXSAVEDPARTICLES   = 20;
  const UInt_t NMAXSAVEDDAUGHTERS   = 20;
  const UInt_t NMAXSAVEDGDAUGHTERS  = 5;
  const UInt_t NMAXSAVEDGGDAUGHTERS = 5;

}

class kaonAnalysis: public baseAnalysis {
 public:
  kaonAnalysis(AnalysisAlgorithm* ana=NULL);
  virtual ~kaonAnalysis(){}

  //---- These are mandatory functions
  void DefineSelections();
  void DefineCorrections();
  void DefineSystematics();
  void DefineConfigurations();
  void DefineMicroTrees(bool addBase=true);
  void DefineTruthTree();

  void FillMicroTrees(bool addBase=true);
  void FillToyVarsInMicroTrees(bool addBase=true);

  bool CheckFillTruthTree(const AnaTrueVertex& vtx);

  using baseAnalysis::FillTruthTree;
  void FillTruthTree(const AnaTrueVertex& vtx);
  //--------------------

  bool Initialize();
  void InitializeBunch();
  void FillCategories();
  void DefineInputConverters();


  /// Get a casted AnaSpillC to AnaSpill from the InputManager
  AnaSpillPD& GetSpill(){return *static_cast<AnaSpillPD*>(&input().GetSpill());}
  
  /// Get a casted AnaBunchBB to AnaBunch from the InputManager
  AnaBunchPD& GetBunch(){return *static_cast<AnaBunchPD*>(&input().GetBunch());}

  /// Get a casted AnaEventC to AnaEvent 
  AnaEventPD& GetEvent(){return *static_cast<AnaEventPD*>(_event);}

  
  /// Create the appropriate event time from an Spill and a Bunch in that spill
  virtual AnaEventC* MakeEvent(){
    return new AnaEventPD(GetSpill(),GetBunch());
  }

  
  /// Returns the ToyBoxPD
  virtual const ToyBoxPD& box(Int_t isel=-1) const {return *static_cast<const ToyBoxPD*>(&boxB(isel));}

  /// Returns the vertex for the ToyBoxPD
  virtual AnaVertexB* GetVertex() const{return box().Vertex;}

  /// Returns the true vertex for the ToyBoxPD
  virtual AnaTrueVertexB* GetTrueVertex() const {return box().TrueVertex;}

  /// Method to add our own categories for color drawing
  void AddCustomCategories();

private:

public:

  // Needed to get the index of the counters from standardPDTree 
  Int_t seltrk_ndau;
  Int_t seltrk_ngdau;
  Int_t ntracks;  

  // Enum with unique indexes for output tree variables
  enum enumStandardMicroTrees_kaonAnalysis{
    seltrk_CNNscore =     standardPDTree::enumStandardMicroTreesLast_standardPDTree,
    seltrk_truedaukaons,
    seltrk_truedaukaon_nmu,
    seltrk_chi2_prot,
    seltrk_chi2_muon,
    seltrk_chi2_ndf,

    seltrk_dau_CNNscore,
    seltrk_dau_chi2_prot,
    seltrk_dau_chi2_muon,
    seltrk_dau_chi2_ndf,

    seltrk_gdau_CNNscore,
    seltrk_gdau_track,
    seltrk_gdau_mom_muon,
    seltrk_gdau_type,
    seltrk_gdau_chi2_prot,
    seltrk_gdau_chi2_muon,
    seltrk_gdau_chi2_ndf,

    seltrk_ggdau_CNNscore0,
    seltrk_ggdau_CNNscore1,
    seltrk_ggdau_CNNscore2,

    enumStandardMicroTreesLast_kaonAnalysis
  };

  enum enumConf_kaonAnalysis{
    detmass_syst=baseAnalysis::enumConfLast_baseAnalysis+1,    
    dedx_syst,
    Lifetime_syst,
    dQdxCalib_syst,
    Recombination_syst,
    dEdxCalib_syst,
    TrackEff_syst,
    enumConfLast_kaonAnalysis
  };

  enum enumSyst_kaonAnalysis{
    kLength=0,
    kLifetime,
    kdQdxCalib,
    kRecombination,
    kdEdxCalib,
    kBeam,
    kTrackEff,
    enumSystLast_kaonAnalysis
  };
  
};

#endif