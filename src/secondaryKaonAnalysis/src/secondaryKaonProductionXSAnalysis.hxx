#ifndef secondaryKaonProductionXSAnalysis_h
#define secondaryKaonProductionXSAnalysis_h

#include <array>
#include <map>
#include "baseAnalysis.hxx"
#include "ToyBoxMuCand.hxx"
#include "standardPDTree.hxx"
#include "pdDataClasses.hxx"

// Análisis de producción de K+ (medida A: sigma de producción de K+ en la primera interacción del haz) con el mu+ de
// K+ -> mu+ nu en reposo (23-nov-2026). Selección en secondaryKaonProductionXSSelection.hxx; variables de los
// candidatos en muCandidateUtils.hxx (las mismas que las macros de test/phase2).
//
// Árbol "ana": una entrada por evento con algún candidato; los candidatos, en vectores cand_* (contador ncand, hasta
// NMAXMUCANDIDATES). cand_pop (0 en la beam slice, 1 fuera) y cand_level (cortes de la cadena de su población que pasa
// seguidos; el primero es el de la población) permiten dibujar cualquier variable en cualquier etapa:
//   ana->Draw("cand_prange", "cand_pop==0 && cand_level>=6")
// accum_level[rama]: el del mejor candidato de cada población (+1 por el corte "hay algún candidato").
// Categorías por candidato: "candmu" (signal muons en capas A / B no A / C no B, pi, p, mu cósmico, EM) y las estándar con prefijo
// "cand" (candparticle, ...).
// MC: los eventos con un K+ en la cadena del haz tienen volcado completo (fulldump = 1: todos los PFParticles
// completos, más de lo que el módulo guarda en un evento normal y en datos). La selección se hace sobre la vista normal
// de esos eventos (muCandUtils::MakeNormalView: lo que el módulo habría guardado), así que la señal y el fondo se cuentan
// en todos los eventos, sin escalar, como en datos. cand_savemode dice cómo guardaría el módulo cada candidato en un
// evento normal (2 completo, 1 resumen). La vista completa no se guarda en cand_* (ver truemu_level_full).
//
// Árbol "truth": una entrada por mu+ de desintegración no cósmico con p true en [230, 241] MeV/c (de K+ en reposo),
// con el origen de su K (truemu_korigin) y accum_level[rama] = 1 + nivel del mejor candidato emparejado en esa rama
// (0 si no hay candidato emparejado), en la vista normal: la eficiencia que se tendría en datos. truemu_level_full[rama]
// es lo mismo con toda la información guardada (vista completa de los eventos con volcado; en los demás, igual que
// accum_level): la diferencia es lo que se pierde porque el módulo guarda en resumen o no guarda.
namespace secondaryKaonProductionXSAnalysisConstants{
  // Máximo de candidatos guardados por evento. No puede pasar de 100: las categorías por objeto de highland
  // (AddStandardObjectCategories) tienen ese tamaño fijo. Con 500 eventos de MC, el percentil 99 es 80 y el máximo
  // 118; si hay más, se guardan primero los que llegan más lejos en la cadena (ver SavedCandidates).
  const UInt_t NMAXMUCANDIDATES = 100;
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

  bool CheckFillTruthTree(const AnaTrueVertex& vtx){(void)vtx; return true;}
  bool CheckFillTruthTree(const AnaTrueParticlePD* part);
  bool FinalizeConfiguration();

  void FillTruthTree();
  void FillTruthTree(const AnaTrueVertex& vtx){(void)vtx; return;}
  void FillTruthTree(const AnaTrueParticlePD& part);
  //--------------------

  bool Initialize();
  void FillCategories();

  /// Guarda en el árbol config el commit de highlandPD con el que se corre (highlandpd_version)
  void FillConfigTree();
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

  /// Returns the ToyBoxMuCand
  virtual const ToyBoxMuCand& box(Int_t isel=-1) const {return *static_cast<const ToyBoxMuCand*>(&boxB(isel));}

  /// Returns the vertex for the ToyBoxMuCand
  virtual AnaVertexB* GetVertex() const{return box().Vertex;}

  /// Returns the true vertex for the ToyBoxMuCand
  virtual AnaTrueVertexB* GetTrueVertex() const {return box().TrueVertex;}

  /// Índices de los candidatos que se guardan en el árbol y en las categorías, en este orden: como mucho
  /// NMAXMUCANDIDATES, primero los de mayor nivel en la cadena (a igual nivel, en el orden original)
  std::vector<int> SavedCandidates() const;

  /// Para el árbol truth: 1 + nivel del mejor candidato emparejado en cada rama, en la vista completa, por ID de la
  /// partícula true (lo rellena FinalizeConfiguration y lo vacía FillTruthTree)
  std::map<int, std::array<int, 2>> _fullLevel;

  /// Categoría true de un candidato (código de la categoría "candmu"): los signal muons en capas excluyentes, 1 señal A,
  /// 2 señal B que no es A, 3 señal C que no es B; el fondo con su PDG
  static int CandidateCategory(const muCandUtils::Candidate& c);

  /// Categoría true de un candidato para una medida (muCandUtils::kMeasC, kMeasB, kMeasA; categorías "candsigC",
  /// "candsigB", "candsigA"): 1 señal de esa medida, 2 signal muon que no es señal de esa medida, y el resto como en candmu
  static int MeasurementCategory(const muCandUtils::Candidate& c, int measurement);

  /// Especie del haz agrupada como la separan los Cherenkov a 6 GeV: 0 e/mu/pi, 1 K, 2 p, 3 ambigua, -1 sin información.
  /// Desde la lista de PDGs compatibles de la instrumentación (datos; en el MC no hay emulación de los Cherenkov y la
  /// lista es siempre {2212}) o desde un PDG true
  static int BeamSpecies(const std::vector<int>& pdgs);
  static int BeamSpecies(int truePDG){ return BeamSpecies(std::vector<int>{truePDG}); }

  /// ¿Es un signal muon para el árbol truth? (muCandUtils::IsSignalMuon)
  bool IsTrueMuFromKaonAtRest(const AnaTrueParticleB* t);

private:

  bool _FillBeamInstrumentationInfo;
  bool _FillBeamParticleInfo;
  bool _FillBeamParticleHitsInfo;
  bool _FillBeamParticleDaughtersInfo;

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

  // Variables de los árboles (el nombre de cada variable en el árbol es el del enum)
  enum enumMicroTrees_secondaryKaonProductionXSAnalysis{
    // evento
    fulldump = standardPDTree::enumStandardMicroTreesLast_standardPDTree,
    sel_fitlevel, sel_countlevel, hasslice, hasvertex, vertex, nmucand_presel, nmucand, hastruevertex, truevertex,
    beam_quality_sce, beam_start_sce, beam_cos_bi, beam_reach, beam_reach_n,
    beam_species, beam_truespecies,
    // candidatos (contador ncand)
    ncand,
    cand_pop, cand_level, cand_savemode, cand_isshower, cand_complete, cand_isprimary, cand_ndau,
    cand_length, cand_prange, cand_chi2mu, cand_chi2mu_other, cand_chi2mu_min, cand_chi2mu_trim, cand_chi2mu_trimcm, cand_chi2mu_trimboth, cand_michel_own, cand_nmichel_own, cand_chi2p, cand_dedx_stop, cand_dedx_start,
    cand_faces_stop, cand_faces_start, cand_axis_min, cand_vtx_min, cand_slice_min,
    cand_cnn_trk, cand_cnn_em, cand_cnn_michel, cand_clearcosmic, cand_cosmiccrit,
    cand_start, cand_stop, cand_cosbeam, cand_diry_stop, cand_trmean, cand_trmean_body, cand_intpeak, cand_intpeak_rr,
    cand_ar_n5, cand_ar_n10, cand_ar_n20, cand_ar_nmichel5, cand_ar_nmichel10, cand_ar_nmichel20,
    cand_ar_michel5, cand_ar_michel10, cand_ar_michel20, cand_ar_charge5, cand_ar_charge10, cand_ar_charge20,
    cand_as_n5, cand_as_n10, cand_as_charge5, cand_as_charge10,
    cand_k_n, cand_k_istrack, cand_k_dedx, cand_k_length,
    cand_ndau_trk, cand_ndau_shw, cand_dau_nearstop, cand_dau_maxmichel,
    cand_vc_dist, cand_vc_cover, cand_vc_maxgap,
    cand_lk_s_dv, cand_lk_s_de, cand_lk_s_len, cand_lk_s_dedx, cand_lk_p_dv, cand_lk_p_de,
    cand_dstart_v, cand_dstop_v,
    cand_lkt_s_dv, cand_lkt_s_de, cand_lkt_p_dv, cand_lkt_p_de, cand_dstart_tv, cand_dstop_tv,
    cand_sh_nsp, cand_sh_lin, cand_sh_lencore, cand_sh_lenbragg, cand_sh_prange_core,
    // truth de los candidatos
    cand_truepdg, cand_trueparentpdg, cand_trueproc, cand_trueorigin, cand_truefrombeam, cand_truesignal, cand_truesignalmuon, cand_truekendmom, cand_truekorigin,
    cand_truepur, cand_trueeff, cand_truelength, cand_truemom, cand_truevtxerr, cand_truelk_s_isk,
    cand_k_truek, cand_k_truepur, cand_k_trueeff, cand_truek_length,
    // árbol truth
    truemu_korigin, truemu_level_full, truemu_mom, truemu_kmomend, truemu_start, truemu_end, truemu_contained,
    truek_length, truek_nreco, truek_besteff, truek_bestpur, truek_bestlength, truek_bestistrack,
    // árbol config
    highlandpd_version,
    enumMicroTreesLast_secondaryKaonProductionXSAnalysis
  };
};

#endif
