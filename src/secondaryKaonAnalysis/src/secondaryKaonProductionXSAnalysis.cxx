#include "secondaryKaonProductionXSAnalysis.hxx"
#include "secondaryKaonProductionXSSelection.hxx"
#include "Parameters.hxx"
#include "CategoriesUtils.hxx"
#include "BasicUtils.hxx"
#include "pdAnalysisUtils.hxx"
#include "baseToyMaker.hxx"
#include "pdDataClasses.hxx"
#include "kaonAnalysisUtils.hxx"
#include "muCandidateUtils.hxx"
#include "TSystem.h"
#include <algorithm>

#include "HighlandMiniTreeConverter.hxx"
#include "PDSPAnalyzerTreeConverter.hxx"

#include "RecombinationVariation.hxx"
#include "dEdxVariation.hxx"
#include "BeamPartIdEffWeight.hxx"
#include "BrokenTrackWeight.hxx"
#include "SCEGeometricVariation.hxx"
#include "NominalBeamMomWeight.hxx"
#include "BeamCompositionWeight.hxx"
#include "ProtonBackgroundWeight.hxx"

#include "ParticlePositionSCECorrection.hxx"
#include "dEdxMCCorrection.hxx"
#include "RecombinationDataCorrection.hxx"

//********************************************************************
secondaryKaonProductionXSAnalysis::secondaryKaonProductionXSAnalysis(AnalysisAlgorithm* ana) : baseAnalysis(ana) {
//********************************************************************

}

namespace {
  // Commit de highlandPD ($HIGHLANDPDPATH), con "+cambios" si hay ficheros seguidos por git modificados sin commitear
  std::string HighlandPDVersion(){
    const char* path = gSystem->Getenv("HIGHLANDPDPATH");
    if(!path)return "desconocida";
    TString hash = gSystem->GetFromPipe(Form("git -C %s rev-parse --short HEAD 2>/dev/null", path));
    if(hash.IsNull())return "desconocida";
    TString dirty = gSystem->GetFromPipe(Form("git -C %s status --porcelain --untracked-files=no 2>/dev/null | head -1", path));
    return std::string(hash.Data()) + (dirty.IsNull() ? "" : "+cambios");
  }
}

//********************************************************************
void secondaryKaonProductionXSAnalysis::FillConfigTree(){
//********************************************************************

  output().AddVar(OutputManager::config, highlandpd_version, "highlandpd_version", "C", "commit de highlandPD con el que se corrió");
  output().FillVar(highlandpd_version, HighlandPDVersion());
}

//********************************************************************
bool secondaryKaonProductionXSAnalysis::Initialize(){
//********************************************************************

  std::cout << "secondaryKaonProductionXSAnalysis: highlandPD " << HighlandPDVersion() << std::endl;

  // Initialize the baseAnalysis (highland/src/highland2/baseAnalysis)
  if(!baseAnalysis::Initialize()) return false;

  // Minimum accum cut level (how many cuts should be passed) to save event into the output tree
  SetMinAccumCutLevelToSave(ND::params().GetParameterI("secondaryKaonAnalysis.xsAnalysis.MinAccumLevelToSave"));

  _FillBeamInstrumentationInfo   = ND::params().GetParameterI("secondaryKaonAnalysis.xsAnalysis.FillBeamInstrumentationInfo");
  _FillBeamParticleInfo          = ND::params().GetParameterI("secondaryKaonAnalysis.xsAnalysis.FillBeamParticleInfo");
  _FillBeamParticleHitsInfo      = ND::params().GetParameterI("secondaryKaonAnalysis.xsAnalysis.FillBeamParticleHitsInfo");
  _FillBeamParticleDaughtersInfo = ND::params().GetParameterI("secondaryKaonAnalysis.xsAnalysis.FillBeamParticleDaughtersInfo");

  _ApplydQdxSystematic = ND::params().GetParameterI("secondaryKaonAnalysis.xsAnalysis.ApplydQdxSystematic");
  _ApplyRecombinationSystematic = ND::params().GetParameterI("secondaryKaonAnalysis.xsAnalysis.ApplyRecombinationSystematic");
  _ApplySCESystematic = ND::params().GetParameterI("secondaryKaonAnalysis.xsAnalysis.ApplySCESystematic");
  _ApplyBrokenTracksSystematic = ND::params().GetParameterI("secondaryKaonAnalysis.xsAnalysis.ApplyBrokenTracksSystematic");
  _ApplyBeamPIDEfficiencySystematic = ND::params().GetParameterI("secondaryKaonAnalysis.xsAnalysis.ApplyBeamPIDEfficiencySystematic");
  _ApplyBeamPartWeightSystematic = ND::params().GetParameterI("secondaryKaonAnalysis.xsAnalysis.ApplyBeamPartWeightSystematic");
  _ApplyBeamMomWeightSystematic = ND::params().GetParameterI("secondaryKaonAnalysis.xsAnalysis.ApplyBeamMomWeightSystematic");
  _ApplyProtonBackgroundWeightSystematic = ND::params().GetParameterI("secondaryKaonAnalysis.xsAnalysis.ApplyProtonBackgroundWeightSystematic");

  _ApplySCECorrection = ND::params().GetParameterI("secondaryKaonAnalysis.xsAnalysis.ApplySCECorrection");
  _ApplydEdxCorrection = ND::params().GetParameterI("secondaryKaonAnalysis.xsAnalysis.ApplydEdxCorrection");
  _ApplyRecombinationCorrection = ND::params().GetParameterI("secondaryKaonAnalysis.xsAnalysis.ApplyRecombinationCorrection");

  // Categories for color drawing
  if(_FillBeamParticleInfo){
    anaUtils::AddStandardCategories();        // primario de la beam slice
    kaonAnaUtils::AddBeamParticleReducedCategory();
    kaonAnaUtils::AddBeamParticleOriginCategory();
  }
  if(_FillBeamInstrumentationInfo)
    anaUtils::AddStandardCategories("beam");  // partícula de la instrumentación del haz
  if(_FillBeamParticleDaughtersInfo)
    anaUtils::AddStandardObjectCategories("dau", standardPDTree::seltrk_ndau, "seltrk_ndau", 1);

  // candidatos a mu: categorías estándar (candparticle, ...) y la propia (candmu)
  anaUtils::AddStandardObjectCategories("cand", ncand, "ncand", 1);
  {
    // signal muons en capas excluyentes (cada candidato un solo código): A; B que no es A; C que no es B.
    // La señal de cada medida es una unión de capas: A = {1}, B = {1, 2}, C = {1, 2, 3}
    std::string types[] = {"signal A", "signal B (not A)", "signal C (not B)",
                           "#pi^{#pm}", "p", "cosmic #mu", "e, #gamma", NAMEOTHER, NAMENOTRUTH};
    int codes[]          = {1, 2, 3, 211, 2212, 13, 11, CATOTHER, CATNOTRUTH};
    int colors[]         = {2, 46, 28, 4, 8, 7, 5, COLOTHER, COLNOTRUTH};
    const int N = sizeof(types)/sizeof(types[0]);
    std::reverse(types, types + N); std::reverse(codes, codes + N); std::reverse(colors, colors + N);
    anaUtils::_categ->AddObjectCategory("candmu", ncand, "ncand", N, types, codes, colors, 1, -(Int_t)secondaryKaonProductionXSAnalysisConstants::NMAXMUCANDIDATES);
  }
  // una categoría por medida (señal en capas: C cualquier signal muon, B de la cascada del haz que dispara, A de la
  // primera interacción): 1 señal de la medida, 2 signal muon de otro origen, y el fondo como en candmu
  {
    const char* sigName[3] = {"candsigC", "candsigB", "candsigA"};
    const char* sigLabel[3] = {"signal C", "signal B", "signal A"};   // en inglés: van a las figuras de la nota
    for(int m = 0; m < 3; m++){
      std::string types[] = {sigLabel[m], "#mu from K^{+} (other origin)", "#pi^{#pm}", "p", "cosmic #mu", "e, #gamma", NAMEOTHER, NAMENOTRUTH};
      int codes[]          = {1, 2, 211, 2212, 13, 11, CATOTHER, CATNOTRUTH};
      int colors[]         = {2, 28, 4, 8, 7, 5, COLOTHER, COLNOTRUTH};
      const int N = sizeof(types)/sizeof(types[0]);
      std::reverse(types, types + N); std::reverse(codes, codes + N); std::reverse(colors, colors + N);
      anaUtils::_categ->AddObjectCategory(sigName[m], ncand, "ncand", N, types, codes, colors, 1, -(Int_t)secondaryKaonProductionXSAnalysisConstants::NMAXMUCANDIDATES);
    }
  }

  return true;
}

//********************************************************************
void secondaryKaonProductionXSAnalysis::DefineInputConverters(){
//********************************************************************

  // add a single converter (a copy of the one in highland/baseAnalysis)
  input().AddConverter("minitree",         new HighlandMiniTreeConverter("highlandana/MiniTree"));
  input().AddConverter("minitreefiltered", new HighlandMiniTreeConverter("MiniTree"));
  input().AddConverter("PDSPAnalyzerTree", new PDSPAnalyzerTreeConverter());
}

//********************************************************************
void secondaryKaonProductionXSAnalysis::DefineSelections(){
//********************************************************************

  sel().AddSelection("XS", "kaon XS selection", new secondaryKaonProductionXSSelection(false));
}

//********************************************************************
void secondaryKaonProductionXSAnalysis::DefineCorrections(){
//********************************************************************

  baseAnalysis::DefineCorrections();
  if(_ApplySCECorrection)
    corr().AddCorrection(0, "sce geometric correction"      , new ParticlePositionSCECorrection());
  if(_ApplydEdxCorrection)
    corr().AddCorrection(1, "dEdx MC correction"            , new dEdxMCCorrection());
  if(_ApplyRecombinationCorrection)
    corr().AddCorrection(2, "Recombination data correction" , new RecombinationDataCorrection());
}

//********************************************************************
void secondaryKaonProductionXSAnalysis::DefineSystematics(){
//********************************************************************

  // Some systematics are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineSystematics();

  if(_ApplySCESystematic)
    evar().AddEventVariation(kSCEGeometric     ,"SCE variation",           new SCEGeometricVariation());
  if(_ApplydQdxSystematic)
    evar().AddEventVariation(kdEdxVariation  ,"dEdx variation",          new dEdxVariation());
  if(_ApplyRecombinationSystematic)
    evar().AddEventVariation(kRecombination    ,"Recombination variation", new RecombinationVariation());
  if(_ApplyBrokenTracksSystematic)
    eweight().AddEventWeight(kBrokenTracks     ,"Broken track weight",     new BrokenTrackWeight());
  if(_ApplyBeamPIDEfficiencySystematic)				          
    eweight().AddEventWeight(kBeamPIDEfficiency,"Beam PID efficiency",     new BeamPartIdEffWeight());
  if(_ApplyBeamMomWeightSystematic)				          
    eweight().AddEventWeight(kBeamMomWeight    ,"Beam Nom Mom Weight",     new NominalBeamMomWeight());
  if(_ApplyBeamPartWeightSystematic)
    eweight().AddEventWeight(kBeamPartWeight   ,"Beam Particle Weight",    new BeamCompositionWeight());
  if(_ApplyProtonBackgroundWeightSystematic)
    eweight().AddEventWeight(kProtonBackgroundWeight   ,"ProtonBackground Weight",    new ProtonBackgroundWeight());
}

//********************************************************************
void secondaryKaonProductionXSAnalysis::DefineConfigurations(){
//********************************************************************

  // Some configurations are defined in baseAnalysis (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineConfigurations();
  // Enable all variation systematics in the all_syst configuration (created in baseAnalysis)
  if(_enableAllSystConfig){
    if(_ApplySCESystematic)
      conf().EnableEventVariation(kSCEGeometric,all_syst);
    if(_ApplydQdxSystematic)
      conf().EnableEventVariation(kdEdxVariation,all_syst);
    if(_ApplyRecombinationSystematic)
      conf().EnableEventVariation(kRecombination,all_syst);
    if(_ApplyBrokenTracksSystematic)
      conf().EnableEventWeight(kBrokenTracks,all_syst);
    if(_ApplyBeamPIDEfficiencySystematic)
      conf().EnableEventWeight(kBeamPIDEfficiency,all_syst);
    if(_ApplyBeamPartWeightSystematic)
      conf().EnableEventWeight(kBeamPartWeight,all_syst);
    if(_ApplyBeamMomWeightSystematic)
      conf().EnableEventWeight(kBeamMomWeight,all_syst);
    if(_ApplyProtonBackgroundWeightSystematic)
      conf().EnableEventWeight(kProtonBackgroundWeight,all_syst);
  }
}

//********************************************************************
void secondaryKaonProductionXSAnalysis::DefineMicroTrees(bool addBase){
//********************************************************************

  // Variables from baseAnalysis (run, event, ...)   (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::DefineMicroTrees(addBase);

  // Add standard sets of variables for ProtoDUNE analysis  (those methods are in highlandPD/src/pdUtils/standardPDTree.cxx)
  standardPDTree::AddStandardVariables_EventInfo(output());
  if(_FillBeamInstrumentationInfo){
    standardPDTree::AddStandardVariables_BeamInstrumentationReco(output());
    standardPDTree::AddStandardVariables_BeamInstrumentationTrue(output());
  }
  if(_FillBeamParticleInfo){
    standardPDTree::AddStandardVariables_BeamParticleTrue(output());
    standardPDTree::AddStandardVariables_BeamParticleReco(output());
    if(_FillBeamParticleHitsInfo)
      standardPDTree::AddStandardVariables_BeamParticleHitsReco(output());
    if(_FillBeamParticleDaughtersInfo){
      standardPDTree::AddStandardVariables_BeamParticleDaughtersTrue(output(),20);
      standardPDTree::AddStandardVariables_BeamParticleDaughtersReco(output(),20);
    }
  }

  // -------- evento --------
  AddVarI  (output(), fulldump,       "MC: volcado completo (evento con un K+ secundario)");
  AddVarI  (output(), hasslice,       "hay beam slice");
  AddVarI  (output(), hasvertex,      "hay vértice (final del primario de la beam slice)");
  AddVar3VF(output(), vertex,         "vértice (cm, sin SCE)");
  AddVarI  (output(), nmucand_presel, "candidatos antes de la preselección geométrica del módulo");
  AddVarI  (output(), nmucand,        "candidatos (puede ser mayor que ncand si se pasa del máximo)");
  AddVarI  (output(), beam_species,     "especie del haz según la instrumentación: 0 e/mu/pi, 1 K, 2 p, 3 ambigua, -1 sin información (en MC siempre 2: sin emulación de los Cherenkov)");
  AddVarI  (output(), beam_truespecies, "MC: especie true del haz agrupada igual (0 e/mu/pi, 1 K, 2 p, -1 otra)");
  AddVarI  (output(), hastruevertex,  "MC: hay vértice perfecto (hits de la partícula true del haz)");
  AddVar3VF(output(), truevertex,     "MC: vértice perfecto, último hit del haz true a lo largo del eje (cm, sistema de la reco)");

  // -------- candidatos (contador ncand) --------
  const UInt_t N = secondaryKaonProductionXSAnalysisConstants::NMAXMUCANDIDATES;
  AddVarMaxSizeVI (output(), cand_pop,           "población: 0 en la beam slice, 1 fuera",                       ncand, N);
  AddVarMaxSizeVI (output(), cand_level,         "cortes de la cadena de su población que pasa seguidos",        ncand, N);
  AddVarMaxSizeVI (output(), cand_isshower,      "sin traza (longitud de la traza forzada)",                     ncand, N);
  AddVarMaxSizeVI (output(), cand_savemode,      "cómo lo guardaría el módulo en un evento normal (2 completo, 1 resumen; en los eventos con volcado, emulado)", ncand, N);
  AddVarMaxSizeVI (output(), cand_complete,      "con hits (no en resumen)",                                     ncand, N);
  AddVarMaxSizeVI (output(), cand_isprimary,     "es el primario de Pandora",                                    ncand, N);
  AddVarMaxSizeVI (output(), cand_ndau,          "hijas de Pandora",                                             ncand, N);
  AddVarMaxSizeVF (output(), cand_length,        "longitud (cm)",                                                ncand, N);
  AddVarMaxSizeVF (output(), cand_prange,        "p de rango, hipótesis mu (MeV/c)",                             ncand, N);
  AddVarMaxSizeVF (output(), cand_chi2mu,        "chi2/ndf de mu con la parada en el extremo de mayor dE/dx",    ncand, N);
  AddVarMaxSizeVF (output(), cand_chi2mu_other,  "chi2/ndf de mu con la orientación contraria",                  ncand, N);
  AddVarMaxSizeVF (output(), cand_chi2mu_min,    "chi2/ndf de mu, el mejor de las dos orientaciones",            ncand, N);
  AddVarMaxSizeVF (output(), cand_chi2mu_trim,   "chi2/ndf de mu recortando 0-10 cm de la parada, el mejor (muCandUtils::Chi2BestTrim)", ncand, N);
  AddVarMaxSizeVF (output(), cand_chi2mu_trimcm, "recorte de cand_chi2mu_trim (cm)",                             ncand, N);
  AddVarMaxSizeVF (output(), cand_chi2mu_trimboth, "como cand_chi2mu_trim, el mejor de las dos orientaciones",   ncand, N);
  AddVarMaxSizeVF (output(), cand_chi2p,         "chi2/ndf de protón, misma orientación",                        ncand, N);
  AddVarMaxSizeVF (output(), cand_dedx_stop,     "dE/dx medio en los 5 cm de la parada (MeV/cm)",                ncand, N);
  AddVarMaxSizeVF (output(), cand_dedx_start,    "dE/dx medio en los 5 cm del inicio (MeV/cm)",                  ncand, N);
  AddVarMaxSizeVF (output(), cand_faces_stop,    "distancia de la parada a las caras (cm)",                      ncand, N);
  AddVarMaxSizeVF (output(), cand_faces_start,   "distancia del inicio a las caras (cm)",                        ncand, N);
  AddVarMaxSizeVF (output(), cand_axis_min,      "distancia mínima de un extremo al eje del haz (cm)",           ncand, N);
  AddVarMaxSizeVF (output(), cand_vtx_min,       "distancia mínima de un extremo al vértice (cm; -1 sin vértice)", ncand, N);
  AddVarMaxSizeVF (output(), cand_slice_min,     "distancia mínima de un extremo a la beam slice (cm; 0 dentro)", ncand, N);
  AddVarMaxSizeVF (output(), cand_cnn_trk,       "CNN traza media de la partícula",                              ncand, N);
  AddVarMaxSizeVF (output(), cand_cnn_em,        "CNN EM media de la partícula",                                 ncand, N);
  AddVarMaxSizeVF (output(), cand_cnn_michel,    "CNN Michel media de la partícula",                             ncand, N);
  AddVarMaxSizeVI (output(), cand_clearcosmic,   "IsClearCosmic de Pandora",                                     ncand, N);
  AddVarMaxSizeVI (output(), cand_cosmiccrit,    "CosmicCriteria del módulo",                                    ncand, N);
  AddVarMaxSize3MF(output(), cand_start,         "posición del inicio (cm, sin SCE)",                            ncand, N);
  AddVarMaxSize3MF(output(), cand_stop,          "posición de la parada (cm, sin SCE)",                          ncand, N);
  AddVarMaxSizeVF (output(), cand_cosbeam,       "coseno inicio -> parada con el eje del haz",                   ncand, N);
  AddVarMaxSizeVF (output(), cand_diry_stop,     "componente y de la dirección hacia la parada",                 ncand, N);
  AddVarMaxSizeVF (output(), cand_trmean,        "media truncada del dE/dx (MeV/cm)",                            ncand, N);
  AddVarMaxSizeVF (output(), cand_trmean_body,   "media truncada sin 5 cm de cada extremo (MeV/cm)",             ncand, N);
  AddVarMaxSizeVF (output(), cand_intpeak,       "pico de dE/dx interior (media de 5 hits; MeV/cm)",             ncand, N);
  AddVarMaxSizeVF (output(), cand_intpeak_rr,    "residual range del pico interior (cm)",                        ncand, N);
  AddVarMaxSizeVI (output(), cand_ar_n5,         "hits (plano 2) a < 5 cm de la parada",                         ncand, N);
  AddVarMaxSizeVI (output(), cand_ar_n10,        "hits (plano 2) a < 10 cm de la parada",                        ncand, N);
  AddVarMaxSizeVI (output(), cand_ar_n20,        "hits (plano 2) a < 20 cm de la parada",                        ncand, N);
  AddVarMaxSizeVI (output(), cand_ar_nmichel5,   "hits con CNN Michel > 0,5 a < 5 cm de la parada",              ncand, N);
  AddVarMaxSizeVI (output(), cand_ar_nmichel10,  "hits con CNN Michel > 0,5 a < 10 cm de la parada",             ncand, N);
  AddVarMaxSizeVI (output(), cand_ar_nmichel20,  "hits con CNN Michel > 0,5 a < 20 cm de la parada",             ncand, N);
  AddVarMaxSizeVF (output(), cand_ar_michel5,    "CNN Michel media de los hits a < 5 cm de la parada",           ncand, N);
  AddVarMaxSizeVF (output(), cand_ar_michel10,   "CNN Michel media de los hits a < 10 cm de la parada",          ncand, N);
  AddVarMaxSizeVF (output(), cand_ar_michel20,   "CNN Michel media de los hits a < 20 cm de la parada",          ncand, N);
  AddVarMaxSizeVF (output(), cand_ar_charge5,    "carga a < 5 cm de la parada (ADC x tick)",                     ncand, N);
  AddVarMaxSizeVF (output(), cand_ar_charge10,   "carga a < 10 cm de la parada (ADC x tick)",                    ncand, N);
  AddVarMaxSizeVF (output(), cand_ar_charge20,   "carga a < 20 cm de la parada (ADC x tick)",                    ncand, N);
  AddVarMaxSizeVI (output(), cand_as_n5,         "hits a < 5 cm del inicio",                                     ncand, N);
  AddVarMaxSizeVI (output(), cand_as_n10,        "hits a < 10 cm del inicio",                                    ncand, N);
  AddVarMaxSizeVF (output(), cand_as_charge5,    "carga a < 5 cm del inicio",                                    ncand, N);
  AddVarMaxSizeVF (output(), cand_as_charge10,   "carga a < 10 cm del inicio",                                   ncand, N);
  AddVarMaxSizeVI (output(), cand_k_n,           "PFParticles con un extremo a < 10 cm del inicio",              ncand, N);
  AddVarMaxSizeVI (output(), cand_k_istrack,     "el de mayor dE/dx de esos es traza",                           ncand, N);
  AddVarMaxSizeVF (output(), cand_k_dedx,        "su dE/dx en ese extremo (MeV/cm)",                             ncand, N);
  AddVarMaxSizeVF (output(), cand_k_length,      "su longitud (cm)",                                             ncand, N);
  AddVarMaxSizeVI (output(), cand_ndau_trk,      "hijas traza (Pandora)",                                        ncand, N);
  AddVarMaxSizeVI (output(), cand_ndau_shw,      "hijas shower (Pandora)",                                       ncand, N);
  AddVarMaxSizeVI (output(), cand_dau_nearstop,  "hijas con un extremo a < 10 cm de la parada",                  ncand, N);
  AddVarMaxSizeVF (output(), cand_dau_maxmichel, "CNN Michel máxima de las hijas",                               ncand, N);
  AddVarMaxSizeVF (output(), cand_vc_dist,       "distancia vértice -> inicio (cm)",                             ncand, N);
  AddVarMaxSizeVF (output(), cand_vc_cover,      "fracción del segmento vértice -> inicio con puntos a < 2 cm",  ncand, N);
  AddVarMaxSizeVF (output(), cand_vc_maxgap,     "tramo vacío más largo de ese segmento (cm)",                   ncand, N);
  AddVarMaxSizeVF (output(), cand_lk_s_dv,       "enlace vértice -> inicio: distancia de su extremo al vértice (cm)", ncand, N);
  AddVarMaxSizeVF (output(), cand_lk_s_de,       "enlace vértice -> inicio: distancia de su otro extremo al inicio (cm)", ncand, N);
  AddVarMaxSizeVF (output(), cand_lk_s_len,      "enlace vértice -> inicio: longitud (cm)",                      ncand, N);
  AddVarMaxSizeVF (output(), cand_lk_s_dedx,     "enlace vértice -> inicio: dE/dx en el extremo del inicio",     ncand, N);
  AddVarMaxSizeVF (output(), cand_lk_p_dv,       "enlace vértice -> parada: distancia al vértice (cm)",          ncand, N);
  AddVarMaxSizeVF (output(), cand_lk_p_de,       "enlace vértice -> parada: distancia a la parada (cm)",         ncand, N);
  AddVarMaxSizeVF (output(), cand_dstart_v,      "distancia inicio - vértice (cm)",                              ncand, N);
  AddVarMaxSizeVF (output(), cand_dstop_v,       "distancia parada - vértice (cm)",                              ncand, N);
  AddVarMaxSizeVF (output(), cand_lkt_s_dv,      "MC, vértice perfecto: enlace con el inicio, extremo - vértice (cm)", ncand, N);
  AddVarMaxSizeVF (output(), cand_lkt_s_de,      "MC, vértice perfecto: enlace con el inicio, otro extremo - inicio (cm)", ncand, N);
  AddVarMaxSizeVF (output(), cand_lkt_p_dv,      "MC, vértice perfecto: enlace con la parada, extremo - vértice (cm)", ncand, N);
  AddVarMaxSizeVF (output(), cand_lkt_p_de,      "MC, vértice perfecto: enlace con la parada, otro extremo - parada (cm)", ncand, N);
  AddVarMaxSizeVF (output(), cand_dstart_tv,     "MC: distancia inicio - vértice perfecto (cm)",                 ncand, N);
  AddVarMaxSizeVF (output(), cand_dstop_tv,      "MC: distancia parada - vértice perfecto (cm)",                 ncand, N);
  AddVarMaxSizeVI (output(), cand_sh_nsp,        "shower: space points",                                         ncand, N);
  AddVarMaxSizeVF (output(), cand_sh_lin,        "shower: linealidad (valor propio mayor / suma)",               ncand, N);
  AddVarMaxSizeVF (output(), cand_sh_lencore,    "shower: longitud del núcleo (cm)",                             ncand, N);
  AddVarMaxSizeVF (output(), cand_sh_lenbragg,   "shower: longitud hasta el pico de Bragg (cm)",                 ncand, N);
  AddVarMaxSizeVF (output(), cand_sh_prange_core,"shower: p de rango con la longitud del núcleo (MeV/c)",        ncand, N);
  // truth
  AddVarMaxSizeVI (output(), cand_truepdg,          "truth: PDG",                                                   ncand, N);
  AddVarMaxSizeVI (output(), cand_trueparentpdg,    "truth: PDG de la madre (0 sin madre guardada)",                ncand, N);
  AddVarMaxSizeVI (output(), cand_trueproc,    "truth: proceso de creación",                                   ncand, N);
  AddVarMaxSizeVI (output(), cand_trueorigin,       "truth: Origin (2 cósmico, 4 haz)",                             ncand, N);
  AddVarMaxSizeVI (output(), cand_truefrombeam,     "truth: desciende del haz que dispara",                         ncand, N);
  AddVarMaxSizeVI (output(), cand_truesignal,       "truth: mu+ de K+ de la cadena del haz (señal antigua)",        ncand, N);
  AddVarMaxSizeVI (output(), cand_truesignalmuon,   "truth: signal muon (mu+ de K+ -> mu+ nu de una K+ en reposo)", ncand, N);
  AddVarMaxSizeVF (output(), cand_truekendmom,      "truth: momento final de su K madre (MeV/c; 0 en reposo; -1 si la madre no es K+)", ncand, N);
  AddVarMaxSizeVI (output(), cand_truekorigin,      "truth: origen de su K (0 1a interacción, 1 con K->K, 2 secundaria, 3 K del haz, 4 cadena cortada, 6 otra raíz, -1 no es mu+ de K+)", ncand, N);
  AddVarMaxSizeVF (output(), cand_truepur,          "truth: pureza",                                                ncand, N);
  AddVarMaxSizeVF (output(), cand_trueeff,          "truth: completitud",                                           ncand, N);
  AddVarMaxSizeVF (output(), cand_truelength,       "truth: distancia inicio-fin (cm)",                             ncand, N);
  AddVarMaxSizeVF (output(), cand_truemom,          "truth: momento inicial (MeV/c)",                               ncand, N);
  AddVarMaxSizeVF (output(), cand_truevtxerr,       "truth: distancia del vértice reco al final true del haz (cm)", ncand, N);
  AddVarMaxSizeVI (output(), cand_truelk_s_isk,     "truth: el enlace vértice -> inicio es la K madre",             ncand, N);
}

//********************************************************************
void secondaryKaonProductionXSAnalysis::DefineTruthTree(){
//********************************************************************

  // Variables from baseAnalysis (run, event, ...)   (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineTruthTree();

  if(_FillBeamInstrumentationInfo)
    standardPDTree::AddStandardVariables_BeamInstrumentationTrue(output());

  AddVarI  (output(), truemu_korigin,   "origen de la K (0 1a interacción, 1 con K->K, 2 secundaria, 3 K del haz, 4 cadena cortada, 6 otra raíz, -1 no es mu+ de K+)");
  AddVarFixVI(output(), truemu_level_full, "como accum_level[rama], con toda la información guardada (vista completa de los eventos con volcado)", 2);
  AddVarF  (output(), truemu_mom,       "momento del mu (MeV/c)");
  AddVarF  (output(), truemu_kmomend,   "momento final de la K madre (MeV/c; -1 sin madre K guardada)");
  AddVar3VF(output(), truemu_start,     "inicio del mu (cm)");
  AddVar3VF(output(), truemu_end,       "final del mu (cm)");
  AddVarI  (output(), truemu_contained, "inicio y final en el volumen activo");
}

//********************************************************************
void secondaryKaonProductionXSAnalysis::FillMicroTrees(bool addBase){
//********************************************************************

  // Variables from baseAnalysis (run, event, ...)  (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::FillMicroTreesBase(addBase);

  // Fill standard variables for the PD analysis (those methods are in highlandPD/src/pdUtils/standardPDTree.cxx)
  standardPDTree::FillStandardVariables_EventInfo(output(), &GetEvent());
  if(_FillBeamInstrumentationInfo){
    standardPDTree::FillStandardVariables_BeamInstrumentationReco(output(), GetSpill().Beam);
    standardPDTree::FillStandardVariables_BeamInstrumentationTrue(output(), GetSpill().Beam);
  }
  // sin beam slice no hay primario (MainTrack nulo): las funciones de standardPDTree no rellenan nada
  if(_FillBeamParticleInfo && box().MainTrack){
    standardPDTree::FillStandardVariables_BeamParticleReco(output(), box().MainTrack);
    standardPDTree::FillStandardVariables_BeamParticleTrue(output(), box().MainTrack);
    if(_FillBeamParticleHitsInfo)
      standardPDTree::FillStandardVariables_BeamParticleHitsReco(output(), box().MainTrack);
    if(_FillBeamParticleDaughtersInfo){
      int ndau = std::min(20,(int)box().MainTrack->Daughters.size());
      for(int i = 0; i < ndau; i++){
        standardPDTree::FillStandardVariables_BeamParticleDaughtersReco(output(), static_cast<AnaParticlePD*>(box().MainTrack->Daughters[i]));
        standardPDTree::FillStandardVariables_BeamParticleDaughtersTrue(output(), static_cast<AnaParticlePD*>(box().MainTrack->Daughters[i]));
        output().IncrementCounter(standardPDTree::seltrk_ndau);
      }
    }
  }

  // -------- evento --------
  const AnaEventInfoPD* info = static_cast<AnaEventInfoPD*>(GetEvent().EventInfo);
  output().FillVar(fulldump, (Int_t)(info && info->FullDump));
  const ToyBoxMuCand& b = box();
  bool slice = false;
  for(const muCandUtils::Candidate& c : b.Candidates)if(c.hasSlice){ slice = true; break; }
  output().FillVar(hasslice, (Int_t)slice);
  output().FillVar(hasvertex, (Int_t)b.HasVertex);
  output().FillVectorVarFromArray(vertex, b.VertexPos, 3);
  output().FillVar(nmucand_presel, b.NPreselected);
  { Int_t n = 0; for(size_t i = 0; i < b.Candidates.size(); i++)if(b.View[i] == 0)n++; output().FillVar(nmucand, n); }
  {
    const AnaBeamPD* beam = static_cast<AnaBeamPD*>(GetSpill().Beam);
    output().FillVar(beam_species, beam ? BeamSpecies(beam->PDGs) : -1);
    const AnaTrueParticleB* bt = (beam && beam->BeamParticle) ? static_cast<AnaTrueParticleB*>(beam->BeamParticle->TrueObject) : nullptr;
    output().FillVar(beam_truespecies, bt ? BeamSpecies(bt->PDG) : -1);
  }
  output().FillVar(hastruevertex, (Int_t)b.HasTrueVertex);
  output().FillVectorVarFromArray(truevertex, b.TrueVertexPos, 3);

  // -------- candidatos --------
  // Uno por posición de los vectores cand_*, en el orden de SavedCandidates (el mismo que en FillCategories). Las
  // variables referidas a la parada (chi2, dE/dx, caras) usan el extremo de mayor dE/dx (ver muCandCuts).
  for(int i : SavedCandidates()){
    const muCandUtils::Candidate& c = b.Candidates[i];
    output().FillVectorVar(cand_pop,           b.Population[i]);
    output().FillVectorVar(cand_level,         b.Level[i]);
    output().FillVectorVar(cand_savemode,      b.SaveMode[i]);
    output().FillVectorVar(cand_isshower,      c.isShower);
    output().FillVectorVar(cand_complete,      c.complete);
    output().FillVectorVar(cand_isprimary,     c.isPrimary);
    output().FillVectorVar(cand_ndau,          c.nDaughters);
    output().FillVectorVar(cand_length,        c.length);
    output().FillVectorVar(cand_prange,        c.pRange);
    output().FillVectorVar(cand_chi2mu,        (Float_t)muCandCuts::Chi2MuStop(c));
    output().FillVectorVar(cand_chi2mu_other,  (Float_t)muCandCuts::Chi2MuOther(c));
    { const double s = muCandCuts::Chi2MuStop(c), o = muCandCuts::Chi2MuOther(c);
      output().FillVectorVar(cand_chi2mu_min, (Float_t)((o >= 0 && (s < 0 || o < s)) ? o : s)); }
    output().FillVectorVar(cand_chi2mu_trim,   c.chi2muTrim);
    output().FillVectorVar(cand_chi2mu_trimcm, c.chi2muTrimCm);
    output().FillVectorVar(cand_chi2mu_trimboth, c.chi2muTrimBoth);
    output().FillVectorVar(cand_chi2p,         (Float_t)muCandCuts::Chi2pStop(c));
    output().FillVectorVar(cand_dedx_stop,     (Float_t)muCandCuts::dEdxStop(c));
    output().FillVectorVar(cand_dedx_start,    (Float_t)muCandCuts::dEdxStart(c));
    output().FillVectorVar(cand_faces_stop,    (Float_t)muCandCuts::FacesStop(c));
    output().FillVectorVar(cand_faces_start,   (Float_t)muCandCuts::FacesStart(c));
    output().FillVectorVar(cand_axis_min,      (Float_t)muCandCuts::AxisMin(c));
    output().FillVectorVar(cand_vtx_min,       (Float_t)muCandCuts::VertexMin(c));
    output().FillVectorVar(cand_slice_min,     (Float_t)std::min(c.sliceA, c.sliceB));
    output().FillVectorVar(cand_cnn_trk,       c.cnnTrk);
    output().FillVectorVar(cand_cnn_em,        c.cnnEM);
    output().FillVectorVar(cand_cnn_michel,    c.cnnMichel);
    output().FillVectorVar(cand_clearcosmic,   c.clearCosmic);
    output().FillVectorVar(cand_cosmiccrit,    c.cosmicCriteria);
    const Float_t start[3] = {c.startX, c.startY, c.startZ}, stop[3] = {c.stopX, c.stopY, c.stopZ};
    output().FillMatrixVarFromArray(cand_start, start, 3);
    output().FillMatrixVarFromArray(cand_stop,  stop,  3);
    output().FillVectorVar(cand_cosbeam,       c.cosBeam);
    output().FillVectorVar(cand_diry_stop,     c.dirYStop);
    output().FillVectorVar(cand_trmean,        c.trMean);
    output().FillVectorVar(cand_trmean_body,   c.trMeanBody);
    output().FillVectorVar(cand_intpeak,       c.intPeak);
    output().FillVectorVar(cand_intpeak_rr,    c.intPeakRR);
    output().FillVectorVar(cand_ar_n5,         c.arN[0]);
    output().FillVectorVar(cand_ar_n10,        c.arN[1]);
    output().FillVectorVar(cand_ar_n20,        c.arN[2]);
    output().FillVectorVar(cand_ar_nmichel5,   c.arNMichel[0]);
    output().FillVectorVar(cand_ar_nmichel10,  c.arNMichel[1]);
    output().FillVectorVar(cand_ar_nmichel20,  c.arNMichel[2]);
    output().FillVectorVar(cand_ar_michel5,    c.arMichel[0]);
    output().FillVectorVar(cand_ar_michel10,   c.arMichel[1]);
    output().FillVectorVar(cand_ar_michel20,   c.arMichel[2]);
    output().FillVectorVar(cand_ar_charge5,    c.arCharge[0]);
    output().FillVectorVar(cand_ar_charge10,   c.arCharge[1]);
    output().FillVectorVar(cand_ar_charge20,   c.arCharge[2]);
    output().FillVectorVar(cand_as_n5,         c.asN5);
    output().FillVectorVar(cand_as_n10,        c.asN10);
    output().FillVectorVar(cand_as_charge5,    c.asCharge5);
    output().FillVectorVar(cand_as_charge10,   c.asCharge10);
    output().FillVectorVar(cand_k_n,           c.kN);
    output().FillVectorVar(cand_k_istrack,     c.kIsTrack);
    output().FillVectorVar(cand_k_dedx,        c.kDEdx);
    output().FillVectorVar(cand_k_length,      c.kLength);
    output().FillVectorVar(cand_ndau_trk,      c.nDauTrk);
    output().FillVectorVar(cand_ndau_shw,      c.nDauShw);
    output().FillVectorVar(cand_dau_nearstop,  c.dauNearStop);
    output().FillVectorVar(cand_dau_maxmichel, c.dauMaxMichel);
    output().FillVectorVar(cand_vc_dist,       c.vcDist);
    output().FillVectorVar(cand_vc_cover,      c.vcCover);
    output().FillVectorVar(cand_vc_maxgap,     c.vcMaxGap);
    output().FillVectorVar(cand_lk_s_dv,       c.lkSdV);
    output().FillVectorVar(cand_lk_s_de,       c.lkSdE);
    output().FillVectorVar(cand_lk_s_len,      c.lkSLen);
    output().FillVectorVar(cand_lk_s_dedx,     c.lkSdEdx);
    output().FillVectorVar(cand_lk_p_dv,       c.lkPdV);
    output().FillVectorVar(cand_lk_p_de,       c.lkPdE);
    output().FillVectorVar(cand_dstart_v,      c.dStartV);
    output().FillVectorVar(cand_dstop_v,       c.dStopV);
    output().FillVectorVar(cand_lkt_s_dv,      c.lktSdV);
    output().FillVectorVar(cand_lkt_s_de,      c.lktSdE);
    output().FillVectorVar(cand_lkt_p_dv,      c.lktPdV);
    output().FillVectorVar(cand_lkt_p_de,      c.lktPdE);
    output().FillVectorVar(cand_dstart_tv,     c.dStartTV);
    output().FillVectorVar(cand_dstop_tv,      c.dStopTV);
    output().FillVectorVar(cand_sh_nsp,        c.shNSP);
    output().FillVectorVar(cand_sh_lin,        c.shLin);
    output().FillVectorVar(cand_sh_lencore,    c.shLenCore);
    output().FillVectorVar(cand_sh_lenbragg,   c.shLenBragg);
    output().FillVectorVar(cand_sh_prange_core,c.shPRangeCore);
    output().FillVectorVar(cand_truepdg,          c.tPDG);
    output().FillVectorVar(cand_trueparentpdg,    c.tParentPDG);
    output().FillVectorVar(cand_trueproc,    c.tProcStart);
    output().FillVectorVar(cand_trueorigin,       c.tOrigin);
    output().FillVectorVar(cand_truefrombeam,     c.tFromBeam);
    output().FillVectorVar(cand_truesignal,       c.tSignal);
    output().FillVectorVar(cand_truesignalmuon,   c.tSignalMuon);
    output().FillVectorVar(cand_truekendmom,      c.tKMomEnd);
    output().FillVectorVar(cand_truekorigin,      c.tKOrigin);
    output().FillVectorVar(cand_truepur,          c.tPur);
    output().FillVectorVar(cand_trueeff,          c.tEff);
    output().FillVectorVar(cand_truelength,       c.tLength);
    output().FillVectorVar(cand_truemom,          c.tMom);
    output().FillVectorVar(cand_truevtxerr,       c.tVtxErr);
    output().FillVectorVar(cand_truelk_s_isk,     c.tLkSIsK);
    output().IncrementCounter(ncand);
  }
}

//********************************************************************
void secondaryKaonProductionXSAnalysis::FillToyVarsInMicroTrees(bool addBase){
//********************************************************************

  // Fill the common variables  (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::FillToyVarsInMicroTreesBase(addBase);
}

//********************************************************************
void secondaryKaonProductionXSAnalysis::FillTruthTree(){
//********************************************************************

  // Fill the "truth" tree

  if(!output().HasTree(OutputManager::truth))return;

  output().SetCurrentTree(OutputManager::truth);

  // Loop over all true particles
  std::vector<AnaTrueParticleB*> trueParts = GetSpill().TrueParticles;
  for(std::vector<AnaTrueParticleB*>::iterator it = trueParts.begin(); it!=trueParts.end(); it++) {
    AnaTrueParticlePD& part = *static_cast<AnaTrueParticlePD*>(*it);

    // Check if this particle needs to be saved in the truth tree
    if (!CheckFillTruthTree(&part)) continue;

    // Initialize the truth tree. We must do that for each saved particle since several entries may correspond to the same spill
    output().InitializeTree(OutputManager::truth);

    // accumulated cut levels to compute efficiencies. This is taken directly from the AnaTrueObject
    Int_t accumLevel=0;  
    for(std::vector<SelectionBase*>::iterator itf = sel().GetSelections().begin(); itf != sel().GetSelections().end(); itf++){
      if(!(*itf)->IsEnabled())continue;

      Int_t isel = (*itf)->GetEnabledIndex();

      for(UInt_t ibranch=0;ibranch<(*itf)->GetNBranches();ibranch++){
        if(part.AccumLevel.size()>0) accumLevel = part.AccumLevel[isel][ibranch];
        if(sel().GetNEnabledSelections()>1){
          if(sel().GetNMaxBranches()>1)
            output().FillMatrixVar(accum_level, accumLevel, isel, ibranch);
          else
            output().FillVectorVar(accum_level, accumLevel, isel);
        }
        else{
          if (sel().GetNMaxBranches()>1)
            output().FillVectorVar(accum_level, accumLevel, ibranch);
          else
            output().FillVar(accum_level, accumLevel);
        }
      }
    }

    // Reset the categories for the current track
    cat().ResetCurrentCategories();

    // Call the derive classes functions, that also fills the categories
    FillTruthTree(part);

    // Fill the truth tree provided the true codes for color drawing
    std::map< std::string, TrackCategoryDefinition* >::iterator its;
    Int_t categ_index = AnalysisAlgorithm::firstCategory;
    for(its=cat().GetCategories().begin();its!=cat().GetCategories().end();its++, categ_index++ ){
      std::string categ_name = its->first;
      TrackCategoryDefinition& categ = *(its->second);
      if(categ.IsMultiType()){
        for(unsigned int i=0;i<categ.GetNTypes();i++)
           output().FillVectorVar(categ_index, (int)cat().CheckCategoryType(categ_name,i),i);
      }
      else if(!categ.IsObject()) output().FillVar(categ_index, cat().GetCode(categ_name));
    }

    // Fill the tree
    output().GetTree(OutputManager::truth)->Fill();
  }
  _fullLevel.clear();   // es de este spill
}

//********************************************************************
int secondaryKaonProductionXSAnalysis::BeamSpecies(const std::vector<int>& pdgs){
//********************************************************************

  bool emupi = false, k = false, p = false;
  for(int pdg : pdgs){
    const int a = std::abs(pdg);
    if(a == 11 || a == 13 || a == 211)emupi = true;
    else if(a == 321)k = true;
    else if(a == 2212)p = true;
  }
  const int n = emupi + k + p;
  if(n == 0)return -1;
  if(n > 1)return 3;
  return emupi ? 0 : (k ? 1 : 2);
}

//********************************************************************
bool secondaryKaonProductionXSAnalysis::IsTrueMuFromKaonAtRest(const AnaTrueParticleB* t){
//********************************************************************

  if(!t || t->PDG != -13 || t->ProcessStart != AnaTrueParticleB::Decay)return false;   // filtro rápido antes del mapa
  std::map<int, AnaTrueParticleB*> byID;
  for(AnaTrueParticleB* x : GetSpill().TrueParticles)byID[x->ID] = x;
  return muCandUtils::IsSignalMuon(t, byID);
}

//********************************************************************
bool secondaryKaonProductionXSAnalysis::CheckFillTruthTree(const AnaTrueParticlePD* part){
//********************************************************************

  return IsTrueMuFromKaonAtRest(part);
}

//********************************************************************
void secondaryKaonProductionXSAnalysis::FillTruthTree(const AnaTrueParticlePD& part){
//********************************************************************

  // Event variables
  output().FillVar(run,    GetSpill().EventInfo->Run);
  output().FillVar(subrun, GetSpill().EventInfo->SubRun);
  output().FillVar(evt,    GetSpill().EventInfo->Event);

  if(_FillBeamInstrumentationInfo)
    standardPDTree::FillStandardVariables_BeamInstrumentationTrue(output(), GetSpill().Beam);

  std::map<int, AnaTrueParticleB*> byID;
  for(AnaTrueParticleB* t : GetSpill().TrueParticles)byID[t->ID] = t;
  const AnaBeamPD* beam = static_cast<AnaBeamPD*>(GetSpill().Beam);
  const AnaTrueParticleB* beamTrue = (beam && beam->BeamParticle) ? static_cast<AnaTrueParticleB*>(beam->BeamParticle->TrueObject) : nullptr;

  output().FillVar(truemu_korigin, muCandUtils::KOrigin(&part, beamTrue, byID));
  { auto it = _fullLevel.find(part.ID);
    for(int ib = 0; ib < 2; ib++)output().FillVectorVar(truemu_level_full, it == _fullLevel.end() ? 0 : it->second[ib], ib); }
  output().FillVar(truemu_mom, (Float_t)(1000.*part.Momentum));
  auto ik = byID.find(part.ParentID);
  output().FillVar(truemu_kmomend, (Float_t)((ik != byID.end() && ik->second->PDG == 321) ? 1000.*ik->second->MomentumEnd : -1.));
  output().FillVectorVarFromArray(truemu_start, part.Position, 3);
  output().FillVectorVarFromArray(truemu_end,   part.PositionEnd, 3);
  output().FillVar(truemu_contained, (Int_t)(muCandUtils::DistToFaces(muCandUtils::V(part.Position)) > 0 &&
                                             muCandUtils::DistToFaces(muCandUtils::V(part.PositionEnd)) > 0));
}

//********************************************************************
std::vector<int> secondaryKaonProductionXSAnalysis::SavedCandidates() const{
//********************************************************************

  // solo la vista normal (la de la selección; ver ToyBoxMuCand::View)
  const ToyBoxMuCand& b = box();
  std::vector<int> idx;
  for(size_t i = 0; i < b.Candidates.size(); i++)if(b.View[i] == 0)idx.push_back(i);
  // solo hace falta ordenar si hay más que el máximo: primero los que llegan más lejos en su cadena
  if(idx.size() > secondaryKaonProductionXSAnalysisConstants::NMAXMUCANDIDATES){
    std::stable_sort(idx.begin(), idx.end(), [&b](int x, int y){ return b.Level[x] > b.Level[y]; });
    idx.resize(secondaryKaonProductionXSAnalysisConstants::NMAXMUCANDIDATES);
  }
  return idx;
}

//********************************************************************
int secondaryKaonProductionXSAnalysis::CandidateCategory(const muCandUtils::Candidate& c){
//********************************************************************

  if(c.tPDG == 0)return CATNOTRUTH;
  if(c.tSignalMuon){
    if(muCandUtils::IsSignalOf(muCandUtils::kMeasA, true, c.tKOrigin))return 1;   // A
    if(muCandUtils::IsSignalOf(muCandUtils::kMeasB, true, c.tKOrigin))return 2;   // B que no es A (secundaria, K del haz)
    return 3;                                                                      // C que no es B (pileup, cadena cortada)
  }
  if(std::abs(c.tPDG) == 211)return 211;
  if(std::abs(c.tPDG) == 2212)return 2212;
  if(std::abs(c.tPDG) == 13 && c.tOrigin == 2)return 13;
  if(std::abs(c.tPDG) == 11 || c.tPDG == 22)return 11;
  return CATOTHER;
}

//********************************************************************
int secondaryKaonProductionXSAnalysis::MeasurementCategory(const muCandUtils::Candidate& c, int measurement){
//********************************************************************

  const int code = CandidateCategory(c);
  if(c.tPDG == 0 || !c.tSignalMuon)return code;
  return muCandUtils::IsSignalOf(measurement, true, c.tKOrigin) ? 1 : 2;
}

//********************************************************************
void secondaryKaonProductionXSAnalysis::FillCategories(){
//********************************************************************

  // primario de la beam slice
  if(_FillBeamParticleInfo && box().MainTrack){
    anaUtils::FillCategories(&GetEvent(), box().MainTrack,""); // method in highland/src/highland2/highlandUtils
    kaonAnaUtils::FillBeamParticleOriginCategory(box().MainTrack);
    kaonAnaUtils::FillBeamParticleReducedCategory(box().MainTrack);
    if(_FillBeamParticleDaughtersInfo){
      int ndau = std::min(20,(int)box().MainTrack->Daughters.size());
      for(int i = 0; i < ndau; i++)
        anaUtils::FillObjectCategories(&GetEvent(), static_cast<AnaParticleB*>(box().MainTrack->Daughters[i]),"dau",1);
    }
  }

  // candidatos, en el mismo orden que los vectores cand_* (SavedCandidates). Las categorías por objeto se van
  // añadiendo en orden: la posición k de candmu corresponde a la posición k de cand_*
  const ToyBoxMuCand& b = box();
  for(int i : SavedCandidates()){
    AnaParticleB* p = const_cast<AnaParticlePD*>(b.Candidates[i].part);
    anaUtils::FillObjectCategories(&GetEvent(), p, "cand", 1);
    const int code = CandidateCategory(b.Candidates[i]);
    anaUtils::_categ->SetObjectCode("candmu", code, code == CATNOTRUTH ? CATNOTRUTH : CATOTHER, -1);
    const char* sigName[3] = {"candsigC", "candsigB", "candsigA"};
    for(int m = 0; m < 3; m++){
      const int cm = MeasurementCategory(b.Candidates[i], m);
      anaUtils::_categ->SetObjectCode(sigName[m], cm, cm == CATNOTRUTH ? CATNOTRUTH : CATOTHER, -1);
    }
  }

  // instrumentación del haz
  if(_FillBeamInstrumentationInfo){
    AnaParticleB* beamPart = static_cast<AnaBeam*>(GetSpill().Beam)->BeamParticle;
    if (beamPart)
      anaUtils::FillCategories(&GetEvent(), beamPart, "beam"); // method in highland/src/highland2/highlandUtils
  }
}

//********************************************************************
bool secondaryKaonProductionXSAnalysis::FinalizeConfiguration(){
//********************************************************************

  // This is called before FillMicroTrees(). Para el árbol truth: accum_level[rama] de cada mu+ de K+ en reposo =
  // 1 + nivel del mejor candidato emparejado con él en esa rama (0 si ninguno; un candidato de la otra población
  // cuenta como nivel 0 en esta rama). Solo la configuración por defecto.
  if(conf().GetCurrentConfigurationIndex() != ConfigurationManager::default_conf)
    return true;

  // best[vista][rama]: vista 0 la normal (la de la selección), 1 la completa (en los eventos normales, la misma)
  const ToyBoxMuCand& b = box();
  _fullLevel.clear();
  for(AnaTrueParticleB* t : GetSpill().TrueParticles){
    if(!IsTrueMuFromKaonAtRest(t))continue;
    int bestv[2][2] = {{-1, -1}, {-1, -1}};
    for(size_t i = 0; i < b.Candidates.size(); i++){
      const AnaParticlePD* p = b.Candidates[i].part;
      if(!p || !p->TrueObject || p->TrueObject->ID != t->ID)continue;
      const int pop = b.Population[i];
      for(int v = 0; v < 2; v++){
        if(b.View[i] != v && !(v == 1 && !b.FullDump))continue;
        bestv[v][pop] = std::max(bestv[v][pop], b.Level[i]);
        bestv[v][1 - pop] = std::max(bestv[v][1 - pop], 0);
      }
    }
    const int* best = bestv[0];
    for(int ib = 0; ib < 2; ib++)_fullLevel[t->ID][ib] = bestv[1][ib] >= 0 ? 1 + bestv[1][ib] : 0;
    if(t->AccumLevel.size() == 0)t->AccumLevel.resize(sel().GetNEnabledSelections());
    for(std::vector<SelectionBase*>::iterator it = sel().GetSelections().begin(); it != sel().GetSelections().end(); it++){
      if(!(*it)->IsEnabled())continue;
      const Int_t isel = (*it)->GetEnabledIndex();
      if(t->AccumLevel[isel].size() == 0)t->AccumLevel[isel].resize((*it)->GetNBranches());
      for(UInt_t ib = 0; ib < (*it)->GetNBranches() && ib < 2; ib++)
        t->AccumLevel[isel][ib] = best[ib] >= 0 ? 1 + best[ib] : 0;
    }
  }
  return true;
}
