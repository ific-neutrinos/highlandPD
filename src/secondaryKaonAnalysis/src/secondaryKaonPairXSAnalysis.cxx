#include "secondaryKaonPairXSAnalysis.hxx"
#include "Parameters.hxx"
#include "CategoriesUtils.hxx"
#include "BasicUtils.hxx"
#include "pdAnalysisUtils.hxx"
#include "baseToyMaker.hxx"
#include "pdDataClasses.hxx"

#include "secondaryKaonPairXSSelection.hxx"
#include "muCandidateUtils.hxx"
#include "secondaryKaonSelectionNoBranches.hxx"
#include "kaonAnalysisUtils.hxx"

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
secondaryKaonPairXSAnalysis::secondaryKaonPairXSAnalysis(AnalysisAlgorithm* ana) : baseAnalysis(ana) {
//********************************************************************

}

//********************************************************************
bool secondaryKaonPairXSAnalysis::Initialize(){
//********************************************************************

  // Initialize the baseAnalysis (highland/src/highland2/baseAnalysis)
  if(!baseAnalysis::Initialize()) return false;

  // Minimum accum cut level (how many cuts should be passed) to save event into the output tree
  SetMinAccumCutLevelToSave(ND::params().GetParameterI("secondaryKaonAnalysis.pairXSAnalysis.MinAccumLevelToSave"));

  _FillBeamInstrumentationInfo = ND::params().GetParameterI("secondaryKaonAnalysis.pairXSAnalysis.FillBeamInstrumentationInfo");
  _FillBeamParticleInfo = ND::params().GetParameterI("secondaryKaonAnalysis.pairXSAnalysis.FillBeamParticleInfo");
  _FillBeamParticleHitsInfo = ND::params().GetParameterI("secondaryKaonAnalysis.pairXSAnalysis.FillBeamParticleHitsInfo");
  _FillBeamParticleDaughtersInfo = ND::params().GetParameterI("secondaryKaonAnalysis.pairXSAnalysis.FillBeamParticleDaughtersInfo");
  _FillKaonCandidatesInfo = ND::params().GetParameterI("secondaryKaonAnalysis.pairXSAnalysis.FillKaonCandidatesInfo");
  _FillKaonCandidatesHitInfo = ND::params().GetParameterI("secondaryKaonAnalysis.pairXSAnalysis.FillKaonCandidatesHitInfo");
  _FillKaonCandidatesPIDInfo = ND::params().GetParameterI("secondaryKaonAnalysis.pairXSAnalysis.FillKaonCandidatePIDInfo");
  _FillBestKaonCandidateInfo = ND::params().GetParameterI("secondaryKaonAnalysis.pairXSAnalysis.FillBestKaonCandidateInfo");
  _FillToyVariablesInfo = ND::params().GetParameterI("secondaryKaonAnalysis.pairXSAnalysis.FillToyVariablesInfo");

  _ApplydQdxSystematic = ND::params().GetParameterI("secondaryKaonAnalysis.pairXSAnalysis.ApplydQdxSystematic");
  _ApplyRecombinationSystematic = ND::params().GetParameterI("secondaryKaonAnalysis.pairXSAnalysis.ApplyRecombinationSystematic");
  _ApplySCESystematic = ND::params().GetParameterI("secondaryKaonAnalysis.pairXSAnalysis.ApplySCESystematic");
  _ApplyBrokenTracksSystematic = ND::params().GetParameterI("secondaryKaonAnalysis.pairXSAnalysis.ApplyBrokenTracksSystematic");
  _ApplyBeamPIDEfficiencySystematic = ND::params().GetParameterI("secondaryKaonAnalysis.pairXSAnalysis.ApplyBeamPIDEfficiencySystematic");
  _ApplyBeamPartWeightSystematic = ND::params().GetParameterI("secondaryKaonAnalysis.pairXSAnalysis.ApplyBeamPartWeightSystematic");
  _ApplyBeamMomWeightSystematic = ND::params().GetParameterI("secondaryKaonAnalysis.pairXSAnalysis.ApplyBeamMomWeightSystematic");
  _ApplyProtonBackgroundWeightSystematic = ND::params().GetParameterI("secondaryKaonAnalysis.pairXSAnalysis.ApplyProtonBackgroundWeightSystematic");

  _ApplySCECorrection = ND::params().GetParameterI("secondaryKaonAnalysis.pairXSAnalysis.ApplySCECorrection");
  _ApplydEdxCorrection = ND::params().GetParameterI("secondaryKaonAnalysis.pairXSAnalysis.ApplydEdxCorrection");
  _ApplyRecombinationCorrection = ND::params().GetParameterI("secondaryKaonAnalysis.pairXSAnalysis.ApplyRecombinationCorrection");
  
  // Define standard categories for color drawing
  if(_FillBeamParticleInfo)
    anaUtils::AddStandardCategories();        // This is for the candidate particle
  if(_FillBeamInstrumentationInfo)
    anaUtils::AddStandardCategories("beam");  // This is for the Beam Instrumentation particle
  if(_FillBestKaonCandidateInfo){
    anaUtils::AddStandardCategories("bestcandidate");  // This is for the best candidate of each event
    anaUtils::AddStandardCategories("bestcandidatedau");  // This is for the best candidate of each event
  }
  
  if(_FillBeamParticleDaughtersInfo)
    anaUtils::AddStandardObjectCategories("dau"   ,standardPDTree::seltrk_ndau ,"seltrk_ndau",1);

  // Add standard categories for the candidates
  if(_FillKaonCandidatesInfo){
    anaUtils::AddStandardObjectCategories("candidate"   ,kaonTree::ncandidates,"ncandidates",1);
    anaUtils::AddStandardObjectCategories("candidatedau",kaonTree::ncandidates,"ncandidates",1);
    // Add our own categories (in secondaryKaonAnalysisUtils.cxx)
  }
  kaonAnaUtils::AddCustomCategories();//to be checked

  return true;
}

//********************************************************************
void secondaryKaonPairXSAnalysis::DefineInputConverters(){
//********************************************************************

  // add a single converter (a copy of the one in highland/baseAnalysis)
  input().AddConverter("minitree",         new HighlandMiniTreeConverter("highlandana/MiniTree"));
  input().AddConverter("minitreefiltered", new HighlandMiniTreeConverter("MiniTree"));
  input().AddConverter("PDSPAnalyzerTree", new PDSPAnalyzerTreeConverter());
}

//********************************************************************
void secondaryKaonPairXSAnalysis::DefineSelections(){
//********************************************************************

  sel().AddSelection("XS", "kaon XS selection", new secondaryKaonPairXSSelection(false));
}

//********************************************************************
void secondaryKaonPairXSAnalysis::DefineCorrections(){
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
void secondaryKaonPairXSAnalysis::DefineSystematics(){
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
void secondaryKaonPairXSAnalysis::DefineConfigurations(){
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
void secondaryKaonPairXSAnalysis::DefineMicroTrees(bool addBase){
//********************************************************************

  // -------- Add variables to the analysis tree ----------------------

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

  // -------- Add candidates variables ----------------------
  if(_FillKaonCandidatesInfo){
    kaonTree::AddKaonVariables_KaonCandidatesReco    (output(),secondaryKaonPairXSAnalysisConstants::NMAXSAVEDCANDIDATES);
    kaonTree::AddKaonVariables_KaonCandidatesTrue    (output(),secondaryKaonPairXSAnalysisConstants::NMAXSAVEDCANDIDATES);
    if(_FillKaonCandidatesHitInfo)
      kaonTree::AddKaonVariables_KaonCandidatesHitsReco(output(),secondaryKaonPairXSAnalysisConstants::NMAXSAVEDCANDIDATES);
    if(_FillKaonCandidatesPIDInfo)
      kaonTree::AddKaonVariables_KaonCandidatesRecoPID (output(),secondaryKaonPairXSAnalysisConstants::NMAXSAVEDCANDIDATES);
    // el contador se pasa por su nombre corto (las macros AddVar* lo convierten en texto y tiene que coincidir con el de kaonTree)
    using kaonTree::ncandidates;
    AddVarMaxSizeVI(output(), candidatedau_truesignalmuon, "hija del candidato: signal muon (mu+ de K+ -> mu+ nu de una K+ en reposo)", ncandidates, secondaryKaonPairXSAnalysisConstants::NMAXSAVEDCANDIDATES);
    AddVarMaxSizeVI(output(), candidatedau_truekorigin,   "hija del candidato: origen de su K (muCandUtils::KOrigin; -1 no es mu+ de K+)", ncandidates, secondaryKaonPairXSAnalysisConstants::NMAXSAVEDCANDIDATES);
  }

  // -------- Add best candidate variables ----------------------
  if(_FillBestKaonCandidateInfo){
    kaonTree::AddKaonVariables_KaonBestCandidateReco    (output());
    kaonTree::AddKaonVariables_KaonBestCandidateHitsReco(output());
    kaonTree::AddKaonVariables_KaonBestCandidateTrue    (output());
  }
}

//********************************************************************
void secondaryKaonPairXSAnalysis::DefineTruthTree(){
//********************************************************************

  // Variables from baseAnalysis (run, event, ...)   (highland/src/highland2/baseAnalysis)
  baseAnalysis::DefineTruthTree();

  // Add standard sets of variables for ProtoDUNE analysis (those methods are in highlandPD/src/pdUtils/standardPDTree.cxx)
  standardPDTree::AddStandardVariables_BeamInstrumentationTrue(output());
  
  // Add specific variables for the kaon candidates
  kaonTree::AddKaonVariables_TrueKaonCandidates(output());
}

//********************************************************************
void secondaryKaonPairXSAnalysis::FillMicroTrees(bool addBase){
//********************************************************************
  
  // Variables from baseAnalysis (run, event, ...)  (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::FillMicroTreesBase(addBase); 
  
  // Fill standard variables for the PD analysis (those methods are in highlandPD/src/pdUtils/standardPDTree.cxx)
  standardPDTree::FillStandardVariables_EventInfo(        output(), &GetEvent());
  if(_FillBeamInstrumentationInfo){
    standardPDTree::FillStandardVariables_BeamInstrumentationReco(         output(), GetSpill().Beam);
    standardPDTree::FillStandardVariables_BeamInstrumentationTrue(         output(), GetSpill().Beam);
  }
  if(_FillBeamParticleInfo){
    standardPDTree::FillStandardVariables_BeamParticleReco(    output(), box().MainTrack);
    standardPDTree::FillStandardVariables_BeamParticleTrue(    output(), box().MainTrack);    
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
  
  // ---------- kaon candidates variables --------------
  if(_FillKaonCandidatesInfo){
    std::map<int, AnaTrueParticleB*> byID;
    for(AnaTrueParticleB* t : GetSpill().TrueParticles)byID[t->ID] = t;
    const AnaBeamPD* beam = static_cast<const AnaBeamPD*>(GetSpill().Beam);
    const AnaTrueParticleB* beamTrue = (beam && beam->BeamParticle) ? static_cast<AnaTrueParticleB*>(beam->BeamParticle->TrueObject) : nullptr;
    if(box().Candidates.size()>0){
      for(int i = 0; i < std::min((int)box().Candidates.size(),(int)secondaryKaonAnalysisConstants::NMAXSAVEDCANDIDATES); i++){
        { const AnaParticlePD* dau = box().Candidates[i]->Daughters.empty() ? nullptr : static_cast<AnaParticlePD*>(box().Candidates[i]->Daughters[0]);
          const AnaTrueParticleB* td = dau ? static_cast<const AnaTrueParticleB*>(dau->TrueObject) : nullptr;
          output().FillVectorVar(candidatedau_truesignalmuon, (Int_t)(td ? muCandUtils::IsSignalMuon(td, byID) : 0));
          output().FillVectorVar(candidatedau_truekorigin,    (Int_t)(td ? muCandUtils::KOrigin(td, beamTrue, byID) : -1)); }
	AnaParticlePD* parent = static_cast<AnaParticlePD*>(anaUtils::GetParticleByID(GetBunch(), box().Candidates[i]->ParentID));
	kaonTree::FillKaonVariables_KaonCandidatesReco    (output(), box().Candidates[i] ,parent);
	kaonTree::FillKaonVariables_KaonCandidatesTrue    (output(), box().Candidates[i]);
	if(_FillKaonCandidatesHitInfo)
	  kaonTree::FillKaonVariables_KaonCandidatesHitsReco(output(), box().Candidates[i]);
	if(_FillKaonCandidatesPIDInfo)
	  kaonTree::FillKaonVariables_KaonCandidatesRecoPID (output(), box().Candidates[i]);
	output().IncrementCounter(kaonTree::ncandidates);
      } 
    }
  }

  // ---------- best kaon candidate variables --------------
  if(_FillBestKaonCandidateInfo){
    //get the kaon selection and get the branch with a larger AccumLevel
    SelectionBase* ksel = sel().GetSelection("XS"); //todo fix this
    int branchmax = 0;
    int max       = 0;  
    for(UInt_t ibranch = 0; ibranch < ksel->GetNBranches(); ibranch++){
      if(ksel->GetAccumCutLevel(ibranch) > max){
	max = ksel->GetAccumCutLevel(ibranch);
	branchmax = ibranch;
      }
    }

    if(box().Candidates.size()>0){
      AnaParticlePD* parent = static_cast<AnaParticlePD*>(anaUtils::GetParticleByID(GetBunch(), box().Candidates[branchmax]->ParentID));
      kaonTree::FillKaonVariables_KaonBestCandidateReco    (output(), box().Candidates[branchmax], parent);
      kaonTree::FillKaonVariables_KaonBestCandidateHitsReco(output(), box().Candidates[branchmax]);
      kaonTree::FillKaonVariables_KaonBestCandidateTrue    (output(), box().Candidates[branchmax], parent);
    }
  }
}

//********************************************************************
void secondaryKaonPairXSAnalysis::FillToyVarsInMicroTrees(bool addBase){
//********************************************************************

  // Fill the common variables  (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::FillToyVarsInMicroTreesBase(addBase);
}

//********************************************************************
void secondaryKaonPairXSAnalysis::FillTruthTree(){
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
}

//********************************************************************
bool secondaryKaonPairXSAnalysis::CheckFillTruthTree(const AnaTrueParticlePD* part){
//********************************************************************

  if(!part)return false;
  else if(abs(part->PDG)==321){
    return true;
  }
  else return false;
}

//********************************************************************
void secondaryKaonPairXSAnalysis::FillTruthTree(const AnaTrueParticlePD& part){
//********************************************************************

  // Fill the common variables (highland/src/highland2/baseAnalysis)

  // Event variables
  output().FillVar(run,    GetSpill().EventInfo->Run);
  output().FillVar(subrun, GetSpill().EventInfo->SubRun);
  output().FillVar(evt,    GetSpill().EventInfo->Event);
  
  // Fill standard variables for the PD analysis (those methods are in highlandPD/src/pdUtils/standardPDTree.cxx)
  standardPDTree::FillStandardVariables_BeamInstrumentationTrue(output(), GetSpill().Beam);

  // Fill true kaons candidate info
  kaonTree::FillKaonVariables_TrueKaonCandidates(output(), &part, NULL);
}

//********************************************************************
void secondaryKaonPairXSAnalysis::FillCategories(){
//********************************************************************

  // For the beam track
  if(_FillBeamParticleInfo){
    if(box().MainTrack){
      anaUtils::FillCategories(&GetEvent(), box().MainTrack,""); // method in highland/src/highland2/highlandUtils
      kaonAnaUtils::FillBeamParticleOriginCategory(box().MainTrack); // method in highland/src/highland2/highlandUtils
      kaonAnaUtils::FillBeamParticleReducedCategory(box().MainTrack);
      if(_FillBeamParticleDaughtersInfo){
	int ndau = std::min(20,(int)box().MainTrack->Daughters.size());
	for(int i = 0; i < ndau; i++)
	  anaUtils::FillObjectCategories(&GetEvent(), static_cast<AnaParticleB*>(box().MainTrack->Daughters[i]),"dau",1);
      }
    }
  }

  // for the candidates
  if(_FillKaonCandidatesInfo){
    if(box().Candidates.size()>0){
      for(int i = 0; i < (int)box().Candidates.size(); i++){
	anaUtils::FillObjectCategories(&GetEvent(), static_cast<AnaParticleB*>(box().Candidates[i]),              "candidate",   1);
	SetIsSignal(box().Candidates[i]);
	kaonAnaUtils::FillCandidateParticleReducedCategory(box().Candidates[i]);
	kaonAnaUtils::FillCandidateDaughterParticleReducedCategory(box().Candidates[i]);
	if(!box().Candidates[i]->Daughters.empty()){
	  anaUtils::FillObjectCategories(&GetEvent(), static_cast<AnaParticleB*>(box().Candidates[i]->Daughters[0]),"candidatedau",1);
	  kaonAnaUtils::FillCandidateDaughterMuonCategory(box().Candidates[i], static_cast<AnaParticlePD*>(box().Candidates[i]->Daughters[0]));
	  kaonAnaUtils::FillCandidateDaughterMuonReducedCategory(box().Candidates[i], static_cast<AnaParticlePD*>(box().Candidates[i]->Daughters[0]));
	}
      }
    }
  }
  
  //get the kaon selection and get the branch with a larger AccumLevel
  if(_FillBestKaonCandidateInfo){
    SelectionBase* ksel = sel().GetSelection("XS");
    int branchmax = 0;
    int max       = 0;  
    for(UInt_t ibranch = 0; ibranch < ksel->GetNBranches(); ibranch++){
      if(ksel->GetAccumCutLevel(ibranch) > max){
	max = ksel->GetAccumCutLevel(ibranch);
	branchmax = ibranch;
      }
    }

    // ---------- best kaon candidate variables --------------
    if(box().Candidates.size()>0){
      anaUtils::FillCategories(&GetEvent(), static_cast<AnaParticleB*>(box().Candidates[branchmax]), "bestcandidate");
      //the signal flag of the true kaon is needed by the category (also when the candidates block is off)
      SetIsSignal(box().Candidates[branchmax]);
      kaonAnaUtils::FillBestCandidateParticleReducedCategory(static_cast<AnaParticlePD*>(box().Candidates[branchmax]));
      if(!box().Candidates[branchmax]->Daughters.empty())
	anaUtils::FillCategories(&GetEvent(), static_cast<AnaParticleB*>(box().Candidates[branchmax]->Daughters[0]), "bestcandidatedau");
    }
  }

  // For the beam
  if(_FillBeamInstrumentationInfo){
    AnaParticleB* beamPart = static_cast<AnaBeam*>(GetSpill().Beam)->BeamParticle;
    if (beamPart)
      anaUtils::FillCategories(&GetEvent(), beamPart, "beam"); // method in highland/src/highland2/highlandUtils
  }
}

//********************************************************************
bool secondaryKaonPairXSAnalysis::FinalizeConfiguration(){
//********************************************************************
  
  // This is called before FillMicroTrees()

  const Int_t CandidatesAccumLevel = 1; //TODO: think about this
  
  // Save the accum level for the true particle associated to the recon one such that efficiencies can be computed from the truth tree
  if(conf().GetCurrentConfigurationIndex() != ConfigurationManager::default_conf)
    return true;

  std::vector<AnaTrueParticleB*> trueParts = GetSpill().TrueParticles;
  // Loop over all true particles
  for(std::vector<AnaTrueParticleB*>::iterator itp = trueParts.begin(); itp!=trueParts.end(); itp++) {
    AnaTrueParticleB* part = *itp;
    if(!part) continue;
    
    // When the AccumLevel has not been already saved for this vertex 
    if(part->AccumLevel.size() == 0)
      part->AccumLevel.resize(sel().GetNEnabledSelections());
    
    for(std::vector<SelectionBase*>::iterator it = sel().GetSelections().begin(); it != sel().GetSelections().end(); it++){
      
      if(!(*it)->IsEnabled())continue;
      
      Int_t isel = (*it)->GetEnabledIndex();

      if (part->AccumLevel[isel].size() == 0)
        part->AccumLevel[isel].resize((*it)->GetNBranches());
      
      // Hotfix (23-nov-2026): esta selección tiene una sola rama. Antes se llamaba a GetAccumCutLevel(i)
      // con i = índice del candidato, que lee fuera del array con más de un candidato (accum_level del
      // árbol truth con basura). Arreglo de fondo en la Fase 4 (highland trabaja por evento, no por candidato).
      const Int_t level = (*it)->GetAccumCutLevel(0);
      bool found = false;
      int max_accum_level = -1;

      for(int i = 0; i < (int)box().Candidates.size(); i++){
	if(!box().Candidates[i]->TrueObject)continue;
	if(box().Candidates[i]->TrueObject->ID == part->ID){
	  if(level>max_accum_level) max_accum_level = level;
	  found = true;
	  //break;
	}
      }

      int min_accum_level = 100;
      for(int i = 0; i < (int)box().Candidates.size(); i++){
	if(level<min_accum_level)min_accum_level=level;
      }

      if(!found){
	part->AccumLevel[0][0] = std::min(min_accum_level, CandidatesAccumLevel);
	std::fill(part->AccumLevel[0].begin()+1, part->AccumLevel[0].end(), -1);
      }
      else{
	part->AccumLevel[0][0] = level;
	std::fill(part->AccumLevel[0].begin()+1, part->AccumLevel[0].end(), -1);
      }
    }
  }
  return true;
}

//********************************************************************
bool secondaryKaonPairXSAnalysis::IsProductionSignal(const AnaTrueParticlePD* kaon){
//********************************************************************

  // Same definition as test/phase2/kaon_origin_minitrees_impl.C
  if(!kaon || kaon->PDG != 321) return false;
  const AnaTrueParticle* beamTrue = pdAnaUtils::GetBeamTrueParticle(GetSpill());
  if(!beamTrue) return false;

  // go up through the K+ -> K+ to the first K+ of the chain
  const AnaTrueParticlePD* first = kaon;
  for(int depth = 0; depth < 1000; depth++){
    if(first->ID == beamTrue->ID) return false;  // the beam K+ itself (also after K+ -> K+): not production
    if(first->ParentID == 0) return false;      // another root (pileup, cosmic...)
    const AnaTrueParticlePD* parent = pdAnaUtils::GetTrueParticle(GetSpill().TrueParticles, first->ParentID);
    if(!parent) return false;                   // broken chain: cannot be classified
    if(parent->PDG != 321)
      return parent->ID == beamTrue->ID;        // signal if the first K+ is a direct daughter of the beam
    first = parent;
  }
  return false;
}

//********************************************************************
void secondaryKaonPairXSAnalysis::SetIsSignal(AnaParticlePD* part){
//********************************************************************

  if(!part) return;
  AnaTrueParticlePD* truepart = static_cast<AnaTrueParticlePD*>(part->TrueObject);
  if(!truepart) return;
  truepart->Signal = IsProductionSignal(truepart);
}
