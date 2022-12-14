#include "cosmicsAnalysis.hxx"
#include "Parameters.hxx"
#include "CategoriesUtils.hxx"
#include "BasicUtils.hxx"
#include "pdAnalysisUtils.hxx"

#include "pdBaseSelection.hxx"

#include "pdDataClasses.hxx"
#include "HighlandMiniTreeConverter.hxx"
#include "PDSPAnalyzerTreeConverter.hxx"

//********************************************************************
cosmicsAnalysis::cosmicsAnalysis(AnalysisAlgorithm* ana) : baseAnalysis(ana) {
//********************************************************************

  // Add the package version (not yet properly done)
  //  ND::versioning().AddPackage("cosmicsAnalysis", anaUtils::GetSoftwareVersionFromPath((std::string)getenv("PDANALYSISROOT")));
}

//********************************************************************
bool cosmicsAnalysis::Initialize(){
//********************************************************************

  // Initialize the baseAnalysis (highland/src/highland2/baseAnalysis)
  if(!baseAnalysis::Initialize()) return false;

  // Minimum accum cut level (how many cuts should be passed) to save event into the output tree
  SetMinAccumCutLevelToSave(ND::params().GetParameterI("pdControlSamples.MinAccumLevelToSave"));

  _cutLength = ND::params().GetParameterD("pdControlSamples.Cosmics.Cuts.Length");
  _cutZmin   = ND::params().GetParameterD("pdControlSamples.Cosmics.Cuts.Zmin");
  _cutZmax   = ND::params().GetParameterD("pdControlSamples.Cosmics.Cuts.Zmax");
  
  // Define standard categories for color drawing
  anaUtils::AddStandardObjectCategories("cosmic",standardPDTree::ntracks,"ntracks",1);  // This is for all the cosmics

  return true;
}

//********************************************************************
void cosmicsAnalysis::DefineInputConverters(){
//********************************************************************
  
  // add a single converter (a copy of the one in highland/baseAnalysis)
  input().AddConverter("minitree",         new HighlandMiniTreeConverter("highlandana/MiniTree"));
  input().AddConverter("minitreefiltered", new HighlandMiniTreeConverter("MiniTree"));
  input().AddConverter("PDSPAnalyzerTree", new PDSPAnalyzerTreeConverter());
}

//********************************************************************
void cosmicsAnalysis::DefineSelections(){
//********************************************************************

  // Dummy selection for the moment
  sel().AddSelection("pdBaseSelection", "dummy selection", new pdBaseSelection());
}

//********************************************************************
void cosmicsAnalysis::DefineMicroTrees(bool addBase){
//********************************************************************

  // Variables from baseAnalysis (run, event, ...)   (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::DefineMicroTrees(addBase);

  // Add standard sets of variables for ProtoDUNE analysis  (those methods are in highlandPD/src/pdUtils/standardPDTree.cxx)
  standardPDTree::AddStandardVariables_AllParticlesReco(output(),cosmicsAnalysisConstants::NMAXSAVEDCOSMICS);
  standardPDTree::AddStandardVariables_AllParticlesTrue(output(),cosmicsAnalysisConstants::NMAXSAVEDCOSMICS);
}

//********************************************************************
void cosmicsAnalysis::FillMicroTrees(bool addBase){
//********************************************************************
  
  // Variables from baseAnalysis (run, event, ...)  (highland/src/highland2/baseAnalysis)
  if (addBase) baseAnalysis::FillMicroTreesBase(addBase); 

  // Fill cosmics info
  AnaParticleB** parts = GetEvent().Particles;
  Int_t ncosmics = 0;
  
  for(size_t nparts=0;nparts < GetEvent().nParticles;nparts++){
    if (ncosmics == (Int_t)cosmicsAnalysisConstants::NMAXSAVEDCOSMICS) break;
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[nparts]);
    if(part->ParentID==-1){

      //      std::cout << nparts << " " << ncosmics << " " << part->Length << " " << part->PositionStart[2] << std::endl;
      if (part->Length>_cutLength && part->PositionStart[2]>_cutZmin && part->PositionStart[2]<_cutZmax){
    
        standardPDTree::FillStandardVariables_AllParticlesReco(output(), part);
        standardPDTree::FillStandardVariables_AllParticlesTrue(output(), part);
        output().IncrementCounter(standardPDTree::ntracks);
        ncosmics++;
      }
    }
  }

}

//********************************************************************
void cosmicsAnalysis::FillTruthTree(const AnaTrueVertex& vtx){
//********************************************************************

  // Fill the common variables (highland/src/highland2/baseAnalysis)
  baseAnalysis::FillTruthTreeBase(vtx);
}

//********************************************************************
void cosmicsAnalysis::FillCategories(){
//********************************************************************

  // For all particles in the event
  AnaParticleB** parts = GetEvent().Particles;
  Int_t nparts   = 0;
  Int_t ncosmics = 0;
  
  while(ncosmics < (Int_t)cosmicsAnalysisConstants::NMAXSAVEDCOSMICS && nparts < GetEvent().nParticles){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(parts[nparts]);
    if(part->ParentID == -1){
      anaUtils::FillObjectCategories(&GetEvent(), static_cast<AnaParticleB*>(part), "cosmic", 1);
      ncosmics++;
    }
    nparts++;
  }
}
