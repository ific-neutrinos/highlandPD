#include "dEdxVariation.hxx"
#include "pdDataFiles.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
dEdxVariation::dEdxVariation():EventVariationBase(),BinnedParams(pd::DataDir("pdSystematics"),"dEdx", BinnedParams::k1D_SYMMETRIC_NOMEAN){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());

 
}

//********************************************************************
void dEdxVariation::Apply(const ToyExperiment& toy, AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event. The SystBox contains vector of objects that are relevant for this systematic.
  // In this way the loop over objects is restricted to the ones that really matter
  SystBoxB* box = GetSystBox(event);

  // A single bin: the same variation applies to every hit of every particle
  // in this toy, so read it once
  BinnedParamsParams params;
  int index;
  if(!GetBinValues(0.5, params, index)){
    std::cerr << "dEdxVariation::Apply: cannot read the dE/dx uncertainty from dEdx.dat. Aborting." << std::endl;
    exit(1);
  }
  const double factor = 1 + params.sigma*toy.GetToyVariations(_index)->Variations[index];

  // Loop over all relevant tracks for this variation
  for(Int_t ipart = 0; ipart < box->nRelevantRecObjects; ipart++){

    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[ipart]);
    if(part->Hits[2].empty())continue;

    // The un-varied particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if(!original)continue;
    //loop over hits
    for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++)
      part->Hits[2][ihit].dEdx = part->Hits[2][ihit].dEdx * factor;
    
    //recompute derived quantities
    std::pair<double, int> chi2pidprot = pdAnaUtils::Chi2PID(*part, 2212);
    part->Chi2Proton  = chi2pidprot.first;
    std::pair<double, int> chi2pidmuon = pdAnaUtils::Chi2PID(*part, 13);
    part->Chi2Muon    = chi2pidmuon.first;
    part->Chi2ndf     = chi2pidmuon.second;
  }
}

//********************************************************************
bool dEdxVariation::UndoSystematic(AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event
  SystBoxB* box = GetSystBox(event);

  for(Int_t ipart = 0; ipart < box->nRelevantRecObjects; ipart++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[ipart]);
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if(!original)continue;

    //loop over hits
    for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++)
      part->Hits[2][ihit].dEdx = original->Hits[2][ihit].dEdx;
    

    // reset the top level quantities
    part->Chi2Proton = original->Chi2Proton;
    part->Chi2ndf    = original->Chi2ndf;
    part->Chi2Muon   = original->Chi2Muon;
  }
  
  // Don't reset the spill to corrected
  return false;
}

//**************************************************
bool dEdxVariation::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const{
//**************************************************

  (void)event;
  
  //the systematic box only includes candidates and daughters, and all of them are relevant
  return true;
}
