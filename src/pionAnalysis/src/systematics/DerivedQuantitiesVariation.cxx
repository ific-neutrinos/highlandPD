#include "DerivedQuantitiesVariation.hxx"
#include "DataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>
#include "EventBoxPD.hxx"

//********************************************************************
DerivedQuantitiesVariation::DerivedQuantitiesVariation(): EventVariationBase(){
//********************************************************************

}

//********************************************************************
void DerivedQuantitiesVariation::Apply(const ToyExperiment& toy, AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event. The SystBox contains vector of objects that are relevant for this systematic.
  // In this way the loop over objects is restricted to the ones that really matter
  SystBoxB* box = GetSystBox(event);


  double calib_factor[3] = {4.81e-3, 4.81e-3, 4.57e-3}; //
  //  double norm_factor[3] = {1.0078, 1.0082, 0.9947};
  double norm_factor[3] = {1.0078, 1.0082, 0.9946};
  double fEventT0=500;
  
  // Loop over all relevant tracks for this variation
  for (Int_t itrk = 0; itrk < box->nRelevantRecObjects; itrk++){

    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[itrk]);

    // The un-corrected particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if (!original)         continue;
        
    //    Fill the dedx_vector for the collection plane
    std::vector<double> dedx_vector;
    for (UInt_t j=0;j<part->Hits[2].size();j++){
      AnaHitPD& hit = part->Hits[2][j];

      // compute the dEdx from the varied dQdx
      hit.dEdx_calib = hit.dEdx = fCalo.dEdx_from_dQdx(hit.dQdx);

      dedx_vector.push_back(hit.dEdx_calib);
    }

    // Recompute the derived quantities. 
    std::pair<double, int> chi2pid = pdAnaUtils::Chi2PID(*part, NULL);
    part->Chi2Proton  = chi2pid.first;
    part->Chi2ndf     = chi2pid.second;
    part->truncLibo_dEdx = pdAnaUtils::ComputeTruncatedMean(0.16,0.16,dedx_vector);
  }
}

//********************************************************************
bool DerivedQuantitiesVariation::UndoSystematic(AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event
  SystBoxB* box = GetSystBox(event);

  for (Int_t itrk=0;itrk<box->nRelevantRecObjects;itrk++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[itrk]);
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if (!original)   continue;

    part->Chi2Proton      = original->Chi2Proton;
    part->Chi2ndf         = original->Chi2ndf;
    part->truncLibo_dEdx  = original->truncLibo_dEdx;
  }

  // Don't reset the spill to corrected
  return false;
}

//**************************************************
bool DerivedQuantitiesVariation::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& track) const{
  //**************************************************

  (void)event;
  (void)track;

  return true;
}

//********************************************************************
Int_t DerivedQuantitiesVariation::GetRelevantRecObjectGroups(const SelectionBase& sel, Int_t* IDs) const{
  //********************************************************************

  Int_t ngroups=0;
  for (UInt_t b=0; b<sel.GetNBranches(); b++){
    IDs[ngroups++] = EventBoxPD::kCandidateAndDaughters;
  }
  return ngroups;
}

