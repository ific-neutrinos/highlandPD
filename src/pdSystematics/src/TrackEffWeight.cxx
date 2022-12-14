#include "TrackEffWeight.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>
#include "EventBoxPD.hxx"
#include "SystematicUtils.hxx"

//********************************************************************
TrackEffWeight::TrackEffWeight():EventWeightBase(), BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","trackEff", BinnedParams::k1D_EFF_ASSYMMETRIC){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());
}

//********************************************************************
Weight_h TrackEffWeight::ComputeWeight(const ToyExperiment& toy, const AnaEventC& event, const ToyBoxB& box){
//********************************************************************
  
  // Initialy the weight is 1
  Weight_h eventWeight = 1;

  // Get the tracking efficiency for bin 0 (there is a single bin);
  BinnedParamsParams params;
  int index;
  if(!GetBinValues(0.5, params, index))	return eventWeight;

  // Get the SystBox for this event, and the appropriate selection and branch
  SystBoxB* SystBox = GetSystBox(event,box.SelectionEnabledIndex,box.SuccessfulBranch);
   
  // Loop over all TrueParticles in the TPC
  for (Int_t itrue=0;itrue< SystBox->nRelevantTrueObjects; itrue++){      
    //    AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(SystBox->RelevantTrueObjects[itrue]);            

    // True-reco association done only once per event in EventWeightBase;
    bool found = (SystBox->RelevantTrueObjectsReco[itrue]!=NULL);

    // compute the weight
    eventWeight *= systUtils::ComputeEffLikeWeight(found, toy, GetIndex(), index, params);
  }

  return eventWeight;
}

//********************************************************************
bool TrackEffWeight::IsRelevantTrueObject(const AnaEventC& event, const AnaTrueObjectC& trueObj) const{
//********************************************************************

  (void)event;
  
  /*const AnaTrueParticleB& truePart = *static_cast<const AnaTrueParticleB*>(&trueObj);

  // consider only charged true particles with momenta below 0.4 (as an example)
  if(truePart.Charge==0) return false;
  if(truePart.Momentum<0.4) return false;*/
    
  return true;
}
