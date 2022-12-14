#include "BeamCompositionWeight.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>
#include "EventBoxPD.hxx"

//********************************************************************
BeamCompositionWeight::BeamCompositionWeight():EventWeightBase(), BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","beamComp", BinnedParams::k2D_SYMMETRIC){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());
}

//********************************************************************
Weight_h BeamCompositionWeight::ComputeWeight(const ToyExperiment& toy, const AnaEventC& event, const ToyBoxB& box){
//********************************************************************

  (void)box;

  // Initialy the weight is 1
  Weight_h eventWeight = 1;

  // Get the true beam particle
  AnaTrueParticlePD* trueBeamPart = pdAnaUtils::GetTrueBeamParticle(event);
  if (!trueBeamPart) return eventWeight;
  
  // We need the errors part of the data file but as well the relative uncertainty for sigma
  Float_t weight_mean, weight_error;
  Int_t weight_index;
    
  // Get the systematics parameters for this particle type and momentum
  if (!GetBinValues(abs(trueBeamPart->PDG), trueBeamPart->Momentum, weight_mean,  weight_error,  weight_index))  return eventWeight;

  // compute the weight
  eventWeight.Systematic = weight_mean +  weight_error*toy.GetToyVariations(_index)->Variations[weight_index];
  eventWeight.Correction = weight_mean;

  return eventWeight;
}
  
