#include "SCEVariation.hxx"
#include "pdServiceManager.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
SCEVariation::SCEVariation():EventVariationBase(),BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","SCE", BinnedParams::k1D_SYMMETRIC_NOMEAN){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());

  //check we have enough bins to run this systematic, otherwise exit
  if(GetNBins() > (int)NMAXSYSTSOURCEBINS){
    std::cout << "SCEVariation::SCEVariation()." << std::endl;
    std::cout << "the source files has more bins than allowed for the systematic propagation!" << std::endl;
    std::cout << "N MAX BINS allowed " << NMAXSYSTSOURCEBINS << std::endl;
    std::cout << "N BINS on source   " << GetNBins() << std::endl;
    std::cout << "if you need more bins, look at psyckeUtils/src/BinnedParams.hxx" << std::endl;
    std::exit(1);
  }
  
  //Toys use the nominal map scaled on the fly, there are no per-toy copies
  _sce = &pd::services().SCE();

  // Shared calorimetry; the toy SCE scaling is passed in the CaloContext
  _cal = &pd::services().Calo();
}

//********************************************************************
SCEVariation::~SCEVariation(){
//********************************************************************

}

//********************************************************************
void SCEVariation::Apply(const ToyExperiment& toy, AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event. The SystBox contains vector of objects that are relevant for this systematic.
  // In this way the loop over objects is restricted to the ones that really matter
  SystBoxB* box = GetSystBox(event);
  const AnaEventPD* eventpd = static_cast<const AnaEventPD*>(&event);
  
  //scaling of the nominal SCE map for this toy
  SCEScaling scale = GetToyScaling(toy);

  //calorimetry with the toy SCE scaling and the appropiate lifetime depending on mc/data and run
  CaloContext ctx = _cal->NominalContext();
  ctx.SCEScale = scale;
  ctx.Lifetime = _cal->GetLifetime(*eventpd);
  
  // Loop over all relevant tracks for this variation
  for(Int_t ipart = 0; ipart < box->nRelevantRecObjects; ipart++){
     AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[ipart]);
  
    // The un-varied particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if(!original)continue;

    if(part->Hits[2].empty())continue;
    
    //important, two things are varied: XYZ position of the hit and pitch
    //loop here to avoid looping more than once
    for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++){
      _sce->ApplyPositionCorrection(part->Hits[2][ihit], scale);
      _cal->ApplySCECorrection(part->Hits[2][ihit], ctx);
      //_cal->ApplyLifetimeCorrection(part->Hits[2][ihit]);
    }
  }
}

//********************************************************************
bool SCEVariation::UndoSystematic(AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event
  SystBoxB* box = GetSystBox(event);

  //_sce->ResetToNominal();
  
  for(Int_t ipart = 0; ipart < box->nRelevantRecObjects; ipart++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[ipart]);
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if(!original)continue;

    //loop over hits
    for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++){
      part->Hits[2][ihit].dQdx = original->Hits[2][ihit].dQdx;
      part->Hits[2][ihit].dQdx_elife = original->Hits[2][ihit].dQdx_elife;
      part->Hits[2][ihit].Position.SetX(original->Hits[2][ihit].Position.X());
      part->Hits[2][ihit].Position.SetY(original->Hits[2][ihit].Position.Y());
      part->Hits[2][ihit].Position.SetZ(original->Hits[2][ihit].Position.Z());
    }
  }
  
  // Don't reset the spill to corrected
  return false;
}

//********************************************************************
SCEScaling SCEVariation::GetToyScaling(const ToyExperiment& toy){
//********************************************************************

  Float_t sigma;
  GetSigmaValueForBin(0, sigma); //only 1 bin

  SCEScaling scale;
  scale.Global = 1 + sigma*toy.GetToyVariations(_index)->Variations[0];
  return scale;
}

//**************************************************
bool SCEVariation::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const{
//**************************************************

  (void)event;
  
  //the systematic box only includes candidates and daughters, and all of them are relevant
  return true;
}
