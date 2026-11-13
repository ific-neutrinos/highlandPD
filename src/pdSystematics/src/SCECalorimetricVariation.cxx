#include "SCECalorimetricVariation.hxx"
#include "pdServiceManager.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
SCECalorimetricVariation::SCECalorimetricVariation():SCEVariationBase(){
//********************************************************************

  // Shared calorimetry; the toy SCE scaling is passed in the CaloContext
  _cal = &pd::services().Calo();
}

//********************************************************************
SCECalorimetricVariation::~SCECalorimetricVariation(){
//********************************************************************

}

//********************************************************************
void SCECalorimetricVariation::Apply(const ToyExperiment& toy, AnaEventC& event){
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
bool SCECalorimetricVariation::UndoSystematic(AnaEventC& event){
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

//**************************************************
bool SCECalorimetricVariation::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const{
//**************************************************

  (void)event;
  
  //the systematic box only includes candidates and daughters, and all of them are relevant
  return true;
}
