#include "LifetimeCorrection.hxx"
#include "pdServiceManager.hxx"
#include "pdDataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>


//********************************************************************
LifetimeCorrection::LifetimeCorrection(){
//********************************************************************

  //always the nominal (MC) lifetime, it does not call GetLifetime(event)
  _cal = &pd::services().Calo();
}

//********************************************************************
void LifetimeCorrection::Apply(AnaSpillC& spillC){
//********************************************************************

  //cast bunch
  AnaSpill& spill = *static_cast<AnaSpill*>(&spillC);
  // spills without bunch (bad DAQ) have nothing to correct
  if(spill.Bunches.empty())return;
  AnaBunch* bunch = static_cast<AnaBunch*>(spill.Bunches[0]);

  // Loop over particles
  for(UInt_t ipart = 0; ipart < bunch->Particles.size(); ipart++){
    
    AnaParticlePD* part = static_cast<AnaParticlePD*>(bunch->Particles[ipart]);
    
    // The un-corrected particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    
    if (!original) continue; //?

    _cal->ApplyLifetimeCorrection(part);
  }
}

