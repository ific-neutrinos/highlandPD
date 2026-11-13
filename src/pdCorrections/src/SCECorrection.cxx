#include "SCECorrection.hxx"
#include "pdDataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>


//********************************************************************
SCECorrection::SCECorrection(){
//********************************************************************

  _sce = new SpaceCharge();
  _sce->Initialize();
  _cal = new Calorimetry();
  _cal->Initialize(_sce);
  
}

//********************************************************************
void SCECorrection::Apply(AnaSpillC& spillC){
//********************************************************************

  //cast bunch
  AnaSpill& spill = *static_cast<AnaSpill*>(&spillC);
  AnaBunch* bunch = static_cast<AnaBunch*>(spill.Bunches[0]);

  // Loop over particles
  for(UInt_t ipart = 0; ipart < bunch->Particles.size(); ipart++){
    
    AnaParticlePD* part = static_cast<AnaParticlePD*>(bunch->Particles[ipart]);
    
    // The un-corrected particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    
    if (!original) continue; //?

    //here position start and end should be modified but tpcid is missing. It could just be estimated. Im lazy right now to do it.
    //I am really sorry but I don't now if the statement above is trueor false, or if it has any relevance in the analysis. 
    //I'm gonna ignore it (Migue, 27/03/2025)

    for(int ihit = 0; ihit < part->Hits[2].size(); ihit++){
      _sce->ApplyPositionCorrection(part->Hits[2][ihit]);
      _cal->ApplySCECorrection(part->Hits[2][ihit]);
    }
  }
}

