#include "BeamMomCorrection.hxx"
#include "pdDataClasses.hxx"
#include <cassert>

//#define DEBUG

//********************************************************************
BeamMomCorrection::BeamMomCorrection(){
//********************************************************************

}

//********************************************************************
void BeamMomCorrection::Apply(AnaSpillC& spillC){
//********************************************************************

  AnaSpill& spill = *static_cast<AnaSpill*>(&spillC);

  // No correction for MC
  if (spill.GetIsMC())
    return;

  // Get particle beam
  AnaBeamPD* beam = static_cast<AnaBeamPD*>(spill.Beam);
  AnaParticle* part = static_cast<AnaParticle*>(beam->BeamParticle);
  if (!part)                        return; //?
    
  // The un-corrected particle
  const AnaParticle* original = static_cast<const AnaParticle*>(part->Original);
 
  if (!original)                    return; //?
  
  float corr = 1.041883726;
      
  part->Momentum = original->Momentum * corr;
  beam->BeamMomentum = beam->BeamMomentum * corr;
  
  //  std::cout << "old momentum " << original->Momentum << " || new momentum " << part->Momentum << std::endl;

}


