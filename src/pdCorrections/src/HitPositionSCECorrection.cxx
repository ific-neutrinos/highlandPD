#include "HitPositionSCECorrection.hxx"
#include "pdDataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>


//********************************************************************
HitPositionSCECorrection::HitPositionSCECorrection(){
//********************************************************************

  sce = new SpaceCharge();
  sce->Initialize();
  
}

//********************************************************************
void HitPositionSCECorrection::Apply(AnaSpillC& spillC){
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

    for (UInt_t j = 0; j < part->Hits[2].size(); j++){
      TVector3 offset = sce->GetCalPosOffsets(part->Hits[2][j].Position_NoSCE, part->Hits[2][j].TPCid);
      //std::cout << "pre " << part->Hits[2][j].Position.X() << " " << part->Hits[2][j].Position.Y() << " " << part->Hits[2][j].Position.Z() << std::endl;
      part->Hits[2][j].Position.SetX(part->Hits[2][j].Position_NoSCE.X() - offset.X());
      part->Hits[2][j].Position.SetY(part->Hits[2][j].Position_NoSCE.Y() + offset.Y());
      part->Hits[2][j].Position.SetZ(part->Hits[2][j].Position_NoSCE.Z() + offset.Z());
      //std::cout << "after " << part->Hits[2][j].Position.X() << " " << part->Hits[2][j].Position.Y() << " " << part->Hits[2][j].Position.Z() << std::endl;
    }
  }
}

