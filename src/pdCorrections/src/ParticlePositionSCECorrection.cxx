#include "ParticlePositionSCECorrection.hxx"
#include "pdDataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include <cassert>


//********************************************************************
ParticlePositionSCECorrection::ParticlePositionSCECorrection(){
//********************************************************************

  _sce = new SpaceCharge();
  _sce->Initialize();
}

//********************************************************************
void ParticlePositionSCECorrection::Apply(AnaSpillC& spillC){
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

    //correct trajectory points
    _sce->ApplyTrjPointPositionCorrection(part);
    _sce->ApplyTrjPointDirectionCorrection(part);
    
    //correct end and start position/direction
    int ntps = part->TrjPoints.size();
    int ifirst = -1;
    for(int itp = 0; itp < ntps; itp++){
      if(part->TrjPoints[itp].IsValid()){
	ifirst = itp;
	break;
      }
    }
    if(ifirst != -1){
      part->PositionStart[0] = part->TrjPoints[ifirst].Position.X();
      part->PositionStart[1] = part->TrjPoints[ifirst].Position.X();
      part->PositionStart[2] = part->TrjPoints[ifirst].Position.X();
      part->DirectionStart[0] = part->TrjPoints[ifirst].Direction.X();
      part->DirectionStart[1] = part->TrjPoints[ifirst].Direction.X();
      part->DirectionStart[2] = part->TrjPoints[ifirst].Direction.X();
    }
    
    int ilast  = -1;
    for(int itp = 1; itp < ntps; itp++){
      if(part->TrjPoints[ntps-itp].IsValid()){
	ilast = itp;
	break;
      }
    }
    if(ilast != -1){
      part->PositionEnd[0] = part->TrjPoints[ilast].Position.X();
      part->PositionEnd[1] = part->TrjPoints[ilast].Position.X();
      part->PositionEnd[2] = part->TrjPoints[ilast].Position.X();
      part->DirectionEnd[0] = part->TrjPoints[ilast].Direction.X();
      part->DirectionEnd[1] = part->TrjPoints[ilast].Direction.X();
      part->DirectionEnd[2] = part->TrjPoints[ilast].Direction.X();
    }
  }
}

