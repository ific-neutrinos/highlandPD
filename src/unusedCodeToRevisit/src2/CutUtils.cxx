#include "CutUtils.hxx"
#include "FiducialVolumeDefinition.hxx"
#include "EventBoxId.hxx"
#include "EventBoxTracker.hxx"
#include <TVector3.h>
#include <stdio.h>
#include <iostream>
#include <math.h>

//**************************************************
bool cutUtils::TrackQualityCut(const AnaParticleB& part){
//**************************************************

  if (part.NHits > 100) return true;
  else return false;
  
}

//**************************************************
bool cutUtils::MuonPIDCut(const AnaParticleB& part){
//**************************************************

  (void)part;
  return true;
/*
    Float_t cut1 = 0.8;
    Float_t cut2 = 0.05;
    Float_t pcut = 500;

    Float_t PIDLikelihood[4];
    anaUtils::GetPIDLikelihood(track, PIDLikelihood, prod5Cut);

    if (((PIDLikelihood[0]+PIDLikelihood[3])/(1-PIDLikelihood[2]) > cut1 || track.Momentum > pcut ) && (PIDLikelihood[0]>cut2)){
        return true; 
    }

    return false;
*/
}
//**************************************************
bool cutUtils::ProtonPIDCut(const AnaParticleB& part){
//**************************************************

  (void)part;
  return true;
/*
    if (anaUtils::GetPIDLikelihood(track,2)>0.9)
        return true;
    else    
        return false;
*/
}

//**************************************************
bool cutUtils::FiducialCut(const AnaParticleB& track, const SubDetId::SubDetEnum det){
//**************************************************
    if (!SubDetId::GetDetectorUsed(track.Detector, det)) return false;
    return FiducialCut(track.PositionStart, det);
}

//**************************************************
bool cutUtils::FiducialCut(const Float_t* pos, const SubDetId::SubDetEnum det){
//**************************************************

    return anaUtils::InFiducialVolume(det,pos);
}

//**************************************************
bool cutUtils::CommonVertexCut(const AnaParticleB& track1, const AnaParticleB& track2, int cutx, int cuty, int cutz){
//**************************************************

    if (fabs(track1.PositionStart[0]-track2.PositionStart[0]) > cutx ) return false; 
    if (fabs(track1.PositionStart[1]-track2.PositionStart[1]) > cuty ) return false; 
    if (fabs(track1.PositionStart[2]-track2.PositionStart[2]) > cutz ) return false;   

    return true;
}


