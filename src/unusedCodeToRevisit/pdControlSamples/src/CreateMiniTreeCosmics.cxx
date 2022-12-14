#include "CreateMiniTreeCosmics.hxx"
#include "Parameters.hxx"

//********************************************************************
CreateMiniTreeCosmics::CreateMiniTreeCosmics(int argc, char *argv[]):CreateMiniTreePD(argc, argv){
//********************************************************************

}

//********************************************************************
bool CreateMiniTreeCosmics::Initialize(){
//********************************************************************

  // call the base class
  CreateMiniTreePD::Initialize();
  
  _cutLength = ND::params().GetParameterD("pdControlSamples.Cosmics.Cuts.Length");
  _cutStartZmin   = ND::params().GetParameterD("pdControlSamples.Cosmics.Cuts.startZmin");
  _cutStartZmax   = ND::params().GetParameterD("pdControlSamples.Cosmics.Cuts.startZmax");
  
  return true;
}


//********************************************************************
bool CreateMiniTreeCosmics::CheckSaveParticle(const AnaParticleB& partB){
//********************************************************************

  const AnaParticlePD& part = *static_cast<const AnaParticlePD*>(&partB);

  if (part.ParentID==-1 &&
      part.Length > _cutLength &&
      part.PositionStart[2]>_cutStartZmin && part.PositionStart[2]<_cutStartZmax){

    return true;
  }
  else
    return false;
}

