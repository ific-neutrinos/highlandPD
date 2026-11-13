#include "SCEGeometricVariation.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdServiceManager.hxx"

//********************************************************************
SCEGeometricVariation::SCEGeometricVariation():EventVariationBase(),BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","SCE", BinnedParams::k1D_SYMMETRIC_NOMEAN){
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
}

//********************************************************************
SCEGeometricVariation::~SCEGeometricVariation(){
//********************************************************************

}

//********************************************************************
void SCEGeometricVariation::Apply(const ToyExperiment& toy, AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event. The SystBox contains vector of objects that are relevant for this systematic.
  // In this way the loop over objects is restricted to the ones that really matter
  SystBoxB* box = GetSystBox(event);

  //scaling of the nominal SCE map for this toy
  SCEScaling scale = GetToyScaling(toy);

  // Loop over all relevant tracks for this variation
  for(Int_t ipart = 0; ipart < box->nRelevantRecObjects; ipart++){
     AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[ipart]);
  
    // The un-varied particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if(!original)continue;
    
    // Apply variation to TrajectoryPoints
    int ntps = part->TrjPoints.size();
    if(ntps > 0){
      _sce->ApplyTrjPointPositionCorrection(part, scale);
      _sce->ApplyTrjPointDirectionCorrection(part);
     
      //recompute length
      part->Length = pdAnaUtils::ComputeTrackLengthFromTrajectoryPoints(part);

      // Apply variation to Position/Direction start/end
      pdAnaUtils::ComputeParticlePositionAndDirection(part);
    }
    else //if no trjpoints, modify at least particle position start/end
      _sce->ApplyParticlePositionCorrection(part, scale);
   
    //apply variation to hits
    //We don't have yet the hit-trjp relationship, this is redundant
    //but there is no other way so far
    if(!part->Hits[2].empty()){
      _sce->ApplyPositionCorrection(part, scale);
      _sce->ApplyPitchCorrection(part, scale);
      pdAnaUtils::ComputeResidualRange(part);
    }
   
    //modify derived quantities
    std::pair<double,int>result = pdAnaUtils::Chi2PID(*part,2212);
    part->Chi2Proton = result.first;
    result = pdAnaUtils::Chi2PID(*part,13);
    part->Chi2Muon = result.first;
    part->Chi2ndf = result.second;
    part->RangeMomentum[0] = pdAnaUtils::ComputeRangeMomentum(part->Length,13);
    part->RangeMomentum[1] = pdAnaUtils::ComputeRangeMomentum(part->Length,2212);
  }
}

//********************************************************************
bool SCEGeometricVariation::UndoSystematic(AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event
  SystBoxB* box = GetSystBox(event);

  //_sce->ResetToNominal();
  
  for(Int_t ipart = 0; ipart < box->nRelevantRecObjects; ipart++){
    AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[ipart]);
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if(!original)continue;

    //loop over trajectorypoints
    for(int itrp = 0; itrp < (int)part->TrjPoints.size(); itrp++){
      part->TrjPoints[itrp].Position = original->TrjPoints[itrp].Position;
      part->TrjPoints[itrp].Direction = original->TrjPoints[itrp].Direction;
    }

    //set lenght and positions/direction
    part->Length = original->Length;
    part->PositionStart[0] = original->PositionStart[0];
    part->PositionStart[1] = original->PositionStart[1];
    part->PositionStart[2] = original->PositionStart[2];
    part->PositionEnd[0] = original->PositionEnd[0];
    part->PositionEnd[1] = original->PositionEnd[1];
    part->PositionEnd[2] = original->PositionEnd[2];
    part->DirectionStart[0] = original->DirectionStart[0];
    part->DirectionStart[1] = original->DirectionStart[1];
    part->DirectionStart[2] = original->DirectionStart[2];
    part->DirectionEnd[0] = original->DirectionEnd[0];
    part->DirectionEnd[1] = original->DirectionEnd[1];
    part->DirectionEnd[2] = original->DirectionEnd[2];

    //loop over hits
    for(int ihit = 0; ihit < (int)part->Hits[2].size(); ihit++){
      part->Hits[2][ihit].Pitch = original->Hits[2][ihit].Pitch;
      part->Hits[2][ihit].ResidualRange = original->Hits[2][ihit].ResidualRange;
      part->Hits[2][ihit].Position.SetX(original->Hits[2][ihit].Position.X());
      part->Hits[2][ihit].Position.SetY(original->Hits[2][ihit].Position.Y());
      part->Hits[2][ihit].Position.SetZ(original->Hits[2][ihit].Position.Z());
    }

    part->Chi2Proton       = original->Chi2Proton;
    part->Chi2Muon         = original->Chi2Muon;
    part->Chi2ndf          = original->Chi2ndf;
    part->RangeMomentum[0] = original->RangeMomentum[0];
    part->RangeMomentum[1] = original->RangeMomentum[1];
  }
  
  // Don't reset the spill to corrected
  return false;
}

//********************************************************************
SCEScaling SCEGeometricVariation::GetToyScaling(const ToyExperiment& toy){
//********************************************************************

  Float_t sigma;
  GetSigmaValueForBin(0, sigma); //only 1 bin

  SCEScaling scale;
  scale.Global = 1 + sigma*toy.GetToyVariations(_index)->Variations[0];
  return scale;
}

//**************************************************
bool SCEGeometricVariation::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const{
//**************************************************

  (void)event;
  
  //the systematic box only includes candidates and daughters, and all of them are relevant
  return true;
}
