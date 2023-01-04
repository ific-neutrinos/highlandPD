#include "SCEResidualRangeVariation.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
SCEResidualRangeVariation::SCEResidualRangeVariation():EventVariationBase(),BinnedParams(std::string(getenv("PDSYSTEMATICSROOT"))+"/data","SCE2", BinnedParams::k3D_SYMMETRIC_NOMEAN){
//********************************************************************

  // Read the systematic source parameters from the data files
  SetNParameters(GetNBins());

  //Initialize SCE object to apply variations
  for(int i = 0; i < 100; i++)_sce[i] = NULL;
}

//********************************************************************
SCEResidualRangeVariation::~SCEResidualRangeVariation(){
//********************************************************************
  
  for(int i = 0; i < 100; i++)delete _sce[i];
}

//********************************************************************
void SCEResidualRangeVariation::Apply(const ToyExperiment& toy, AnaEventC& event){
//********************************************************************

  // Get the SystBox for this event. The SystBox contains vector of objects that are relevant for this systematic.
  // In this way the loop over objects is restricted to the ones that really matter
  SystBoxB* box = GetSystBox(event);

  //initialize the SCE for this toy if it has not been initialized yet 
  int toy_index = toy.GetToyIndex();
  if(!_sce[toy_index]){
    _sce[toy_index] = new SpaceCharge();
    _sce[toy_index]->Initialize();
  }
  
  //Vary SCE map
  if(!_sce[toy_index]->IsVaried())
    VarySCEMap(toy);

  // Loop over all relevant tracks for this variation
  for(Int_t ipart = 0; ipart < box->nRelevantRecObjects; ipart++){
     AnaParticlePD* part = static_cast<AnaParticlePD*>(box->RelevantRecObjects[ipart]);
  
    // The un-varied particle
    const AnaParticlePD* original = static_cast<const AnaParticlePD*>(part->Original);
    if(!original)continue;

    if(part->Hits[2].empty())continue;
    
    _sce[toy_index]->ApplyPositionCorrection(part);

    pdAnaUtils::ComputeResidualRange(part);
  }
}

//********************************************************************
bool SCEResidualRangeVariation::UndoSystematic(AnaEventC& event){
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
void SCEResidualRangeVariation::VarySCEMap(const ToyExperiment& toy){
//********************************************************************
  
  int toy_index = toy.GetToyIndex();

  std::cout << "varying toy map " << toy_index << std::endl;

  //get map binning
  int nbinsx = _sce[toy_index]->GetNbinsX();
  int nbinsy = _sce[toy_index]->GetNbinsY();
  int nbinsz = _sce[toy_index]->GetNbinsZ();

  //loop over voxels
  Float_t xcenter,ycenter,zcenter;
  for(int ix = 0; ix < nbinsx; ix++){
    xcenter = _sce[toy_index]->GetBinCenterX(ix+1);
    for(int iy = 0; iy < nbinsy; iy++){
      ycenter = _sce[toy_index]->GetBinCenterY(iy+1);
      for(int iz = 0; iz < nbinsz; iz++){
	zcenter = _sce[toy_index]->GetBinCenterZ(iz+1);
	Int_t index = 0;
	Float_t sigma = 0;
	if(!GetBinSigmaValue(xcenter,ycenter,zcenter,sigma,index))
	  continue;
	_sce[toy_index]->ApplyVoxelVariation(ix+1,iy+1,iz+1,
					     1 + sigma*toy.GetToyVariations(_index)->Variations[index],
					     false); //do not reset splines now, only at the end
      }
    }
  }
  _sce[toy_index]->ResetSplines(); //now recompute splines;
  
}

//**************************************************
bool SCEResidualRangeVariation::IsRelevantRecObject(const AnaEventC& event, const AnaRecObjectC& part) const{
//**************************************************

  (void)event;
  
  //the systematic box only includes candidates and daughters, and all of them are relevant
  return true;
}