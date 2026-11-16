#include "SCEVariationBase.hxx"
#include "pdDataFiles.hxx"
#include "pdServiceManager.hxx"
#include "Parameters.hxx"

#include <fstream>
#include <sstream>
#include <cstdlib>

//********************************************************************
SCEVariationBase::SCEVariationBase():EventVariationBase(),BinnedParams(pd::DataDir("pdSystematics"),"SCE", BinnedParams::k1D_SYMMETRIC_NOMEAN){
//********************************************************************

  //Toys use the nominal map scaled on the fly, there are no per-toy copies
  _sce = &pd::services().SCE();

  _voxelMode = (ND::params().GetParameterI("pdSystematics.SCE.VariationMode") == 1);

  if(!_voxelMode){
    // Global: the parameters come from data/SCE.dat
    SetNParameters(GetNBins());

    //check we have enough bins to run this systematic, otherwise exit
    if(GetNBins() > (int)NMAXSYSTSOURCEBINS){
      std::cerr << "SCEVariationBase::SCEVariationBase()." << std::endl;
      std::cerr << "the source files has more bins than allowed for the systematic propagation!" << std::endl;
      std::cerr << "N MAX BINS allowed " << NMAXSYSTSOURCEBINS << std::endl;
      std::cerr << "N BINS on source   " << GetNBins() << std::endl;
      std::cout << "if you need more bins, look at psyckeUtils/src/BinnedParams.hxx" << std::endl;
      std::exit(1);
    }
  }
  else{
    // Voxel: our own reader, BinnedParams is limited to NMAXSYSTSOURCEBINS bins
    InitializeVoxels(pd::DataFile("pdSystematics", 
                     ND::params().GetParameterS("pdSystematics.SCE.VoxelFile")));
  }
}

//********************************************************************
void SCEVariationBase::InitializeVoxels(const std::string& filename){
//********************************************************************

  struct Box{double min[3], max[3]; float sigma;};
  std::vector<Box> boxes;

  std::ifstream in(filename.c_str());
  if(!in){
    std::cerr << "SCEVariationBase::InitializeVoxels(). Can't open file " << filename << std::endl;
    std::exit(1);
  }
  std::string line;
  while(std::getline(in, line)){
    std::istringstream ss(line);
    Box b;
    if(ss >> b.min[0] >> b.max[0] >> b.min[1] >> b.max[1] >> b.min[2] >> b.max[2] >> b.sigma)
      boxes.push_back(b);
  }
  if(boxes.empty()){
    std::cerr << "SCEVariationBase::InitializeVoxels(). No boxes in " << filename << std::endl;
    std::exit(1);
  }

  // one systematic parameter per box
  SetNParameters(boxes.size());

  // assign each node of the calibration maps of each drift volume to its box.
  // All the maps of a drift volume share binning: use Dx.
  // CalSCEhistograms: 6 = Dx x<0 (volume 0), 0 = Dx x>0 (volume 1)
  int nvaried[2] = {0, 0};
  for(int dv = 0; dv < 2; dv++){
    const TH3F* h = _sce->GetCalHistogram(dv == 0 ? 6 : 0);
    const TAxis* axes[3] = {h->GetXaxis(), h->GetYaxis(), h->GetZaxis()};
    int nbins = (h->GetNbinsX()+2)*(h->GetNbinsY()+2)*(h->GetNbinsZ()+2);
    _voxelSigmas.Index[dv].assign(nbins, -1);
    _voxelSigmas.Sigma[dv].assign(nbins, 0.);

    for(int ib = 0; ib < (int)boxes.size(); ib++){
      // range of nodes whose centers are inside the box, for each axis
      int first[3], last[3];
      for(int a = 0; a < 3; a++){
        first[a] = axes[a]->GetNbins()+1; last[a] = 0;
        for(int i = 1; i <= axes[a]->GetNbins(); i++){
          double c = axes[a]->GetBinCenter(i);
          if(c >= boxes[ib].min[a] && c < boxes[ib].max[a]){
            if(i < first[a]) first[a] = i;
            last[a] = i;
          }
        }
      }
      for(int ix = first[0]; ix <= last[0]; ix++)
        for(int iy = first[1]; iy <= last[1]; iy++)
          for(int iz = first[2]; iz <= last[2]; iz++){
            int bin = h->GetBin(ix, iy, iz);
            if(_voxelSigmas.Index[dv][bin] >= 0) continue; // first box wins
            _voxelSigmas.Index[dv][bin] = ib;
            _voxelSigmas.Sigma[dv][bin] = boxes[ib].sigma;
            nvaried[dv]++;
          }
    }
  }

  std::cout << "SCEVariationBase::InitializeVoxels(). " << boxes.size() << " boxes from " << filename
            << "; varied nodes: " << nvaried[0] << " (x<0), " << nvaried[1] << " (x>0)" << std::endl;
}

//********************************************************************
SCEScaling SCEVariationBase::GetToyScaling(const ToyExperiment& toy){
//********************************************************************

  SCEScaling scale;

  if(_voxelMode){
    scale.Voxels        = &_voxelSigmas;
    scale.ToyVariations = toy.GetToyVariations(_index)->Variations;
  }
  else{
    Float_t sigma;
    GetSigmaValueForBin(0, sigma); //only 1 bin
    scale.Global = 1 + sigma*toy.GetToyVariations(_index)->Variations[0];
  }
  return scale;
}
