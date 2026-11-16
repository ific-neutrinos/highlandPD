#ifndef SCEVariationBase_h
#define SCEVariationBase_h

#include "EventVariationBase.hxx"
#include "BinnedParams.hxx"
#include "SpaceCharge.hxx"

/// Common part of the space charge systematics: how the nominal SCE map
/// (pd::services().SCE()) is varied in each toy. The derived classes
/// decide what the varied map is applied to:
///  - SCEGeometricVariation: track geometry (trajectory points, length,
///    hit positions and pitch), not the calorimetry, which is covered
///    by dEdxVariation.
///  - SCECalorimetricVariation: hit positions and the pitch correction
///    of dQ/dx.
///
/// Two variation modes, chosen with pdSystematics.SCE.VariationMode:
///  - 0, global: the whole map is scaled by 1 + sigma*v, with one
///    parameter whose sigma comes from data/SCE.dat.
///  - 1, voxel: each node of the map is scaled by 1 + sigma_i*v_i. The
///    boxes and their sigmas come from data/<pdSystematics.SCE.VoxelFile>,
///    one line per box: xmin xmax ymin ymax zmin zmax sigma. Each box is
///    one systematic parameter; a node belongs to the first box that
///    contains its center (min <= value < max), and nodes in no box are
///    not varied.
/// In both modes the map is never copied: the nominal map is scaled on the
/// fly (see SCEScaling in SpaceCharge.hxx).
class SCEVariationBase: public EventVariationBase, public BinnedParams{
public:

  SCEVariationBase();
  virtual ~SCEVariationBase(){}

  /// Scaling of the nominal SCE map for this toy.
  /// Use it with the nominal map, pd::services().SCE().
  SCEScaling GetToyScaling(const ToyExperiment& toy);

  bool IsVoxelMode() const {return _voxelMode;}

protected:

  /// Read the voxel boxes and assign each map node to its box
  void InitializeVoxels(const std::string& filename);

  const SpaceCharge* _sce; ///< nominal map, owned by pd::services()

  bool _voxelMode;
  SCEVoxelSigmas _voxelSigmas;
};

#endif
