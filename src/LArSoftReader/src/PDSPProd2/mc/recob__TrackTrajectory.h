//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Fri Apr 24 09:37:06 2020 by ROOT version 6.18/04)
//      from the StreamerInfo in file /pnfs/dune/archive/sam_managed_users/calcuttj/data/5/c/5/d/PDSPProd2_protoDUNE_sp_reco_35ms_sce_off_23503706_0_124d8226-6d53-4dbc-af94-340521eab19b_pma.root
//////////////////////////////////////////////////////////


#ifndef recob__TrackTrajectory_h
#define recob__TrackTrajectory_h
namespace recob {
class TrackTrajectory;
} // end of namespace.

#include "Rtypes.h"
#include "recob__Trajectory.h"
#include "Riostream.h"
#include <vector>
#include "recob__TrajectoryPointFlags.h"

namespace recob {
class TrackTrajectory : public recob::Trajectory {

public:
// Nested classes declaration.

public:
// Data Members.
   vector<recob::TrajectoryPointFlags> fFlags;      //Flags of each of the points in trajectory

   TrackTrajectory();
   TrackTrajectory(TrackTrajectory && ) = default;
   TrackTrajectory(const TrackTrajectory & );
   virtual ~TrackTrajectory();

   ClassDef(TrackTrajectory,12); // Generated by MakeProject.
};
} // namespace
#endif
