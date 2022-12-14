//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Tue Jun  4 02:53:16 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /pnfs/dune/scratch/users/tjyang/v08_21_00/reco/mcc12_protodune_beam_p1GeV_cosmics_12ms_test/19285494_84/mcc12_protodune_beam_p1GeV_cosmics_12ms_test_84_20190601T145951_reco0.root
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
   TrackTrajectory(const TrackTrajectory & );
   virtual ~TrackTrajectory();

   ClassDef(TrackTrajectory,12); // Generated by MakeProject.
};
} // namespace
#endif
