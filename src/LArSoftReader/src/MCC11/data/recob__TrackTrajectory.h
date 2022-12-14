//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Sat Mar  9 10:19:31 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file ../../../protoDuneExampleAnalysis/v0r0/cmt/np04_raw_run005250_0014_dl9_reco_13061912_0_20181110T220718.root
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

   ClassDef(TrackTrajectory,11); // Generated by MakeProject.
};
} // namespace
#endif
