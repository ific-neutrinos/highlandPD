//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Thu Feb 14 10:40:59 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/MCC11/mcc11_protoDUNE_sp_reco_12231045_0_1ee4e6cf-25fb-482c-a231-fbd91e816aa4.root
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
