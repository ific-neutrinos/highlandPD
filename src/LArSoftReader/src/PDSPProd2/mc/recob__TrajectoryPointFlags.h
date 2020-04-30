//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Fri Apr 24 09:37:06 2020 by ROOT version 6.18/04)
//      from the StreamerInfo in file /pnfs/dune/archive/sam_managed_users/calcuttj/data/5/c/5/d/PDSPProd2_protoDUNE_sp_reco_35ms_sce_off_23503706_0_124d8226-6d53-4dbc-af94-340521eab19b_pma.root
//////////////////////////////////////////////////////////


#ifndef recob__TrajectoryPointFlags_h
#define recob__TrajectoryPointFlags_h
namespace recob {
class TrajectoryPointFlags;
} // end of namespace.

#include "Rtypes.h"
#include "util__flags__FlagSet_32_unsigned_int_.h"

namespace recob {
class TrajectoryPointFlags {

public:
// Nested classes declaration.

public:
// Data Members.
   unsigned int fFromHit;    //Index of the original hit.
   util::flags::FlagSet<32,unsigned int> fFlags;      //Set of flags

   TrajectoryPointFlags();
   TrajectoryPointFlags(TrajectoryPointFlags && ) = default;
   TrajectoryPointFlags(const TrajectoryPointFlags & );
   virtual ~TrajectoryPointFlags();

};
} // namespace
#endif
