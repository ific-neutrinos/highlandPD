//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Fri Apr 24 09:37:05 2020 by ROOT version 6.18/04)
//      from the StreamerInfo in file /pnfs/dune/archive/sam_managed_users/calcuttj/data/5/c/5/d/PDSPProd2_protoDUNE_sp_reco_35ms_sce_off_23503706_0_124d8226-6d53-4dbc-af94-340521eab19b_pma.root
//////////////////////////////////////////////////////////


#ifndef art__RunID_h
#define art__RunID_h
namespace art {
class RunID;
} // end of namespace.

#include "Rtypes.h"

namespace art {
class RunID {

public:
// Nested classes declaration.

public:
// Data Members.
   unsigned int run_;        //

   RunID();
   RunID(RunID && ) = default;
   RunID(const RunID & );
   virtual ~RunID();

   ClassDef(RunID,11); // Generated by MakeProject.
};
} // namespace
#endif
