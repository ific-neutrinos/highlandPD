//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Fri Apr 24 09:37:05 2020 by ROOT version 6.18/04)
//      from the StreamerInfo in file /pnfs/dune/archive/sam_managed_users/calcuttj/data/5/c/5/d/PDSPProd2_protoDUNE_sp_reco_35ms_sce_off_23503706_0_124d8226-6d53-4dbc-af94-340521eab19b_pma.root
//////////////////////////////////////////////////////////


#ifndef art__HLTPathStatus_h
#define art__HLTPathStatus_h
namespace art {
class HLTPathStatus;
} // end of namespace.

#include "Rtypes.h"

namespace art {
class HLTPathStatus {

public:
// Nested classes declaration.

public:
// Data Members.
   unsigned short status_;     //

   HLTPathStatus();
   HLTPathStatus(HLTPathStatus && ) = default;
   HLTPathStatus(const HLTPathStatus & );
   virtual ~HLTPathStatus();

   ClassDef(HLTPathStatus,11); // Generated by MakeProject.
};
} // namespace
#endif
