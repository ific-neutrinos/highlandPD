//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Fri Apr 24 09:37:05 2020 by ROOT version 6.18/04)
//      from the StreamerInfo in file /pnfs/dune/archive/sam_managed_users/calcuttj/data/5/c/5/d/PDSPProd2_protoDUNE_sp_reco_35ms_sce_off_23503706_0_124d8226-6d53-4dbc-af94-340521eab19b_pma.root
//////////////////////////////////////////////////////////


#ifndef CRT__Trigger_h
#define CRT__Trigger_h
namespace CRT {
class Trigger;
} // end of namespace.

#include "Rtypes.h"
#include "Riostream.h"
#include <vector>
#include "CRT__Hit.h"

namespace CRT {
class Trigger {

public:
// Nested classes declaration.

public:
// Data Members.
   unsigned short fChannel;    //Mapping to CRT module that was Triggered.  Index into the Geometry service's array of AuxDetGeos
   ULong64_t      fTimestamp;    //Timestamp when this Trigger occurred.  First 32 bits are Linux time since epoch; last 32 bits are time in
   vector<CRT::Hit> fHits;         //All activity in CRT strips within this module when it was read out

   Trigger();
   Trigger(Trigger && ) = default;
   Trigger(const Trigger & );
   virtual ~Trigger();

   ClassDef(Trigger,11); // Generated by MakeProject.
};
} // namespace
#endif
