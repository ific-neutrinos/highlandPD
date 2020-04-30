//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Fri Apr 24 09:37:05 2020 by ROOT version 6.18/04)
//      from the StreamerInfo in file /pnfs/dune/archive/sam_managed_users/calcuttj/data/5/c/5/d/PDSPProd2_protoDUNE_sp_reco_35ms_sce_off_23503706_0_124d8226-6d53-4dbc-af94-340521eab19b_pma.root
//////////////////////////////////////////////////////////


#ifndef sim__ProtoDUNEbeamsim_h
#define sim__ProtoDUNEbeamsim_h
namespace sim {
class ProtoDUNEbeamsim;
} // end of namespace.

#include "Rtypes.h"
#include "Riostream.h"
#include <vector>
#include "sim__ProtoDUNEBeamInstrument.h"

namespace sim {
class ProtoDUNEbeamsim {

public:
// Nested classes declaration.

public:
// Data Members.
   vector<sim::ProtoDUNEBeamInstrument> fAllInstruments;    //

   ProtoDUNEbeamsim();
   ProtoDUNEbeamsim(ProtoDUNEbeamsim && ) = default;
   ProtoDUNEbeamsim(const ProtoDUNEbeamsim & );
   virtual ~ProtoDUNEbeamsim();

   ClassDef(ProtoDUNEbeamsim,14); // Generated by MakeProject.
};
} // namespace
#endif
