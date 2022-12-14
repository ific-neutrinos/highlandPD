//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Fri Apr 24 09:37:05 2020 by ROOT version 6.18/04)
//      from the StreamerInfo in file /pnfs/dune/archive/sam_managed_users/calcuttj/data/5/c/5/d/PDSPProd2_protoDUNE_sp_reco_35ms_sce_off_23503706_0_124d8226-6d53-4dbc-af94-340521eab19b_pma.root
//////////////////////////////////////////////////////////


#ifndef art__ProcessHistory_h
#define art__ProcessHistory_h
namespace art {
class ProcessHistory;
} // end of namespace.

#include "Rtypes.h"
#include "Riostream.h"
#include <vector>
#include "art__ProcessConfiguration.h"
#include "art__Transient_art__ProcessHistory__Transients_.h"
namespace art {
#ifndef art__ProcessHistory_h
#define art__ProcessHistory_h
enum ProcessHistory::Transients { kDefault_ProcessHistory::Transients };
#endif
} // end of namespace.

namespace art {
class ProcessHistory {

public:
// Nested classes forward declaration.
enum Transients { kDefault_Transients };

public:
// Nested classes declaration.

public:
// Data Members.
   vector<art::ProcessConfiguration> data_;       //
   art::Transient<art::ProcessHistory::Transients> transients_;    //

   ProcessHistory();
   ProcessHistory(ProcessHistory && ) = default;
   ProcessHistory(const ProcessHistory & );
   virtual ~ProcessHistory();

   ClassDef(ProcessHistory,11); // Generated by MakeProject.
};
} // namespace
#endif
