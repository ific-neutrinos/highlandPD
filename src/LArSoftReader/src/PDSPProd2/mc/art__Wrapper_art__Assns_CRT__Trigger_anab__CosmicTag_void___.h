//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Fri Apr 24 09:37:05 2020 by ROOT version 6.18/04)
//      from the StreamerInfo in file /pnfs/dune/archive/sam_managed_users/calcuttj/data/5/c/5/d/PDSPProd2_protoDUNE_sp_reco_35ms_sce_off_23503706_0_124d8226-6d53-4dbc-af94-340521eab19b_pma.root
//////////////////////////////////////////////////////////


#ifndef art__Wrapper_art__Assns_CRT__Trigger_anab__CosmicTag_void____h
#define art__Wrapper_art__Assns_CRT__Trigger_anab__CosmicTag_void____h
namespace CRT {
class Trigger;
} // end of namespace.
namespace anab {
class CosmicTag;
} // end of namespace.
namespace art {
template <typename T0, typename T1, typename T2> class Assns;
template <> class Assns<CRT::Trigger,anab::CosmicTag,void>;
} // end of namespace.
namespace art {
template <typename T0> class Wrapper;
template <> class Wrapper<art::Assns<CRT::Trigger,anab::CosmicTag,void> >;
} // end of namespace.

#include "Rtypes.h"
#include "art__Assns_CRT__Trigger_anab__CosmicTag_void_.h"
#include "art__EDProduct.h"
#include "CRT__Trigger.h"
#include "anab__CosmicTag.h"

namespace art {
template <typename T0> class Wrapper;
template <> class Wrapper<art::Assns<CRT::Trigger,anab::CosmicTag,void> > : public art::EDProduct {

public:
// Nested classes declaration.

public:
// Data Members.
   bool        present;     //
   unsigned int rangeSetID;    //
   art::Assns<CRT::Trigger,anab::CosmicTag,void> obj;           //

   Wrapper();
   Wrapper(Wrapper && ) = default;
   Wrapper(const Wrapper & );
   virtual ~Wrapper();

   ClassDef(Wrapper,12); // Generated by MakeProject.
};
} // namespace
#endif
