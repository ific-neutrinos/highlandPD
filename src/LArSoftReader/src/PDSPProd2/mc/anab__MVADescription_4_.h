//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Fri Apr 24 09:37:06 2020 by ROOT version 6.18/04)
//      from the StreamerInfo in file /pnfs/dune/archive/sam_managed_users/calcuttj/data/5/c/5/d/PDSPProd2_protoDUNE_sp_reco_35ms_sce_off_23503706_0_124d8226-6d53-4dbc-af94-340521eab19b_pma.root
//////////////////////////////////////////////////////////


#ifndef anab__MVADescription_4__h
#define anab__MVADescription_4__h
namespace anab {
template <int T0> class MVADescription;
template <> class MVADescription<4>;
} // end of namespace.

#include "Rtypes.h"
#include "Riostream.h"
#include <string>

namespace anab {
template <int T0> class MVADescription;
template <> class MVADescription<4> {

public:
// Nested classes declaration.

public:
// Data Members.
   string      fDataTag;    //Tag of the reco data products (art::InputTag format)
   string      fOutputInstance;    //Instance name of the feature vector collection
   string      fOutputNames[4];    //Feature vector entries names/meaning

   MVADescription();
   MVADescription(MVADescription && ) = default;
   MVADescription(const MVADescription & );
   virtual ~MVADescription();

   ClassDef(MVADescription,11); // Generated by MakeProject.
};
} // namespace
#endif
