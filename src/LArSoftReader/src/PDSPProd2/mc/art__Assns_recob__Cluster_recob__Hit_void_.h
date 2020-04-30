//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Fri Apr 24 09:37:05 2020 by ROOT version 6.18/04)
//      from the StreamerInfo in file /pnfs/dune/archive/sam_managed_users/calcuttj/data/5/c/5/d/PDSPProd2_protoDUNE_sp_reco_35ms_sce_off_23503706_0_124d8226-6d53-4dbc-af94-340521eab19b_pma.root
//////////////////////////////////////////////////////////


#ifndef art__Assns_recob__Cluster_recob__Hit_void__h
#define art__Assns_recob__Cluster_recob__Hit_void__h
namespace recob {
class Cluster;
} // end of namespace.
namespace recob {
class Hit;
} // end of namespace.
namespace art {
template <typename T0, typename T1, typename T2> class Assns;
template <> class Assns<recob::Cluster,recob::Hit,void>;
} // end of namespace.

#include "Rtypes.h"
#include "recob__Cluster.h"
#include "recob__Hit.h"
#include "art__detail__AssnsBase.h"
#include "Riostream.h"
#include <vector>
#include <utility>
#include "art__RefCore.h"

namespace art {
template <typename T0, typename T1, typename T2> class Assns;
template <> class Assns<recob::Cluster,recob::Hit,void> : public art::detail::AssnsBase {

public:
// Nested classes declaration.

public:
// Data Members.
   vector<pair<art::RefCore,unsigned long> > ptr_data_1_;    //
   vector<pair<art::RefCore,unsigned long> > ptr_data_2_;    //

   Assns();
   Assns(Assns && ) = default;
   Assns(const Assns & );
   virtual ~Assns();

   ClassDef(Assns,12); // Generated by MakeProject.
};
} // namespace
#ifdef __MAKECINT__
#pragma link C++ class pair<art::RefCore,unsigned long>+;
#endif
#endif
