//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Tue Jun  4 02:53:16 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /pnfs/dune/scratch/users/tjyang/v08_21_00/reco/mcc12_protodune_beam_p1GeV_cosmics_12ms_test/19285494_84/mcc12_protodune_beam_p1GeV_cosmics_12ms_test_84_20190601T145951_reco0.root
//////////////////////////////////////////////////////////


#ifndef art__Wrapper_art__Assns_recob__Cluster_recob__Hit_void____h
#define art__Wrapper_art__Assns_recob__Cluster_recob__Hit_void____h
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
namespace art {
template <typename T0> class Wrapper;
template <> class Wrapper<art::Assns<recob::Cluster,recob::Hit,void> >;
} // end of namespace.

#include "Rtypes.h"
#include "art__Assns_recob__Cluster_recob__Hit_void_.h"
#include "art__EDProduct.h"
#include "recob__Cluster.h"
#include "recob__Hit.h"

namespace art {
template <typename T0> class Wrapper;
template <> class Wrapper<art::Assns<recob::Cluster,recob::Hit,void> > : public art::EDProduct {

public:
// Nested classes declaration.

public:
// Data Members.
   bool        present;     //
   unsigned int rangeSetID;    //
   art::Assns<recob::Cluster,recob::Hit,void> obj;           //

   Wrapper();
   Wrapper(const Wrapper & );
   virtual ~Wrapper();

   ClassDef(Wrapper,12); // Generated by MakeProject.
};
} // namespace
#endif
