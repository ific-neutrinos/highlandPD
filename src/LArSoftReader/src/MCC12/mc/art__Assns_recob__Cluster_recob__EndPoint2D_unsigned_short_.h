//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Tue Jun  4 02:53:16 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /pnfs/dune/scratch/users/tjyang/v08_21_00/reco/mcc12_protodune_beam_p1GeV_cosmics_12ms_test/19285494_84/mcc12_protodune_beam_p1GeV_cosmics_12ms_test_84_20190601T145951_reco0.root
//////////////////////////////////////////////////////////


#ifndef art__Assns_recob__Cluster_recob__EndPoint2D_unsigned_short__h
#define art__Assns_recob__Cluster_recob__EndPoint2D_unsigned_short__h
namespace recob {
class Cluster;
} // end of namespace.
namespace recob {
class EndPoint2D;
} // end of namespace.
namespace art {
template <typename T0, typename T1, typename T2> class Assns;
template <> class Assns<recob::Cluster,recob::EndPoint2D,unsigned short>;
} // end of namespace.

#include "Rtypes.h"
#include "recob__Cluster.h"
#include "recob__EndPoint2D.h"
#include "art__Assns_recob__Cluster_recob__EndPoint2D_void_.h"
#include "Riostream.h"
#include <vector>

namespace art {
template <typename T0, typename T1, typename T2> class Assns;
template <> class Assns<recob::Cluster,recob::EndPoint2D,unsigned short> : public art::Assns<recob::Cluster,recob::EndPoint2D,void> {

public:
// Nested classes declaration.

public:
// Data Members.
   vector<unsigned short> data_;       //

   Assns();
   Assns(const Assns & );
   virtual ~Assns();

   ClassDef(Assns,12); // Generated by MakeProject.
};
} // namespace
#endif
