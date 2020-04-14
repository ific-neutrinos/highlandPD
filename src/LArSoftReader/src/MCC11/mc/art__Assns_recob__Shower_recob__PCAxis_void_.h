//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Thu Feb 14 10:40:59 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/MCC11/mcc11_protoDUNE_sp_reco_12231045_0_1ee4e6cf-25fb-482c-a231-fbd91e816aa4.root
//////////////////////////////////////////////////////////


#ifndef art__Assns_recob__Shower_recob__PCAxis_void__h
#define art__Assns_recob__Shower_recob__PCAxis_void__h
namespace recob {
class Shower;
} // end of namespace.
namespace recob {
class PCAxis;
} // end of namespace.
namespace art {
template <typename T0, typename T1, typename T2> class Assns;
template <> class Assns<recob::Shower,recob::PCAxis,void>;
} // end of namespace.

#include "Rtypes.h"
#include "recob__Shower.h"
#include "recob__PCAxis.h"
#include "art__detail__AssnsBase.h"
#include "Riostream.h"
#include <vector>
#include <utility>
#include "art__RefCore.h"

namespace art {
template <typename T0, typename T1, typename T2> class Assns;
template <> class Assns<recob::Shower,recob::PCAxis,void> : public art::detail::AssnsBase {

public:
// Nested classes declaration.

public:
// Data Members.
   vector<pair<art::RefCore,unsigned long> > ptr_data_1_;    //
   vector<pair<art::RefCore,unsigned long> > ptr_data_2_;    //

   Assns();
   Assns(const Assns & );
   virtual ~Assns();

   ClassDef(Assns,12); // Generated by MakeProject.
};
} // namespace
#ifdef __MAKECINT__
#pragma link C++ class pair<art::RefCore,unsigned long>+;
#endif
#endif
