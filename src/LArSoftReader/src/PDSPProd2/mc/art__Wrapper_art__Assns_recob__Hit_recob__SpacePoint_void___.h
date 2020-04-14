//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Thu Feb 13 18:14:53 2020 by ROOT version 6.16/00)
//      from the StreamerInfo in file ../PDSPProd2_protoDUNE_sp_reco_35ms_sce_off_23473772_0_452d9f89-a2a1-4680-ab72-853a3261da5d.root
//////////////////////////////////////////////////////////


#ifndef art__Wrapper_art__Assns_recob__Hit_recob__SpacePoint_void____h
#define art__Wrapper_art__Assns_recob__Hit_recob__SpacePoint_void____h
namespace recob {
class Hit;
} // end of namespace.
namespace recob {
class SpacePoint;
} // end of namespace.
namespace art {
template <typename T0, typename T1, typename T2> class Assns;
template <> class Assns<recob::Hit,recob::SpacePoint,void>;
} // end of namespace.
namespace art {
template <typename T0> class Wrapper;
template <> class Wrapper<art::Assns<recob::Hit,recob::SpacePoint,void> >;
} // end of namespace.

#include "Rtypes.h"
#include "art__Assns_recob__Hit_recob__SpacePoint_void_.h"
#include "art__EDProduct.h"
#include "recob__Hit.h"
#include "recob__SpacePoint.h"

namespace art {
template <typename T0> class Wrapper;
template <> class Wrapper<art::Assns<recob::Hit,recob::SpacePoint,void> > : public art::EDProduct {

public:
// Nested classes declaration.

public:
// Data Members.
   bool        present;     //
   unsigned int rangeSetID;    //
   art::Assns<recob::Hit,recob::SpacePoint,void> obj;           //

   Wrapper();
   Wrapper(const Wrapper & );
   virtual ~Wrapper();

   ClassDef(Wrapper,12); // Generated by MakeProject.
};
} // namespace
#endif
