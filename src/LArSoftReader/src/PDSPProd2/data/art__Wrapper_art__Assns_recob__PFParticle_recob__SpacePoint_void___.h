//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Thu Feb 13 18:08:35 2020 by ROOT version 6.16/00)
//      from the StreamerInfo in file ../np04_raw_run005834_0126_dl7_reco_22481932_0_20190828T190321.root
//////////////////////////////////////////////////////////


#ifndef art__Wrapper_art__Assns_recob__PFParticle_recob__SpacePoint_void____h
#define art__Wrapper_art__Assns_recob__PFParticle_recob__SpacePoint_void____h
namespace recob {
class PFParticle;
} // end of namespace.
namespace recob {
class SpacePoint;
} // end of namespace.
namespace art {
template <typename T0, typename T1, typename T2> class Assns;
template <> class Assns<recob::PFParticle,recob::SpacePoint,void>;
} // end of namespace.
namespace art {
template <typename T0> class Wrapper;
template <> class Wrapper<art::Assns<recob::PFParticle,recob::SpacePoint,void> >;
} // end of namespace.

#include "Rtypes.h"
#include "art__Assns_recob__PFParticle_recob__SpacePoint_void_.h"
#include "art__EDProduct.h"
#include "recob__PFParticle.h"
#include "recob__SpacePoint.h"

namespace art {
template <typename T0> class Wrapper;
template <> class Wrapper<art::Assns<recob::PFParticle,recob::SpacePoint,void> > : public art::EDProduct {

public:
// Nested classes declaration.

public:
// Data Members.
   bool        present;     //
   unsigned int rangeSetID;    //
   art::Assns<recob::PFParticle,recob::SpacePoint,void> obj;           //

   Wrapper();
   Wrapper(const Wrapper & );
   virtual ~Wrapper();

   ClassDef(Wrapper,12); // Generated by MakeProject.
};
} // namespace
#endif
