//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Thu Feb 14 10:40:59 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/MCC11/mcc11_protoDUNE_sp_reco_12231045_0_1ee4e6cf-25fb-482c-a231-fbd91e816aa4.root
//////////////////////////////////////////////////////////


#ifndef art__Wrapper_art__Assns_recob__Track_anab__T0_void____h
#define art__Wrapper_art__Assns_recob__Track_anab__T0_void____h
namespace recob {
class Track;
} // end of namespace.
namespace anab {
class T0;
} // end of namespace.
namespace art {
template <typename T0, typename T1, typename T2> class Assns;
template <> class Assns<recob::Track,anab::T0,void>;
} // end of namespace.
namespace art {
template <typename T0> class Wrapper;
template <> class Wrapper<art::Assns<recob::Track,anab::T0,void> >;
} // end of namespace.

#include "Rtypes.h"
#include "art__Assns_recob__Track_anab__T0_void_.h"
#include "art__EDProduct.h"
#include "recob__Track.h"
#include "anab__T0.h"

namespace art {
template <typename T0> class Wrapper;
template <> class Wrapper<art::Assns<recob::Track,anab::T0,void> > : public art::EDProduct {

public:
// Nested classes declaration.

public:
// Data Members.
   bool        present;     //
   unsigned int rangeSetID;    //
   art::Assns<recob::Track,anab::T0,void> obj;           //

   Wrapper();
   Wrapper(const Wrapper & );
   virtual ~Wrapper();

   ClassDef(Wrapper,12); // Generated by MakeProject.
};
} // namespace
#endif
