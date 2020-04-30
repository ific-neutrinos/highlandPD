//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Fri Apr 24 09:37:06 2020 by ROOT version 6.18/04)
//      from the StreamerInfo in file /pnfs/dune/archive/sam_managed_users/calcuttj/data/5/c/5/d/PDSPProd2_protoDUNE_sp_reco_35ms_sce_off_23503706_0_124d8226-6d53-4dbc-af94-340521eab19b_pma.root
//////////////////////////////////////////////////////////


#ifndef art__Wrapper_vector_anab__FeatureVector_4______h
#define art__Wrapper_vector_anab__FeatureVector_4______h
#include <vector>
namespace std {} using namespace std;
namespace anab {
template <int T0> class FeatureVector;
template <> class FeatureVector<4>;
} // end of namespace.
namespace art {
template <typename T0> class Wrapper;
template <> class Wrapper<vector<anab::FeatureVector<4> > >;
} // end of namespace.

namespace std {} using namespace std;
#include "anab__FeatureVector_4_.h"
#include "art__EDProduct.h"
#include "Riostream.h"

namespace art {
template <typename T0> class Wrapper;
template <> class Wrapper<vector<anab::FeatureVector<4> > > : public art::EDProduct {

public:
// Nested classes declaration.

public:
// Data Members.
   bool        present;     //
   unsigned int rangeSetID;    //
   vector<anab::FeatureVector<4> > obj;           // (anab::FeatureVector<4>)

   Wrapper();
   Wrapper(Wrapper && ) = default;
   Wrapper(const Wrapper & );
   virtual ~Wrapper();

   ClassDef(Wrapper,12); // Generated by MakeProject.
};
} // namespace
#endif
