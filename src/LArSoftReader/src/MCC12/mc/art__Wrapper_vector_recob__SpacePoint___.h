//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Tue Jun  4 02:53:16 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /pnfs/dune/scratch/users/tjyang/v08_21_00/reco/mcc12_protodune_beam_p1GeV_cosmics_12ms_test/19285494_84/mcc12_protodune_beam_p1GeV_cosmics_12ms_test_84_20190601T145951_reco0.root
//////////////////////////////////////////////////////////


#ifndef art__Wrapper_vector_recob__SpacePoint____h
#define art__Wrapper_vector_recob__SpacePoint____h
#include <vector>
namespace std {} using namespace std;
namespace recob {
class SpacePoint;
} // end of namespace.
namespace art {
template <typename T0> class Wrapper;
template <> class Wrapper<vector<recob::SpacePoint> >;
} // end of namespace.

namespace std {} using namespace std;
#include "recob__SpacePoint.h"
#include "art__EDProduct.h"
#include "Riostream.h"

namespace art {
template <typename T0> class Wrapper;
template <> class Wrapper<vector<recob::SpacePoint> > : public art::EDProduct {

public:
// Nested classes declaration.

public:
// Data Members.
   bool        present;     //
   unsigned int rangeSetID;    //
   vector<recob::SpacePoint> obj;           // (recob::SpacePoint)

   Wrapper();
   Wrapper(const Wrapper & );
   virtual ~Wrapper();

   ClassDef(Wrapper,12); // Generated by MakeProject.
};
} // namespace
#endif
