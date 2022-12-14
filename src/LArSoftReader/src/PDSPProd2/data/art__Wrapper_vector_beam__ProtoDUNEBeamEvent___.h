//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Thu Feb 13 18:08:35 2020 by ROOT version 6.16/00)
//      from the StreamerInfo in file ../np04_raw_run005834_0126_dl7_reco_22481932_0_20190828T190321.root
//////////////////////////////////////////////////////////


#ifndef art__Wrapper_vector_beam__ProtoDUNEBeamEvent____h
#define art__Wrapper_vector_beam__ProtoDUNEBeamEvent____h
#include <vector>
namespace std {} using namespace std;
namespace beam {
class ProtoDUNEBeamEvent;
} // end of namespace.
namespace art {
template <typename T0> class Wrapper;
template <> class Wrapper<vector<beam::ProtoDUNEBeamEvent> >;
} // end of namespace.

namespace std {} using namespace std;
#include "beam__ProtoDUNEBeamEvent.h"
#include "art__EDProduct.h"
#include "Riostream.h"

namespace art {
template <typename T0> class Wrapper;
template <> class Wrapper<vector<beam::ProtoDUNEBeamEvent> > : public art::EDProduct {

public:
// Nested classes declaration.

public:
// Data Members.
   bool        present;     //
   unsigned int rangeSetID;    //
   vector<beam::ProtoDUNEBeamEvent> obj;           // (beam::ProtoDUNEBeamEvent)

   Wrapper();
   Wrapper(const Wrapper & );
   virtual ~Wrapper();

   ClassDef(Wrapper,12); // Generated by MakeProject.
};
} // namespace
#endif
