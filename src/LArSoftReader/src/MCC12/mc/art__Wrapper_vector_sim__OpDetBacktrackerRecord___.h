//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Tue Jun  4 02:53:16 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /pnfs/dune/scratch/users/tjyang/v08_21_00/reco/mcc12_protodune_beam_p1GeV_cosmics_12ms_test/19285494_84/mcc12_protodune_beam_p1GeV_cosmics_12ms_test_84_20190601T145951_reco0.root
//////////////////////////////////////////////////////////


#ifndef art__Wrapper_vector_sim__OpDetBacktrackerRecord____h
#define art__Wrapper_vector_sim__OpDetBacktrackerRecord____h
#include <vector>
namespace std {} using namespace std;
namespace sim {
class OpDetBacktrackerRecord;
} // end of namespace.
namespace art {
template <typename T0> class Wrapper;
template <> class Wrapper<vector<sim::OpDetBacktrackerRecord> >;
} // end of namespace.

namespace std {} using namespace std;
#include "sim__OpDetBacktrackerRecord.h"
#include "art__EDProduct.h"
#include "Riostream.h"

namespace art {
template <typename T0> class Wrapper;
template <> class Wrapper<vector<sim::OpDetBacktrackerRecord> > : public art::EDProduct {

public:
// Nested classes declaration.

public:
// Data Members.
   bool        present;     //
   unsigned int rangeSetID;    //
   vector<sim::OpDetBacktrackerRecord> obj;           // (sim::OpDetBacktrackerRecord)

   Wrapper();
   Wrapper(const Wrapper & );
   virtual ~Wrapper();

   ClassDef(Wrapper,12); // Generated by MakeProject.
};
} // namespace
#endif
