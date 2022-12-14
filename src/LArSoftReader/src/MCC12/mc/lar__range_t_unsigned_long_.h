//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Tue Jun  4 02:53:17 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /pnfs/dune/scratch/users/tjyang/v08_21_00/reco/mcc12_protodune_beam_p1GeV_cosmics_12ms_test/19285494_84/mcc12_protodune_beam_p1GeV_cosmics_12ms_test_84_20190601T145951_reco0.root
//////////////////////////////////////////////////////////


#ifndef lar__range_t_unsigned_long__h
#define lar__range_t_unsigned_long__h
namespace lar {
template <typename T0> class range_t;
template <> class range_t<unsigned long>;
} // end of namespace.

#include "Rtypes.h"

namespace lar {
template <typename T0> class range_t;
template <> class range_t<unsigned long> {

public:
// Nested classes declaration.

public:
// Data Members.
   unsigned long offset;      //offset (absolute index) of the first element
   unsigned long last;        //offset (absolute index) after the last element

   range_t();
   range_t(const range_t & );
   virtual ~range_t();

};
} // namespace
#endif
