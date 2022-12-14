//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Thu Feb 14 10:40:59 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/MCC11/mcc11_protoDUNE_sp_reco_12231045_0_1ee4e6cf-25fb-482c-a231-fbd91e816aa4.root
//////////////////////////////////////////////////////////


#ifndef util__flags__BitMask_unsigned_int__h
#define util__flags__BitMask_unsigned_int__h
namespace util {
namespace flags {
template <typename T0> class BitMask;
template <> class BitMask<unsigned int>;
} // end of namespace.
} // end of namespace.

#include "Rtypes.h"
#include "util__flags__Bits_t_unsigned_int_.h"

namespace util {
namespace flags {
template <typename T0> class BitMask;
template <> class BitMask<unsigned int> {

public:
// Nested classes declaration.

public:
// Data Members.
   util::flags::Bits_t<unsigned int> values;      //Storage of value bits.
   util::flags::Bits_t<unsigned int> presence;    //Storage of definition bits.

   BitMask();
   BitMask(const BitMask & );
   virtual ~BitMask();

};
} // namespace
} // namespace
#endif
