//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Thu Feb 14 10:40:59 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/MCC11/mcc11_protoDUNE_sp_reco_12231045_0_1ee4e6cf-25fb-482c-a231-fbd91e816aa4.root
//////////////////////////////////////////////////////////


#ifndef util__flags__FlagSet_32_unsigned_int__h
#define util__flags__FlagSet_32_unsigned_int__h
namespace util {
namespace flags {
template <int T0, typename T1> class FlagSet;
template <> class FlagSet<32,unsigned int>;
} // end of namespace.
} // end of namespace.

#include "Rtypes.h"
#include "util__flags__BitMask_unsigned_int_.h"

namespace util {
namespace flags {
template <int T0, typename T1> class FlagSet;
template <> class FlagSet<32,unsigned int> : public util::flags::BitMask<unsigned int> {

public:
// Nested classes declaration.

public:
// Data Members.

   FlagSet();
   FlagSet(const FlagSet & );
   virtual ~FlagSet();

};
} // namespace
} // namespace
#endif
