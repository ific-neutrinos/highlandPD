//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Mon Jun  3 14:12:50 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/ana_-043d807f-0b76-401e-9889-c1184b96a71c.root
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
