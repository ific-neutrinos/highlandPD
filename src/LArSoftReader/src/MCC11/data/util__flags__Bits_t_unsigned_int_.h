//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Sat Mar  9 10:19:31 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file ../../../protoDuneExampleAnalysis/v0r0/cmt/np04_raw_run005250_0014_dl9_reco_13061912_0_20181110T220718.root
//////////////////////////////////////////////////////////


#ifndef util__flags__Bits_t_unsigned_int__h
#define util__flags__Bits_t_unsigned_int__h
namespace util {
namespace flags {
template <typename T0> class Bits_t;
template <> class Bits_t<unsigned int>;
} // end of namespace.
} // end of namespace.

#include "Rtypes.h"

namespace util {
namespace flags {
template <typename T0> class Bits_t;
template <> class Bits_t<unsigned int> {

public:
// Nested classes declaration.

public:
// Data Members.
   unsigned int data;        //The bits representing all set bits

   Bits_t();
   Bits_t(const Bits_t & );
   virtual ~Bits_t();

};
} // namespace
} // namespace
#endif
