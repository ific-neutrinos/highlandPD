//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Thu Feb 13 18:08:35 2020 by ROOT version 6.16/00)
//      from the StreamerInfo in file ../np04_raw_run005834_0126_dl7_reco_22481932_0_20190828T190321.root
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
