//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Sat Mar  9 10:19:31 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file ../../../protoDuneExampleAnalysis/v0r0/cmt/np04_raw_run005250_0014_dl9_reco_13061912_0_20181110T220718.root
//////////////////////////////////////////////////////////


#ifndef art__Transient_art__ProductProvenance__Transients__h
#define art__Transient_art__ProductProvenance__Transients__h
namespace art {
#ifndef art__ProductProvenance_h
#define art__ProductProvenance_h
enum ProductProvenance::Transients { kDefault_ProductProvenance::Transients };
#endif
} // end of namespace.
namespace art {
template <typename T0> class Transient;
} // end of namespace.

#include "Rtypes.h"
namespace art {
#ifndef art__ProductProvenance_h
#define art__ProductProvenance_h
enum ProductProvenance::Transients { kDefault_ProductProvenance::Transients };
#endif
} // end of namespace.

#ifndef template_art__Transient_h
#define template_art__Transient_h
namespace art {
template <typename T0> class Transient {

public:
// Nested classes declaration.

public:
// Data Members.

   Transient() {
   }
   Transient(const Transient & rhs )
{
   // This is NOT a copy constructor. This is actually a move constructor (for stl container's sake).
   // Use at your own risk!
   (void)rhs; // avoid warning about unused parameter
   }
   virtual ~Transient() {
   }

};
} // namespace
#endif // generic template declaration
#endif
