//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Thu Feb 14 10:40:59 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/MCC11/mcc11_protoDUNE_sp_reco_12231045_0_1ee4e6cf-25fb-482c-a231-fbd91e816aa4.root
//////////////////////////////////////////////////////////


#ifndef art__Hash_2__h
#define art__Hash_2__h
namespace art {
template <int T0> class Hash;
template <> class Hash<2>;
} // end of namespace.

#include "Rtypes.h"
#include "Riostream.h"
#include <string>

namespace art {
template <int T0> class Hash;
template <> class Hash<2> {

public:
// Nested classes declaration.

public:
// Data Members.
   string      hash_;       //

   Hash();
   Hash(const Hash & );
   virtual ~Hash();

   ClassDef(Hash,11); // Generated by MakeProject.
};
} // namespace
#endif
