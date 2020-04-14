//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Mon Jun  3 14:12:49 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/ana_-043d807f-0b76-401e-9889-c1184b96a71c.root
//////////////////////////////////////////////////////////


#ifndef art__ProductProvenance_h
#define art__ProductProvenance_h
namespace art {
class ProductProvenance;
} // end of namespace.

#include "Rtypes.h"
#include "art__ProductID.h"
#include "art__Hash_5_.h"
#include "art__Transient_art__ProductProvenance__Transients_.h"
namespace art {
#ifndef art__ProductProvenance_h
#define art__ProductProvenance_h
enum ProductProvenance::Transients { kDefault_ProductProvenance::Transients };
#endif
} // end of namespace.

namespace art {
class ProductProvenance {

public:
// Nested classes forward declaration.
enum Transients { kDefault_Transients };

public:
// Nested classes declaration.

public:
// Data Members.
   art::ProductID productID_;    //
   unsigned char  productStatus_;    //
   art::Hash<5>   parentageID_;      //
   art::Transient<art::ProductProvenance::Transients> transients_;       //

   ProductProvenance();
   ProductProvenance(const ProductProvenance & );
   virtual ~ProductProvenance();

   ClassDef(ProductProvenance,12); // Generated by MakeProject.
};
} // namespace
#endif
