//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Thu Feb 13 18:08:35 2020 by ROOT version 6.16/00)
//      from the StreamerInfo in file ../np04_raw_run005834_0126_dl7_reco_22481932_0_20190828T190321.root
//////////////////////////////////////////////////////////


#ifndef anab__MVADescription_4__h
#define anab__MVADescription_4__h
namespace anab {
template <int T0> class MVADescription;
template <> class MVADescription<4>;
} // end of namespace.

#include "Rtypes.h"
#include "Riostream.h"
#include <string>

namespace anab {
template <int T0> class MVADescription;
template <> class MVADescription<4> {

public:
// Nested classes declaration.

public:
// Data Members.
   string      fDataTag;    //Tag of the reco data products (art::InputTag format)
   string      fOutputInstance;    //Instance name of the feature vector collection
   string      fOutputNames[4];    //Feature vector entries names/meaning

   MVADescription();
   MVADescription(const MVADescription & );
   virtual ~MVADescription();

   ClassDef(MVADescription,11); // Generated by MakeProject.
};
} // namespace
#endif
