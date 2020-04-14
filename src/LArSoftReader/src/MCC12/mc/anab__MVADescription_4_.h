//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Tue Jun  4 02:53:16 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /pnfs/dune/scratch/users/tjyang/v08_21_00/reco/mcc12_protodune_beam_p1GeV_cosmics_12ms_test/19285494_84/mcc12_protodune_beam_p1GeV_cosmics_12ms_test_84_20190601T145951_reco0.root
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
