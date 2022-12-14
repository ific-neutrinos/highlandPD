//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Thu Feb 14 10:40:59 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/MCC11/mcc11_protoDUNE_sp_reco_12231045_0_1ee4e6cf-25fb-482c-a231-fbd91e816aa4.root
//////////////////////////////////////////////////////////


#ifndef recob__Wire_h
#define recob__Wire_h
namespace recob {
class Wire;
} // end of namespace.

#include "Rtypes.h"
#include "lar__sparse_vector_float_.h"

namespace recob {
class Wire {

public:
// Nested classes declaration.

public:
// Data Members.
   unsigned int fChannel;    //ID of the associated channel
   Int_t        fView;       //View corresponding to the plane of this wire
   lar::sparse_vector<float> fSignalROI;    //Signal on the channel

   Wire();
   Wire(const Wire & );
   virtual ~Wire();

   ClassDef(Wire,20); // Generated by MakeProject.
};
} // namespace
#endif
