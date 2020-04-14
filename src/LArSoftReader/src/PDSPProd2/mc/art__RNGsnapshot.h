//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Thu Feb 13 18:14:53 2020 by ROOT version 6.16/00)
//      from the StreamerInfo in file ../PDSPProd2_protoDUNE_sp_reco_35ms_sce_off_23473772_0_452d9f89-a2a1-4680-ab72-853a3261da5d.root
//////////////////////////////////////////////////////////


#ifndef art__RNGsnapshot_h
#define art__RNGsnapshot_h
namespace art {
class RNGsnapshot;
} // end of namespace.

#include "Rtypes.h"
#include "Riostream.h"
#include <string>
#include <vector>

namespace art {
class RNGsnapshot {

public:
// Nested classes declaration.

public:
// Data Members.
   string      engine_kind_;    //
   string      label_;          //
   vector<unsigned int> state_;          //

   RNGsnapshot();
   RNGsnapshot(const RNGsnapshot & );
   virtual ~RNGsnapshot();

   ClassDef(RNGsnapshot,11); // Generated by MakeProject.
};
} // namespace
#endif
