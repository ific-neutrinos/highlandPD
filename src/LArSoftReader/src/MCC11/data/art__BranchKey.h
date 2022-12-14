//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Sat Mar  9 10:19:31 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file ../../../protoDuneExampleAnalysis/v0r0/cmt/np04_raw_run005250_0014_dl9_reco_13061912_0_20181110T220718.root
//////////////////////////////////////////////////////////


#ifndef art__BranchKey_h
#define art__BranchKey_h
namespace art {
class BranchKey;
} // end of namespace.

#include "Rtypes.h"
#include "Riostream.h"
#include <string>

namespace art {
class BranchKey {

public:
// Nested classes declaration.

public:
// Data Members.
   string      friendlyClassName_;    //
   string      moduleLabel_;          //
   string      productInstanceName_;    //
   string      processName_;            //
   int         branchType_;             //

   BranchKey();
   BranchKey(const BranchKey & );
   virtual ~BranchKey();

   ClassDef(BranchKey,12); // Generated by MakeProject.
};
} // namespace
#endif
