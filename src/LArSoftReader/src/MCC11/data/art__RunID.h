//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Sat Mar  9 10:19:31 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file ../../../protoDuneExampleAnalysis/v0r0/cmt/np04_raw_run005250_0014_dl9_reco_13061912_0_20181110T220718.root
//////////////////////////////////////////////////////////


#ifndef art__RunID_h
#define art__RunID_h
namespace art {
class RunID;
} // end of namespace.

#include "Rtypes.h"

namespace art {
class RunID {

public:
// Nested classes declaration.

public:
// Data Members.
   unsigned int run_;        //

   RunID();
   RunID(const RunID & );
   virtual ~RunID();

   ClassDef(RunID,11); // Generated by MakeProject.
};
} // namespace
#endif
