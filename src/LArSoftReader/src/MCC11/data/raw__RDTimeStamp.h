//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Sat Mar  9 10:19:31 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file ../../../protoDuneExampleAnalysis/v0r0/cmt/np04_raw_run005250_0014_dl9_reco_13061912_0_20181110T220718.root
//////////////////////////////////////////////////////////


#ifndef raw__RDTimeStamp_h
#define raw__RDTimeStamp_h
namespace raw {
class RDTimeStamp;
} // end of namespace.

#include "Rtypes.h"

namespace raw {
class RDTimeStamp {

public:
// Nested classes declaration.

public:
// Data Members.
   ULong64_t   fTimeStamp;    //
   unsigned short fFlags;        //

   RDTimeStamp();
   RDTimeStamp(const RDTimeStamp & );
   virtual ~RDTimeStamp();

   ClassDef(RDTimeStamp,13); // Generated by MakeProject.
};
} // namespace
#endif
