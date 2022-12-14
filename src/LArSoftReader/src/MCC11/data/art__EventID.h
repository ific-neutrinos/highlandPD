//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Sat Mar  9 10:19:31 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file ../../../protoDuneExampleAnalysis/v0r0/cmt/np04_raw_run005250_0014_dl9_reco_13061912_0_20181110T220718.root
//////////////////////////////////////////////////////////


#ifndef art__EventID_h
#define art__EventID_h
namespace art {
class EventID;
} // end of namespace.

#include "Rtypes.h"
#include "art__SubRunID.h"

namespace art {
class EventID {

public:
// Nested classes declaration.

public:
// Data Members.
   art::SubRunID subRun_;     //
   unsigned int  event_;      //

   EventID();
   EventID(const EventID & );
   virtual ~EventID();

   ClassDef(EventID,13); // Generated by MakeProject.
};
} // namespace
#endif
