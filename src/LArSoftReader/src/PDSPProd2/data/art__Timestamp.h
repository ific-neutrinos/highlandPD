//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Thu Feb 13 18:08:35 2020 by ROOT version 6.16/00)
//      from the StreamerInfo in file ../np04_raw_run005834_0126_dl7_reco_22481932_0_20190828T190321.root
//////////////////////////////////////////////////////////


#ifndef art__Timestamp_h
#define art__Timestamp_h
namespace art {
class Timestamp;
} // end of namespace.

#include "Rtypes.h"

namespace art {
class Timestamp {

public:
// Nested classes declaration.

public:
// Data Members.
   unsigned int timeLow_;    //
   unsigned int timeHigh_;    //

   Timestamp();
   Timestamp(const Timestamp & );
   virtual ~Timestamp();

   ClassDef(Timestamp,11); // Generated by MakeProject.
};
} // namespace
#endif
