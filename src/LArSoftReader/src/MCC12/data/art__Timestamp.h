//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Mon Jun  3 14:12:49 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/ana_-043d807f-0b76-401e-9889-c1184b96a71c.root
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
