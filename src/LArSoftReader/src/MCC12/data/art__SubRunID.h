//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Mon Jun  3 14:12:49 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/ana_-043d807f-0b76-401e-9889-c1184b96a71c.root
//////////////////////////////////////////////////////////


#ifndef art__SubRunID_h
#define art__SubRunID_h
namespace art {
class SubRunID;
} // end of namespace.

#include "Rtypes.h"
#include "art__RunID.h"

namespace art {
class SubRunID {

public:
// Nested classes declaration.

public:
// Data Members.
   art::RunID  run_;        //
   unsigned int subRun_;     //

   SubRunID();
   SubRunID(const SubRunID & );
   virtual ~SubRunID();

   ClassDef(SubRunID,11); // Generated by MakeProject.
};
} // namespace
#endif
