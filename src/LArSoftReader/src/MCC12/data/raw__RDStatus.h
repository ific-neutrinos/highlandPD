//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Mon Jun  3 14:12:50 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/ana_-043d807f-0b76-401e-9889-c1184b96a71c.root
//////////////////////////////////////////////////////////


#ifndef raw__RDStatus_h
#define raw__RDStatus_h
namespace raw {
class RDStatus;
} // end of namespace.

#include "Rtypes.h"

namespace raw {
class RDStatus {

public:
// Nested classes declaration.

public:
// Data Members.
   bool        fCorruptDataDropped;    //True if some data suspected to be corrupt were not decoded
   bool        fCorruptDataKept;       //True if some data suspected to be corrupt were decoded
   unsigned int fStatWord;              //

   RDStatus();
   RDStatus(const RDStatus & );
   virtual ~RDStatus();

   ClassDef(RDStatus,11); // Generated by MakeProject.
};
} // namespace
#endif
