//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Mon Jun  3 14:12:49 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/ana_-043d807f-0b76-401e-9889-c1184b96a71c.root
//////////////////////////////////////////////////////////


#ifndef geo__TPCID_h
#define geo__TPCID_h
namespace geo {
class TPCID;
} // end of namespace.

#include "Rtypes.h"
#include "geo__CryostatID.h"

namespace geo {
class TPCID : public geo::CryostatID {

public:
// Nested classes declaration.

public:
// Data Members.
   unsigned int TPC;         //Index of the TPC within its cryostat.

   TPCID();
   TPCID(const TPCID & );
   virtual ~TPCID();

   ClassDef(TPCID,19); // Generated by MakeProject.
};
} // namespace
#endif
