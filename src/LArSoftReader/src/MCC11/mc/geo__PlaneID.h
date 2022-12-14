//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Thu Feb 14 10:40:59 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/MCC11/mcc11_protoDUNE_sp_reco_12231045_0_1ee4e6cf-25fb-482c-a231-fbd91e816aa4.root
//////////////////////////////////////////////////////////


#ifndef geo__PlaneID_h
#define geo__PlaneID_h
namespace geo {
class PlaneID;
} // end of namespace.

#include "Rtypes.h"
#include "geo__TPCID.h"

namespace geo {
class PlaneID : public geo::TPCID {

public:
// Nested classes declaration.

public:
// Data Members.
   unsigned int Plane;       //Index of the plane within its TPC.

   PlaneID();
   PlaneID(const PlaneID & );
   virtual ~PlaneID();

   ClassDef(PlaneID,4); // Generated by MakeProject.
};
} // namespace
#endif
