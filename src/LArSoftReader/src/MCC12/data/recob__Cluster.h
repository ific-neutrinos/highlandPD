//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Mon Jun  3 14:12:50 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/ana_-043d807f-0b76-401e-9889-c1184b96a71c.root
//////////////////////////////////////////////////////////


#ifndef recob__Cluster_h
#define recob__Cluster_h
namespace recob {
class Cluster;
} // end of namespace.

#include "Rtypes.h"
#include "geo__PlaneID.h"

namespace recob {
class Cluster {

public:
// Nested classes declaration.

public:
// Data Members.
   unsigned int fNHits;      //Number of hits in the cluster
   float        fEndWires[2];    //
   float        fSigmaEndWires[2];    //
   float        fEndTicks[2];         //
   float        fSigmaEndTicks[2];    //
   float        fEndCharges[2];       //
   float        fAngles[2];           //
   float        fOpeningAngles[2];    //
   float        fChargeSum[2];        //
   float        fChargeStdDev[2];     //
   float        fChargeAverage[2];    //
   float        fMultipleHitDensity;    //
   float        fWidth;                 //
   int          fID;                    //
   Int_t        fView;                  //View for this cluster
   geo::PlaneID fPlaneID;               //Location of the start of the cluster

   Cluster();
   Cluster(const Cluster & );
   virtual ~Cluster();

   ClassDef(Cluster,17); // Generated by MakeProject.
};
} // namespace
#endif
