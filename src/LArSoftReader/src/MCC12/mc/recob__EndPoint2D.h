//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Tue Jun  4 02:53:16 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /pnfs/dune/scratch/users/tjyang/v08_21_00/reco/mcc12_protodune_beam_p1GeV_cosmics_12ms_test/19285494_84/mcc12_protodune_beam_p1GeV_cosmics_12ms_test_84_20190601T145951_reco0.root
//////////////////////////////////////////////////////////


#ifndef recob__EndPoint2D_h
#define recob__EndPoint2D_h
namespace recob {
class EndPoint2D;
} // end of namespace.

#include "Rtypes.h"
#include "geo__WireID.h"

namespace recob {
class EndPoint2D {

public:
// Nested classes declaration.

public:
// Data Members.
   double      fDriftTime;    //vertex's drift time
   geo::WireID fWireID;       //vertex's wireID
   int         fID;           //vertex's ID
   double      fStrength;     //vertex's strength
   Int_t       fView;         //view for this end point
   double      fTotalCharge;    //total charge of hits associated with end point

   EndPoint2D();
   EndPoint2D(const EndPoint2D & );
   virtual ~EndPoint2D();

   ClassDef(EndPoint2D,15); // Generated by MakeProject.
};
} // namespace
#endif
