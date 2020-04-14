//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Thu Feb 13 18:14:53 2020 by ROOT version 6.16/00)
//      from the StreamerInfo in file ../PDSPProd2_protoDUNE_sp_reco_35ms_sce_off_23473772_0_452d9f89-a2a1-4680-ab72-853a3261da5d.root
//////////////////////////////////////////////////////////


#ifndef sim__OpDet_Time_Chans_h
#define sim__OpDet_Time_Chans_h
namespace sim {
class OpDet_Time_Chans;
} // end of namespace.

#include "Rtypes.h"
#include "Riostream.h"
#include <vector>
#include "sim__Chan_Phot.h"

namespace sim {
class OpDet_Time_Chans {

public:
// Nested classes declaration.

public:
// Data Members.
   double      time;        //
   vector<sim::Chan_Phot> phots;       //

   OpDet_Time_Chans();
   OpDet_Time_Chans(const OpDet_Time_Chans & );
   virtual ~OpDet_Time_Chans();

};
} // namespace
#endif
