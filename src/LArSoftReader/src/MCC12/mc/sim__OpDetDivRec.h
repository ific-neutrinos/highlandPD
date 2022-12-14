//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Tue Jun  4 02:53:16 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /pnfs/dune/scratch/users/tjyang/v08_21_00/reco/mcc12_protodune_beam_p1GeV_cosmics_12ms_test/19285494_84/mcc12_protodune_beam_p1GeV_cosmics_12ms_test_84_20190601T145951_reco0.root
//////////////////////////////////////////////////////////


#ifndef sim__OpDetDivRec_h
#define sim__OpDetDivRec_h
namespace sim {
class OpDetDivRec;
} // end of namespace.

#include "Rtypes.h"
#include "Riostream.h"
#include <vector>
#include "sim__OpDet_Time_Chans.h"

namespace sim {
class OpDetDivRec {

public:
// Nested classes declaration.

public:
// Data Members.
   int         fOpDetNum;    //Move this to be private.
   vector<sim::OpDet_Time_Chans> time_chans;    //Move this to private //This must be filled with emplace to keep it sorted.

   OpDetDivRec();
   OpDetDivRec(const OpDetDivRec & );
   virtual ~OpDetDivRec();

   ClassDef(OpDetDivRec,14); // Generated by MakeProject.
};
} // namespace
#endif
