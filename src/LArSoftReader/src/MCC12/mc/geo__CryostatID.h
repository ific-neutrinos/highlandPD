//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Tue Jun  4 02:53:16 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /pnfs/dune/scratch/users/tjyang/v08_21_00/reco/mcc12_protodune_beam_p1GeV_cosmics_12ms_test/19285494_84/mcc12_protodune_beam_p1GeV_cosmics_12ms_test_84_20190601T145951_reco0.root
//////////////////////////////////////////////////////////


#ifndef geo__CryostatID_h
#define geo__CryostatID_h
namespace geo {
class CryostatID;
} // end of namespace.

#include "Rtypes.h"

namespace geo {
class CryostatID {

public:
// Nested classes declaration.

public:
// Data Members.
   bool        isValid;     //Whether this ID points to a valid element.
   unsigned int Cryostat;    //Index of cryostat.

   CryostatID();
   CryostatID(const CryostatID & );
   virtual ~CryostatID();

   ClassDef(CryostatID,19); // Generated by MakeProject.
};
} // namespace
#endif
