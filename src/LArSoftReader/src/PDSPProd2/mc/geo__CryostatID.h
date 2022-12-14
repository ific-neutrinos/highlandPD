//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Fri Apr 24 09:37:06 2020 by ROOT version 6.18/04)
//      from the StreamerInfo in file /pnfs/dune/archive/sam_managed_users/calcuttj/data/5/c/5/d/PDSPProd2_protoDUNE_sp_reco_35ms_sce_off_23503706_0_124d8226-6d53-4dbc-af94-340521eab19b_pma.root
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
   CryostatID(CryostatID && ) = default;
   CryostatID(const CryostatID & );
   virtual ~CryostatID();

   ClassDef(CryostatID,19); // Generated by MakeProject.
};
} // namespace
#endif
