//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Fri Apr 24 09:37:05 2020 by ROOT version 6.18/04)
//      from the StreamerInfo in file /pnfs/dune/archive/sam_managed_users/calcuttj/data/5/c/5/d/PDSPProd2_protoDUNE_sp_reco_35ms_sce_off_23503706_0_124d8226-6d53-4dbc-af94-340521eab19b_pma.root
//////////////////////////////////////////////////////////


#ifndef fhicl__ParameterSetID_h
#define fhicl__ParameterSetID_h
namespace fhicl {
class ParameterSetID;
} // end of namespace.

#include "Rtypes.h"

namespace fhicl {
class ParameterSetID {

public:
// Nested classes declaration.

public:
// Data Members.
   bool        valid_;      //
   unsigned char id_[20];     //

   ParameterSetID();
   ParameterSetID(ParameterSetID && ) = default;
   ParameterSetID(const ParameterSetID & );
   virtual ~ParameterSetID();

   ClassDef(ParameterSetID,12); // Generated by MakeProject.
};
} // namespace
#endif
