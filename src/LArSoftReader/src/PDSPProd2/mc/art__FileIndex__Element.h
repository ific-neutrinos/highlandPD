//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Fri Apr 24 09:37:05 2020 by ROOT version 6.18/04)
//      from the StreamerInfo in file /pnfs/dune/archive/sam_managed_users/calcuttj/data/5/c/5/d/PDSPProd2_protoDUNE_sp_reco_35ms_sce_off_23503706_0_124d8226-6d53-4dbc-af94-340521eab19b_pma.root
//////////////////////////////////////////////////////////


#ifndef art__FileIndex__Element_h
#define art__FileIndex__Element_h
namespace art {
namespace FileIndex {
class Element;
} // end of namespace.
} // end of namespace.

#include "Rtypes.h"
#include "art__EventID.h"

namespace art {
namespace FileIndex {
class Element {

public:
// Nested classes declaration.

public:
// Data Members.
   art::EventID eventID_;    //
   Long64_t     entry_;      //

   Element();
   Element(Element && ) = default;
   Element(const Element & );
   virtual ~Element();

   ClassDef(Element,12); // Generated by MakeProject.
};
} // namespace
} // namespace
#endif
