//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Fri Apr 24 09:37:05 2020 by ROOT version 6.18/04)
//      from the StreamerInfo in file /pnfs/dune/archive/sam_managed_users/calcuttj/data/5/c/5/d/PDSPProd2_protoDUNE_sp_reco_35ms_sce_off_23503706_0_124d8226-6d53-4dbc-af94-340521eab19b_pma.root
//////////////////////////////////////////////////////////


#ifndef art__ResultsAuxiliary_h
#define art__ResultsAuxiliary_h
namespace art {
class ResultsAuxiliary;
} // end of namespace.

#include "Rtypes.h"
#include "art__Hash_2_.h"
#include "Riostream.h"
#include <set>

namespace art {
class ResultsAuxiliary {

public:
// Nested classes declaration.

public:
// Data Members.
   art::Hash<2> processHistoryID_;    //
   std::vector<art::Hash<2> > allEventsProcessHistories_;    //

   ResultsAuxiliary();
   ResultsAuxiliary(ResultsAuxiliary && ) = default;
   ResultsAuxiliary(const ResultsAuxiliary & );
   virtual ~ResultsAuxiliary();

   ClassDef(ResultsAuxiliary,11); // Generated by MakeProject.
};
} // namespace
#endif
