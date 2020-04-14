//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Tue Jun  4 02:53:16 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /pnfs/dune/scratch/users/tjyang/v08_21_00/reco/mcc12_protodune_beam_p1GeV_cosmics_12ms_test/19285494_84/mcc12_protodune_beam_p1GeV_cosmics_12ms_test_84_20190601T145951_reco0.root
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
   ResultsAuxiliary(const ResultsAuxiliary & );
   virtual ~ResultsAuxiliary();

   ClassDef(ResultsAuxiliary,11); // Generated by MakeProject.
};
} // namespace
#endif
