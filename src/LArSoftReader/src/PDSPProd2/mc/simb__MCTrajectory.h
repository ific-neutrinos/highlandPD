//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Fri Apr 24 09:37:05 2020 by ROOT version 6.18/04)
//      from the StreamerInfo in file /pnfs/dune/archive/sam_managed_users/calcuttj/data/5/c/5/d/PDSPProd2_protoDUNE_sp_reco_35ms_sce_off_23503706_0_124d8226-6d53-4dbc-af94-340521eab19b_pma.root
//////////////////////////////////////////////////////////


#ifndef simb__MCTrajectory_h
#define simb__MCTrajectory_h
namespace simb {
class MCTrajectory;
} // end of namespace.

#include "Rtypes.h"
#include "Riostream.h"
#include <vector>
#include <utility>
#include "TLorentzVector.h"

namespace simb {
class MCTrajectory {

public:
// Nested classes declaration.

public:
// Data Members.
   vector<pair<TLorentzVector,TLorentzVector> > ftrajectory;    //The list of trajectory points
   vector<pair<unsigned long,unsigned char> >   fTrajectoryProcess;    //map of the scattering process to index

   MCTrajectory();
   MCTrajectory(MCTrajectory && ) = default;
   MCTrajectory(const MCTrajectory & );
   virtual ~MCTrajectory();

   ClassDef(MCTrajectory,13); // Generated by MakeProject.
};
} // namespace
#endif
