//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Thu Feb 14 10:40:59 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/MCC11/mcc11_protoDUNE_sp_reco_12231045_0_1ee4e6cf-25fb-482c-a231-fbd91e816aa4.root
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
   MCTrajectory(const MCTrajectory & );
   virtual ~MCTrajectory();

   ClassDef(MCTrajectory,13); // Generated by MakeProject.
};
} // namespace
#endif
