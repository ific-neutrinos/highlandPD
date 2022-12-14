//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Thu Feb 14 10:40:59 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/MCC11/mcc11_protoDUNE_sp_reco_12231045_0_1ee4e6cf-25fb-482c-a231-fbd91e816aa4.root
//////////////////////////////////////////////////////////


#ifndef simb__MCTruth_h
#define simb__MCTruth_h
namespace simb {
class MCTruth;
} // end of namespace.

#include "Rtypes.h"
#include "Riostream.h"
#include <vector>
#include "simb__MCParticle.h"
#include "simb__MCNeutrino.h"

namespace simb {
class MCTruth {

public:
// Nested classes declaration.

public:
// Data Members.
   vector<simb::MCParticle> fPartList;    //list of particles in this event
   simb::MCNeutrino         fMCNeutrino;    //reference to neutrino info - null if not a neutrino
   Int_t                    fOrigin;        //origin for this event
   bool                     fNeutrinoSet;    //flag for whether the neutrino information has been set

   MCTruth();
   MCTruth(const MCTruth & );
   virtual ~MCTruth();

   ClassDef(MCTruth,12); // Generated by MakeProject.
};
} // namespace
#endif
