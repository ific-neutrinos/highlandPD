//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Tue Jun  4 02:53:16 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /pnfs/dune/scratch/users/tjyang/v08_21_00/reco/mcc12_protodune_beam_p1GeV_cosmics_12ms_test/19285494_84/mcc12_protodune_beam_p1GeV_cosmics_12ms_test_84_20190601T145951_reco0.root
//////////////////////////////////////////////////////////


#ifndef art__Wrapper_art__Assns_simb__MCTruth_simb__MCParticle_sim__GeneratedParticleInfo____h
#define art__Wrapper_art__Assns_simb__MCTruth_simb__MCParticle_sim__GeneratedParticleInfo____h
namespace simb {
class MCTruth;
} // end of namespace.
namespace simb {
class MCParticle;
} // end of namespace.
namespace sim {
class GeneratedParticleInfo;
} // end of namespace.
namespace art {
template <typename T0, typename T1, typename T2> class Assns;
template <> class Assns<simb::MCTruth,simb::MCParticle,sim::GeneratedParticleInfo>;
} // end of namespace.
namespace art {
template <typename T0> class Wrapper;
template <> class Wrapper<art::Assns<simb::MCTruth,simb::MCParticle,sim::GeneratedParticleInfo> >;
} // end of namespace.

#include "Rtypes.h"
#include "art__Assns_simb__MCTruth_simb__MCParticle_sim__GeneratedParticleInfo_.h"
#include "art__EDProduct.h"
#include "simb__MCTruth.h"
#include "simb__MCParticle.h"
#include "sim__GeneratedParticleInfo.h"

namespace art {
template <typename T0> class Wrapper;
template <> class Wrapper<art::Assns<simb::MCTruth,simb::MCParticle,sim::GeneratedParticleInfo> > : public art::EDProduct {

public:
// Nested classes declaration.

public:
// Data Members.
   bool        present;     //
   unsigned int rangeSetID;    //
   art::Assns<simb::MCTruth,simb::MCParticle,sim::GeneratedParticleInfo> obj;           //

   Wrapper();
   Wrapper(const Wrapper & );
   virtual ~Wrapper();

   ClassDef(Wrapper,12); // Generated by MakeProject.
};
} // namespace
#endif
