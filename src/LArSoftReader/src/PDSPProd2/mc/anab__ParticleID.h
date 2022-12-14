//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Fri Apr 24 09:37:06 2020 by ROOT version 6.18/04)
//      from the StreamerInfo in file /pnfs/dune/archive/sam_managed_users/calcuttj/data/5/c/5/d/PDSPProd2_protoDUNE_sp_reco_35ms_sce_off_23503706_0_124d8226-6d53-4dbc-af94-340521eab19b_pma.root
//////////////////////////////////////////////////////////


#ifndef anab__ParticleID_h
#define anab__ParticleID_h
namespace anab {
class ParticleID;
} // end of namespace.

#include "Rtypes.h"
#include "geo__PlaneID.h"

namespace anab {
class ParticleID {

public:
// Nested classes declaration.

public:
// Data Members.
   int         fPdg;        //determined particle ID
   int         fNdf;        //ndf for chi2 test
   double      fMinChi2;    //Minimum reduced chi2
   double      fDeltaChi2;    //difference between two lowest reduced chi2's
   double      fChi2Proton;    //reduced chi2 using proton template
   double      fChi2Kaon;      //reduced chi2 using kaon template
   double      fChi2Pion;      //reduced chi2 using pion template
   double      fChi2Muon;      //reduced chi2 using muon template
   double      fMissingE;      //missing energy from dead wires for contained particle
   double      fMissingEavg;    //missing energy from dead wires using average dEdx
   double      fPIDA;           //PID developed by Bruce Baller
   geo::PlaneID fPlaneID;        //

   ParticleID();
   ParticleID(ParticleID && ) = default;
   ParticleID(const ParticleID & );
   virtual ~ParticleID();

   ClassDef(ParticleID,14); // Generated by MakeProject.
};
} // namespace
#endif
