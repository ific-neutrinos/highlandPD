//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Mon Jun  3 14:12:49 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/ana_-043d807f-0b76-401e-9889-c1184b96a71c.root
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
   ParticleID(const ParticleID & );
   virtual ~ParticleID();

   ClassDef(ParticleID,14); // Generated by MakeProject.
};
} // namespace
#endif
