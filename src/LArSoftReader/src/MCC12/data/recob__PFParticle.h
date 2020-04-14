//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Mon Jun  3 14:12:50 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/ana_-043d807f-0b76-401e-9889-c1184b96a71c.root
//////////////////////////////////////////////////////////


#ifndef recob__PFParticle_h
#define recob__PFParticle_h
namespace recob {
class PFParticle;
} // end of namespace.

#include "Rtypes.h"
#include "Riostream.h"
#include <vector>

namespace recob {
class PFParticle {

public:
// Nested classes declaration.

public:
// Data Members.
   int         fPdgCode;    //A preliminary estimate of the PFParticle type using the PDG code
   unsigned long fSelf;       //Self reference
   unsigned long fParent;     //Index into PFParticle collection for parent
   vector<unsigned long> fDaughters;    //Vector of indices into PFParticle Collection for daughters

   PFParticle();
   PFParticle(const PFParticle & );
   virtual ~PFParticle();

   ClassDef(PFParticle,14); // Generated by MakeProject.
};
} // namespace
#endif
