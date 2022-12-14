//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Thu Feb 13 18:08:35 2020 by ROOT version 6.16/00)
//      from the StreamerInfo in file ../np04_raw_run005834_0126_dl7_reco_22481932_0_20190828T190321.root
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
