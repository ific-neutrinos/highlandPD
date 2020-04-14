//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Tue Jun  4 02:53:16 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /pnfs/dune/scratch/users/tjyang/v08_21_00/reco/mcc12_protodune_beam_p1GeV_cosmics_12ms_test/19285494_84/mcc12_protodune_beam_p1GeV_cosmics_12ms_test_84_20190601T145951_reco0.root
//////////////////////////////////////////////////////////


#ifndef art__BranchDescription_h
#define art__BranchDescription_h
namespace art {
class BranchDescription;
} // end of namespace.

#include "Rtypes.h"
#include "Riostream.h"
#include <string>
#include "art__ProductID.h"
#include <set>
#include "fhicl__ParameterSetID.h"
#include "art__Hash_3_.h"
#include "art__Transient_art__BranchDescription__Transients_.h"
namespace art {
#ifndef art__BranchDescription_h
#define art__BranchDescription_h
enum BranchDescription::Transients { kDefault_BranchDescription::Transients };
#endif
} // end of namespace.

namespace art {
class BranchDescription {

public:
// Nested classes forward declaration.
enum Transients { kDefault_Transients };

public:
// Nested classes declaration.

public:
// Data Members.
   Int_t       branchType_;    //
   string      moduleLabel_;    //
   string      processName_;    //
   art::ProductID productID_;      //
   string         producedClassName_;    //
   string         friendlyClassName_;    //
   string         productInstanceName_;    //
   bool           supportsView_;           //
   std::vector<fhicl::ParameterSetID> psetIDs_;                //
   std::vector<art::Hash<3> >         processConfigurationIDs_;    //
   art::Transient<art::BranchDescription::Transients> transients_;                 //

   BranchDescription();
   BranchDescription(const BranchDescription & );
   virtual ~BranchDescription();

   ClassDef(BranchDescription,16); // Generated by MakeProject.
};
} // namespace
#endif
