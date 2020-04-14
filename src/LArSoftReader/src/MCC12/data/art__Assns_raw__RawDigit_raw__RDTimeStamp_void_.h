//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Mon Jun  3 14:12:49 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/ana_-043d807f-0b76-401e-9889-c1184b96a71c.root
//////////////////////////////////////////////////////////


#ifndef art__Assns_raw__RawDigit_raw__RDTimeStamp_void__h
#define art__Assns_raw__RawDigit_raw__RDTimeStamp_void__h
namespace raw {
#ifndef raw__RawDigit_h
#define raw__RawDigit_h
enum RawDigit { kDefault_RawDigit };
#endif
} // end of namespace.
namespace raw {
class RDTimeStamp;
} // end of namespace.
namespace art {
template <typename T0, typename T1, typename T2> class Assns;
template <> class Assns<raw::RawDigit,raw::RDTimeStamp,void>;
} // end of namespace.

#include "Rtypes.h"
namespace raw {
#ifndef raw__RawDigit_h
#define raw__RawDigit_h
enum RawDigit { kDefault_RawDigit };
#endif
} // end of namespace.
#include "raw__RDTimeStamp.h"
#include "art__detail__AssnsBase.h"
#include "Riostream.h"
#include <vector>
#include <utility>
#include "art__RefCore.h"

namespace art {
template <typename T0, typename T1, typename T2> class Assns;
template <> class Assns<raw::RawDigit,raw::RDTimeStamp,void> : public art::detail::AssnsBase {

public:
// Nested classes declaration.

public:
// Data Members.
   vector<pair<art::RefCore,unsigned long> > ptr_data_1_;    //
   vector<pair<art::RefCore,unsigned long> > ptr_data_2_;    //

   Assns();
   Assns(const Assns & );
   virtual ~Assns();

   ClassDef(Assns,12); // Generated by MakeProject.
};
} // namespace
#ifdef __MAKECINT__
#pragma link C++ class pair<art::RefCore,unsigned long>+;
#endif
#endif
