//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Thu Feb 13 18:08:35 2020 by ROOT version 6.16/00)
//      from the StreamerInfo in file ../np04_raw_run005834_0126_dl7_reco_22481932_0_20190828T190321.root
//////////////////////////////////////////////////////////


#ifndef art__EventAuxiliary_h
#define art__EventAuxiliary_h
namespace art {
class EventAuxiliary;
} // end of namespace.

#include "Rtypes.h"
#include "art__EventID.h"
#include "art__Timestamp.h"

namespace art {
class EventAuxiliary {

public:
// Nested classes forward declaration.
enum ExperimentType { kDefault_ExperimentType };

public:
// Nested classes declaration.

public:
// Data Members.
   art::EventID id_;         //
   art::Timestamp time_;       //
   bool           isRealData_;    //
   Int_t          experimentType_;    //

   EventAuxiliary();
   EventAuxiliary(const EventAuxiliary & );
   virtual ~EventAuxiliary();

   ClassDef(EventAuxiliary,12); // Generated by MakeProject.
};
} // namespace
#endif
