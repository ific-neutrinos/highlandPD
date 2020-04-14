//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Mon Jun  3 14:12:50 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/ana_-043d807f-0b76-401e-9889-c1184b96a71c.root
//////////////////////////////////////////////////////////


#ifndef art__RefCore_h
#define art__RefCore_h
namespace art {
class RefCore;
} // end of namespace.

#include "Rtypes.h"
#include "art__ProductID.h"
#include "art__RefCore.h"

namespace art {
class RefCore {

public:
// Nested classes forward declaration.
class RefCoreTransients;

public:
// Nested classes declaration.
class RefCoreTransients {

public:
// Nested classes declaration.

public:
// Data Members.

   RefCoreTransients();
   RefCoreTransients(const RefCoreTransients & );
   virtual ~RefCoreTransients();

   ClassDef(RefCoreTransients,11); // Generated by MakeProject.
};

public:
// Data Members.
   art::ProductID id_;         //
   art::RefCore::RefCoreTransients transients_;    //

   RefCore();
   RefCore(const RefCore & );
   virtual ~RefCore();

   ClassDef(RefCore,11); // Generated by MakeProject.
};
} // namespace
#endif
