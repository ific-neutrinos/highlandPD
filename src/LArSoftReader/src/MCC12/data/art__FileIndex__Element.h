//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Mon Jun  3 14:12:49 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/ana_-043d807f-0b76-401e-9889-c1184b96a71c.root
//////////////////////////////////////////////////////////


#ifndef art__FileIndex__Element_h
#define art__FileIndex__Element_h
namespace art {
namespace FileIndex {
class Element;
} // end of namespace.
} // end of namespace.

#include "Rtypes.h"
#include "art__EventID.h"

namespace art {
namespace FileIndex {
class Element {

public:
// Nested classes declaration.

public:
// Data Members.
   art::EventID eventID_;    //
   Long64_t     entry_;      //

   Element();
   Element(const Element & );
   virtual ~Element();

   ClassDef(Element,12); // Generated by MakeProject.
};
} // namespace
} // namespace
#endif
