//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Mon Jun  3 14:12:49 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/ana_-043d807f-0b76-401e-9889-c1184b96a71c.root
//////////////////////////////////////////////////////////


#ifndef art__FileFormatVersion_h
#define art__FileFormatVersion_h
namespace art {
class FileFormatVersion;
} // end of namespace.

#include "Rtypes.h"
#include "Riostream.h"
#include <string>

namespace art {
class FileFormatVersion {

public:
// Nested classes declaration.

public:
// Data Members.
   int         value_;      //
   string      era_;        //

   FileFormatVersion();
   FileFormatVersion(const FileFormatVersion & );
   virtual ~FileFormatVersion();

   ClassDef(FileFormatVersion,11); // Generated by MakeProject.
};
} // namespace
#endif
