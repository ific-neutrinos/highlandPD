//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Thu Feb 13 18:08:35 2020 by ROOT version 6.16/00)
//      from the StreamerInfo in file ../np04_raw_run005834_0126_dl7_reco_22481932_0_20190828T190321.root
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
