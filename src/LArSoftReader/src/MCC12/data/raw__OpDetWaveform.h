//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Mon Jun  3 14:12:50 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/ana_-043d807f-0b76-401e-9889-c1184b96a71c.root
//////////////////////////////////////////////////////////


#ifndef raw__OpDetWaveform_h
#define raw__OpDetWaveform_h
namespace raw {
class OpDetWaveform;
} // end of namespace.

#include "Rtypes.h"
#include "Riostream.h"
#include <vector>

namespace raw {
class OpDetWaveform : public vector<short> {

public:
// Nested classes declaration.

public:
// Data Members.
   unsigned int fChannel;    //
   double       fTimeStamp;    //

   OpDetWaveform();
   OpDetWaveform(const OpDetWaveform & );
   virtual ~OpDetWaveform();

   ClassDef(OpDetWaveform,14); // Generated by MakeProject.
};
} // namespace
#endif
