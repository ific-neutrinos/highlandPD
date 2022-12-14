//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Sat Mar  9 10:19:31 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file ../../../protoDuneExampleAnalysis/v0r0/cmt/np04_raw_run005250_0014_dl9_reco_13061912_0_20181110T220718.root
//////////////////////////////////////////////////////////


#ifndef raw__RawDigit_h
#define raw__RawDigit_h
namespace raw {
class RawDigit;
} // end of namespace.

#include "Rtypes.h"
#include "Riostream.h"
#include <vector>

namespace raw {
class RawDigit {

public:
// Nested classes declaration.

public:
// Data Members.
   vector<short> fADC;        //ADC readout per tick, before pedestal subtraction
   unsigned int  fChannel;    //channel number in the readout
   unsigned short fSamples;    //number of ticks of the clock
   float          fPedestal;    //pedestal for this channel
   float          fSigma;       //sigma of the pedestal counts for this channel
   Int_t          fCompression;    //compression scheme used for the ADC vector

   RawDigit();
   RawDigit(const RawDigit & );
   virtual ~RawDigit();

   ClassDef(RawDigit,15); // Generated by MakeProject.
};
} // namespace
#endif
