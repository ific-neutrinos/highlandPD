//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Thu Feb 14 10:40:59 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/MCC11/mcc11_protoDUNE_sp_reco_12231045_0_1ee4e6cf-25fb-482c-a231-fbd91e816aa4.root
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
