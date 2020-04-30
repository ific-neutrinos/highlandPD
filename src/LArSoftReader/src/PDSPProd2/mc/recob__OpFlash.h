//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Fri Apr 24 09:37:06 2020 by ROOT version 6.18/04)
//      from the StreamerInfo in file /pnfs/dune/archive/sam_managed_users/calcuttj/data/5/c/5/d/PDSPProd2_protoDUNE_sp_reco_35ms_sce_off_23503706_0_124d8226-6d53-4dbc-af94-340521eab19b_pma.root
//////////////////////////////////////////////////////////


#ifndef recob__OpFlash_h
#define recob__OpFlash_h
namespace recob {
class OpFlash;
} // end of namespace.

#include "Rtypes.h"
#include "Riostream.h"
#include <vector>

namespace recob {
class OpFlash {

public:
// Nested classes declaration.

public:
// Data Members.
   double      fTime;       //Time relative to trigger
   double      fTimeWidth;    //Width of the flash in time
   double      fAbsTime;      //Time by PMT readout clock
   unsigned int fFrame;        //Frame number
   vector<double> fPEperOpDet;    //Number of PE on each PMT
   vector<double> fWireCenters;    //Geometric center in each view
   vector<double> fWireWidths;     //Geometric width in each view
   double         fYCenter;        //Geometric center in y
   double         fYWidth;         //Geometric width in y
   double         fZCenter;        //Geometric center in z
   double         fZWidth;         //Geometric width in z
   double         fFastToTotal;    //Fast to total light ratio
   bool           fInBeamFrame;    //Is this in the beam frame?
   int            fOnBeamTime;     //Is this in time with beam?

   OpFlash();
   OpFlash(OpFlash && ) = default;
   OpFlash(const OpFlash & );
   virtual ~OpFlash();

   ClassDef(OpFlash,19); // Generated by MakeProject.
};
} // namespace
#endif
