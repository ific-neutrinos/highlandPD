//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Tue Jun  4 02:53:16 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /pnfs/dune/scratch/users/tjyang/v08_21_00/reco/mcc12_protodune_beam_p1GeV_cosmics_12ms_test/19285494_84/mcc12_protodune_beam_p1GeV_cosmics_12ms_test_84_20190601T145951_reco0.root
//////////////////////////////////////////////////////////


#ifndef recob__PCAxis_h
#define recob__PCAxis_h
namespace recob {
class PCAxis;
} // end of namespace.

#include "Rtypes.h"
#include "Riostream.h"
#include <vector>
namespace std {} using namespace std;

namespace recob {
class PCAxis {

public:
// Nested classes declaration.

public:
// Data Members.
   bool        fSvdOK;      //SVD Decomposition was successful
   int         fNumHitsUsed;    //Number of hits in the decomposition
   double      fEigenValues[3];    //Eigen values from SVD decomposition
   vector<vector<double> > fEigenVectors;      //The three principle axes
   double                  fAvePosition[3];    //Average position of hits fed to PCA
   double                  fAveHitDoca;        //Average doca of hits used in PCA
   unsigned long           fID;                //axis ID

   PCAxis();
   PCAxis(const PCAxis & );
   virtual ~PCAxis();

   ClassDef(PCAxis,13); // Generated by MakeProject.
};
} // namespace
#endif
