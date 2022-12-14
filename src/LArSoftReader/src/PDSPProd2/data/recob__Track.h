//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Thu Feb 13 18:08:35 2020 by ROOT version 6.16/00)
//      from the StreamerInfo in file ../np04_raw_run005834_0126_dl7_reco_22481932_0_20190828T190321.root
//////////////////////////////////////////////////////////


#ifndef recob__Track_h
#define recob__Track_h
namespace recob {
class Track;
} // end of namespace.

#include "Rtypes.h"
#include "recob__TrackTrajectory.h"
#include "Math/SMatrix.h"
#include "Math/MatrixRepresentationsStatic.h"

namespace recob {
class Track {

public:
// Nested classes declaration.

public:
// Data Members.
   recob::TrackTrajectory fTraj;       //Stored trajectory data member
   int                    fPId;        //Particle ID hypothesis used in the fit (if any)
   float                  fChi2;       //Fit chi2
   int                    fNdof;       //Number of degrees of freedom of the fit
   ROOT::Math::SMatrix<Double32_t,5,5,ROOT::Math::MatRepSym<Double32_t,5> > fCovVertex;    //Covariance matrix (local 5D representation) at start point (vertex)
   ROOT::Math::SMatrix<Double32_t,5,5,ROOT::Math::MatRepSym<Double32_t,5> > fCovEnd;       //Covariance matrix (local 5D representation) at end point
   int                                                                      fID;           //track's ID

   Track();
   Track(const Track & );
   virtual ~Track();

   ClassDef(Track,18); // Generated by MakeProject.
};
} // namespace
#endif
