//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Tue Jun  4 02:53:16 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /pnfs/dune/scratch/users/tjyang/v08_21_00/reco/mcc12_protodune_beam_p1GeV_cosmics_12ms_test/19285494_84/mcc12_protodune_beam_p1GeV_cosmics_12ms_test_84_20190601T145951_reco0.root
//////////////////////////////////////////////////////////


#ifndef recob__Vertex_h
#define recob__Vertex_h
namespace recob {
class Vertex;
} // end of namespace.

#include "Rtypes.h"
#include "Math/GenVector/PositionVector3D.h"
#include "Math/GenVector/Cartesian3D.h"
#include "Math/GenVector/CoordinateSystemTags.h"
#include "Math/SMatrix.h"
#include "Math/MatrixRepresentationsStatic.h"

namespace recob {
class Vertex {

public:
// Nested classes forward declaration.
enum Status { kDefault_Status };

public:
// Nested classes declaration.

public:
// Data Members.
   ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<Double32_t>,ROOT::Math::GlobalCoordinateSystemTag> pos_;        //Vertex 3D position
   ROOT::Math::SMatrix<Double32_t,3,3,ROOT::Math::MatRepSym<Double32_t,3> >                                cov_;        //Vertex covariance matrix 3x3
   double                                                                                                  chi2_;       //Vertex fit chi2
   int                                                                                                     ndof_;       //Vertex fit degrees of freedom
   Int_t                                                                                                   status_;     //Vertex status, as define in Vertex::Status enum
   int                                                                                                     id_;         //id number for vertex

   Vertex();
   Vertex(const Vertex & );
   virtual ~Vertex();

   ClassDef(Vertex,16); // Generated by MakeProject.
};
} // namespace
#endif
