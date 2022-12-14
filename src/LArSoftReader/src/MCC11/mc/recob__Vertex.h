//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Thu Feb 14 10:40:59 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/MCC11/mcc11_protoDUNE_sp_reco_12231045_0_1ee4e6cf-25fb-482c-a231-fbd91e816aa4.root
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
   ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double>,ROOT::Math::GlobalCoordinateSystemTag> pos_;        //Vertex 3D position
   ROOT::Math::SMatrix<double,3,3,ROOT::Math::MatRepSym<double,3> >                                    cov_;        //Vertex covariance matrix 3x3
   double                                                                                              chi2_;       //Vertex fit chi2
   int                                                                                                 ndof_;       //Vertex fit degrees of freedom
   Int_t                                                                                               status_;     //Vertex status, as define in Vertex::Status enum
   int                                                                                                 id_;         //id number for vertex

   Vertex();
   Vertex(const Vertex & );
   virtual ~Vertex();

   ClassDef(Vertex,15); // Generated by MakeProject.
};
} // namespace
#endif
