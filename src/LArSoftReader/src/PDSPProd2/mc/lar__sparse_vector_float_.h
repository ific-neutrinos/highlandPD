//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Fri Apr 24 09:37:06 2020 by ROOT version 6.18/04)
//      from the StreamerInfo in file /pnfs/dune/archive/sam_managed_users/calcuttj/data/5/c/5/d/PDSPProd2_protoDUNE_sp_reco_35ms_sce_off_23503706_0_124d8226-6d53-4dbc-af94-340521eab19b_pma.root
//////////////////////////////////////////////////////////


#ifndef lar__sparse_vector_float__h
#define lar__sparse_vector_float__h
namespace lar {
template <typename T0> class sparse_vector;
template <> class sparse_vector<float>;
} // end of namespace.

#include "Rtypes.h"
#include "Riostream.h"
#include <vector>
#include "lar__sparse_vector_float_.h"
#include "lar__range_t_unsigned_long_.h"

namespace lar {
template <typename T0> class sparse_vector;
template <> class sparse_vector<float> {

public:
// Nested classes forward declaration.
class datarange_t;

public:
// Nested classes declaration.
class datarange_t : public lar::range_t<unsigned long> {

public:
// Nested classes declaration.

public:
// Data Members.
   vector<float> values;      //data in the range

   datarange_t();
   datarange_t(datarange_t && ) = default;
   datarange_t(const datarange_t & );
   virtual ~datarange_t();

};

public:
// Data Members.
   unsigned long nominal_size;    //current size
   vector<lar::sparse_vector<float>::datarange_t> ranges;          //list of ranges

   sparse_vector();
   sparse_vector(sparse_vector && ) = default;
   sparse_vector(const sparse_vector & );
   virtual ~sparse_vector();

};
} // namespace
#endif
