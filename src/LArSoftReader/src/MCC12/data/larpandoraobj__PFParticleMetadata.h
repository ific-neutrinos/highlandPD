//////////////////////////////////////////////////////////
//   This class has been generated by TFile::MakeProject
//     (Mon Jun  3 14:12:50 2019 by ROOT version 6.12/06)
//      from the StreamerInfo in file /hep/DUNE/DataDir/ana_-043d807f-0b76-401e-9889-c1184b96a71c.root
//////////////////////////////////////////////////////////


#ifndef larpandoraobj__PFParticleMetadata_h
#define larpandoraobj__PFParticleMetadata_h
namespace larpandoraobj {
class PFParticleMetadata;
} // end of namespace.

#include "Rtypes.h"
#include "Riostream.h"
#include <map>

namespace larpandoraobj {
class PFParticleMetadata {

public:
// Nested classes declaration.

public:
// Data Members.
//   std::vector<std::pair<std::string,float> > m_propertiesMap;    //The properties map
   std::map<std::string,float > m_propertiesMap;    //The properties map
   
   PFParticleMetadata();
   PFParticleMetadata(const PFParticleMetadata & );
   virtual ~PFParticleMetadata();

   ClassDef(PFParticleMetadata,13); // Generated by MakeProject.
};
} // namespace
#endif
