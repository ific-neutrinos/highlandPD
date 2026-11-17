// Implementación de compare_module_minitrees.C (se compila con ACLiC).
#include "pdDataClasses.hxx"
#include "TBufferFile.h"
#include "TClass.h"
#include "TFile.h"
#include "TTree.h"
#include <cmath>
#include <iostream>
#include <memory>
#include <string>

namespace {

  TTree* GetMiniTree(TFile* f){
    TTree* t = dynamic_cast<TTree*>(f->Get("highlandana/MiniTree"));
    if(!t)t = dynamic_cast<TTree*>(f->Get("MiniTree"));
    return t;
  }

  // bytes del spill serializado (incluye todo lo que cuelga de él)
  std::string Serialize(AnaSpillPD* spill){
    TBufferFile buf(TBuffer::kWrite);
    buf.WriteObjectAny(spill, TClass::GetClass("AnaSpillPD"));
    return std::string(buf.Buffer(), buf.Length());
  }

  bool Same(float a, float b){ return a == b || (std::isnan(a) && std::isnan(b)); }

  // diferencias de las partículas de una entrada, para localizar el cambio
  void PrintParticleDiffs(AnaSpillPD* a, AnaSpillPD* b, Long64_t entry){
    std::cout << "  entrada " << entry << ": " << a->Bunches.size() << " / " << b->Bunches.size() << " bunches, "
	      << a->TrueParticles.size() << " / " << b->TrueParticles.size() << " partículas true" << std::endl;
    if(a->Bunches.size() != b->Bunches.size() || a->Bunches.empty())return;
    AnaBunch* ba = static_cast<AnaBunch*>(a->Bunches[0]);
    AnaBunch* bb = static_cast<AnaBunch*>(b->Bunches[0]);
    if(ba->Particles.size() != bb->Particles.size()){
      std::cout << "    " << ba->Particles.size() << " / " << bb->Particles.size() << " partículas reco" << std::endl;
      return;
    }
    int nprint = 0;
    for(size_t i = 0; i < ba->Particles.size() && nprint < 20; i++){
      AnaParticlePD* pa = static_cast<AnaParticlePD*>(ba->Particles[i]);
      AnaParticlePD* pb = static_cast<AnaParticlePD*>(bb->Particles[i]);
      std::string what;
      if(pa->UniqueID != pb->UniqueID)                         what += " UniqueID";
      if(pa->NHits != pb->NHits)                               what += " NHits";
      if(pa->SliceID != pb->SliceID)                           what += " SliceID";
      if(pa->IsClearCosmic != pb->IsClearCosmic)               what += " IsClearCosmic";
      if(!Same(pa->PandoraTrackScore, pb->PandoraTrackScore))  what += " PandoraTrackScore";
      if(!Same(pa->T0, pb->T0))                                what += " T0";
      if(!Same(pa->BeamCosmicScore, pb->BeamCosmicScore))      what += " BeamCosmicScore";
      for(int k = 0; k < 3; k++)
	if(!Same(pa->PandoraVertex[k], pb->PandoraVertex[k])){ what += " PandoraVertex"; break; }
      if(pa->Type != pb->Type)                                 what += " Type";
      if(!Same(pa->TruePur, pb->TruePur))                      what += " TruePur";
      if(!Same(pa->TrueEff, pb->TrueEff))                      what += " TrueEff";
      const int ta = pa->TrueObject ? pa->TrueObject->ID : -1;
      const int tb = pb->TrueObject ? pb->TrueObject->ID : -1;
      if(ta != tb)                                             what += " TrueObject";
      if(!Same(pa->Chi2Proton, pb->Chi2Proton))                what += " Chi2Proton";
      for(int k = 0; k < 3; k++)
	if(!Same(pa->CNNscore[k], pb->CNNscore[k])){ what += " CNNscore"; break; }
      for(int k = 0; k < 3; k++)
	if(pa->Hits[k].size() != pb->Hits[k].size() || pa->HitsNoCalo[k].size() != pb->HitsNoCalo[k].size()){ what += " nHits"; break; }
      if(pa->ForcedTrack.size() != pb->ForcedTrack.size())     what += " ForcedTrack";
      if(!what.empty()){
	std::cout << "    partícula " << i << " (PFParticle " << pa->UniqueID << "):" << what << std::endl;
	nprint++;
      }
    }
    if(nprint == 0)std::cout << "    (las partículas coinciden en los campos principales: la diferencia está en otro miembro)" << std::endl;
  }
}

int compare_module_minitrees_impl(const char* fileA, const char* fileB, int maxPrint){
  std::unique_ptr<TFile> fa(TFile::Open(fileA)), fb(TFile::Open(fileB));
  if(!fa || fa->IsZombie() || !fb || fb->IsZombie()){ std::cerr << "No se pueden abrir los ficheros" << std::endl; return -1; }
  TTree* ta = GetMiniTree(fa.get());
  TTree* tb = GetMiniTree(fb.get());
  if(!ta || !tb){ std::cerr << "No se encuentra highlandana/MiniTree" << std::endl; return -1; }

  const Long64_t na = ta->GetEntries(), nb = tb->GetEntries();
  if(na != nb)std::cout << "Número de entradas distinto: " << na << " / " << nb << std::endl;
  const Long64_t n = std::min(na, nb);

  AnaSpillPD* sa = nullptr;
  AnaSpillPD* sb = nullptr;
  ta->SetBranchAddress("Spill", &sa);
  tb->SetBranchAddress("Spill", &sb);

  int ndiff = 0;
  for(Long64_t i = 0; i < n; i++){
    ta->GetEntry(i);
    tb->GetEntry(i);
    if(Serialize(sa) == Serialize(sb))continue;
    ndiff++;
    if(ndiff <= maxPrint){
      sa->CopyVectorsIntoArrays(); sa->RedoLinks();
      sb->CopyVectorsIntoArrays(); sb->RedoLinks();
      PrintParticleDiffs(sa, sb, i);
    }
  }

  std::cout << n << " entradas comparadas: " << (ndiff == 0 && na == nb ? "IDÉNTICAS" : Form("%d distintas", ndiff)) << std::endl;
  return (ndiff == 0 && na == nb) ? 0 : 1;
}
