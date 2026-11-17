// Implementación de size_breakdown_minitree.C (se compila con ACLiC).
#include "pdDataClasses.hxx"
#include "TBufferFile.h"
#include "TClass.h"
#include "TFile.h"
#include "TTree.h"
#include <cstdio>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace {

  // bytes serializados (sin comprimir) de un objeto de la clase dada
  double Bytes(const void* obj, const char* cl){
    TClass* c = TClass::GetClass(cl);
    if(!obj || !c)return 0;
    TBufferFile buf(TBuffer::kWrite);
    buf.WriteObjectAny(obj, c);
    return buf.Length();
  }

  struct Counter{
    double bytes = 0;
    double n = 0;   // número de objetos (hits, puntos, partículas...)
  };

  std::map<std::string, Counter> gC;
  std::vector<std::string> gOrder;

  void Add(const std::string& key, double bytes, double n){
    if(!gC.count(key))gOrder.push_back(key);
    gC[key].bytes += bytes;
    gC[key].n += n;
  }

  // Reparte una partícula en su parte fija y sus vectores, por diferencia: se serializa, se vacía
  // un vector y se vuelve a serializar (así las partes suman exactamente el total). Vacía la
  // partícula: llamar después de haber usado su contenido. Devuelve el tamaño total.
  template <class V>
  double Drop(void* obj, const char* cl, V& v, double& current){
    v.clear();
    const double now = Bytes(obj, cl);
    const double d = current - now;
    current = now;
    return d;
  }

  double AddParticle(AnaParticlePD* p, const std::string& who){
    const char* cl = "AnaParticlePD";
    const double total = Bytes(p, cl);
    double cur = total;
    const char* plane[3] = {"U", "V", "W (colección)"};
    for(int k = 0; k < 3; k++){
      const double nh = p->Hits[k].size();
      Add(who + " / Hits " + plane[k], Drop(p, cl, p->Hits[k], cur), nh);
    }
    double bnc = 0, nnc = 0;
    for(int k = 0; k < 3; k++){ nnc += p->HitsNoCalo[k].size(); bnc += Drop(p, cl, p->HitsNoCalo[k], cur); }
    Add(who + " / HitsNoCalo (3 planos)", bnc, nnc);
    double n = p->TrjPoints.size();
    Add(who + " / TrjPoints", Drop(p, cl, p->TrjPoints, cur), n);
    n = p->SpacePoints.size();
    Add(who + " / SpacePoints", Drop(p, cl, p->SpacePoints, cur), n);
    n = p->ForcedTrack.size();
    Add(who + " / ForcedTrack", Drop(p, cl, p->ForcedTrack, cur), n);
    Add(who + " / resto de la partícula", cur, 1);
    return total;
  }
}

int size_breakdown_minitree_impl(const char* file, Long64_t maxEntries){
  std::unique_ptr<TFile> f(TFile::Open(file));
  if(!f || f->IsZombie()){ std::cerr << "No se puede abrir " << file << std::endl; return -1; }
  TTree* t = dynamic_cast<TTree*>(f->Get("highlandana/MiniTree"));
  if(!t)t = dynamic_cast<TTree*>(f->Get("MiniTree"));
  if(!t){ std::cerr << "No se encuentra highlandana/MiniTree" << std::endl; return -1; }

  AnaSpillPD* spill = nullptr;
  t->SetBranchAddress("Spill", &spill);
  Long64_t n = t->GetEntries();
  if(maxEntries > 0 && maxEntries < n)n = maxEntries;

  double totalSpill = 0;
  int nGoodDaq = 0, nBeamSlice = 0;
  for(Long64_t i = 0; i < n; i++){
    t->GetEntry(i);
    const double sb = Bytes(spill, "AnaSpillPD");
    totalSpill += sb;
    double parts = 0;

    if(!spill->Bunches.empty())nGoodDaq++;
    for(AnaBunchC* bc : spill->Bunches){
      AnaBunchPD* b = static_cast<AnaBunchPD*>(bc);
      bool hasBeam = false;
      for(AnaParticleB* pb : b->Particles){
	AnaParticlePD* p = static_cast<AnaParticlePD*>(pb);
	if(p->BeamOrigin)hasBeam = true;
	parts += AddParticle(p, p->BeamOrigin ? "Beam slice" : "Cósmicos");
      }
      if(hasBeam)nBeamSlice++;
      const double bu = Bytes(&b->UnassociatedHits, "vector<AnaHitPD>");
      Add("Hits no asociados", bu, b->UnassociatedHits.size());
      parts += bu;
    }

    for(AnaTrueParticleB* tb : spill->TrueParticles){
      AnaTrueParticlePD* tp = static_cast<AnaTrueParticlePD*>(tb);
      const double tot = Bytes(tp, "AnaTrueParticlePD");
      double cur = tot;
      const double ntj = tp->TrjPoints.size();
      Add("Partículas true / TrjPoints", Drop(tp, "AnaTrueParticlePD", tp->TrjPoints, cur), ntj);
      Add("Partículas true / resto", cur, 1);
      parts += tot;
    }
    Add("Resto del spill (EventInfo, Beam, DataQuality, cabeceras)", sb > parts ? sb - parts : 0, 1);
  }

  // factor de compresión global del árbol
  const double zip = t->GetZipBytes(), tot = t->GetTotBytes();
  const double cf = tot > 0 ? zip / tot : 1;
  const double sumC = [&]{ double s = 0; for(auto& kv : gC)s += kv.second.bytes; return s; }();

  printf("\n%s\n", file);
  printf("%lld entradas (%d con bunch, %d con beam slice). En disco: %.1f MB (%.3f MB/entrada); compresión %.2f\n",
	 n, nGoodDaq, nBeamSlice, zip/1e6*n/t->GetEntries(), zip/1e6/t->GetEntries(), cf);
  printf("Sin comprimir: %.1f MB los componentes, %.1f MB los spills serializados\n", sumC/1e6, totalSpill/1e6);
  printf("\n%-58s %10s %8s %12s %10s\n", "componente", "MB disco*", "%", "objetos", "B/objeto");
  for(const std::string& k : gOrder){
    const Counter& c = gC[k];
    if(c.bytes <= 0)continue;
    printf("%-58s %10.2f %7.1f%% %12.0f %10.0f\n", k.c_str(), c.bytes*cf/1e6, 100*c.bytes/sumC, c.n, c.n > 0 ? c.bytes/c.n : 0);
  }
  printf("* estimado con el factor de compresión global\n");
  return 0;
}
