// Comprobación del contenido de los minitrees de una producción nueva (tag d01, 25-nov-2026): hits de la traza forzada
// (solo con longitud en ForcedTrackHitsLength, en los planos de SaveHitsAllPlanes) y tiempo de las partículas true.
#include "minitree_loop.h"
#include <map>
#include <cmath>
#include <fstream>
#include <string>

int check_production_minitree_impl(const char* pattern, double lmin, double lmax){
  TChain chain("highlandana/MiniTree");
  // lista .txt (una URL xrootd por línea, como en la receta 4) o patrón de ficheros
  int nFiles = 0;
  if(TString(pattern).EndsWith(".txt")){
    std::ifstream in(pattern); std::string l;
    while(std::getline(in, l))if(!l.empty() && l[0] != '#')nFiles += chain.Add(l.c_str());
  }
  else nFiles = chain.Add(pattern);
  AnaSpillPD* spill = nullptr;
  chain.SetBranchAddress("Spill", &spill);
  long nFT = 0, nFTin = 0, nFTinHits = 0, nFToutHits = 0, nFTinNoHits = 0, nFTplane[3] = {0, 0, 0}, nFTtrj = 0, nFTchi2 = 0;
  long nTrue = 0, nTime = 0, nTimeEnd = 0, nBeamT = 0, nEv = 0; double tMin = 1e30, tMax = -1e30;
  std::map<int, long> tPileup;   // primarias del haz distintas de la que dispara, por intervalo de tiempo (ms)
  LoopMiniTree(chain, spill, -1, nFiles, [&](AnaSpillPD* s){
    nEv++;
    for(AnaBunchC* bc : s->Bunches)for(AnaParticleB* pb : static_cast<AnaBunchB*>(bc)->Particles){
      const AnaParticlePD* p = static_cast<const AnaParticlePD*>(pb);
      for(const AnaForcedTrackPD& ft : p->ForcedTrack){
        nFT++;
        const bool in = ft.Length > lmin && ft.Length < lmax;
        size_t nh = 0; for(int pl = 0; pl < 3; pl++){ nh += ft.Hits[pl].size(); if(in && !ft.Hits[pl].empty())nFTplane[pl]++; }
        if(in){ nFTin++; if(nh)nFTinHits++; else nFTinNoHits++; if(!ft.TrjPoints.empty())nFTtrj++; }
        else if(nh)nFToutHits++;
        if(ft.Chi2Muon > 0)nFTchi2++;
      }
    }
    const AnaBeamPD* beam = static_cast<const AnaBeamPD*>(s->Beam);
    const AnaTrueParticleB* bt = (beam && beam->BeamParticle) ? static_cast<const AnaTrueParticleB*>(beam->BeamParticle->TrueObject) : nullptr;
    for(AnaTrueParticleB* t : s->TrueParticles){
      nTrue++;
      const double t0 = t->Position[3], t1 = t->PositionEnd[3];
      if(t0 != 0 && t0 > -900){ nTime++; tMin = std::min(tMin, t0); tMax = std::max(tMax, t0); }
      if(t1 != 0 && t1 > -900)nTimeEnd++;
      if(bt && t->ID == bt->ID && t0 > -900)nBeamT++;
      if(bt && t->ParentID == 0 && static_cast<const AnaTrueParticlePD*>(t)->Origin == 4 && t->ID != bt->ID && t0 > -900)tPileup[(int)std::floor((t0 - bt->Position[3])/1e6)]++;
    }
  });
  printf("Ficheros %d, entradas %ld\n", nFiles, nEv);
  printf("Trazas forzadas: %ld (con chi2 de mu %ld)\n", nFT, nFTchi2);
  printf("  con longitud en (%.2f, %.2f) cm: %ld; con hits %ld, sin hits %ld; con puntos de trayectoria %ld\n", lmin, lmax, nFTin, nFTinHits, nFTinNoHits, nFTtrj);
  printf("  con hits por plano (dentro de la ventana): U %ld, V %ld, W %ld\n", nFTplane[0], nFTplane[1], nFTplane[2]);
  printf("  fuera de la ventana con hits (debería ser 0): %ld\n", nFToutHits);
  printf("Partículas true: %ld; con tiempo inicial (!= 0 y > -900) %ld (%.1f %%), final %ld; rango del inicial [%.3g, %.3g] ns; haz con tiempo: %ld de %ld eventos\n",
         nTrue, nTime, nTrue ? 100.*nTime/nTrue : 0., nTimeEnd, tMin, tMax, nBeamT, nEv);
  printf("Primarias de pileup frente al haz que dispara, por ms: ");
  for(auto& kv : tPileup)printf(" [%d,%d) %ld", kv.first, kv.first + 1, kv.second);
  printf("\n");
  return 0;
}
