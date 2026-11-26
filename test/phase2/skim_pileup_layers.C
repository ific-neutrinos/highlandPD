// Capas A' y B' para el pileup del haz (26-nov-2026): los signal muons con KOrigin 6 (la cadena de madres acaba en una
// primaria del haz que no es la que dispara) clasificados como KOrigin pero respecto a su raíz en lugar del haz que
// dispara. A' = 0-1 (primera interacción de la raíz), B' = 0-3 (cascada de la raíz). Es lo que ganarían A y B si el BI
// identificara también las partículas de pileup (pendiente de ver en datos: BI con varias trazas, d02).
//
//   root -l -b -q 'test/phase2/skim_pileup_layers.C+("truth_skim_v100d01.root")'
//
// Necesita el tiempo true (minitrees desde d01). Por capa: PDG y momento de la raíz, dt respecto al haz que dispara (ms),
// si la raíz llega al TPC y si el mu está reconstruido y en la beam slice.
#include "skim_common.h"
#include <algorithm>
#include <map>
#include <string>

namespace {
  // Como muCandUtils::KOrigin, con la raíz r en el papel del haz que dispara. -1: no aplica.
  int KOriginRoot(const skim::SkimReader& r, int iMu, int& iRoot){
    const auto& ev = r.ev;
    iRoot = -1;
    const int iK = r.Parent(iMu);
    if(iK < 0 || (*ev.t_pdg)[iK] != 321)return -1;
    std::vector<int> chain{iK};
    int cur = iK;
    for(int depth = 0; depth < 1000 && (*ev.t_parent)[cur] != 0; depth++){
      const int m = r.Parent(cur);
      if(m < 0)return 4;
      cur = m; chain.push_back(cur);
    }
    iRoot = cur;
    size_t i0 = 0;
    while(i0 + 1 < chain.size() && (*ev.t_pdg)[chain[i0 + 1]] == 321)i0++;
    if(i0 + 1 == chain.size())return 3;
    if(i0 + 2 == chain.size())return i0 == 0 ? 0 : 1;
    return 2;
  }
  std::string PdgName(int p){
    switch(p){ case 211: return "pi+"; case -211: return "pi-"; case 321: return "K+"; case 2212: return "p"; case 130: return "K0L";
      case 2112: return "n"; case -13: return "mu+"; case 13: return "mu-"; default: return std::to_string(p); }
  }
  std::string MsBin(double dt){
    if(std::abs(dt) < 1e4)return "|dt| < 10 us";
    const long k = (long)std::floor(dt/1e6);
    return "[" + std::to_string(k) + "," + std::to_string(k+1) + ") ms";
  }
  std::string PBin(double p){ return p < 1 ? "< 1" : p < 3 ? "1-3" : p < 5 ? "3-5" : p < 5.7 ? "5-5,7" : p < 6.3 ? "5,7-6,3" : "> 6,3"; }
  void Print(const char* title, const std::map<std::string, long>& m, long tot){
    std::vector<std::pair<long, std::string>> v; for(auto& kv : m)v.push_back({kv.second, kv.first});
    std::sort(v.rbegin(), v.rend());
    printf("    %s:", title); for(auto& x : v)printf("  %s %ld (%.0f %%)", x.second.c_str(), x.first, tot ? 100.*x.first/tot : 0.); printf("\n");
  }
  struct Layer {
    long n = 0, toTPC = 0, inSlice = 0, outSlice = 0, noReco = 0;
    std::map<std::string, long> pdg, dt, p;
  };
}

int skim_pileup_layers(const char* input){
  skim::SkimReader r(input);
  long nSig = 0, nSigByK[8] = {0};
  std::map<int, long> kOr6;          // KOrigin respecto a la raíz, para los KOrigin 6
  Layer L[3];                        // 0: A' (0-1), 1: B' no A' (2-3), 2: resto (4 o -1)
  long noBeamT = 0;
  for(Long64_t e = 0; r.Next(e); e++){
    const auto& ev = r.ev;
    if(!ev.isMC)continue;
    const int iB = ev.beamIdx;
    const double tB = iB >= 0 ? (*ev.t_t)[iB] : 0;
    for(size_t i = 0; i < ev.t_id->size(); i++){
      if((*ev.t_sigMu)[i] != 1)continue;
      nSig++;
      const int ko = (*ev.t_kOrigin)[i];
      if(ko >= 0 && ko < 8)nSigByK[ko]++;
      if(ko != 6)continue;
      int iR = -1;
      const int kr = KOriginRoot(r, i, iR);
      kOr6[kr]++;
      Layer& l = L[(kr == 0 || kr == 1) ? 0 : (kr == 2 || kr == 3) ? 1 : 2];
      l.n++;
      if(iR >= 0){
        l.pdg[PdgName((*ev.t_pdg)[iR]) + ((*ev.t_origin)[iR] == 4 ? "" : " (Origin " + std::to_string((*ev.t_origin)[iR]) + ")")]++;
        l.p[PBin((*ev.t_p)[iR])]++;
        if(iB >= 0)l.dt[MsBin((*ev.t_t)[iR] - tB)]++; else { l.dt["sin haz"]++; noBeamT++; }
        if((*ev.t_ze)[iR] > 0)l.toTPC++;
      }
      if((*ev.t_recoN)[i] <= 0)l.noReco++; else if((*ev.t_recoInSlice)[i] > 0)l.inSlice++; else l.outSlice++;
    }
  }
  printf("Signal muons %ld; por KOrigin:", nSig); for(int k = 0; k < 8; k++)if(nSigByK[k])printf("  %d: %ld", k, nSigByK[k]); printf("\n");
  printf("  A (0-1) %ld, B (0-3) %ld\n", nSigByK[0] + nSigByK[1], nSigByK[0] + nSigByK[1] + nSigByK[2] + nSigByK[3]);
  printf("KOrigin 6 (%ld) clasificados respecto a su raíz:", nSigByK[6]); for(auto& kv : kOr6)printf("  %d: %ld", kv.first, kv.second); printf("\n");
  const char* names[3] = {"A' (primera interacción de la raíz, 0-1)", "B' no A' (2-3)", "resto (cadena cortada o sin K+ madre)"};
  for(int k = 0; k < 3; k++){
    const Layer& l = L[k];
    printf("\n== %s: %ld ==\n", names[k], l.n);
    if(!l.n)continue;
    Print("raíz", l.pdg, l.n);
    Print("p de la raíz (GeV/c)", l.p, l.n);
    Print("dt de la raíz respecto al haz que dispara", l.dt, l.n);
    printf("    la raíz llega al TPC (final en z > 0): %ld (%.0f %%)\n", l.toTPC, 100.*l.toTPC/l.n);
    printf("    mu: en la beam slice %ld, fuera %ld, sin reconstruir %ld\n", l.inSlice, l.outSlice, l.noReco);
  }
  const long A = nSigByK[0] + nSigByK[1], B = A + nSigByK[2] + nSigByK[3];
  printf("\nGanancia si el BI identificara el pileup: A %ld -> %ld (x%.2f), B %ld -> %ld (x%.2f)\n",
         A, A + L[0].n, A ? 1. + (double)L[0].n/A : 0., B, B + L[0].n + L[1].n, B ? 1. + (double)(L[0].n + L[1].n)/B : 0.);
  return 0;
}
