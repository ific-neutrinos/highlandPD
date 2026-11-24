// Lectura del truth skim (truth_skim_minitrees.C) para las macros de test/phase2/skim_*.C: no necesitan highland,
// solo ROOT, y corren en local sobre el fichero de ~3 GB.
//
//   SkimReader r("truth_skim_v100s.root");
//   for(Long64_t e = 0; r.Next(e); e++){ ... r.ev.t_pdg[i] ..., r.Index(id) ... }
//
// Los vectores t_* son paralelos (una posición por partícula true no cósmica); Index(id) da la posición de la
// partícula con ese ID de Geant4 (-1 si no está guardada). Unidades del minitree: GeV, GeV/c, cm, ns.
#ifndef SKIM_COMMON_H
#define SKIM_COMMON_H

#include "TChain.h"
#include "TString.h"
#include <cmath>
#include <cstdio>
#include <unordered_map>
#include <vector>

namespace skim {
  // AnaTrueParticleB::ProcessEnum (highland)
  enum Proc { kUnknown = 0, kPrimary = 1, kDecay = 2, kKaonPlusInelastic = 3 };

  struct Event {
    Int_t run = 0, subrun = 0, event = 0, isMC = 0, fullDump = 0, beamIdx = -1, nTrue = 0, nSlice = 0, hasTVtx = 0;
    Float_t bi_p = 0;
    std::vector<int> *t_id = nullptr, *t_parent = nullptr, *t_pdg = nullptr, *t_origin = nullptr, *t_gen = nullptr, *t_procStart = nullptr,
                     *t_procEnd = nullptr, *t_nDau = nullptr, *t_nHits2 = nullptr, *t_sigMu = nullptr, *t_kOrigin = nullptr, *t_recoN = nullptr,
                     *t_recoType = nullptr, *t_recoInSlice = nullptr, *bi_pdg = nullptr;
    std::vector<float> *t_p = nullptr, *t_pEnd = nullptr, *t_x = nullptr, *t_y = nullptr, *t_z = nullptr, *t_t = nullptr, *t_xe = nullptr,
                       *t_ye = nullptr, *t_ze = nullptr, *t_te = nullptr, *t_dx = nullptr, *t_dy = nullptr, *t_dz = nullptr, *t_dxe = nullptr,
                       *t_dye = nullptr, *t_dze = nullptr, *t_eDep = nullptr, *t_recoLength = nullptr;
  };

  class SkimReader {
  public:
    TChain chain;
    Event ev;
    explicit SkimReader(const char* input): chain("skim"){
      chain.Add(input);
      chain.SetBranchStatus("*", 0);
#define SKB(x) chain.SetBranchStatus(#x, 1); chain.SetBranchAddress(#x, &ev.x);
      SKB(run) SKB(subrun) SKB(event) SKB(isMC) SKB(fullDump) SKB(beamIdx) SKB(nTrue) SKB(nSlice) SKB(hasTVtx) SKB(bi_p) SKB(bi_pdg)
      SKB(t_id) SKB(t_parent) SKB(t_pdg) SKB(t_origin) SKB(t_gen) SKB(t_procStart) SKB(t_procEnd) SKB(t_nDau) SKB(t_nHits2) SKB(t_sigMu)
      SKB(t_kOrigin) SKB(t_recoN) SKB(t_recoType) SKB(t_recoInSlice) SKB(t_p) SKB(t_pEnd) SKB(t_x) SKB(t_y) SKB(t_z) SKB(t_t) SKB(t_xe)
      SKB(t_ye) SKB(t_ze) SKB(t_te) SKB(t_dx) SKB(t_dy) SKB(t_dz) SKB(t_dxe) SKB(t_dye) SKB(t_dze) SKB(t_eDep) SKB(t_recoLength)
#undef SKB
    }
    /// Lee la entrada e (false al acabar) y rehace el índice por ID
    bool Next(Long64_t e){
      if(chain.LoadTree(e) < 0 || chain.GetEntry(e) <= 0)return false;
      idx.clear();
      for(size_t i = 0; i < ev.t_id->size(); i++)idx[(*ev.t_id)[i]] = i;
      if(e > 0 && e % 10000 == 0)fprintf(stderr, "[avance] %lld entradas\n", e);
      return true;
    }
    int Index(int id) const { auto it = idx.find(id); return it == idx.end() ? -1 : it->second; }
    int Parent(int i) const { return Index((*ev.t_parent)[i]); }
    /// Raíz de la cadena de madres de la partícula i (la última guardada si se corta; cut = true en ese caso)
    int Root(int i, bool& cut) const {
      cut = false;
      for(int depth = 0; depth < 1000; depth++){
        if((*ev.t_parent)[i] == 0)return i;
        const int m = Parent(i);
        if(m < 0){ cut = true; return i; }
        i = m;
      }
      return i;
    }
  private:
    std::unordered_map<int, int> idx;
  };

  inline double Dist(double x0, double y0, double z0, double x1, double y1, double z1){ return std::sqrt(std::pow(x1-x0,2) + std::pow(y1-y0,2) + std::pow(z1-z0,2)); }
}

#endif
