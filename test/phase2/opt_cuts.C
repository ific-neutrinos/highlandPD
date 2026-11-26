// Optimización conjunta de los umbrales de secondaryKaonProductionXS sobre el árbol plano de opt_flat.C (26-nov-2026).
//
//   root -l -b -q 'test/phase2/opt_cuts.C+("opt_flat_v100to105d01.root")'
//
// Figura de mérito: error relativo de un recuento de la señal (B = capas A y B no A; o A) en la ventana de p de rango,
// con la estadística de los datos: sqrt(S+B)/S/sqrt(kHalf), kHalf = 2 x 1,196 (cada mitad es la mitad del MC y el MC
// entero es 1/1,196 de los datos a 6 GeV). C no B es fondo. Búsqueda al azar (nTry configuraciones) alrededor de los
// umbrales actuales, por población, con la mitad 0 (half == 0); la mejor se evalúa en la mitad 1, y al revés.
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace optc {
  struct Cand { int pop, cat, half; float prange, faces, vtx, axis, clear, chi2mu, chi2p, dstop, dstart, gapmin, nmich, cnnm; };
  struct Th { double faces, vtx, axis, chi2mu, chi2p, dstop, dstart, gap, lo, hi; int nmich, clear; };
  const double kHalf = 2*1.196;
  bool Pass(const Cand& c, const Th& t){
    if(c.prange < t.lo || c.prange > t.hi)return false;
    if(c.chi2mu >= t.chi2mu || c.dstop <= t.dstop || c.gapmin < t.gap)return false;
    if(t.nmich > 0 && (c.nmich < 0 ? c.cnnm <= 0.02 : c.nmich < t.nmich))return false;
    if(c.pop == 0){ if(c.faces <= t.faces)return false; if(!(c.vtx < 0 || c.vtx > t.vtx))return false; }
    else{ if(t.clear && c.clear)return false; if(c.axis >= t.axis || c.chi2p <= t.chi2p || c.dstart <= t.dstart)return false; }
    return true;
  }
  // error relativo (%) de A (which 0) o de B (which 1) en la mitad h
  double Err(const std::vector<Cand>& v, const Th& t, int h, int which, double* sOut = nullptr, double* nOut = nullptr){
    double s = 0, n = 0;
    for(const Cand& c : v){ if(c.half != h || !Pass(c, t))continue; n++; if(which == 0 ? c.cat == 0 : c.cat <= 1)s++; }
    if(sOut)*sOut = s; if(nOut)*nOut = n;
    return s > 0 ? 100.*std::sqrt(n)/s/std::sqrt(kHalf) : 1e9;
  }
  void Print(const char* tag, const Th& t, int pop){
    if(pop == 0)printf("    %s: faces > %.0f, vtx > %.1f, chi2mu < %.2f, dE/dx stop > %.2f, gap >= %.1f, Michel >= %d, ventana [%.0f, %.0f]\n",
                       tag, t.faces, t.vtx, t.chi2mu, t.dstop, t.gap, t.nmich, t.lo, t.hi);
    else printf("    %s: %sclear, eje < %.0f, chi2mu < %.2f, dE/dx stop > %.2f, chi2p > %.0f, dE/dx start > %.2f, gap >= %.1f, Michel >= %d, ventana [%.0f, %.0f]\n",
                tag, t.clear ? "no " : "con ", t.axis, t.chi2mu, t.dstop, t.chi2p, t.dstart, t.gap, t.nmich, t.lo, t.hi);
  }
}

int opt_cuts(const char* in, int nTry = 20000, int seed = 1){
  using namespace optc;
  TFile* f = TFile::Open(in); TTree* t = (TTree*)f->Get("flat");
  Float_t pop, prange, faces, vtx, axis, clear, chi2mu, chi2p, dstop, dstart, zs, ze, nmich, cnnm; Int_t cat, half;
  t->SetBranchAddress("pop", &pop); t->SetBranchAddress("prange", &prange); t->SetBranchAddress("faces", &faces); t->SetBranchAddress("vtx", &vtx);
  t->SetBranchAddress("axis", &axis); t->SetBranchAddress("clear", &clear); t->SetBranchAddress("chi2mu", &chi2mu); t->SetBranchAddress("chi2p", &chi2p);
  t->SetBranchAddress("dedx_stop", &dstop); t->SetBranchAddress("dedx_start", &dstart); t->SetBranchAddress("zstart", &zs); t->SetBranchAddress("zstop", &ze);
  t->SetBranchAddress("nmich10", &nmich); t->SetBranchAddress("cnn_michel", &cnnm); t->SetBranchAddress("cat", &cat); t->SetBranchAddress("half", &half);
  auto gap = [](double z){ return std::min(std::abs(z - 230.), std::abs(z - 462.)); };
  std::vector<Cand> v[2];
  for(Long64_t e = 0; e < t->GetEntries(); e++){
    t->GetEntry(e);
    const int p = (int)pop;
    // preselección suelta (los umbrales buscados están por encima)
    if(chi2mu > 20 || (p == 0 && faces < 10) || (p == 1 && axis > 150))continue;
    v[p].push_back({p, cat, half, prange, faces, vtx, axis, clear, chi2mu, chi2p, dstop, dstart, (float)std::min(gap(zs), gap(ze)), nmich, cnnm});
  }
  printf("Candidatos tras la preselección suelta: dentro %zu, fuera %zu\n", v[0].size(), v[1].size());
  const Th cur[2] = {{35, 15, 0, 6, 0, 3, 0, 5, 205, 240, 2, 0}, {0, 0, 80, 6, 80, 3, 1.8, 5, 205, 240, 2, 1}};
  TRandom3 r(seed);
  const char* popName[2] = {"dentro de la beam slice", "fuera"};
  for(int p = 0; p < 2; p++){
    printf("\n== %s ==\n", popName[p]);
    for(int which = 1; which >= 0; which--){
      printf("  Figura de mérito: error de %s\n", which ? "B" : "A");
      for(int h = 0; h < 2; h++){
        double s, n;
        const double e0 = Err(v[p], cur[p], h, which, &s, &n);
        printf("   mitad %d, cortes actuales: S %.0f, N %.0f, error %.2f %%\n", h, s, n, e0);
      }
      for(int hTrain = 0; hTrain < 2; hTrain++){
        Th best = cur[p]; double eBest = Err(v[p], cur[p], hTrain, which);
        for(int i = 0; i < nTry; i++){
          Th x = cur[p];
          x.faces = r.Uniform(0, 60); x.vtx = r.Uniform(0, 30); x.chi2mu = r.Uniform(2, 15); x.dstop = r.Uniform(0, 5); x.gap = r.Uniform(0, 12);
          x.nmich = (int)r.Uniform(0, 5); x.lo = r.Uniform(195, 215); x.hi = r.Uniform(230, 250);
          if(p == 1){ x.axis = r.Uniform(40, 150); x.chi2p = r.Uniform(0, 150); x.dstart = r.Uniform(0, 3); x.clear = r.Rndm() < 0.8; }
          const double e = Err(v[p], x, hTrain, which);
          if(e < eBest){ eBest = e; best = x; }
        }
        double s1, n1; const int hTest = 1 - hTrain;
        const double eTest = Err(v[p], best, hTest, which, &s1, &n1);
        printf("   entrenada en la mitad %d: error %.2f %% ahí; en la mitad %d: S %.0f, N %.0f, error %.2f %% (actuales %.2f %%)\n",
               hTrain, eBest, hTest, s1, n1, eTest, Err(v[p], cur[p], hTest, which));
        Print("umbrales", best, p);
      }
    }
  }
  return 0;
}
