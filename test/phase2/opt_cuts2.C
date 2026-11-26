// Optimización de cortes con las variables nuevas que el BDT señala como útiles (26-nov-2026): carga media por hit
// alrededor del inicio del mu (as_mean5), K al inicio (k_n, k_dedx) y Michel (nmich10, michel_own, cnn_michel), además
// de los cortes actuales. Árbol plano de opt_flat.C (opt_flat2, con las variables del inicio).
//
//   root -l -b -q 'test/phase2/opt_cuts2.C+("opt_flat2_v100to105d01.root", 0, 0)'   // población 0, modo Y
//   root -l -b -q 'test/phase2/opt_cuts2.C+("opt_flat2_v100to105d01.root", 0, 1)'   // población 0, modo O
//
// Modo Y (0): todos los cortes a la vez; cada uno puede quedar apagado si su umbral sale en el mínimo del barrido.
// Modo O (1): los cortes "duros" de antes más dos indicadores, Michel (M: nmich10 >= n o michel_own > m) y K al inicio
// (K: as_mean5 > a o k_dedx > d), cada uno con un umbral flojo y uno estricto; pasa si (M estricto y K flojo) o
// (M flojo y K estricto): basta con que uno de los dos sea claro si el otro no falla del todo.
// Figura de mérito, mitades y escala como opt_cuts.C.
#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace optc2 {
  struct Cand { int cat, half; float prange, faces, vtx, axis, clear, chi2mu, chi2p, dstop, dstart, gapmin, nmich, cnnm, mown, asm5, kdedx, kn, cosb; };
  struct Th {
    double faces, vtx, axis, chi2mu, chi2p, dstop, dstart, gap, lo, hi; int clear;
    // indicadores: [0] flojo, [1] estricto (en modo Y solo se usa [0])
    int nmich[2]; double mown[2], asm5[2], kdedx[2], cnnm, cosb;
  };
  const double kHalf = 2*1.196;
  bool M(const Cand& c, const Th& t, int s){ return (c.nmich < 0 ? c.cnnm > 0.02 : c.nmich >= t.nmich[s]) || c.mown > t.mown[s]; }
  bool K(const Cand& c, const Th& t, int s){ return (c.asm5 < 0 || c.asm5 > t.asm5[s]) || (c.kn > 0 && c.kdedx > t.kdedx[s]); }
  bool Pass(const Cand& c, const Th& t, int pop, int mode){
    if(c.prange < t.lo || c.prange > t.hi)return false;
    if(c.chi2mu >= t.chi2mu || c.dstop <= t.dstop || c.gapmin < t.gap || c.cnnm <= t.cnnm || c.cosb >= t.cosb)return false;
    if(pop == 0){ if(c.faces <= t.faces)return false; if(!(c.vtx < 0 || c.vtx > t.vtx))return false; }
    else{ if(t.clear && c.clear)return false; if(c.axis >= t.axis || c.chi2p <= t.chi2p || c.dstart <= t.dstart)return false; }
    if(mode == 0){
      // Y: Michel por hits de alrededor (como ahora) y, aparte, los cortes nuevos
      if(t.nmich[0] > 0 && (c.nmich < 0 ? c.cnnm <= 0.02 : c.nmich < t.nmich[0]))return false;
      if(c.mown <= t.mown[0] && t.mown[0] > -0.5)return false;
      if(!(c.asm5 < 0 || c.asm5 > t.asm5[0]))return false;
      if(t.kdedx[0] > 0 && !(c.kn > 0 && c.kdedx > t.kdedx[0]))return false;
      return true;
    }
    return (M(c, t, 1) && K(c, t, 0)) || (M(c, t, 0) && K(c, t, 1));
  }
  double Err(const std::vector<Cand>& v, const Th& t, int pop, int mode, int h, int which, double* sOut = nullptr, double* nOut = nullptr){
    double s = 0, n = 0;
    for(const Cand& c : v){ if(c.half != h || !Pass(c, t, pop, mode))continue; n++; if(which == 0 ? c.cat == 0 : c.cat <= 1)s++; }
    if(sOut)*sOut = s; if(nOut)*nOut = n;
    return s > 0 ? 100.*std::sqrt(n)/s/std::sqrt(kHalf) : 1e9;
  }
  void Print(const Th& t, int pop, int mode){
    printf("      comunes: chi2mu < %.2f, dE/dx stop > %.2f, gap >= %.1f, cnn_michel > %.2f, cosbeam < %.2f, ventana [%.0f, %.0f]", t.chi2mu, t.dstop, t.gap, t.cnnm, t.cosb, t.lo, t.hi);
    if(pop == 0)printf(", faces > %.0f, vtx > %.1f\n", t.faces, t.vtx);
    else printf(", %sclear, eje < %.0f, chi2p > %.0f, dE/dx start > %.2f\n", t.clear ? "no " : "con ", t.axis, t.chi2p, t.dstart);
    if(mode == 0)printf("      Y: Michel >= %d, michel_own > %.2f, as_mean5 > %.0f, k_dedx > %.2f\n", t.nmich[0], t.mown[0], t.asm5[0], t.kdedx[0]);
    else printf("      O: M flojo (Michel >= %d o own > %.2f), M estricto (>= %d o own > %.2f); K flojo (as_mean5 > %.0f o k_dedx > %.2f), K estricto (> %.0f o > %.2f)\n",
                t.nmich[0], t.mown[0], t.nmich[1], t.mown[1], t.asm5[0], t.kdedx[0], t.asm5[1], t.kdedx[1]);
  }
}

int opt_cuts2(const char* in, int pop, int mode, int nTry = 20000, int seed = 2){
  using namespace optc2;
  TFile* f = TFile::Open(in); TTree* t = (TTree*)f->Get("flat");
  Float_t fpop, prange, faces, vtx, axis, clear, chi2mu, chi2p, dstop, dstart, zs, ze, nmich, cnnm, mown, asm5, kdedx, kn, cosb; Int_t cat, half;
  auto B = [&](const char* n, void* a){ t->SetBranchAddress(n, a); };
  B("pop", &fpop); B("prange", &prange); B("faces", &faces); B("vtx", &vtx); B("axis", &axis); B("clear", &clear); B("chi2mu", &chi2mu);
  B("chi2p", &chi2p); B("dedx_stop", &dstop); B("dedx_start", &dstart); B("zstart", &zs); B("zstop", &ze); B("nmich10", &nmich);
  B("cnn_michel", &cnnm); B("michel_own", &mown); B("as_mean5", &asm5); B("k_dedx", &kdedx); B("k_n", &kn); B("cosbeam", &cosb);
  B("cat", &cat); B("half", &half);
  auto gap = [](double z){ return std::min(std::abs(z - 230.), std::abs(z - 462.)); };
  std::vector<Cand> v;
  for(Long64_t e = 0; e < t->GetEntries(); e++){
    t->GetEntry(e);
    if((int)fpop != pop || prange < 195 || prange > 250 || chi2mu > 15 || (pop == 1 && axis > 150))continue;
    v.push_back({cat, half, prange, faces, vtx, axis, clear, chi2mu, chi2p, dstop, dstart, (float)std::min(gap(zs), gap(ze)), nmich, cnnm, mown, asm5, kdedx, kn, cosb});
  }
  printf("Población %d (%s), modo %s: %zu candidatos tras la preselección suelta\n", pop, pop ? "fuera" : "dentro", mode ? "O" : "Y", v.size());
  // punto de partida: los cortes actuales, con los nuevos apagados
  Th cur;
  cur.faces = 35; cur.vtx = 15; cur.axis = 80; cur.chi2mu = 6; cur.chi2p = pop ? 80 : -1e9; cur.dstop = 3; cur.dstart = pop ? 1.8 : -1e9; cur.gap = 5;
  cur.lo = 205; cur.hi = 240; cur.clear = pop; cur.cnnm = -1; cur.cosb = 2;
  cur.nmich[0] = cur.nmich[1] = 2; cur.mown[0] = cur.mown[1] = -1; cur.asm5[0] = cur.asm5[1] = -1; cur.kdedx[0] = cur.kdedx[1] = 0;
  if(mode == 1){ cur.mown[0] = cur.mown[1] = 2; }   // O: michel_own > 2 = apagado; K trivial (as_mean5 > -1): M estricto o M flojo = Michel >= 2, la cadena actual
  TRandom3 r(seed);
  for(int which = 1; which >= 0; which--){
    printf("\n  Figura de mérito: error de %s. Cortes actuales: mitad 0 %.2f %%, mitad 1 %.2f %%\n", which ? "B" : "A",
           Err(v, cur, pop, mode, 0, which), Err(v, cur, pop, mode, 1, which));
    for(int hTrain = 0; hTrain < 2; hTrain++){
      Th best = cur; double eBest = Err(v, cur, pop, mode, hTrain, which);
      for(int i = 0; i < nTry; i++){
        Th x = cur;
        x.faces = r.Uniform(0, 60); x.vtx = r.Uniform(0, 30); x.chi2mu = r.Uniform(2, 15); x.dstop = r.Uniform(0, 5); x.gap = r.Uniform(0, 12);
        x.lo = r.Uniform(195, 215); x.hi = r.Uniform(230, 250);
        x.cnnm = r.Rndm() < 0.5 ? -1 : r.Uniform(0, 0.3); x.cosb = r.Rndm() < 0.5 ? 2 : r.Uniform(0.3, 1);
        if(pop == 1){ x.axis = r.Uniform(40, 150); x.chi2p = r.Uniform(0, 150); x.dstart = r.Uniform(0, 3); x.clear = r.Rndm() < 0.8; }
        for(int s = 0; s < 2; s++){
          x.nmich[s] = (int)r.Uniform(0, 5); x.mown[s] = r.Rndm() < 0.4 ? (mode == 1 ? 2. : -1.) : r.Uniform(0, 0.8);   // apagado: -1 en Y (sin corte), 2 en O (no da M)
          x.asm5[s] = r.Rndm() < 0.3 ? -1 : r.Uniform(30, 110); x.kdedx[s] = r.Rndm() < 0.5 ? 0 : r.Uniform(2, 8);
        }
        if(mode == 1){   // estricto por encima del flojo
          if(x.nmich[1] < x.nmich[0])std::swap(x.nmich[0], x.nmich[1]);
          if(x.mown[1] < x.mown[0])std::swap(x.mown[0], x.mown[1]);
          if(x.asm5[1] < x.asm5[0])std::swap(x.asm5[0], x.asm5[1]);
          if(x.kdedx[1] < x.kdedx[0])std::swap(x.kdedx[0], x.kdedx[1]);
        }
        const double e = Err(v, x, pop, mode, hTrain, which);
        if(e < eBest){ eBest = e; best = x; }
      }
      double s1, n1; const int hTest = 1 - hTrain;
      const double eTest = Err(v, best, pop, mode, hTest, which, &s1, &n1);
      printf("    elegidos en la mitad %d: %.2f %% ahí; en la mitad %d: S %.0f, N %.0f, error %.2f %%\n", hTrain, eBest, hTest, s1, n1, eTest);
      Print(best, pop, mode);
    }
  }
  return 0;
}
