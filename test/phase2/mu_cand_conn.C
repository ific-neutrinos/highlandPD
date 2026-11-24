// Conexión de los candidatos a mu con el vértice primario y origen de su K (23-nov-2026), sobre el árbol cand de
// mu_candidates_minitrees.C con las ramas del 23-nov (mu_cand_v100b.root o posteriores).
//
//   root -l -b -q 'test/phase2/mu_cand_conn.C+("mu_cand_v100b.root")'
//
// Clases de los candidatos:
//   - mu+ de K+ en reposo (mu+ de desintegración, no cósmico, p true en [230, 241] MeV/c, madre K+ o no guardada), por
//     el origen de su K (tKOrigin): A (primera interacción del haz que dispara, 0 y 1), secundaria (2), K+ del haz (3),
//     pileup (otra raíz, 6), cadena cortada (4, 5). De todos los eventos (casi todos con volcado completo).
//   - fondo que no es de K: pi, p, mu cósmico, EM, resto. Solo de eventos normales, escalado a todos (como en
//     mu_cand_cutflow.C).
// Conectado (definición provisional): un PFParticle con un extremo a < d cm del vértice y el otro a < d cm de un
// extremo del mu, o un extremo del mu a < d cm del vértice (d = 10 cm; también 5 y 15). Sin vértice (sin beam slice),
// no conectado.
// Etapas: "suelta" (cortes del módulo y p de rango en [205, 240]) y "cadena" (la cadena corta de mu_cand_common.h, con
// las opciones por defecto).
#include "mu_cand_common.h"
#include <map>

namespace mconn {
  using mucut::C;
  const double kSc = 73070./65977.;
  enum Cls { kA = 0, kSec, kBeamK, kPileup, kCut, kPi, kP, kMuCos, kEM, kRest, kNCls };
  const char* kClsName[kNCls] = {"mu de K: A", "mu de K: secundaria", "mu de K: K+ del haz", "mu de K: pileup", "mu de K: cadena cortada",
                                 "pi", "p", "mu cósmico", "EM", "resto"};

  bool IsMuK(const C& c){
    return c.tPDG == -13 && c.tProcStart == 2 && c.tOrigin != 2 && c.tMom > 230 && c.tMom < 241 && c.tParentPDG == 321;
  }
  // clase; -1 si no se usa (fondo de eventos con volcado)
  int Class(const C& c){
    if(IsMuK(c)){
      switch(c.tKOrigin){ case 0: case 1: return kA; case 2: return kSec; case 3: return kBeamK; case 6: return kPileup; default: return kCut; }
    }
    if(c.fullDump)return -1;
    if(std::abs(c.tPDG) == 211)return kPi;
    if(std::abs(c.tPDG) == 2212)return kP;
    if(std::abs(c.tPDG) == 13 && c.tOrigin == 2)return kMuCos;
    if(std::abs(c.tPDG) == 11 || c.tPDG == 22)return kEM;
    return kRest;
  }
  double Weight(int cls){ return cls >= kPi ? kSc : 1.; }
  bool Connected(const C& c, double d){
    if(c.dStartV < 0)return false;   // sin vértice
    return (c.lkSdV >= 0 && c.lkSdV < d && c.lkSdE < d) || (c.lkPdV >= 0 && c.lkPdV < d && c.lkPdE < d) || c.dStartV < d || c.dStopV < d;
  }
}

int mu_cand_conn(const char* in){
  using namespace mucut; using namespace mconn;
  gKeepFullDump = true;
  std::vector<C> cands = Load(in);
  if(cands.empty())return 1;
  const std::vector<Step> steps[2] = {OutSliceSteps(), InSliceSteps()};
  const char* popName[2] = {"fuera de la beam slice", "en la beam slice"};

  for(int pop = 1; pop >= 0; pop--){
    for(int st = 0; st < 2; st++){
      double n[kNCls] = {0}, conn[3][kNCls] = {{0}}, vtx[kNCls] = {0};
      const double dd[3] = {5., 10., 15.};
      for(const C& c : cands){
        if(c.inSlice != (bool)pop)continue;
        if(st == 0 ? !(c.pRange >= 205 && c.pRange <= 240) : !Pass(steps[pop], c))continue;
        const int k = Class(c); if(k < 0)continue;
        const double w = Weight(k);
        n[k] += w; if(c.dStartV >= 0)vtx[k] += w;
        for(int j = 0; j < 3; j++)if(Connected(c, dd[j]))conn[j][k] += w;
      }
      std::cout << std::endl << "== " << popName[pop] << ", " << (st == 0 ? "suelta (módulo + [205, 240])" : "cadena corta") << " ==" << std::endl
                << Form("  %-26s %9s %9s   %-26s", "clase", "candidatos", "con vért.", "conectados (5 / 10 / 15 cm)") << std::endl;
      double nA = 0, cA = 0, nAll = 0, cAll = 0;
      for(int k = 0; k < kNCls; k++){
        std::cout << Form("  %-26s %9.1f %8.0f %%   %5.0f %% / %5.0f %% / %5.0f %%", kClsName[k], n[k], n[k] ? 100.*vtx[k]/n[k] : 0.,
                          n[k] ? 100.*conn[0][k]/n[k] : 0., n[k] ? 100.*conn[1][k]/n[k] : 0., n[k] ? 100.*conn[2][k]/n[k] : 0.) << std::endl;
        nAll += n[k]; cAll += conn[1][k];
        if(k == kA){ nA = n[k]; cA = conn[1][k]; }
      }
      std::cout << Form("  pureza para A: sin conexión %.3f (%.0f de %.1f); exigiendo conexión (10 cm) %.3f (%.0f de %.1f)",
                        nAll ? nA/nAll : 0., nA, nAll, cAll ? cA/cAll : 0., cA, cAll) << std::endl;
    }
  }

  // Señal A: conexión según la calidad del vértice
  std::cout << std::endl << "Señal A (mu de K en reposo de la primera interacción), cadena corta: conexión (10 cm) según el error del vértice reco (tVtxErr):" << std::endl;
  for(int pop = 1; pop >= 0; pop--){
    const double edges[] = {0, 5, 10, 20, 50, 1e9};
    for(int b = 0; b < 5; b++){
      double n = 0, c10 = 0;
      for(const C& c : cands){
        if(c.inSlice != (bool)pop || Class(c) != kA || !Pass(steps[pop], c))continue;
        if(c.tVtxErr < edges[b] || c.tVtxErr >= edges[b+1])continue;
        n++; if(Connected(c, 10))c10++;
      }
      std::cout << Form("  %-24s tVtxErr %4.0f-%-4s cm: %4.0f, conectados %3.0f %%", popName[pop], edges[b], b == 4 ? "" : Form("%.0f", edges[b+1]), n, n ? 100.*c10/n : 0.) << std::endl;
    }
    double n = 0, nNo = 0;
    for(const C& c : cands)if(c.inSlice == (bool)pop && Class(c) == kA && Pass(steps[pop], c)){ n++; if(c.tVtxErr < 0)nNo++; }
    std::cout << Form("  %-24s sin vértice o sin truth del haz: %.0f de %.0f", popName[pop], nNo, n) << std::endl;
  }

  // Pileup: dónde está y qué ve Pandora
  std::cout << std::endl << "Pileup (mu de K cuya cadena acaba en otra partícula del haz), cortes del módulo:" << std::endl;
  for(int st = 0; st < 2; st++){
    double in = 0, out = 0, clear = 0, outNoSlice = 0;
    for(const C& c : cands){
      if(Class(c) != kPileup)continue;
      if(st == 1 && !Pass(c.inSlice ? steps[1] : steps[0], c))continue;
      (c.inSlice ? in : out)++;
      if(c.clearCosmic)clear++;
      if(!c.inSlice && !c.hasSlice)outNoSlice++;
    }
    std::cout << Form("  %-14s en la beam slice %4.0f, fuera %4.0f (de ellos en eventos sin beam slice %4.0f); IsClearCosmic %4.0f",
                      st == 0 ? "módulo:" : "cadena corta:", in, out, outNoSlice, clear) << std::endl;
  }
  return 0;
}
