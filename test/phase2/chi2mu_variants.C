// Variantes del corte de chi2 de mu, al final de la cadena: para los candidatos que llegan al corte de chi2 de mu
// (cand_level >= 3 en las dos ramas), se sustituye ese corte por "variante < umbral" y se rehacen los cortes de
// después con las variables del microtree (los mismos umbrales que muCandCuts). Para cada variante y umbral: señal
// (capas A, B, C de candmu), fondo y S/sqrt(S+B) de B al final de la cadena. Sobre el árbol "ana" de
// secondaryKaonProductionXS (microtrees desde d764c9d para las variantes nuevas; con microtrees anteriores solo
// cand_chi2mu). La señal y el fondo se cuentan en todos los eventos (vista normal desde c0c0a92; dentro de la beam
// slice da igual).
//
//   root -l -b -q 'test/phase2/chi2mu_variants.C+("xs_v100s.root", 0)'      // 0 dentro de la beam slice, 1 fuera
//
// Variantes: cand_chi2mu (la de siempre: parada en el extremo de mayor dE/dx), cand_chi2mu_min (mejor de las dos
// orientaciones), cand_chi2mu_trim (mejor desplazando la ventana de rango residual 0-10 cm), cand_chi2mu_trimboth
// (las dos cosas). Comprobación: con cand_chi2mu < 6 tiene que salir lo mismo que cand_level al final de la cadena.
#include "TChain.h"
#include "TMath.h"
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace c2v {
  const int NMAX = 100;
  const double kGap[2] = {230., 462.};
  double GapDist(double z){ return z < -900 ? 1e9 : std::min(std::abs(z - kGap[0]), std::abs(z - kGap[1])); }
}

int chi2mu_variants(const char* input, int pop = 0){
  using namespace c2v;
  TChain t("ana"); t.Add(input);
  if(t.GetEntries() <= 0){ printf("No hay entradas en %s\n", input); return -1; }
  const std::vector<std::string> vars = {"cand_chi2mu", "cand_chi2mu_min", "cand_chi2mu_trim", "cand_chi2mu_trimboth"};
  std::vector<std::string> have;
  for(const std::string& v : vars)if(t.GetBranch(v.c_str()))have.push_back(v);

  Int_t ncand, cpop[NMAX], level[NMAX], mu[NMAX];
  Int_t nmichel[NMAX]; Float_t prange[NMAX], dedxStop[NMAX], dedxStart[NMAX], chi2p[NMAX], cnnMichel[NMAX], start[NMAX][3], stop[NMAX][3];
  std::vector<std::vector<Float_t>> var(have.size(), std::vector<Float_t>(NMAX));
  t.SetBranchStatus("*", 0);
  auto B = [&t](const char* n, void* a){ t.SetBranchStatus(n, 1); t.SetBranchAddress(n, a); };
  B("ncand", &ncand); B("cand_pop", cpop); B("cand_level", level); B("candmu", mu); B("cand_prange", prange);
  B("cand_dedx_stop", dedxStop); B("cand_dedx_start", dedxStart); B("cand_chi2p", chi2p); B("cand_ar_nmichel10", nmichel);
  B("cand_cnn_michel", cnnMichel); B("cand_start", start); B("cand_stop", stop);
  for(size_t k = 0; k < have.size(); k++)B(have[k].c_str(), var[k].data());

  // Corte de Michel sin hit en la parada (arNMichel -1): desde el 24-nov se usa la CNN media de la partícula (microtrees
  // con cand_nmichel_own); antes solo con -999, que no se daba
  const double noHit = t.GetBranch("cand_nmichel_own") ? 0. : -900.;
  // Cortes de después del chi2 de mu, en el orden de secondaryKaonProductionXSSelection
  auto after = [pop, noHit](int i, const Float_t* ds, const Float_t* dst, const Float_t* c2p, const Float_t* pr, const Int_t* nm,
                     const Float_t* cm, const Float_t (*st)[3], const Float_t (*sp)[3]) -> bool {
    if(!(ds[i] > 3.))return false;
    if(pop == 1 && !(c2p[i] > 80.))return false;
    if(pop == 1 && !(dst[i] > 1.8))return false;
    if(GapDist(st[i][2]) < 5. || GapDist(sp[i][2]) < 5.)return false;
    if(pr[i] < 205. || pr[i] > 240.)return false;
    if(nm[i] < noHit ? !(cm[i] > 0.02) : !(nm[i] >= 2))return false;
    return true;
  };
  const int lastLevel = pop == 0 ? 8 : 10;
  const std::vector<double> thr = {2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 30, 1e9};
  // n[variante][umbral][0 A, 1 B no A, 2 C no B, 3 fondo]
  std::vector<std::vector<std::vector<double>>> n(have.size(), std::vector<std::vector<double>>(thr.size(), std::vector<double>(4, 0)));
  double ref[4] = {0, 0, 0, 0}, reach[4] = {0, 0, 0, 0};
  auto cls = [](int code){ return code == 1 ? 0 : code == 2 ? 1 : code == 3 ? 2 : 3; };

  for(Long64_t e = 0; e < t.GetEntries(); e++){
    t.GetEntry(e);
    for(int i = 0; i < std::min(ncand, NMAX); i++){
      if(cpop[i] != pop || level[i] < 3)continue;
      const int c = cls(mu[i]);
      reach[c]++;
      if(level[i] >= lastLevel)ref[c]++;
      if(!after(i, dedxStop, dedxStart, chi2p, prange, nmichel, cnnMichel, start, stop))continue;
      for(size_t k = 0; k < have.size(); k++)
        for(size_t j = 0; j < thr.size(); j++)
          if(var[k][i] >= 0 && var[k][i] < thr[j])n[k][j][c]++;
    }
  }

  printf("%s, %s la beam slice. Llegan al corte de chi2 de mu: A %.0f, B no A %.0f, C no B %.0f, fondo %.0f\n", input, pop == 0 ? "dentro de" : "fuera de",
         reach[0], reach[1], reach[2], reach[3]);
  printf("Final de la cadena según cand_level (referencia): A %.0f, B no A %.0f, C no B %.0f, fondo %.0f\n", ref[0], ref[1], ref[2], ref[3]);
  for(size_t k = 0; k < have.size(); k++){
    printf("\n== %s < umbral, y el resto de la cadena ==\n   umbral      A   B-A   C-B   fondo   effB(%%)   purB   S/sqrt(S+B) (B)\n", have[k].c_str());
    for(size_t j = 0; j < thr.size(); j++){
      const auto& x = n[k][j];
      const double s = x[0] + x[1], b = x[2] + x[3], sB0 = reach[0] + reach[1];
      printf("   %6s  %5.0f %5.0f %5.0f %7.0f   %6.1f   %5.3f   %6.2f%s\n", thr[j] > 1e8 ? "sin" : Form("%.0f", thr[j]), x[0], x[1], x[2], x[3],
             sB0 ? 100.*s/sB0 : 0., s + b ? s/(s + b) : 0., s + b ? s/std::sqrt(s + b) : 0.,
             (k == 0 && thr[j] == 6 && (x[0] != ref[0] || x[1] != ref[1] || x[2] != ref[2] || x[3] != ref[3])) ? "   <-- NO COINCIDE CON cand_level" : "");
    }
  }
  printf("\nB = A + (B no A); el fondo de B incluye la capa C no B (signal muons que no son de la cascada del haz).\n");
  return 0;
}
