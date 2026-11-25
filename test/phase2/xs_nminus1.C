// N-1 de cada corte de secondaryKaonProductionXS sobre el microtree "ana": para cada población y cada corte, los
// candidatos que pasan todos los demás cortes de su cadena, y qué quita el corte entre ellos: signal muons por capa
// (A, B no A, C no B) y fondo por especie true. Para el corte de vértice (dentro), además, qué es el fondo que quita.
// Los cortes se recalculan con las variables del microtree (como en secondaryKaonProductionXSSelection; 25-nov-2026).
//
//   root -l -b -q 'test/phase2/xs_nminus1.C("xs_v100s_fix.root")'
#include "TFile.h"
#include "TTree.h"
#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace xsn1 {
  const int N = 100;
  struct C {
    Int_t pop, level, clear, sig, kor, pdg, ppdg, fromBeam, shower;
    Float_t faces_stop, faces_start, vtx, chi2mu, dedx_stop, dedx_start, chi2p, axis, nmich, cnnmich, prange, start[3], stop[3];
  };
  double GapDist(double z){ return z < -900 ? 1e9 : std::min(std::abs(z - 230.), std::abs(z - 462.)); }
  struct Cut { std::string name; std::function<bool(const C&)> pass; };
  // categoría: 0 A, 1 B no A, 2 C no B, 3 pi, 4 p, 5 mu (no señal), 6 e/gamma, 7 otro
  int Cat(const C& c){
    if(c.sig == 1){ if(c.kor == 0 || c.kor == 1)return 0; if(c.kor == 2 || c.kor == 3)return 1; return 2; }
    const int a = std::abs(c.pdg);
    return a == 211 ? 3 : a == 2212 ? 4 : a == 13 ? 5 : (a == 11 || a == 22) ? 6 : 7;
  }
  const char* kCat[8] = {"A", "B no A", "C no B", "pi", "p", "mu", "e/gamma", "otro"};
}

int xs_nminus1(const char* in){
  using namespace xsn1;
  TFile* f = TFile::Open(in);
  if(!f || f->IsZombie())return -1;
  TTree* t = (TTree*)f->Get("ana");
  Int_t nc = 0;
  Int_t pop[N], level[N], clear[N], sig[N], kor[N], pdg[N], ppdg[N], fromBeam[N], shower[N];
  Float_t faces_stop[N], faces_start[N], vtx[N], chi2mu[N], dedx_stop[N], dedx_start[N], chi2p[N], axis[N], nmich[N], cnnmich[N], prange[N], start[N][3], stop[N][3];
  Int_t inmich[N];
  t->SetBranchStatus("*", 0);
  auto B = [&](const char* n, void* a){ t->SetBranchStatus(n, 1); t->SetBranchAddress(n, a); };
  B("ncand", &nc); B("cand_pop", pop); B("cand_level", level); B("cand_clearcosmic", clear); B("cand_truesignalmuon", sig);
  B("cand_truekorigin", kor); B("cand_truepdg", pdg); B("cand_trueparentpdg", ppdg); B("cand_truefrombeam", fromBeam); B("cand_isshower", shower);
  B("cand_faces_stop", faces_stop); B("cand_faces_start", faces_start); B("cand_vtx_min", vtx); B("cand_chi2mu", chi2mu);
  B("cand_dedx_stop", dedx_stop); B("cand_dedx_start", dedx_start); B("cand_chi2p", chi2p); B("cand_axis_min", axis);
  B("cand_ar_nmichel10", inmich); B("cand_cnn_michel", cnnmich); B("cand_prange", prange); B("cand_start", start); B("cand_stop", stop);

  // cadenas de secondaryKaonProductionXSSelection::DefineSteps, sin el corte de población
  auto michel = [](const C& c){ return c.nmich < 0 ? c.cnnmich > 0.02 : c.nmich >= 2; };
  auto gaps = [](const C& c){ return GapDist(c.start[2]) >= 5 && GapDist(c.stop[2]) >= 5; };
  auto chi2 = [](const C& c){ return c.chi2mu >= 0 && c.chi2mu < 6; };
  std::vector<Cut> chain[2] = {
    {{"faces > 35", [](const C& c){ return std::min(c.faces_stop, c.faces_start) > 35; }},
     {"vertex > 15", [](const C& c){ return c.vtx < 0 || c.vtx > 15; }},
     {"chi2mu < 6", chi2}, {"dEdx stop > 3", [](const C& c){ return c.dedx_stop > 3; }}, {"gaps", gaps}, {"Michel >= 2", michel},
     {"p range [205, 240]", [](const C& c){ return c.prange >= 205 && c.prange <= 240; }}},
    {{"no IsClearCosmic", [](const C& c){ return !c.clear; }}, {"axis < 60", [](const C& c){ return c.axis < 60; }},
     {"chi2mu < 6", chi2}, {"dEdx stop > 3", [](const C& c){ return c.dedx_stop > 3; }},
     {"chi2p > 80", [](const C& c){ return c.chi2p > 80; }}, {"dEdx start > 1.8", [](const C& c){ return c.dedx_start > 1.8; }},
     {"gaps", gaps}, {"Michel >= 2", michel}, {"p range [205, 240]", [](const C& c){ return c.prange >= 205 && c.prange <= 240; }}}
  };
  // nN1[pop][cut][cat]: pasan todos los demás; nRem: y este no
  std::vector<std::vector<std::vector<double>>> nN1(2), nRem(2);
  for(int p = 0; p < 2; p++){ nN1[p].assign(chain[p].size(), std::vector<double>(8, 0)); nRem[p].assign(chain[p].size(), std::vector<double>(8, 0)); }
  double nAll[2][8] = {{0}}, check[2][8] = {{0}};
  std::map<std::string, double> vtxRemBkg;   // fondo que quita el corte de vértice (N-1): PDG / madre / del haz
  std::vector<double> vtxRemSig[3];          // distancia de la señal quitada, por capa
  std::vector<double> vtxN1[8];              // distancia al vértice de la muestra N-1 del corte de vértice, por categoría
  std::vector<double> axisN1[8];             // distancia al eje de la muestra N-1 del corte del eje (fuera)
  for(Long64_t e = 0; e < t->GetEntries(); e++){
    t->GetEntry(e);
    for(int j = 0; j < nc && j < N; j++){
      C c{pop[j], level[j], clear[j], sig[j], kor[j], pdg[j], ppdg[j], fromBeam[j], shower[j], faces_stop[j], faces_start[j], vtx[j], chi2mu[j],
          dedx_stop[j], dedx_start[j], chi2p[j], axis[j], (Float_t)inmich[j], cnnmich[j], prange[j], {start[j][0], start[j][1], start[j][2]}, {stop[j][0], stop[j][1], stop[j][2]}};
      if(c.pop < 0 || c.pop > 1 || c.level < 1)continue;   // tiene que pasar el corte de población
      const int k = Cat(c), p = c.pop;
      auto& ch = chain[p];
      std::vector<bool> ok(ch.size()); int nfail = 0;
      for(size_t i = 0; i < ch.size(); i++){ ok[i] = ch[i].pass(c); if(!ok[i])nfail++; }
      nAll[p][k]++;
      if(nfail == 0){ check[p][k]++; }
      for(size_t i = 0; i < ch.size(); i++){
        if(nfail - (ok[i] ? 0 : 1) != 0)continue;   // falla algún otro
        nN1[p][i][k]++;
        if(p == 0 && i == 1)vtxN1[k].push_back(c.vtx);
        if(p == 1 && i == 1)axisN1[k].push_back(c.axis);
        if(!ok[i]){
          nRem[p][i][k]++;
          if(p == 0 && i == 1){
            if(k >= 3)vtxRemBkg[Form("%d (madre %d)%s", c.pdg, c.ppdg, c.fromBeam ? ", desciende del haz" : "")]++;
            else vtxRemSig[k].push_back(c.vtx);
          }
        }
      }
    }
  }
  const char* popName[2] = {"dentro de la beam slice", "fuera"};
  for(int p = 0; p < 2; p++){
    printf("\n== %s: candidatos que pasan todos los cortes (comprobación frente a cand_level):", popName[p]);
    for(int k = 0; k < 8; k++)printf("  %s %.0f", kCat[k], check[p][k]); printf("\n");
    printf("  corte                N-1: señal A / B no A / C no B | fondo      quita: A / B no A / C no B | fondo (pi, p, mu, e/gamma, otro)   S/sqrt(S+B) de B con / sin el corte\n");
    for(size_t i = 0; i < chain[p].size(); i++){
      const auto& n = nN1[p][i]; const auto& r = nRem[p][i];
      double bN = 0, bR = 0; for(int k = 3; k < 8; k++){ bN += n[k]; bR += r[k]; }
      const double sN = n[0] + n[1], sR = r[0] + r[1];
      const double sWith = sN - sR, bWith = bN - bR;
      printf("  %-20s %5.0f / %5.0f / %5.0f | %6.0f      %4.0f (%4.1f%%) / %4.0f (%4.1f%%) / %4.0f | %6.0f (%4.1f%%: %.0f, %.0f, %.0f, %.0f, %.0f)   %.2f / %.2f\n",
             chain[p][i].name.c_str(), n[0], n[1], n[2], bN, r[0], n[0] ? 100*r[0]/n[0] : 0., r[1], n[1] ? 100*r[1]/n[1] : 0., r[2], bR, bN ? 100*bR/bN : 0.,
             r[3], r[4], r[5], r[6], r[7], sWith + bWith > 0 ? sWith/std::sqrt(sWith + bWith) : 0., sN + bN > 0 ? sN/std::sqrt(sN + bN) : 0.);
    }
  }
  printf("\n== Corte de vértice (dentro), N-1: fondo que quita, por PDG true (madre):\n");
  std::vector<std::pair<double, std::string>> v; for(auto& kv : vtxRemBkg)v.push_back({kv.second, kv.first});
  std::sort(v.rbegin(), v.rend()); for(size_t i = 0; i < v.size() && i < 12; i++)printf("    %-45s %5.0f\n", v[i].second.c_str(), v[i].first);
  for(int k = 0; k < 3; k++){
    auto& x = vtxRemSig[k]; std::sort(x.begin(), x.end());
    printf("  señal %-7s quitada: %zu, distancia al vértice (cm):", kCat[k], x.size()); for(double d : x)printf(" %.1f", d); printf("\n");
  }
  printf("\n  Barrido del umbral del corte de vértice (muestra N-1; sin vértice pasa): umbral (cm): A / B no A / fondo, S/sqrt(S+B) de B, eff x pur de A y de B\n");
  for(double th : {0., 2., 4., 6., 8., 10., 12., 15., 20., 25.}){
    double n[8] = {0}; for(int k = 0; k < 8; k++)for(double d : vtxN1[k])if(d < 0 || d > th)n[k]++;
    double b = 0; for(int k = 3; k < 8; k++)b += n[k];
    const double sB = n[0] + n[1], tot = n[0] + n[1] + n[2] + b;
    printf("    %4.0f: %4.0f / %4.0f / %5.0f   %.2f   %.3f  %.3f\n", th, n[0], n[1], b, sB/std::sqrt(sB + b + n[2]), n[0]/vtxN1[0].size()*n[0]/tot,
           sB/(vtxN1[0].size() + vtxN1[1].size())*sB/tot);
  }
  printf("\n  Barrido del umbral del corte del eje (fuera, muestra N-1): umbral (cm): A / B no A / C no B / fondo, S/sqrt(S+B) de B (C no B es fondo), eff x pur de B\n");
  for(double th : {40., 50., 60., 80., 100., 120., 150., 200., 1e9}){
    double n[8] = {0}; for(int k = 0; k < 8; k++)for(double d : axisN1[k])if(d < th)n[k]++;
    double b = n[2]; for(int k = 3; k < 8; k++)b += n[k];
    const double sB = n[0] + n[1];
    printf("    %5.0f: %4.0f / %4.0f / %4.0f / %5.0f   %.2f   %.4f\n", th > 1e8 ? -1. : th, n[0], n[1], n[2], b - n[2], sB/std::sqrt(sB + b),
           sB/(axisN1[0].size() + axisN1[1].size())*sB/(sB + b));
  }
  return 0;
}
