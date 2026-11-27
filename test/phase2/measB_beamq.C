// Medida B (26-nov-2026): cortes de calidad del haz, composición del grupo e/mu/pi y separación del pileup fuera de la
// beam slice, para los mu de K seleccionados por la cadena de cortes (en la ventana).
//
//   root -l -b -q 'test/phase2/measB_beamq.C+("opt_flat2_v100to105d01.root", "xs_v100to105d01.root")'
//
// Calidad del haz (como pdAnaUtils::GetBeamQuality): el inicio del primario de la beam slice (seltrk_pos) a menos de 3
// sigma de la media en x, y, z y cos(dirección del primario, dirección del BI al final, beam_enddir) > 0,95. Las medias y
// sigmas de pdUtils.parameters.dat para MC a 6 GeV (z = 0,0 +- 0,2 cm) no valen para estos minitrees (el primario empieza
// en z ~ 30 cm): se miden aquí con mediana y 1,4826 MAD de los eventos con beam slice cuyo primario es la partícula del haz.
// Sin beam slice no hay calidad (falla). Clases de los mu de K seleccionados: A o secundaria con el haz acabando en el TPC
// (z > 0) o antes (z < 0), K del haz, pileup.
// sce = 1 (27-nov-2026, por defecto): la calidad es beam_quality_sce de los xs_*d01b (inicio corregido por SCE, cortes de
// pdUtils.parameters.dat); las medias sin SCE se siguen midiendo e imprimiendo. sce = 0: la provisional.
#include "TFile.h"
#include "TTree.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <vector>

namespace bq {
  double Median(std::vector<double> v){ if(v.empty())return 0; std::sort(v.begin(), v.end()); return v[v.size()/2]; }
  double Mad(const std::vector<double>& v, double m){ std::vector<double> d; for(double x : v)d.push_back(std::abs(x - m)); return 1.4826*Median(d); }
}

int measB_beamq(const char* flat, const char* xs, int sce = 1){
  using namespace bq;
  TFile fx(xs); TTree* ana = (TTree*)fx.Get("ana");
  Float_t pos[4], dir[3], bdir[3], tend[4]; Int_t hasslice, sp, tpdg, qsce = -1;
  ana->SetBranchStatus("*", 0);
  for(const char* b : {"seltrk_pos", "seltrk_dir", "beam_enddir", "beam_trueendpos", "hasslice", "beam_truespecies", "beam_truepdg"})ana->SetBranchStatus(b, 1);
  ana->SetBranchAddress("seltrk_pos", pos); ana->SetBranchAddress("seltrk_dir", dir); ana->SetBranchAddress("beam_enddir", bdir);
  ana->SetBranchAddress("beam_trueendpos", tend); ana->SetBranchAddress("hasslice", &hasslice);
  ana->SetBranchAddress("beam_truespecies", &sp); ana->SetBranchAddress("beam_truepdg", &tpdg);
  if(sce){ ana->SetBranchStatus("beam_quality_sce", 1); ana->SetBranchAddress("beam_quality_sce", &qsce); }
  const Long64_t nEv = ana->GetEntries();
  struct Ev { float x, y, z, cos, tz; int slice, sp, pdg, qsce; };
  std::vector<Ev> ev(nEv);
  std::vector<double> vx, vy, vz;
  for(Long64_t e = 0; e < nEv; e++){
    ana->GetEntry(e);
    double c = 0; for(int i = 0; i < 3; i++)c += dir[i]*bdir[i];
    ev[e] = {pos[0], pos[1], pos[2], (float)c, tend[2], hasslice, sp, tpdg, qsce};
    // para medir medias y sigmas: con beam slice, primario hacia delante y el haz true entrando en el TPC
    if(hasslice && pos[2] > -900 && c > 0.95 && tend[2] > 50){ vx.push_back(pos[0]); vy.push_back(pos[1]); vz.push_back(pos[2]); }
  }
  const double mx = Median(vx), my = Median(vy), mz = Median(vz), sx = Mad(vx, mx), sy = Mad(vy, my), sz = Mad(vz, mz);
  printf("Inicio del primario (eventos con beam slice, cos > 0,95, haz true acabando a z > 50): x %.1f +- %.2f, y %.1f +- %.2f, z %.1f +- %.2f cm (n = %zu)\n",
         mx, sx, my, sy, mz, sz, vx.size());
  printf("Calidad del haz: %s\n", sce ? "beam_quality_sce (inicio corregido por SCE)" : "provisional (sin SCE, medias de arriba)");
  auto quality = [&](const Ev& v){ if(sce)return v.slice && v.qsce == 1; return v.slice && v.z > -900 && std::abs(v.x - mx) < 3*sx && std::abs(v.y - my) < 3*sy && std::abs(v.z - mz) < 3*sz && v.cos > 0.95; };
  // denominador: partículas del haz por especie, con y sin calidad, y cuántas acaban antes del TPC
  double nB[3][2] = {{0}}, nUp[3][2] = {{0}};
  for(auto& v : ev){ if(v.sp < 0 || v.sp > 2)continue; const bool q = quality(v); nB[v.sp][0]++; if(v.tz < 0)nUp[v.sp][0]++; if(q){ nB[v.sp][1]++; if(v.tz < 0)nUp[v.sp][1]++; } }
  const char* sn[3] = {"e/mu/pi", "K", "p"};
  printf("\nPartículas del haz (todas -> con calidad), y las que acaban antes del TPC (z true < 0):\n");
  for(int s = 0; s < 3; s++)printf("  %-8s %7.0f -> %7.0f (%4.1f %%);  antes del TPC %5.0f (%4.1f %%) -> %4.0f (%4.1f %% de las de calidad)\n", sn[s], nB[s][0], nB[s][1],
                                   100*nB[s][1]/nB[s][0], nUp[s][0], 100*nUp[s][0]/nB[s][0], nUp[s][1], 100*nUp[s][1]/nB[s][1]);
  // mu de K seleccionados (cadena entera, ventana), por clase, con y sin calidad
  TFile ff(flat); TTree* t = (TTree*)ff.Get("flat");
  Float_t fpop, pr, kor, smin; Int_t cat, chain, chain80; Long64_t entry;
  t->SetBranchStatus("*", 0);
  for(const char* b : {"pop", "prange", "kor", "cat", "chain", "chain80", "entry", "slice_min"})t->SetBranchStatus(b, 1);
  t->SetBranchAddress("pop", &fpop); t->SetBranchAddress("prange", &pr); t->SetBranchAddress("kor", &kor); t->SetBranchAddress("cat", &cat);
  t->SetBranchAddress("chain", &chain); t->SetBranchAddress("chain80", &chain80); t->SetBranchAddress("entry", &entry); t->SetBranchAddress("slice_min", &smin);
  const char* cn[6] = {"cascada del haz, haz acaba en el TPC", "cascada del haz, haz acaba antes", "K del haz (KOrigin 3)", "pileup (KOrigin 6)", "otro", "fondo (no mu de K)"};
  double n[2][3][6][2] = {{{{0}}}};   // [pob][especie][clase][calidad no/sí]
  std::map<int, double> pdgSig;        // e/mu/pi, dentro, señal de la cascada: PDG true del haz
  std::vector<double> sminPU, sminOwn, sminBkg;   // fuera: distancia a la beam slice
  for(Long64_t i = 0; i < t->GetEntries(); i++){
    t->GetEntry(i);
    const int p = (int)fpop; if(p < 0 || p > 1 || !(p == 0 ? chain : chain80) || pr < 205 || pr > 240)continue;
    const Ev& v = ev[entry]; if(v.sp < 0 || v.sp > 2)continue;
    const int k = (int)kor; int c;
    if(cat > 2)c = 5; else if(k >= 0 && k <= 2)c = v.tz > 0 ? 0 : 1; else if(k == 3)c = 2; else if(k == 6)c = 3; else c = 4;
    const bool q = quality(v);
    n[p][v.sp][c][0]++; if(q)n[p][v.sp][c][1]++;
    if(p == 0 && v.sp == 0 && c == 0)pdgSig[v.pdg]++;
    if(p == 1 && v.sp == 0 && v.slice){ if(c == 3)sminPU.push_back(smin); else if(c == 0 || c == 1)sminOwn.push_back(smin); else if(c == 5)sminBkg.push_back(smin); }
  }
  for(int p = 0; p < 2; p++){
    printf("\n== %s: candidatos seleccionados en la ventana (todos -> con calidad del haz)\n", p ? "fuera" : "dentro");
    for(int s = 0; s < 3; s++){
      printf("  haz %s\n", sn[s]);
      for(int c = 0; c < 6; c++)if(n[p][s][c][0] > 0)printf("    %-40s %5.0f -> %5.0f (%3.0f %%)\n", cn[c], n[p][s][c][0], n[p][s][c][1], 100*n[p][s][c][1]/n[p][s][c][0]);
    }
  }
  printf("\nComposición (e/mu/pi, dentro, mu de la cascada con el haz en el TPC): PDG true del haz:");
  double tot = 0; for(auto& kv : pdgSig)tot += kv.second;
  for(auto& kv : pdgSig)printf("  %d: %.0f (%.1f %%)", kv.first, kv.second, 100*kv.second/tot);
  printf("\n");
  // fracción del grupo e/mu/pi por PDG (todas las partículas y con calidad)
  std::map<int, double> pdgAll[2];
  for(auto& v : ev)if(v.sp == 0){ pdgAll[0][v.pdg]++; if(quality(v))pdgAll[1][v.pdg]++; }
  for(int q = 0; q < 2; q++){ printf("Grupo e/mu/pi (%s), PDG true:", q ? "con calidad" : "todas"); double tt = 0; for(auto& kv : pdgAll[q])tt += kv.second;
    for(auto& kv : pdgAll[q])if(kv.second/tt > 0.001)printf("  %d: %.1f %%", kv.first, 100*kv.second/tt); printf("\n"); }
  // slice_min fuera (eventos con beam slice): pileup frente a la cascada propia
  auto quant = [](std::vector<double> v){ std::sort(v.begin(), v.end()); return v.empty() ? std::vector<double>{} : std::vector<double>{v[v.size()/10], v[v.size()/4], v[v.size()/2], v[3*v.size()/4], v[9*v.size()/10]}; };
  printf("\nFuera, e/mu/pi, eventos con beam slice: distancia a la beam slice (slice_min, cm), cuantiles 10/25/50/75/90\n");
  for(int k = 0; k < 3; k++){
    const std::vector<double>& v = k == 0 ? sminOwn : k == 1 ? sminPU : sminBkg;
    const auto q = quant(v); if(q.empty())continue;
    printf("  %-22s n %4zu: %6.1f %6.1f %6.1f %6.1f %6.1f\n", k == 0 ? "cascada del haz" : k == 1 ? "pileup" : "fondo", v.size(), q[0], q[1], q[2], q[3], q[4]);
  }
  for(double th : {20., 40., 60., 100.}){
    double a = 0, b = 0, c = 0; for(double x : sminOwn)if(x < th)a++; for(double x : sminPU)if(x < th)b++; for(double x : sminBkg)if(x < th)c++;
    printf("  slice_min < %3.0f: cascada %3.0f %%, pileup %3.0f %%, fondo %3.0f %%\n", th, 100*a/std::max<size_t>(1, sminOwn.size()), 100*b/std::max<size_t>(1, sminPU.size()), 100*c/std::max<size_t>(1, sminBkg.size()));
  }
  return 0;
}
