// Medida B (27-nov-2026): denominador por pi+ con la fracción de pi del grupo e/mu/pi medida en el TPC (opción b del
// usuario). Separación sencilla del primario de la beam slice, solo con variables reconstruidas:
//   e  : puntuación EM de la CNN del primario (seltrk_CNNscore[1]) > emCut (0,8: con 0,5 el 11,8 % de los pi pasa por e;
//        con 0,8 el 2,2 %, con el 86 % de los e);
//   mu : no e y la traza llega lejos sin interaccionar: z del final "empalmado" > zCut. Pandora corta las trazas largas
//        en los huecos entre APAs (z ~ 230 y ~ 462 cm): se empalma con la hija del primario que empieza a < dMax del final
//        (cualquiera de sus extremos) y es colineal (|cos| > cosMin con la dirección final del primario), y se toma la z
//        de su otro extremo (un nivel de hijas: es lo que hay en el microtree). Si el microtree tiene beam_reach
//        (27-nov, muCandUtils::BeamReach: empalme con cualquier partícula del evento, varios trozos), se usa ese;
//   pi : el resto.
// Matriz de confusión frente a la partícula true del haz (MC), con la calidad del haz con SCE y beam slice, y barrido
// de zCut. Estimación de la fracción de pi corrigiendo con la matriz (inversión 3x3) y su sensibilidad a la
// supervivencia de los pi (escala de la fila pi->mu, que depende de la sección eficaz pi-Ar).
//
//   root -l -b -q 'test/phase2/measB_beampid.C+("xs_v100to105d01b.root")'
#include "TFile.h"
#include "TTree.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include <cmath>
#include <cstdio>
#include <vector>

int measB_beampid(const char* xs, double emCut = 0.8, double dMax = 30, double cosMin = 0.95){
  TFile fx(xs); TTree* t = (TTree*)fx.Get("ana");
  const int kMaxDau = 200;
  Int_t sp, hs, q, tpdg, nd;
  Float_t end[4], edir[3], cnn[3], dpos[kMaxDau][4], dend[kMaxDau][4], ddir[kMaxDau][3], reach[3] = {-999, -999, -999};
  t->SetBranchStatus("*", 0);
  for(const char* b : {"beam_truespecies", "hasslice", "beam_quality_sce", "beam_truepdg", "seltrk_ndau", "seltrk_endpos", "seltrk_enddir",
                       "seltrk_CNNscore", "seltrk_dau_pos", "seltrk_dau_endpos", "seltrk_dau_dir"})t->SetBranchStatus(b, 1);
  t->SetBranchAddress("beam_truespecies", &sp); t->SetBranchAddress("hasslice", &hs); t->SetBranchAddress("beam_quality_sce", &q);
  t->SetBranchAddress("beam_truepdg", &tpdg); t->SetBranchAddress("seltrk_ndau", &nd); t->SetBranchAddress("seltrk_endpos", end);
  t->SetBranchAddress("seltrk_enddir", edir); t->SetBranchAddress("seltrk_CNNscore", cnn);
  t->SetBranchAddress("seltrk_dau_pos", dpos); t->SetBranchAddress("seltrk_dau_endpos", dend); t->SetBranchAddress("seltrk_dau_dir", ddir);
  const bool hasReach = t->GetBranch("beam_reach");
  if(hasReach){ t->SetBranchStatus("beam_reach", 1); t->SetBranchAddress("beam_reach", reach); }
  printf("Empalme: %s\n", hasReach ? "beam_reach (todas las partículas del evento)" : "hijas del primario (no hay beam_reach)");
  if(t->GetMaximum("seltrk_ndau") > kMaxDau){ printf("ERROR: seltrk_ndau > %d\n", kMaxDau); return 1; }
  struct Ev { int tr; float em, z0, z1; };   // tr: 0 pi, 1 mu, 2 e; z0 fin del primario, z1 fin empalmado
  std::vector<Ev> ev;
  long nStitch[3] = {0}, nTr[3] = {0};
  for(Long64_t e = 0; e < t->GetEntries(); e++){
    t->GetEntry(e);
    if(sp != 0 || !hs || q != 1)continue;
    const int a = std::abs(tpdg); const int tr = a == 211 ? 0 : a == 13 ? 1 : a == 11 ? 2 : -1; if(tr < 0)continue;
    float z1 = end[2]; bool st = false;
    for(int i = 0; i < nd; i++){
      double d0 = 0, d1 = 0; for(int k = 0; k < 3; k++){ d0 += std::pow(dpos[i][k] - end[k], 2); d1 += std::pow(dend[i][k] - end[k], 2); }
      d0 = std::sqrt(d0); d1 = std::sqrt(d1);
      if(std::min(d0, d1) > dMax || dpos[i][2] < -900 || dend[i][2] < -900)continue;
      double c = 0; for(int k = 0; k < 3; k++)c += ddir[i][k]*edir[k];
      if(std::abs(c) < cosMin)continue;
      const float zf = d0 < d1 ? dend[i][2] : dpos[i][2];
      if(zf > z1){ z1 = zf; st = true; }
    }
    if(hasReach){ st = reach[2] > end[2] + 1; z1 = reach[2] > -900 ? reach[2] : end[2]; }
    nTr[tr]++; if(st)nStitch[tr]++;
    ev.push_back({tr, cnn[1], end[2], z1});
  }
  const char* nm[3] = {"pi+", "mu+", "e+"};
  printf("Grupo e/mu/pi, beam slice, calidad del haz con SCE: pi %ld, mu %ld, e %ld\n", nTr[0], nTr[1], nTr[2]);
  printf("Empalmados (trozo colineal a < %.0f cm del final que llega más lejos): pi %.1f %%, mu %.1f %%, e %.1f %%\n", dMax,
         100.*nStitch[0]/nTr[0], 100.*nStitch[1]/nTr[1], 100.*nStitch[2]/nTr[2]);
  // fin del primario y fin empalmado, por bins de 50 cm
  for(int w = 0; w < 2; w++){
    printf("%s (bins de 50 cm, 0-700, %%):\n", w ? "fin empalmado" : "fin del primario");
    for(int tr = 0; tr < 3; tr++){
      std::vector<double> h(15, 0);
      for(auto& x : ev)if(x.tr == tr){ const float z = w ? x.z1 : x.z0; int b = (int)(z/50); if(b < 0)b = 0; if(b > 14)b = 14; h[b]++; }
      printf("  %-4s", nm[tr]); for(int b = 0; b < 15; b++)printf(" %5.1f", 100*h[b]/nTr[tr]); printf("\n");
    }
  }
  // matriz de confusión y fracción de pi corregida, para varios zCut
  printf("\nClasificación: e si CNN EM > %.2f; mu si no e y fin empalmado > zCut; pi el resto. Filas: true; columnas: reco (%%).\n", emCut);
  for(double zCut : {300., 400., 480., 550., 600., 650.}){
    double m[3][3] = {{0}};
    for(auto& x : ev){ const int r = x.em > emCut ? 2 : x.z1 > zCut ? 1 : 0; m[x.tr][r]++; }
    printf("zCut %.0f cm\n", zCut);
    double nr[3] = {0}; for(int r = 0; r < 3; r++)for(int tr = 0; tr < 3; tr++)nr[r] += m[tr][r];
    for(int tr = 0; tr < 3; tr++)printf("  %-4s -> pi %5.1f  mu %5.1f  e %5.1f\n", nm[tr], 100*m[tr][0]/nTr[tr], 100*m[tr][1]/nTr[tr], 100*m[tr][2]/nTr[tr]);
    printf("  pureza de lo clasificado: pi %5.1f %%, mu %5.1f %%, e %5.1f %%\n", 100*m[0][0]/nr[0], 100*m[1][1]/nr[1], 100*m[2][2]/nr[2]);
    // corrección: reco = M^T true, con M[tr][r] = P(r|tr); fracción de pi = solución / total
    TMatrixD M(3, 3); for(int tr = 0; tr < 3; tr++)for(int r = 0; r < 3; r++)M(r, tr) = m[tr][r]/nTr[tr];
    TVectorD R(3); for(int r = 0; r < 3; r++)R(r) = nr[r];
    const double tot = nr[0] + nr[1] + nr[2];
    printf("  fracción de pi: true %.4f; sin corregir (reco) %.4f; corregida con la matriz %.4f\n", nTr[0]/tot, nr[0]/tot, (TMatrixD(TMatrixD::kInverted, M)*R)(0)/tot);
    // sensibilidad: si en los datos los pi sobreviven un 10 % más (fila pi->mu x 1,1), con la matriz del MC
    for(double s : {0.9, 1.1}){
      double mm[3][3]; for(int a = 0; a < 3; a++)for(int b = 0; b < 3; b++)mm[a][b] = m[a][b];
      const double d = mm[0][1]*(s - 1); mm[0][1] += d; mm[0][0] -= d;
      TVectorD R2(3); for(int r = 0; r < 3; r++){ R2(r) = 0; for(int tr = 0; tr < 3; tr++)R2(r) += mm[tr][r]; }
      printf("    pi que llegan como mu x%.1f en los datos: fracción corregida %.4f (true %.4f, %+.2f %%)\n", s,
             (TMatrixD(TMatrixD::kInverted, M)*R2)(0)/tot, nTr[0]/tot, 100*((TMatrixD(TMatrixD::kInverted, M)*R2)(0)/nTr[0] - 1));
    }
  }
  return 0;
}
