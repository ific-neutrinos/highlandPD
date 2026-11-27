// Medida C (27-nov-2026): parejas K -> mu con la forced mother (cand_fmc_*, cand_fmk_*; highlandPD 310615b). Como C no
// usa el ajuste del momento, se corta en el momento del mu y, de los que sobreviven, se buscan los que tienen madre.
//   Opción 2: mu del BDT (score > cut) en la ventana de p, y madre (la más cercana, fmc, o la de menor chi2 de K, fmk)
//             que cumpla requisitos de K (distancia, chi2 de K frente a p, longitud).
//   Opción 1: madre preseleccionada como K (chi2) y mu con un umbral del BDT más suelto.
// El BDT (10 variables, dentro de la beam slice) se evalúa leyendo el microtree con SetBranchAddress (con TTreeFormula
// salía mal); con oddOnly = 1 solo los eventos impares (la mitad par lo entrenó), x2. Pareja correcta: el mu es signal
// muon y la madre es su madre true. Proyección a los datos de 6 GeV: x 1/1,196.
//
//   root -l -b -q 'test/phase2/measC_fm.C+("xs_v100to105d01c.root", "opt_bdt_final/dataset/weights/bdt0_BDT.weights.xml", 1)'
#include "TFile.h"
#include "TTree.h"
#include "TMVA/Reader.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

int measC_fm(const char* xs, const char* weights, int oddOnly = 1, double cut = 0.19, double pMin = 205, double pMax = 240){
  TFile fx(xs); TTree* t = (TTree*)fx.Get("ana");
  const int N = 1000;
  if(t->GetMaximum("ncand") > N){ printf("ERROR: ncand > %d\n", N); return 1; }
  const bool hasFM = t->GetBranch("cand_fmc_dist");
  Int_t nc, pop[N], lev[N], sig[N], asn5[N], kn[N], nm10[N], ktrue[N];
  Float_t pr[N], chi2mu[N], fs[N], fa[N], asc5[N], kd[N], cnnm[N], mown[N], cosb[N], dstop[N], arc10[N], ipk[N];
  // forced mother: [0] fmc, [1] fmk
  Float_t mDist[2][N], mLen[2][N], mChi2K[2][N], mChi2P[2][N], mDedx[2][N]; Int_t mTrueK[2][N], mTrk[2][N];
  t->SetBranchStatus("*", 0);
  auto I = [&](const char* b, void* a){ t->SetBranchStatus(b, 1); t->SetBranchAddress(b, a); };
  I("ncand", &nc); I("cand_pop", pop); I("cand_level", lev); I("cand_truesignalmuon", sig); I("cand_prange", pr); I("cand_chi2mu", chi2mu);
  I("cand_faces_stop", fs); I("cand_faces_start", fa); I("cand_as_n5", asn5); I("cand_as_charge5", asc5); I("cand_k_dedx", kd); I("cand_k_n", kn);
  I("cand_cnn_michel", cnnm); I("cand_michel_own", mown); I("cand_ar_nmichel10", nm10); I("cand_cosbeam", cosb); I("cand_dedx_stop", dstop);
  I("cand_ar_charge10", arc10); I("cand_intpeak", ipk); I("cand_k_truek", ktrue);
  const char* pre[2] = {"fmc", "fmk"};
  if(hasFM)for(int m = 0; m < 2; m++){
    I(Form("cand_%s_dist", pre[m]), mDist[m]); I(Form("cand_%s_length", pre[m]), mLen[m]); I(Form("cand_%s_chi2k", pre[m]), mChi2K[m]);
    I(Form("cand_%s_chi2p", pre[m]), mChi2P[m]); I(Form("cand_%s_dedx", pre[m]), mDedx[m]); I(Form("cand_%s_truek", pre[m]), mTrueK[m]); I(Form("cand_%s_istrack", pre[m]), mTrk[m]);
  }
  TMVA::Reader rd("!Color:Silent"); float buf[10];
  const char* vs[10] = {"as_mean5", "k_dedx", "k_n", "cnn_michel", "michel_own", "nmich10", "cosbeam", "dedx_stop", "ar_charge10", "intpeak"};
  for(int i = 0; i < 10; i++)rd.AddVariable(vs[i], &buf[i]);
  rd.BookMVA("BDT", weights);
  // candidatos dentro, preseleccionados, en la ventana: score y lo que hace falta para las parejas
  struct C { float score; int sig, ktrue; float d[2], len[2], c2k[2], c2p[2], dedx[2]; int tk[2], trk[2]; };
  std::vector<C> cs;
  for(Long64_t e = oddOnly ? 1 : 0; e < t->GetEntries(); e += oddOnly ? 2 : 1){
    t->GetEntry(e);
    for(int j = 0; j < nc; j++){
      if(pop[j] != 0 || lev[j] < 1 || chi2mu[j] < 0 || chi2mu[j] >= 20 || std::min(fs[j], fa[j]) <= 10 || pr[j] < pMin || pr[j] > pMax)continue;
      buf[0] = asn5[j] > 0 ? asc5[j]/asn5[j] : -1; buf[1] = kd[j]; buf[2] = kn[j]; buf[3] = cnnm[j]; buf[4] = mown[j]; buf[5] = nm10[j];
      buf[6] = cosb[j]; buf[7] = dstop[j]; buf[8] = arc10[j]; buf[9] = ipk[j];
      C c; c.score = rd.EvaluateMVA("BDT"); c.sig = sig[j]; c.ktrue = ktrue[j];
      for(int m = 0; m < 2; m++){ c.d[m] = hasFM ? mDist[m][j] : -1; c.len[m] = hasFM ? mLen[m][j] : -1; c.c2k[m] = hasFM ? mChi2K[m][j] : -1;
        c.c2p[m] = hasFM ? mChi2P[m][j] : -1; c.dedx[m] = hasFM ? mDedx[m][j] : -1; c.tk[m] = hasFM ? mTrueK[m][j] : -1; c.trk[m] = hasFM ? mTrk[m][j] : -1; }
      cs.push_back(c);
    }
  }
  const double w = (oddOnly ? 2. : 1.)/1.196;
  auto report = [&](const char* name, std::function<int(const C&)> pass){   // pass: 0 no, 1 sí (mu); para parejas devuelve el índice de la madre + 1
    double n = 0, ok = 0, nsig = 0;
    for(auto& c : cs){ const int p = pass(c); if(!p)continue; n++; if(c.sig == 1)nsig++;
      if(c.sig == 1 && (p == 1 ? true : p == 2 ? c.tk[0] == 1 : p == 3 ? c.tk[1] == 1 : c.ktrue == 1))ok++; }
    printf("  %-62s %6.0f %6.0f %5.1f %% | datos: %5.0f %5.0f  (mu signal %5.0f)\n", name, n, ok, n > 0 ? 100*ok/n : 0., n*w, ok*w, nsig*w);
  };
  printf("Candidatos dentro, preseleccionados, p de rango en [%.0f, %.0f]: %zu (%s)\n", pMin, pMax, cs.size(), oddOnly ? "eventos impares" : "todos");
  printf("  %-62s %6s %6s %7s | datos: total, correctas\n", "selección", "total", "correc", "pureza");
  report(Form("BDT > %.2f: mu (correcto = signal muon)", cut), [&](const C& c){ return c.score > cut ? 1 : 0; });
  report(Form("BDT > %.2f + candidato a K al inicio (cand_k_*)", cut), [&](const C& c){ return c.score > cut && c.ktrue >= 0 ? 4 : 0; });
  if(!hasFM){ printf("(sin variables de forced mother en este microtree)\n"); return 0; }
  printf("Opción 2: BDT > %.2f y madre con requisitos\n", cut);
  for(int m = 0; m < 2; m++)for(double dmax : {3., 5., 10.}){
    report(Form("%s a < %.0f cm", pre[m], dmax), [&, m, dmax](const C& c){ return c.score > cut && c.d[m] >= 0 && c.d[m] < dmax ? 2 + m : 0; });
    report(Form("%s a < %.0f cm, chi2 K < chi2 p", pre[m], dmax), [&, m, dmax](const C& c){ return c.score > cut && c.d[m] >= 0 && c.d[m] < dmax && c.c2k[m] >= 0 && c.c2k[m] < c.c2p[m] ? 2 + m : 0; });
    report(Form("%s a < %.0f cm, chi2 K < chi2 p, dE/dx > 5", pre[m], dmax), [&, m, dmax](const C& c){ return c.score > cut && c.d[m] >= 0 && c.d[m] < dmax && c.c2k[m] >= 0 && c.c2k[m] < c.c2p[m] && c.dedx[m] > 5 ? 2 + m : 0; });
  }
  printf("Opción 1: madre K preseleccionada (fmk a < 5 cm, chi2 K < chi2 p, traza) y BDT del mu más suelto\n");
  for(double c0 : {-0.2, -0.1, 0., 0.1, cut})
    report(Form("BDT > %.2f", c0), [&, c0](const C& c){ return c.score > c0 && c.d[1] >= 0 && c.d[1] < 5 && c.c2k[1] >= 0 && c.c2k[1] < c.c2p[1] && c.trk[1] == 1 ? 3 : 0; });
  return 0;
}
