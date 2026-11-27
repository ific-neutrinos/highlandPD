// Medida B (27-nov-2026): incertidumbres sistemáticas del ajuste, con la medida simulada de measB_mock.C (MC =
// entrada%4 == 1, pseudo-datos = entrada%4 == 3 de v1_00-v1_05, beam slice, calidad del haz con SCE, BDT). Se repite la
// medida cambiando una cosa cada vez y se da R y su cambio frente a la nominal (por partícula del grupo; el cambio
// relativo es el mismo por pi+, la fracción de pi no depende del ajuste):
//   fondo: exp(pol2) (nominal), exp(pol3), (p-p0)^k exp(-(p-p0)/tau), (1-exp(-(p-p0)/r)) exp(-(p-p0)/tau)
//   rango del ajuste: 150-400 (nominal), 160-400, 150-350, 170-380 MeV/c
//   binning: 5 (nominal), 2,5 y 10 MeV/c
//   pico de K->pi pi0: libre (nominal), fijo a la fracción del MC respecto al mu de K, sin él
//   forma del mu de K: escala de p de rango +-1 % (-+2,3 MeV/c en el pico) y resolución +2 MeV/c (suavizado gaussiano)
//   umbral del BDT: 0,15, 0,19 (nominal), 0,23 (r y las formas se rehacen con el umbral)
// Los cambios incluyen fluctuaciones estadísticas (los pseudo-datos son los mismos, pero el MC y el ajuste no son
// iguales): se comparan con el error estadístico.
//
//   root -l -b -q 'test/phase2/measB_syst.C+("xs_v100to105d01b.root", "opt_bdt_final/dataset/weights/bdt0_BDT.weights.xml")'
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TMVA/Reader.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace mbs {
  struct Cand { int role, sp, cls; float p, score; };   // cls: 0 mu de K, 1 pi de K->pipi, 2 resto
  struct Var { std::string name; int bkg = 0; double lo = 150, hi = 400, bw = 5; int kpi = 0; double scale = 1, smear = 0, cut = 0.19; };
  double Bkg(int m, double x, const double* p){
    const double u = (x - 150)/100;
    switch(m){
      case 0: return std::exp(p[0] + p[1]*u + p[2]*u*u);
      case 1: return std::exp(p[0] + p[1]*u + p[2]*u*u + p[3]*u*u*u);
      case 2: { const double d = x - p[3]; return d > 0 ? std::exp(p[0])*std::pow(d/100, p[1])*std::exp(-d/p[2]) : 0; }
      default: { const double d = x - p[3]; return d > 0 ? std::exp(p[0])*(1 - std::exp(-d/p[1]))*std::exp(-d/p[2]) : 0; }
    }
  }
  const int kNB[4] = {3, 4, 4, 4};
}

int measB_syst(const char* xs, const char* weights){
  using namespace mbs;
  TFile fx(xs); TTree* t = (TTree*)fx.Get("ana"); TTree* tru = (TTree*)fx.Get("truth");
  using Key = std::tuple<int, int, int>;
  std::map<Key, int> nBtrue;
  { Int_t run, sub, evt, kor;
    tru->SetBranchStatus("*", 0); for(const char* b : {"run", "subrun", "evt", "truemu_korigin"})tru->SetBranchStatus(b, 1);
    tru->SetBranchAddress("run", &run); tru->SetBranchAddress("subrun", &sub); tru->SetBranchAddress("evt", &evt); tru->SetBranchAddress("truemu_korigin", &kor);
    for(Long64_t e = 0; e < tru->GetEntries(); e++){ tru->GetEntry(e); if(kor >= 0 && kor <= 3)nBtrue[Key{run, sub, evt}]++; } }
  const int N = 1000;
  Int_t run, sub, evt, sp, hs, q, nc, pop[N], lev[N], sig[N], asn5[N], kn[N], nm10[N], tp[N], tpp[N], tpr[N];
  Float_t pr[N], chi2mu[N], fs[N], fa[N], asc5[N], kd[N], cnnm[N], mown[N], cosb[N], dstop[N], arc10[N], ipk[N];
  t->SetBranchStatus("*", 0);
  auto I = [&](const char* b, void* a){ t->SetBranchStatus(b, 1); t->SetBranchAddress(b, a); };
  I("run", &run); I("subrun", &sub); I("evt", &evt); I("beam_truespecies", &sp); I("hasslice", &hs); I("beam_quality_sce", &q);
  I("ncand", &nc); I("cand_pop", pop); I("cand_level", lev); I("cand_truesignalmuon", sig); I("cand_prange", pr); I("cand_chi2mu", chi2mu);
  I("cand_faces_stop", fs); I("cand_faces_start", fa); I("cand_as_n5", asn5); I("cand_as_charge5", asc5); I("cand_k_dedx", kd); I("cand_k_n", kn);
  I("cand_cnn_michel", cnnm); I("cand_michel_own", mown); I("cand_ar_nmichel10", nm10); I("cand_cosbeam", cosb); I("cand_dedx_stop", dstop);
  I("cand_ar_charge10", arc10); I("cand_intpeak", ipk); I("cand_truepdg", tp); I("cand_trueparentpdg", tpp); I("cand_trueproc", tpr);
  TMVA::Reader rd("!Color:Silent"); float buf[10];
  const char* vs[10] = {"as_mean5", "k_dedx", "k_n", "cnn_michel", "michel_own", "nmich10", "cosbeam", "dedx_stop", "ar_charge10", "intpeak"};
  for(int i = 0; i < 10; i++)rd.AddVariable(vs[i], &buf[i]);
  rd.BookMVA("BDT", weights);
  std::vector<Cand> cs; double nBeam[2][2] = {{0}}, nB[2][2] = {{0}};
  for(Long64_t e = 1; e < t->GetEntries(); e += 2){
    const int role = (e % 4 == 1) ? 0 : 1;
    t->GetEntry(e);
    if(sp < 0 || sp > 1 || !hs || q != 1)continue;
    nBeam[role][sp]++; { auto it = nBtrue.find(Key{run, sub, evt}); if(it != nBtrue.end())nB[role][sp] += it->second; }
    for(int j = 0; j < nc; j++){
      if(pop[j] != 0 || lev[j] < 1 || chi2mu[j] < 0 || chi2mu[j] >= 20 || std::min(fs[j], fa[j]) <= 10 || pr[j] < 140 || pr[j] >= 410)continue;
      buf[0] = asn5[j] > 0 ? asc5[j]/asn5[j] : -1; buf[1] = kd[j]; buf[2] = kn[j]; buf[3] = cnnm[j]; buf[4] = mown[j]; buf[5] = nm10[j];
      buf[6] = cosb[j]; buf[7] = dstop[j]; buf[8] = arc10[j]; buf[9] = ipk[j];
      const float sc = rd.EvaluateMVA("BDT"); if(sc <= 0.10)continue;
      const int cls = sig[j] == 1 ? 0 : (tp[j] == 211 && tpp[j] == 321 && tpr[j] == 2) ? 1 : 2;
      cs.push_back({role, sp, cls, pr[j], sc});
    }
  }
  // una variante: R por partícula del grupo y su error, para la especie s
  auto measure = [&](const Var& v, int s, double& R, double& eR, double& nK, double& eK) -> bool {
    const int nbins = int((v.hi - v.lo)/v.bw + 0.5);
    TH1D hD("hD", "", nbins, v.lo, v.hi), tK("tK", "", nbins, v.lo, v.hi), tP("tP", "", nbins, v.lo, v.hi);
    TRandom3 rnd(12345);
    double nKspec = 0, nPspec = 0;
    for(auto& c : cs){
      if(c.sp != s || c.score <= v.cut)continue;
      if(c.role == 1){ hD.Fill(c.p); continue; }
      if(c.cls == 0){ const double p = c.p*v.scale + (v.smear > 0 ? rnd.Gaus(0, v.smear) : 0.); tK.Fill(p); if(c.p >= 150 && c.p < 400)nKspec++; }
      else if(c.cls == 1){ tP.Fill(c.p); if(c.p >= 150 && c.p < 400)nPspec++; }
    }
    // r con el espectro nominal (150-400) del MC: N_K ajustado se refiere a los mu de K de la forma, que se normaliza a
    // la fracción de mu de K del MC dentro del rango de esta variante
    const double fIn = tK.Integral()/std::max(1., nKspec);   // mu de K de la forma dentro del rango / en 150-400
    const double r = nB[0][s]/nKspec;
    const double kpiRatio = nPspec/std::max(1., nKspec);
    if(tK.Integral() <= 0)return false;
    tK.Scale(1./tK.Integral()); if(tP.Integral() > 0)tP.Scale(1./tP.Integral());
    const int nbk = kNB[v.bkg];
    std::function<double(double*, double*)> fn = [&](double* x, double* p) -> double {
      const int b = std::min(nbins, std::max(1, int((x[0] - v.lo)/v.bw) + 1));
      const double npi = v.kpi == 0 ? p[1] : v.kpi == 1 ? kpiRatio*p[0]/std::max(1e-9, fIn) : 0.;
      return p[0]*tK.GetBinContent(b) + npi*tP.GetBinContent(b) + Bkg(v.bkg, x[0], p + 2); };
    TF1 f("fs", fn, v.lo, v.hi, 2 + nbk);
    const double tot = std::max(1., hD.Integral());
    f.SetParameter(0, 0.2*tot); f.SetParLimits(0, 0, 2*tot);
    f.SetParameter(1, 0.02*tot); if(v.kpi == 0)f.SetParLimits(1, 0, 2*tot); else f.FixParameter(1, 0);
    const double l = std::log(std::max(1., hD.GetBinContent(1)));
    if(v.bkg == 0){ f.SetParameter(2, l); f.SetParameter(3, -1); f.SetParameter(4, 0); }
    else if(v.bkg == 1){ f.SetParameter(2, l); f.SetParameter(3, -1); f.SetParameter(4, 0); f.SetParameter(5, 0); }
    else if(v.bkg == 2){ f.SetParameter(2, l + 1); f.SetParameter(3, 0.5); f.SetParameter(4, 80); f.SetParameter(5, 100); f.SetParLimits(3, 0, 5); f.SetParLimits(4, 10, 1000); f.SetParLimits(5, 0, v.lo - 1); }
    else { f.SetParameter(2, l + 0.3); f.SetParameter(3, 20); f.SetParameter(4, 120); f.SetParameter(5, 120); f.SetParLimits(3, 1, 200); f.SetParLimits(4, 10, 1000); f.SetParLimits(5, 0, v.lo - 1); }
    const int st = hD.Fit(&f, "LQNR0");
    nK = f.GetParameter(0)/fIn; eK = f.GetParError(0)/fIn;   // mu de K en 150-400
    R = nK*r/nBeam[1][s]; eR = eK*r/nBeam[1][s];
    return st == 0;
  };
  std::vector<Var> vars;
  auto add = [&](Var v){ vars.push_back(v); };
  { Var v; v.name = "nominal"; add(v); }
  { Var v; v.name = "fondo exp(pol3)"; v.bkg = 1; add(v); }
  { Var v; v.name = "fondo (p-p0)^k exp"; v.bkg = 2; add(v); }
  { Var v; v.name = "fondo (1-exp) exp"; v.bkg = 3; add(v); }
  { Var v; v.name = "rango 160-400"; v.lo = 160; add(v); }
  { Var v; v.name = "rango 150-350"; v.hi = 350; add(v); }
  { Var v; v.name = "rango 170-380"; v.lo = 170; v.hi = 380; add(v); }
  { Var v; v.name = "bins de 2,5"; v.bw = 2.5; add(v); }
  { Var v; v.name = "bins de 10"; v.bw = 10; add(v); }
  { Var v; v.name = "K->pipi fijo a la fracción del MC"; v.kpi = 1; add(v); }
  { Var v; v.name = "sin K->pipi"; v.kpi = 2; add(v); }
  { Var v; v.name = "escala del mu +1 %"; v.scale = 1.01; add(v); }
  { Var v; v.name = "escala del mu -1 %"; v.scale = 0.99; add(v); }
  { Var v; v.name = "resolución del mu +2 MeV/c"; v.smear = 2; add(v); }
  { Var v; v.name = "BDT > 0,15"; v.cut = 0.15; add(v); }
  { Var v; v.name = "BDT > 0,23"; v.cut = 0.23; add(v); }
  const char* sn[2] = {"e/mu/pi", "K"};
  printf("Sistemáticos del ajuste de la medida B (medida simulada: MC = entrada%%4 == 1, PD = entrada%%4 == 3; beam slice, calidad SCE)\n");
  for(int s = 0; s < 2; s++){
    const double Rt = nB[1][s]/nBeam[1][s];
    double R0 = 0, eR0 = 0;
    printf("\n== haz %s: verdad R = %.5f (por partícula del grupo)\n", sn[s], Rt);
    printf("  %-36s %9s %9s %8s %9s %9s\n", "variante", "R", "error", "N_K", "cambio", "cierre");
    for(size_t i = 0; i < vars.size(); i++){
      double R, eR, nK, eK; const bool ok = measure(vars[i], s, R, eR, nK, eK);
      if(i == 0){ R0 = R; eR0 = eR; }
      printf("  %-36s %9.5f %9.5f %8.1f %+8.1f %% %+8.1f %%%s\n", vars[i].name.c_str(), R, eR, nK, 100*(R/R0 - 1), 100*(R/Rt - 1), ok ? "" : "  (no convergió)");
    }
    printf("  (error estadístico nominal: %.1f %%)\n", 100*eR0/R0);
  }
  return 0;
}
