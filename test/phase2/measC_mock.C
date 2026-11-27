// Medida C simulada (27-nov-2026): sección eficaz de dispersión de K+ (por encima de un ángulo theta) en K+ que se paran
// y se etiquetan por su mu (K+ -> mu+ nu en reposo). Necesita las variables de highlandPD d6a3539 (forced mother,
// dispersión reconstruida, dispersiones true): xs_*d01c.
//
// Método (blanco fino, una dispersión por K como mucho): sigma(> theta) = N_disp / (n_Ar * suma de L), con N_disp el número
// de K con una dispersión > theta y L la longitud de cada K. n_Ar = 1,396 g/cm3 * N_A / 39,95 g/mol = 2,104e22 cm^-3.
//   Reco (pseudo-datos): parejas = mu del BDT (dentro de la beam slice, p de rango en [pMin, pMax]) con forced mother
//     más cercana (fmc) que es traza, a < dMax cm, con chi2 de K < chi2 de p y longitud > lMin; dispersión si el mayor
//     cambio de dirección a lo largo de la K (angmax, a más de rrMin cm de la parada) o el ángulo con el trozo anterior
//     (prevangle, si hay trozo a < 3 cm) supera theta.
//   Corrección (MC): c(theta) = sigma_true / sigma_reco en el MC, con sigma_true de las K madre true de los signal muons de
//     las mismas parejas correctas... (aquí: de todos los signal muons con K madre y trayectoria, árbol truth, con su
//     longitud true y si tienen algún paso > theta: n10, n20, n30).
//   sigma_medida = sigma_reco(PD) * c(theta); verdad = sigma_true(PD).
// Papeles: con oddOnly, MC = entrada%4 == 1 y PD = entrada%4 == 3; sin él (bloques independientes del entrenamiento del
// BDT), MC = entrada par y PD = impar. También se imprime la estadística esperada a la de los datos de 6 GeV.
//
//   root -l -b -q 'test/phase2/measC_mock.C+("xs_v100to105d01c.root", "opt_bdt_final/dataset/weights/bdt0_BDT.weights.xml", 1)'
#include "TFile.h"
#include "TTree.h"
#include "TMVA/Reader.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <tuple>
#include <vector>

int measC_mock(const char* xs, const char* weights, int oddOnly = 1, double cut = 0.19, double pMin = 205, double pMax = 240,
               double dMax = 5, double lMin = 5, double rrMin = 2, double nEvData = 306000){
  const double nAr = 1.396*6.02214e23/39.95;   // cm^-3
  const double kTh[3] = {10, 20, 30};
  TFile fx(xs); TTree* t = (TTree*)fx.Get("ana"); TTree* tru = (TTree*)fx.Get("truth");
  if(!t->GetBranch("cand_fmc_angmax") || !tru->GetBranch("truek_n20")){ printf("ERROR: faltan las variables de d6a3539 (forced mother, dispersiones)\n"); return 1; }
  auto roleOf = [&](Long64_t e){ return oddOnly ? (e % 4 == 1 ? 0 : e % 4 == 3 ? 1 : -1) : int(e % 2); };
  // ---- truth: K madre de los signal muons (todas las capas): longitud y dispersiones; papel por evento
  using Key = std::tuple<int, int, int>;
  std::map<Key, int> role;
  { Int_t run, sub, evt; t->SetBranchStatus("*", 0); for(const char* b : {"run", "subrun", "evt"})t->SetBranchStatus(b, 1);
    t->SetBranchAddress("run", &run); t->SetBranchAddress("subrun", &sub); t->SetBranchAddress("evt", &evt);
    for(Long64_t e = 0; e < t->GetEntries(); e++){ t->GetEntry(e); role[Key{run, sub, evt}] = roleOf(e); } }
  double tL[2] = {0}, tN[2][3] = {{0}}, tK[2] = {0};   // [papel]: suma de longitudes, K con paso > theta, K
  { Int_t run, sub, evt, kor, n10, n20, n30, ntrj; Float_t len;
    tru->SetBranchStatus("*", 0);
    for(const char* b : {"run", "subrun", "evt", "truemu_korigin", "truek_length", "truek_n10", "truek_n20", "truek_n30", "truek_ntrj"})tru->SetBranchStatus(b, 1);
    tru->SetBranchAddress("run", &run); tru->SetBranchAddress("subrun", &sub); tru->SetBranchAddress("evt", &evt); tru->SetBranchAddress("truemu_korigin", &kor);
    tru->SetBranchAddress("truek_length", &len); tru->SetBranchAddress("truek_n10", &n10); tru->SetBranchAddress("truek_n20", &n20); tru->SetBranchAddress("truek_n30", &n30);
    tru->SetBranchAddress("truek_ntrj", &ntrj);
    for(Long64_t e = 0; e < tru->GetEntries(); e++){
      tru->GetEntry(e); auto it = role.find(Key{run, sub, evt}); if(it == role.end() || it->second < 0)continue;
      if(kor < 0 || len <= 0 || ntrj < 2)continue;
      const int r = it->second; tK[r]++; tL[r] += len; const int n[3] = {n10, n20, n30};
      for(int i = 0; i < 3; i++)if(n[i] > 0)tN[r][i]++;
    } }
  // ---- reco: parejas
  const int N = 1000;
  Int_t nc, pop[N], lev[N], sig[N], asn5[N], kn[N], nm10[N], fTrk[N], fTK[N];
  Float_t pr[N], chi2mu[N], fs[N], fa[N], asc5[N], kd[N], cnnm[N], mown[N], cosb[N], dstop[N], arc10[N], ipk[N];
  Float_t fD[N], fL[N], fC2K[N], fC2P[N], fAng[N], fAngRR[N], fPD[N], fPA[N];
  t->SetBranchStatus("*", 0);
  auto I = [&](const char* b, void* a){ t->SetBranchStatus(b, 1); t->SetBranchAddress(b, a); };
  I("ncand", &nc); I("cand_pop", pop); I("cand_level", lev); I("cand_truesignalmuon", sig); I("cand_prange", pr); I("cand_chi2mu", chi2mu);
  I("cand_faces_stop", fs); I("cand_faces_start", fa); I("cand_as_n5", asn5); I("cand_as_charge5", asc5); I("cand_k_dedx", kd); I("cand_k_n", kn);
  I("cand_cnn_michel", cnnm); I("cand_michel_own", mown); I("cand_ar_nmichel10", nm10); I("cand_cosbeam", cosb); I("cand_dedx_stop", dstop);
  I("cand_ar_charge10", arc10); I("cand_intpeak", ipk);
  I("cand_fmc_dist", fD); I("cand_fmc_length", fL); I("cand_fmc_chi2k", fC2K); I("cand_fmc_chi2p", fC2P); I("cand_fmc_istrack", fTrk);
  I("cand_fmc_truek", fTK); I("cand_fmc_angmax", fAng); I("cand_fmc_angmax_rr", fAngRR); I("cand_fmc_prevdist", fPD); I("cand_fmc_prevangle", fPA);
  TMVA::Reader rd("!Color:Silent"); float buf[10];
  const char* vs[10] = {"as_mean5", "k_dedx", "k_n", "cnn_michel", "michel_own", "nmich10", "cosbeam", "dedx_stop", "ar_charge10", "intpeak"};
  for(int i = 0; i < 10; i++)rd.AddVariable(vs[i], &buf[i]);
  rd.BookMVA("BDT", weights);
  double rN[2] = {0}, rOK[2] = {0}, rL[2] = {0}, rS[2][3] = {{0}}, rSok[2][3] = {{0}}, nEv[2] = {0};
  for(Long64_t e = 0; e < t->GetEntries(); e++){
    const int r = roleOf(e); if(r < 0)continue;
    t->GetEntry(e); nEv[r]++;
    for(int j = 0; j < nc; j++){
      if(pop[j] != 0 || lev[j] < 1 || chi2mu[j] < 0 || chi2mu[j] >= 20 || std::min(fs[j], fa[j]) <= 10 || pr[j] < pMin || pr[j] > pMax)continue;
      buf[0] = asn5[j] > 0 ? asc5[j]/asn5[j] : -1; buf[1] = kd[j]; buf[2] = kn[j]; buf[3] = cnnm[j]; buf[4] = mown[j]; buf[5] = nm10[j];
      buf[6] = cosb[j]; buf[7] = dstop[j]; buf[8] = arc10[j]; buf[9] = ipk[j];
      if(rd.EvaluateMVA("BDT") <= cut)continue;
      if(fD[j] < 0 || fD[j] > dMax || fTrk[j] != 1 || fL[j] < lMin || fC2K[j] < 0 || fC2K[j] >= fC2P[j])continue;
      const bool ok = sig[j] == 1 && fTK[j] == 1;
      rN[r]++; if(ok)rOK[r]++; rL[r] += fL[j];
      const double ang = std::max(fAngRR[j] > rrMin ? fAng[j] : -1.f, (fPD[j] >= 0 && fPD[j] < 3) ? fPA[j] : -1.f);
      for(int i = 0; i < 3; i++)if(ang > kTh[i]){ rS[r][i]++; if(ok)rSok[r][i]++; }
    }
  }
  const char* rn[2] = {"MC", "PD"};
  printf("Medida C simulada (%s): parejas = BDT > %.2f, p en [%.0f, %.0f], fmc traza a < %.0f cm, chi2 K < chi2 p, L > %.0f cm\n",
         oddOnly ? "MC = entrada%4 == 1, PD = entrada%4 == 3" : "MC = par, PD = impar", cut, pMin, pMax, dMax, lMin);
  for(int r = 0; r < 2; r++)printf("  %s: %.0f eventos; truth: %.0f K madre, suma L %.0f cm | reco: %.0f parejas (%.0f correctas, %.0f %%), suma L %.0f cm\n",
                                   rn[r], nEv[r], tK[r], tL[r], rN[r], rOK[r], 100*rOK[r]/std::max(1., rN[r]), rL[r]);
  printf("  theta | sigma true MC, PD (mb)        | reco MC, PD: K dispersadas (correctas) y sigma_reco (mb) | c = true/reco (MC) | sigma medida PD (mb) | cierre | parejas dispersadas en datos\n");
  for(int i = 0; i < 3; i++){
    double sT[2], sR[2];
    for(int r = 0; r < 2; r++){ sT[r] = tL[r] > 0 ? tN[r][i]/(nAr*tL[r])*1e27 : 0; sR[r] = rL[r] > 0 ? rS[r][i]/(nAr*rL[r])*1e27 : 0; }
    const double c = sR[0] > 0 ? sT[0]/sR[0] : 0, sM = sR[1]*c, eM = rS[1][i] > 0 ? sM/std::sqrt(rS[1][i]) : 0;
    printf("  %4.0f  | %7.1f %7.1f | %4.0f (%3.0f) %7.1f ; %4.0f (%3.0f) %7.1f | %5.2f | %7.1f +- %5.1f | %+6.1f %% | %5.0f\n", kTh[i], sT[0], sT[1],
           rS[0][i], rSok[0][i], sR[0], rS[1][i], rSok[1][i], sR[1], c, sM, eM, sT[1] > 0 ? 100*(sM/sT[1] - 1) : 0., rS[1][i]*nEvData/std::max(1., nEv[1]));
  }
  return 0;
}
