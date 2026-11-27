// Medida B simulada (27-nov-2026), con lo decidido hasta hoy: beam slice, calidad del haz con SCE (beam_quality_sce),
// selección BDT final (10 variables, > cut), ajuste del p de rango con la forma del mu de K + la forma del pi de
// K+ -> pi+ pi0 en reposo + exp(pol2), y resultado por pi+ con la fracción de pi del grupo e/mu/pi medida en el TPC
// (opción b: e por la CNN EM del primario, mu por la traza que llega lejos; matriz de confusión del MC). Sustituye a
// mockB.C (26-nov: opt_flat2, BDT anterior, sin calidad del haz, por partícula del grupo).
//
// Papeles, por la entrada del evento en el microtree (v1_00-v1_05): la mitad par entrenó el BDT y no se usa;
// entrada%4 == 1 = "MC" (formas, correcciones, matriz), entrada%4 == 3 = pseudo-datos (solo su espectro, su número de
// partículas del haz y sus clases pi/mu/e reconstruidas; su truth solo para comprobar). Pulls: 10 trozos de los
// pseudo-datos ((entrada/4) % 10). Con v1_06-v1_09 los pseudo-datos serán esos bloques enteros.
//
// R = N_K (ajuste) x r / N_haz, con r = (mu de B true en los eventos con calidad) / (mu de K del espectro seleccionado),
// del MC. N_haz: partículas del grupo con calidad (por partícula del grupo) o por la fracción de pi medida (por pi+).
// Verdad: mu de B true (árbol truth, KOrigin 0-3) en los eventos con calidad del grupo, sobre las partículas del grupo
// o sobre los pi+ true. Error: el del ajuste (el de la fracción de pi, aparte).
//
//   root -l -b -q 'test/phase2/measB_mock.C+("xs_v100to105d01b.root", "opt_bdt_final/dataset/weights/bdt0_BDT.weights.xml")'
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMVA/Reader.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <tuple>
#include <vector>

namespace mbm {
  const double x0 = 150, x1 = 400; const int nb = 50;
  const int kNSub = 10;
  // ajuste: N_K y su error; tK, tP formas normalizadas (MC)
  bool Fit(const TH1D& hD, const TH1D& tK, const TH1D& tP, double& nK, double& eK){
    std::function<double(double*, double*)> fn = [&tK, &tP](double* x, double* p) -> double {
      const int b = std::min(nb, std::max(1, int((x[0] - x0)/5) + 1)); const double u = (x[0] - x0)/100;
      return p[0]*tK.GetBinContent(b) + p[1]*tP.GetBinContent(b) + std::exp(p[2] + p[3]*u + p[4]*u*u); };
    TF1 f("fmock", fn, x0, x1, 5);
    const double tot = std::max(1., hD.Integral());
    f.SetParameters(0.2*tot, 0.02*tot, std::log(hD.GetBinContent(1) + 1), -1, 0);
    f.SetParLimits(0, 0, 2*tot); f.SetParLimits(1, 0, 2*tot);
    TH1D h(hD); const int st = h.Fit(&f, "LQNR0");
    nK = f.GetParameter(0); eK = f.GetParError(0);
    return st == 0;
  }
}

int measB_mock(const char* xs, const char* weights, double cut = 0.19, double emCut = 0.8, double zMu = 480){
  using namespace mbm;
  TFile fx(xs); TTree* t = (TTree*)fx.Get("ana"); TTree* tru = (TTree*)fx.Get("truth");
  // ---- truth: mu de B por evento
  using Key = std::tuple<int, int, int>;
  std::map<Key, int> nBtrue;
  { Int_t run, sub, evt, kor;
    tru->SetBranchStatus("*", 0); for(const char* b : {"run", "subrun", "evt", "truemu_korigin"})tru->SetBranchStatus(b, 1);
    tru->SetBranchAddress("run", &run); tru->SetBranchAddress("subrun", &sub); tru->SetBranchAddress("evt", &evt); tru->SetBranchAddress("truemu_korigin", &kor);
    for(Long64_t e = 0; e < tru->GetEntries(); e++){ tru->GetEntry(e); if(kor >= 0 && kor <= 3)nBtrue[Key{run, sub, evt}]++; } }
  // ---- microtree
  const int N = 1000, ND = 200;
  if(t->GetMaximum("ncand") > N || t->GetMaximum("seltrk_ndau") > ND){ printf("ERROR: tamaño de los arrays\n"); return 1; }
  Int_t run, sub, evt, sp, hs, q, tpdg, nd, nc, pop[N], lev[N], sig[N], asn5[N], kn[N], nm10[N], tp[N], tpp[N], tpr[N];
  Float_t end[4], edir[3], cnn[3], dpos[ND][4], dend[ND][4], ddir[ND][3];
  Float_t pr[N], chi2mu[N], fs[N], fa[N], asc5[N], kd[N], cnnm[N], mown[N], cosb[N], dstop[N], arc10[N], ipk[N];
  t->SetBranchStatus("*", 0);
  auto I = [&](const char* b, void* a){ t->SetBranchStatus(b, 1); t->SetBranchAddress(b, a); };
  I("run", &run); I("subrun", &sub); I("evt", &evt); I("beam_truespecies", &sp); I("hasslice", &hs); I("beam_quality_sce", &q); I("beam_truepdg", &tpdg);
  I("seltrk_ndau", &nd); I("seltrk_endpos", end); I("seltrk_enddir", edir); I("seltrk_CNNscore", cnn);
  I("seltrk_dau_pos", dpos); I("seltrk_dau_endpos", dend); I("seltrk_dau_dir", ddir);
  I("ncand", &nc); I("cand_pop", pop); I("cand_level", lev); I("cand_truesignalmuon", sig); I("cand_prange", pr); I("cand_chi2mu", chi2mu);
  I("cand_faces_stop", fs); I("cand_faces_start", fa); I("cand_as_n5", asn5); I("cand_as_charge5", asc5); I("cand_k_dedx", kd); I("cand_k_n", kn);
  I("cand_cnn_michel", cnnm); I("cand_michel_own", mown); I("cand_ar_nmichel10", nm10); I("cand_cosbeam", cosb); I("cand_dedx_stop", dstop);
  I("cand_ar_charge10", arc10); I("cand_intpeak", ipk); I("cand_truepdg", tp); I("cand_trueparentpdg", tpp); I("cand_trueproc", tpr);
  TMVA::Reader rd("!Color:Silent"); float buf[10];
  const char* vs[10] = {"as_mean5", "k_dedx", "k_n", "cnn_michel", "michel_own", "nmich10", "cosbeam", "dedx_stop", "ar_charge10", "intpeak"};
  for(int i = 0; i < 10; i++)rd.AddVariable(vs[i], &buf[i]);
  rd.BookMVA("BDT", weights);
  // histogramas: [papel 0 MC, 1 PD][especie 0 e/mu/pi, 1 K]; PD también por trozo
  TH1D *hAll[2][2], *hK[2][2], *hP[2][2], *hSub[kNSub][2];
  for(int r = 0; r < 2; r++)for(int s = 0; s < 2; s++){
    hAll[r][s] = new TH1D(Form("a%d%d", r, s), "", nb, x0, x1); hK[r][s] = new TH1D(Form("k%d%d", r, s), "", nb, x0, x1); hP[r][s] = new TH1D(Form("p%d%d", r, s), "", nb, x0, x1); }
  for(int k = 0; k < kNSub; k++)for(int s = 0; s < 2; s++)hSub[k][s] = new TH1D(Form("s%d%d", k, s), "", nb, x0, x1);
  // recuentos: partículas del haz, pi+ true, mu de B true; por papel, especie y trozo (trozo solo PD)
  double nBeam[2][2] = {{0}}, nPiTrue[2] = {0}, nB[2][2] = {{0}}, nBeamSub[kNSub][2] = {{0}}, nPiSub[kNSub] = {0}, nBSub[kNSub][2] = {{0}};
  // clases pi/mu/e: matriz del MC [true][reco] y reco en PD (total y por trozo)
  double mat[3][3] = {{0}}, recoPD[3] = {0}, recoSub[kNSub][3] = {{0}};
  for(Long64_t e = 1; e < t->GetEntries(); e += 2){
    const int role = (e % 4 == 1) ? 0 : 1, k = (e/4) % kNSub;
    t->GetEntry(e);
    if(sp < 0 || sp > 1 || !hs || q != 1)continue;
    nBeam[role][sp]++; if(role)nBeamSub[k][sp]++;
    const int nbt = nBtrue.count(Key{run, sub, evt}) ? nBtrue[Key{run, sub, evt}] : 0;
    nB[role][sp] += nbt; if(role)nBSub[k][sp] += nbt;
    if(sp == 0){
      const int a = std::abs(tpdg), tr = a == 211 ? 0 : a == 13 ? 1 : a == 11 ? 2 : -1;
      if(tr == 0){ nPiTrue[role]++; if(role)nPiSub[k]++; }
      // clase reconstruida (como measB_beampid.C, empalme con las hijas)
      float z1 = end[2];
      for(int i = 0; i < nd; i++){
        double d0 = 0, d1 = 0; for(int c = 0; c < 3; c++){ d0 += std::pow(dpos[i][c] - end[c], 2); d1 += std::pow(dend[i][c] - end[c], 2); }
        d0 = std::sqrt(d0); d1 = std::sqrt(d1);
        if(std::min(d0, d1) > 30 || dpos[i][2] < -900 || dend[i][2] < -900)continue;
        double c = 0; for(int j = 0; j < 3; j++)c += ddir[i][j]*edir[j];
        if(std::abs(c) < 0.95)continue;
        z1 = std::max(z1, d0 < d1 ? dend[i][2] : dpos[i][2]);
      }
      const int rc = cnn[1] > emCut ? 2 : z1 > zMu ? 1 : 0;
      if(role == 0){ if(tr >= 0)mat[tr][rc]++; } else { recoPD[rc]++; recoSub[k][rc]++; }
    }
    for(int j = 0; j < nc; j++){
      if(pop[j] != 0 || lev[j] < 1 || chi2mu[j] < 0 || chi2mu[j] >= 20 || std::min(fs[j], fa[j]) <= 10 || pr[j] < x0 || pr[j] >= x1)continue;
      buf[0] = asn5[j] > 0 ? asc5[j]/asn5[j] : -1; buf[1] = kd[j]; buf[2] = kn[j]; buf[3] = cnnm[j]; buf[4] = mown[j]; buf[5] = nm10[j];
      buf[6] = cosb[j]; buf[7] = dstop[j]; buf[8] = arc10[j]; buf[9] = ipk[j];
      if(rd.EvaluateMVA("BDT") <= cut)continue;
      hAll[role][sp]->Fill(pr[j]); if(role)hSub[k][sp]->Fill(pr[j]);
      if(sig[j] == 1)hK[role][sp]->Fill(pr[j]);
      else if(tp[j] == 211 && tpp[j] == 321 && tpr[j] == 2)hP[role][sp]->Fill(pr[j]);
    }
  }
  // fracción de pi: matriz del MC (P(reco|true)) e inversión con las clases reconstruidas
  TMatrixD M(3, 3); double nTr[3] = {0}; for(int a = 0; a < 3; a++)for(int b = 0; b < 3; b++)nTr[a] += mat[a][b];
  for(int a = 0; a < 3; a++)for(int b = 0; b < 3; b++)M(b, a) = mat[a][b]/nTr[a];
  TMatrixD Minv(TMatrixD::kInverted, M);
  auto fPi = [&](const double* reco){ TVectorD R(3); double tot = 0; for(int b = 0; b < 3; b++){ R(b) = reco[b]; tot += reco[b]; } const TVectorD T = Minv*R; return tot > 0 ? T(0)/tot : 0.; };
  printf("Medida B simulada: BDT > %.2f dentro de la beam slice, calidad del haz con SCE, ajuste mu de K + pi de K->pipi + exp(pol2)\n", cut);
  printf("MC = entrada%%4 == 1, pseudo-datos (PD) = entrada%%4 == 3 (v1_00-v1_05; el BDT se entrenó con la mitad par)\n");
  printf("Fracción de pi del grupo (e: CNN EM > %.2f; mu: final empalmado > %.0f cm): MC true %.4f; PD true %.4f, reconstruida sin corregir %.4f, corregida %.4f\n",
         emCut, zMu, nPiTrue[0]/nBeam[0][0], nPiTrue[1]/nBeam[1][0], recoPD[0]/(recoPD[0] + recoPD[1] + recoPD[2]), fPi(recoPD));
  const char* sn[2] = {"e/mu/pi", "K"};
  for(int s = 0; s < 2; s++){
    TH1D tK(*hK[0][s]), tP(*hP[0][s]);
    const double nKspecMC = tK.Integral(), r = nB[0][s]/nKspecMC;
    tK.Scale(1./nKspecMC); if(tP.Integral() > 0)tP.Scale(1./tP.Integral());
    printf("\n== haz %s: partículas con calidad MC %.0f / PD %.0f; mu de B true MC %.0f / PD %.0f; mu de K en el espectro MC %.0f / PD %.0f; r = %.3f\n",
           sn[s], nBeam[0][s], nBeam[1][s], nB[0][s], nB[1][s], nKspecMC, hK[1][s]->Integral(), r);
    double nK, eK; const bool ok = Fit(*hAll[1][s], tK, tP, nK, eK);
    printf("   ajuste PD: N_K = %.1f +- %.1f (true %.0f: %+.1f %%)%s\n", nK, eK, hK[1][s]->Integral(), 100*(nK/hK[1][s]->Integral() - 1), ok ? "" : " (ajuste no convergió)");
    // por partícula del grupo y por pi+ (solo e/mu/pi)
    for(int per = 0; per < (s == 0 ? 2 : 1); per++){
      const double fpi = per ? fPi(recoPD) : 1.;
      const double den = nBeam[1][s]*fpi, R = nK*r/den, eR = eK*r/den;
      const double Rt = per ? nB[1][s]/nPiTrue[1] : nB[1][s]/nBeam[1][s];
      printf("   %s: R = %.5f +- %.5f (%.1f %%); verdad %.5f; cierre %+.1f %% = %+.2f sigma\n", per ? "por pi+ (fracción de pi medida)" : (s == 0 ? "por partícula del grupo" : "por K del haz"),
             R, eR, 100*eR/R, Rt, 100*(R/Rt - 1), (R - Rt)/eR);
      // pulls
      std::vector<double> pull;
      for(int k = 0; k < kNSub; k++){
        double nk, ek; if(!Fit(*hSub[k][s], tK, tP, nk, ek) || ek <= 0)continue;
        const double fk = per ? fPi(recoSub[k]) : 1.;
        const double Rk = nk*r/(nBeamSub[k][s]*fk), eRk = ek*r/(nBeamSub[k][s]*fk);
        const double Rtk = per ? nBSub[k][s]/nPiSub[k] : nBSub[k][s]/nBeamSub[k][s];
        pull.push_back((Rk - Rtk)/eRk);
      }
      double m = 0, v = 0; for(double x : pull)m += x; m /= std::max<size_t>(1, pull.size()); for(double x : pull)v += (x - m)*(x - m); v = std::sqrt(v/std::max<size_t>(1, pull.size() - 1));
      printf("     pulls en %zu trozos: media %+.2f, anchura %.2f\n", pull.size(), m, v);
    }
  }
  return 0;
}
