// Medida B simulada (26-nov-2026): mu+ de K+ en reposo (capas A + B no A: cadena del haz que dispara) por partícula
// del haz de cada especie, medida sobre pseudo-datos de MC con todo lo demás sacado de otra muestra de MC.
//
//   root -l -b -q 'test/phase2/mockB.C+("opt_flat2_v100to105d01.root", "xs_v100to105d01.root", "<pesos BDT dentro>", "<pesos BDT fuera>", "<dir de figuras>")'
//
// Papeles (por ahora, hasta tener v1_06-v1_09), por cuartos de la entrada del evento en v1_00-v1_05: entrada%4 = 0 y 2
// (la mitad 0) entrenan el BDT y no se usan aquí (26-nov: usar la mitad de entrenamiento como MC infla la eficiencia por
// sobreentrenamiento: -34 % fuera); entrada%4 == 1 = "MC": plantillas del espectro de p de rango de la señal (todos los mu de K en
// reposo seleccionados, capas A, B no A y C no B) y del fondo, eficiencia por especie y población (mu de B seleccionados
// sobre mu de B true de esa especie) y fracción de pileup (C no B) en el pico seleccionado. entrada%4 == 3 = pseudo-datos: solo
// se usa su espectro y su número de partículas del haz por especie; su truth solo para comprobar.
// Selección: BDT de opt_bdt*.C (variables y umbrales abajo) con la preselección suelta de opt_bdt.C, p de rango 150-400.
// Especie del haz: la true (en el MC no hay emulación de los Cherenkov; en datos, el PID del BI con su pureza).
// Ajustes del espectro de la mitad 1 (verosimilitud de Poisson, bins de 5 MeV/c): (T) plantillas de la mitad 0 con las
// dos normalizaciones libres; (G) gaussiana + exp(pol2), sin plantillas. N_B = N_pico (1 - f_pileup) / eficiencia.
// Cierre: N_B medido frente a los mu de B true de la mitad 1; pulls en 10 submuestras independientes de la mitad 1.
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TMath.h"
#include "TMVA/Reader.h"
#include <cmath>
#include <functional>
#include <algorithm>
#include <cstdio>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace mockb {
  const char* kSp[3] = {"e/mu/pi", "K", "p"};
  const char* kPop[2] = {"dentro", "fuera"};
  const double x0 = 150, x1 = 400, bw = 5;
  const int nb = int((x1 - x0)/bw + 0.5);
  const int kNSub = 10;   // submuestras de los pseudo-datos para los pulls ((entrada/4) % 10)
  // papel de un evento: 0 MC (entrada%4 == 1), 1 pseudo-datos (entrada%4 == 3), -1 entrenamiento del BDT
  inline int Role(Long64_t e){ return e % 4 == 1 ? 0 : e % 4 == 3 ? 1 : -1; }
  struct FitOut { double n, en; bool ok; };
  // (T) plantillas: f(x) = Ns*Ts(x) + Nb*Tb(x), Ts y Tb normalizadas a 1
  FitOut FitTemplates(TH1D* h, TH1D* ts, TH1D* tb){
    std::function<double(double*, double*)> fn = [ts, tb](double* x, double* p) -> double {
      const int b = std::min(nb, std::max(1, int((x[0] - x0)/bw) + 1));
      return p[0]*ts->GetBinContent(b) + p[1]*tb->GetBinContent(b); };
    TF1 f(Form("ft_%s", h->GetName()), fn, x0, x1, 2);
    const double tot = h->Integral();
    f.SetParameters(0.1*tot, 0.9*tot); f.SetParLimits(0, 0, 3*tot + 10); f.SetParLimits(1, 0, 3*tot + 10);
    const int st = h->Fit(&f, "LQNR0");
    return {f.GetParameter(0), f.GetParError(0), st == 0};
  }
  // (G) gaussiana + exp(pol2) (como mu_cand_fit.C)
  FitOut FitGaus(TH1D* h){
    TF1 f(Form("fg_%s", h->GetName()), Form("[0]*%g*TMath::Gaus(x,[1],[2],true) + exp([3] + [4]*(x-%g)/100 + [5]*((x-%g)/100)^2)", bw, x0, x0), x0, x1);
    const double tot = h->Integral();
    f.SetParameters(0.1*tot, 228, 9, std::log(std::max(1., h->GetBinContent(1))), -1, 0);
    f.SetParLimits(0, 0, 3*tot + 10); f.SetParLimits(1, 215, 240); f.SetParLimits(2, 4, 20);
    const int st = h->Fit(&f, "LQNR0");
    return {f.GetParameter(0), f.GetParError(0), st == 0};
  }
}

int mockB(const char* flat, const char* xs, const char* wIn, const char* wOut, const char* figdir,
          double cutIn = 0.20, double cutOut = 0.29,
          const char* varsIn = "as_mean5,k_dedx,k_n,cnn_michel,michel_own,nmich10,cosbeam,dedx_stop",
          const char* varsOut = "as_mean5,k_dedx,k_n,cnn_michel,michel_own,nmich10,cosbeam,dedx_stop,axis,clear"){
  using namespace mockb;
  gStyle->SetOptStat(0); gSystem->mkdir(figdir, true);
  // ---- eventos: especie true y mitad, por entrada de ana; y mapa (run, subrun, evt) -> entrada
  TFile fx(xs); TTree* ana = (TTree*)fx.Get("ana"); TTree* tru = (TTree*)fx.Get("truth");
  Int_t sp, run, sub, evt;
  ana->SetBranchStatus("*", 0);
  for(const char* b : {"beam_truespecies", "run", "subrun", "evt"})ana->SetBranchStatus(b, 1);
  ana->SetBranchAddress("beam_truespecies", &sp); ana->SetBranchAddress("run", &run); ana->SetBranchAddress("subrun", &sub); ana->SetBranchAddress("evt", &evt);
  const Long64_t nEv = ana->GetEntries();
  std::vector<int> evSp(nEv);
  std::map<std::tuple<int, int, int>, Long64_t> evIdx;
  double nBeam[3][2] = {{0}}, nBeamSub[3][kNSub] = {{0}};
  for(Long64_t e = 0; e < nEv; e++){
    ana->GetEntry(e); evSp[e] = sp; evIdx[{run, sub, evt}] = e;
    if(sp < 0 || sp > 2)continue;
    const int r = Role(e); if(r < 0)continue;
    nBeam[sp][r]++; if(r == 1)nBeamSub[sp][(e/4) % kNSub]++;
  }
  // ---- truth: mu de B (KOrigin 0-3) por especie y mitad (y submuestra)
  Int_t kor;
  tru->SetBranchStatus("*", 0);
  for(const char* b : {"truemu_korigin", "run", "subrun", "evt"})tru->SetBranchStatus(b, 1);
  tru->SetBranchAddress("truemu_korigin", &kor); tru->SetBranchAddress("run", &run); tru->SetBranchAddress("subrun", &sub); tru->SetBranchAddress("evt", &evt);
  double nTrueB[3][2] = {{0}}, nTrueBSub[3][kNSub] = {{0}};
  for(Long64_t i = 0; i < tru->GetEntries(); i++){
    tru->GetEntry(i);
    if(kor < 0 || kor > 3)continue;
    auto it = evIdx.find({run, sub, evt}); if(it == evIdx.end())continue;
    const Long64_t e = it->second; const int s = evSp[e]; if(s < 0 || s > 2)continue;
    const int r = Role(e); if(r < 0)continue;
    nTrueB[s][r]++; if(r == 1)nTrueBSub[s][(e/4) % kNSub]++;
  }
  // ---- candidatos seleccionados por el BDT
  TFile ff(flat); TTree* t = (TTree*)ff.Get("flat");
  auto split = [](const char* l){ std::vector<std::string> v; std::string s = std::string(l) + ","; for(size_t a = 0, b; (b = s.find(',', a)) != std::string::npos; a = b + 1)if(b > a)v.push_back(s.substr(a, b - a)); return v; };
  const std::vector<std::string> vn[2] = {split(varsIn), split(varsOut)};
  std::map<std::string, Float_t> val;
  for(int p = 0; p < 2; p++)for(auto& x : vn[p])val[x] = 0;
  for(const char* x : {"pop", "prange", "chi2mu", "faces", "axis"})val[x] = 0;
  Int_t cat; Long64_t entry;
  t->SetBranchStatus("*", 0);
  for(auto& kv : val){ t->SetBranchStatus(kv.first.c_str(), 1); t->SetBranchAddress(kv.first.c_str(), &kv.second); }
  t->SetBranchStatus("cat", 1); t->SetBranchAddress("cat", &cat); t->SetBranchStatus("entry", 1); t->SetBranchAddress("entry", &entry);
  std::unique_ptr<TMVA::Reader> rd[2];
  std::vector<float> buf[2];
  for(int p = 0; p < 2; p++){
    rd[p].reset(new TMVA::Reader("!Color:Silent")); buf[p].resize(vn[p].size());
    for(size_t i = 0; i < vn[p].size(); i++)rd[p]->AddVariable(vn[p][i].c_str(), &buf[p][i]);
    rd[p]->BookMVA("BDT", p == 0 ? wIn : wOut);
  }
  // histogramas: [especie][población]; mitad 0: señal (mu de K) y fondo; mitad 1: total (pseudo-datos) y por submuestra
  auto H = [](const char* n){ TH1D* h = new TH1D(n, "", nb, x0, x1); h->Sumw2(); return h; };
  TH1D *hS0[3][2], *hB0[3][2], *hD1[3][2], *hSub[3][2][kNSub];
  double selB0[3][2] = {{0}}, selK0[3][2] = {{0}}, selC0[3][2] = {{0}}, selK1[3][2] = {{0}}, selB1[3][2] = {{0}};
  for(int s = 0; s < 3; s++)for(int p = 0; p < 2; p++){
    hS0[s][p] = H(Form("hS0_%d_%d", s, p)); hB0[s][p] = H(Form("hB0_%d_%d", s, p)); hD1[s][p] = H(Form("hD1_%d_%d", s, p));
    for(int k = 0; k < kNSub; k++)hSub[s][p][k] = H(Form("hSub_%d_%d_%d", s, p, k));
  }
  for(Long64_t i = 0; i < t->GetEntries(); i++){
    t->GetEntry(i);
    const int p = (int)val["pop"]; const double pr = val["prange"];
    if(p < 0 || p > 1 || pr < x0 || pr >= x1)continue;
    if(val["chi2mu"] >= 20 || (p == 0 && val["faces"] <= 10) || (p == 1 && val["axis"] >= 150))continue;
    const int s = evSp[entry]; if(s < 0 || s > 2)continue;
    for(size_t k = 0; k < vn[p].size(); k++)buf[p][k] = val[vn[p][k]];
    if(rd[p]->EvaluateMVA("BDT") <= (p == 0 ? cutIn : cutOut))continue;
    const bool isK = cat <= 2, isB = cat <= 1;
    const int role = Role(entry); if(role < 0)continue;
    if(role == 0){ (isK ? hS0 : hB0)[s][p]->Fill(pr); if(isB)selB0[s][p]++; if(isK)selK0[s][p]++; if(cat == 2)selC0[s][p]++; }
    else{ hD1[s][p]->Fill(pr); hSub[s][p][(entry/4) % kNSub]->Fill(pr); if(isK)selK1[s][p]++; if(isB)selB1[s][p]++; }
  }
  // ---- medida
  printf("Medida B simulada. Selección: BDT > %.2f dentro, > %.2f fuera. MC = entrada%%4 == 1, pseudo-datos = entrada%%4 == 3 (el BDT se entrenó con la mitad par).\n", cutIn, cutOut);
  printf("Partículas del haz (MC / pseudo-datos): e/mu/pi %.0f / %.0f, K %.0f / %.0f, p %.0f / %.0f\n\n",
         nBeam[0][0], nBeam[0][1], nBeam[1][0], nBeam[1][1], nBeam[2][0], nBeam[2][1]);
  TCanvas cv("cv", "", 1400, 900); cv.Divide(3, 2);
  for(int s = 0; s < 3; s++){
    const double rTrue = nTrueB[s][1]/nBeam[s][1];
    printf("== haz %s: verdad (pseudo-datos) %.0f mu de B en %.0f partículas: R = %.5f (MC: %.5f)\n", kSp[s], nTrueB[s][1], nBeam[s][1], rTrue, nTrueB[s][0]/nBeam[s][0]);
    double comb[2][2] = {{0}};   // [ajuste][suma de pesos, suma de pesos*valor]
    for(int p = 0; p < 2; p++){
      const double eff = nTrueB[s][0] > 0 ? selB0[s][p]/nTrueB[s][0] : 0, fpu = selK0[s][p] > 0 ? selC0[s][p]/selK0[s][p] : 0;
      if(selK0[s][p] < 5 || eff <= 0){ printf("   %s: sin estadística en la mitad 0 (%.0f mu de K seleccionados)\n", kPop[p], selK0[s][p]); continue; }
      TH1D ts(*hS0[s][p]), tb(*hB0[s][p]);
      ts.Scale(1./ts.Integral()); if(tb.Integral() > 0)tb.Scale(1./tb.Integral());
      const FitOut fT = FitTemplates(hD1[s][p], &ts, &tb), fG = FitGaus(hD1[s][p]);
      printf("   %s: eficiencia %.4f, f_pileup %.3f | pico seleccionado true (pseudo-datos) %.0f\n", kPop[p], eff, fpu, selK1[s][p]);
      const FitOut* fr[2] = {&fT, &fG}; const char* fn[2] = {"plantillas", "gauss+exp"};
      for(int a = 0; a < 2; a++){
        const double nB = fr[a]->n*(1 - fpu)/eff, enB = fr[a]->en*(1 - fpu)/eff, R = nB/nBeam[s][1], eR = enB/nBeam[s][1];
        printf("     (%s) N_pico %.0f ± %.0f (%s) -> R = %.5f ± %.5f (%.1f %%): cierre %+.1f %% = %+.2f sigma\n", fn[a], fr[a]->n, fr[a]->en,
               fr[a]->ok ? "ok" : "FALLO", R, eR, 100*eR/R, 100*(R/rTrue - 1), (R - rTrue)/eR);
        if(fr[a]->ok && eR > 0){ comb[a][0] += 1/(eR*eR); comb[a][1] += R/(eR*eR); }
      }
      // pulls en las submuestras (plantillas y gauss+exp)
      for(int a = 0; a < 2; a++){
        double sum = 0, sum2 = 0; int n = 0;
        for(int k = 0; k < kNSub; k++){
          if(nTrueBSub[s][k] < 1)continue;
          const FitOut f = a == 0 ? FitTemplates(hSub[s][p][k], &ts, &tb) : FitGaus(hSub[s][p][k]);
          if(!f.ok || f.en <= 0)continue;
          const double R = f.n*(1 - fpu)/eff/nBeamSub[s][k], eR = f.en*(1 - fpu)/eff/nBeamSub[s][k], Rt = nTrueBSub[s][k]/nBeamSub[s][k];
          const double pull = (R - Rt)/eR; sum += pull; sum2 += pull*pull; n++;
        }
        if(n > 1)printf("     pulls (%s, %d submuestras): media %+.2f, anchura %.2f\n", fn[a], n, sum/n, std::sqrt(std::max(0., sum2/n - (sum/n)*(sum/n))));
      }
      // figura: pseudo-datos y ajuste de plantillas (señal abajo)
      cv.cd(1 + s + 3*p);
      TH1D* hs = (TH1D*)ts.Clone(Form("fs_%d_%d", s, p)); hs->Scale(fT.n); hs->SetFillColor(kRed); hs->SetLineColor(kRed);
      TH1D* hb = (TH1D*)tb.Clone(Form("fb_%d_%d", s, p)); hb->Scale(hD1[s][p]->Integral() > 0 ? (hD1[s][p]->Integral() - fT.n) : 0); hb->Add(hs); hb->SetFillColor(kAzure + 1); hb->SetLineColor(kAzure + 1);
      hD1[s][p]->SetTitle(Form("haz %s, %s;p de rango (MeV/c);candidatos / 5 MeV/c", kSp[s], kPop[p]));
      hD1[s][p]->SetMarkerStyle(20); hD1[s][p]->SetMarkerSize(0.6);
      hD1[s][p]->SetMaximum(1.3*hD1[s][p]->GetMaximum()); hD1[s][p]->Draw("E");
      hb->Draw("hist same"); hs->Draw("hist same"); hD1[s][p]->Draw("E same");
    }
    for(int a = 0; a < 2; a++)if(comb[a][0] > 0){
      const double R = comb[a][1]/comb[a][0], eR = 1/std::sqrt(comb[a][0]);
      printf("   combinado (%s): R = %.5f ± %.5f (%.1f %%), cierre %+.1f %% = %+.2f sigma\n", a == 0 ? "plantillas" : "gauss+exp", R, eR, 100*eR/R, 100*(R/rTrue - 1), (R - rTrue)/eR);
    }
    printf("\n");
  }
  cv.SaveAs(Form("%s/mockB_fits.pdf", figdir)); cv.SaveAs(Form("%s/mockB_fits.png", figdir));
  printf("Nota: errores estadísticos de un cuarto del MC (~0,2 x los datos a 6 GeV). Sin incertidumbre de las plantillas, de la eficiencia ni del PID del haz.\n");
  return 0;
}
