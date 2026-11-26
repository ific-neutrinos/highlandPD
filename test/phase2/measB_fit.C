// Medida B, paso 1 (26-nov-2026): modelo del ajuste del pico de p de rango. Ajusta el espectro de la figura de la nota
// (cadena de cortes sin la ventana, en la beam slice, calidad del haz provisional, una especie del haz) con tres modelos
// y compara el número de mu de K ajustado con el true. Verosimilitud de Poisson, bins de 5 MeV/c, 150-400 MeV/c.
//   (a) gaussiana + exp(pol2);
//   (b) forma de la señal del MC (todos los mu de K: capas A, B no A, C no B) + fondo exp(pol2);
//   (c) formas de señal y fondo del MC, normalizaciones libres.
// Aquí el "dato" es el mismo MC del que salen las formas (paso 1: ¿el modelo describe el espectro?); la medida simulada
// con muestras independientes es el paso siguiente.
//
//   root -l -b -q 'test/phase2/measB_fit.C+("opt_flat2_v100to105d01.root", "xs_v100to105d01.root", 0, "figs/measB_fit_pi")'
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include <cmath>
#include <cstdio>
#include <functional>
#include <vector>

int measB_fit(const char* flat, const char* xs, int species, const char* out){
  gStyle->SetOptStat(0);
  const double x0 = 150, x1 = 400, bw = 5; const int nb = 50;
  TFile fx(xs); TTree* a = (TTree*)fx.Get("ana"); Int_t sp, hs; Float_t pos[4], dir[3], bdir[3];
  a->SetBranchStatus("*", 0);
  for(const char* b : {"beam_truespecies", "hasslice", "seltrk_pos", "seltrk_dir", "beam_enddir"})a->SetBranchStatus(b, 1);
  a->SetBranchAddress("beam_truespecies", &sp); a->SetBranchAddress("hasslice", &hs); a->SetBranchAddress("seltrk_pos", pos);
  a->SetBranchAddress("seltrk_dir", dir); a->SetBranchAddress("beam_enddir", bdir);
  std::vector<char> ok(a->GetEntries());
  double nBeam = 0;
  for(Long64_t e = 0; e < a->GetEntries(); e++){
    a->GetEntry(e); double c = 0; for(int i = 0; i < 3; i++)c += dir[i]*bdir[i];
    const bool q = hs && std::abs(pos[0] + 28.7) < 3*4.17 && std::abs(pos[1] - 422.4) < 3*3.91 && std::abs(pos[2] - 29.7) < 3*0.59 && c > 0.95;
    ok[e] = (sp == species && q); if(ok[e])nBeam++;
  }
  TFile ff(flat); TTree* t = (TTree*)ff.Get("flat"); Float_t pop, pr; Int_t cat, chain; Long64_t entry;
  t->SetBranchStatus("*", 0); for(const char* b : {"pop", "prange", "cat", "chain", "entry"})t->SetBranchStatus(b, 1);
  t->SetBranchAddress("pop", &pop); t->SetBranchAddress("prange", &pr); t->SetBranchAddress("cat", &cat); t->SetBranchAddress("chain", &chain); t->SetBranchAddress("entry", &entry);
  TH1D hD("hD", "", nb, x0, x1), hS("hS", "", nb, x0, x1), hB("hB", "", nb, x0, x1);
  double nSigB = 0;   // capas A + B no A (señal de la medida)
  for(Long64_t i = 0; i < t->GetEntries(); i++){
    t->GetEntry(i);
    if((int)pop != 0 || !chain || pr < x0 || pr >= x1 || !ok[entry])continue;
    hD.Fill(pr); (cat <= 2 ? hS : hB).Fill(pr); if(cat <= 1)nSigB++;
  }
  const double nK = hS.Integral();   // todos los mu de K del espectro (lo que mide el pico)
  TH1D ts(hS), tb(hB); ts.Scale(1./ts.Integral()); tb.Scale(1./tb.Integral());
  const double tot = hD.Integral();
  // (a)
  TF1 fa("fa", Form("[0]*%g*TMath::Gaus(x,[1],[2],true) + exp([3] + [4]*(x-%g)/100 + [5]*((x-%g)/100)^2)", bw, x0, x0), x0, x1);
  fa.SetParameters(0.1*tot, 228, 8, std::log(hD.GetBinContent(1) + 1), -1, 0); fa.SetParLimits(0, 0, tot); fa.SetParLimits(1, 215, 240); fa.SetParLimits(2, 3, 20);
  const int sa = hD.Fit(&fa, "LQNR0");
  // (b) forma de señal del MC + exp(pol2)
  std::function<double(double*, double*)> fbf = [&ts, x0, bw, nb](double* x, double* p) -> double {
    const int b = std::min(nb, std::max(1, int((x[0] - x0)/bw) + 1));
    return p[0]*ts.GetBinContent(b) + std::exp(p[1] + p[2]*(x[0] - x0)/100 + p[3]*std::pow((x[0] - x0)/100, 2)); };
  TF1 fb("fb", fbf, x0, x1, 4);
  fb.SetParameters(0.1*tot, std::log(hD.GetBinContent(1) + 1), -1, 0); fb.SetParLimits(0, 0, tot);
  const int sb = hD.Fit(&fb, "LQNR0");
  // (c) formas del MC
  std::function<double(double*, double*)> fcf = [&ts, &tb, x0, bw, nb](double* x, double* p) -> double {
    const int b = std::min(nb, std::max(1, int((x[0] - x0)/bw) + 1));
    return p[0]*ts.GetBinContent(b) + p[1]*tb.GetBinContent(b); };
  TF1 fc("fc", fcf, x0, x1, 2);
  fc.SetParameters(0.1*tot, 0.9*tot); fc.SetParLimits(0, 0, tot); fc.SetParLimits(1, 0, 2*tot);
  const int sc = hD.Fit(&fc, "LQNR0");
  const char* spn[3] = {"e/mu/pi", "K", "p"};
  printf("Haz %s, en la beam slice, calidad provisional: %.0f partículas, %.0f candidatos; mu de K true en el espectro %.0f (señal A+B %.0f)\n",
         spn[species], nBeam, tot, nK, nSigB);
  struct R { const char* n; double v, e; int st; double chi2; int ndf; };
  const R res[3] = {{"(a) gaussiana + exp(pol2)", fa.GetParameter(0), fa.GetParError(0), sa, fa.GetChisquare(), fa.GetNDF()},
                    {"(b) señal del MC + exp(pol2)", fb.GetParameter(0), fb.GetParError(0), sb, fb.GetChisquare(), fb.GetNDF()},
                    {"(c) señal y fondo del MC", fc.GetParameter(0), fc.GetParError(0), sc, fc.GetChisquare(), fc.GetNDF()}};
  const char* resEn[3] = {"(a) Gaussian + exp(pol2)", "(b) MC signal shape + exp(pol2)", "(c) MC signal and background shapes"};
  for(const R& r : res)printf("  %-30s N_K = %6.1f +- %5.1f (%4.1f %%), frente al true %+5.1f %% (%+5.2f sigma); chi2/ndf %.1f/%d%s\n", r.n, r.v, r.e, 100*r.e/r.v,
                              100*(r.v/nK - 1), (r.v - nK)/r.e, r.chi2, r.ndf, r.st ? "  FALLO" : "");
  if(species == 0)printf("  (a) media %.1f, sigma %.1f MeV/c\n", fa.GetParameter(1), fa.GetParameter(2));
  // figura: los tres ajustes sobre el espectro, con la señal true abajo
  TCanvas cv("cv", "", 1500, 500); cv.Divide(3, 1);
  TF1* fs[3] = {&fa, &fb, &fc};
  for(int m = 0; m < 3; m++){
    cv.cd(m + 1); gPad->SetLeftMargin(0.13);
    TH1D* sig = (TH1D*)hS.Clone(Form("sig%d", m)); sig->SetFillColor(kRed - 7); sig->SetLineColor(kRed - 7);
    TH1D* d = (TH1D*)hD.Clone(Form("d%d", m));
    d->SetTitle(Form("%s beam, %s;Momentum from range, #mu hypothesis (MeV/c);Candidates / 5 MeV/c", species == 0 ? "e/#mu/#pi" : spn[species], resEn[m]));
    d->SetMarkerStyle(20); d->SetMarkerSize(0.6); d->SetMinimum(0); d->SetMaximum(1.25*hD.GetMaximum());
    d->Draw("E"); sig->Draw("hist same"); d->Draw("E same");
    fs[m]->SetLineColor(kBlue); fs[m]->SetNpx(500); fs[m]->Draw("same");
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.04);
    tx.DrawLatex(0.45, 0.82, Form("N_{K} fit %.0f #pm %.0f", res[m].v, res[m].e));
    tx.DrawLatex(0.45, 0.77, Form("N_{K} true %.0f", nK));
    tx.DrawLatex(0.45, 0.72, Form("#chi^{2}/ndf %.0f/%d", res[m].chi2, res[m].ndf));
  }
  cv.SaveAs(Form("%s.pdf", out)); cv.SaveAs(Form("%s.png", out));
  return 0;
}
