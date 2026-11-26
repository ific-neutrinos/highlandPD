// Medida B (26-nov-2026): función para el fondo del ajuste del pico. Con la misma muestra que measB_fit.C (en la beam
// slice, calidad del haz provisional, cadena sin la ventana, 150-400 MeV/c): (1) ajusta solo el fondo true del MC con
// varias funciones (chi2 de Pearson con el modelo, verosimilitud de Poisson); (2) ajusta el espectro completo con la
// forma de la señal del MC más cada función, y compara N_K con el true. u = (p - 150)/100 salvo donde se indica.
//   F1 exp(pol2)                                   F2 exp(pol3)
//   F3 A (p - p0)^k exp(-(p - p0)/tau), p0 < 150   (sube como una potencia, cae como una exponencial)
//   F4 A (1 - exp(-(p - p0)/r)) exp(-(p - p0)/tau) (sube y satura, cae como una exponencial)
//
//   root -l -b -q 'test/phase2/measB_bkgshape.C+("opt_flat2_v100to105d01.root", "xs_v100to105d01.root", 0, "figs/measB_bkg_pi")'
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

namespace bks {
  const double x0 = 150, x1 = 400, bw = 5; const int nb = 50;
  // funciones de fondo: valor por bin (densidad por 5 MeV/c)
  double F(int m, double x, const double* p){
    const double u = (x - x0)/100;
    switch(m){
      case 0: return std::exp(p[0] + p[1]*u + p[2]*u*u);
      case 1: return std::exp(p[0] + p[1]*u + p[2]*u*u + p[3]*u*u*u);
      case 2: { const double d = x - p[3]; return d > 0 ? std::exp(p[0])*std::pow(d/100, p[1])*std::exp(-d/p[2]) : 0; }
      default: { const double d = x - p[3]; return d > 0 ? std::exp(p[0])*(1 - std::exp(-d/p[1]))*std::exp(-d/p[2]) : 0; }
    }
  }
  const int kNPar[4] = {3, 4, 4, 4};
  const char* kName[4] = {"exp(pol2)", "exp(pol3)", "(p-p0)^k exp(-(p-p0)/tau)", "(1-exp(-(p-p0)/r)) exp(-(p-p0)/tau)"};
  void Init(int m, TF1& f, double y0){
    const double l = std::log(std::max(1., y0));
    if(m == 0)f.SetParameters(l, -1, 0);
    else if(m == 1)f.SetParameters(l, -1, 0, 0);
    else if(m == 2){ f.SetParameters(l + 1, 0.5, 80, 100); f.SetParLimits(1, 0, 5); f.SetParLimits(2, 10, 1000); f.SetParLimits(3, 0, 149); }
    else{ f.SetParameters(l + 0.3, 20, 120, 120); f.SetParLimits(1, 1, 200); f.SetParLimits(2, 10, 1000); f.SetParLimits(3, 0, 149); }
  }
  double Pearson(const TH1D& h, TF1& f){ double c = 0; for(int b = 1; b <= h.GetNbinsX(); b++){ const double mu = f.Eval(h.GetBinCenter(b)); if(mu > 0)c += std::pow(h.GetBinContent(b) - mu, 2)/mu; } return c; }
}

int measB_bkgshape(const char* flat, const char* xs, int species, const char* out){
  using namespace bks;
  gStyle->SetOptStat(0);
  TFile fx(xs); TTree* a = (TTree*)fx.Get("ana"); Int_t sp, hs; Float_t pos[4], dir[3], bdir[3];
  a->SetBranchStatus("*", 0);
  for(const char* b : {"beam_truespecies", "hasslice", "seltrk_pos", "seltrk_dir", "beam_enddir"})a->SetBranchStatus(b, 1);
  a->SetBranchAddress("beam_truespecies", &sp); a->SetBranchAddress("hasslice", &hs); a->SetBranchAddress("seltrk_pos", pos);
  a->SetBranchAddress("seltrk_dir", dir); a->SetBranchAddress("beam_enddir", bdir);
  std::vector<char> ok(a->GetEntries());
  for(Long64_t e = 0; e < a->GetEntries(); e++){
    a->GetEntry(e); double c = 0; for(int i = 0; i < 3; i++)c += dir[i]*bdir[i];
    ok[e] = sp == species && hs && std::abs(pos[0] + 28.7) < 3*4.17 && std::abs(pos[1] - 422.4) < 3*3.91 && std::abs(pos[2] - 29.7) < 3*0.59 && c > 0.95;
  }
  TFile ff(flat); TTree* t = (TTree*)ff.Get("flat"); Float_t pop, pr; Int_t cat, chain; Long64_t entry;
  t->SetBranchStatus("*", 0); for(const char* b : {"pop", "prange", "cat", "chain", "entry"})t->SetBranchStatus(b, 1);
  t->SetBranchAddress("pop", &pop); t->SetBranchAddress("prange", &pr); t->SetBranchAddress("cat", &cat); t->SetBranchAddress("chain", &chain); t->SetBranchAddress("entry", &entry);
  TH1D hD("hD", "", nb, x0, x1), hS("hS", "", nb, x0, x1), hB("hB", "", nb, x0, x1);
  for(Long64_t i = 0; i < t->GetEntries(); i++){
    t->GetEntry(i);
    if((int)pop != 0 || !chain || pr < x0 || pr >= x1 || !ok[entry])continue;
    hD.Fill(pr); (cat <= 2 ? hS : hB).Fill(pr);
  }
  const double nK = hS.Integral();
  TH1D ts(hS); ts.Scale(1./ts.Integral());
  printf("Haz %s: fondo true %.0f, mu de K true %.0f\n", species == 0 ? "e/mu/pi" : species == 1 ? "K" : "p", hB.Integral(), nK);
  TCanvas cv("cv", "", 1400, 550); cv.Divide(2, 1);
  cv.cd(1); hB.SetTitle(Form("%s beam: MC background only;Momentum from range, #mu hypothesis (MeV/c);Candidates / 5 MeV/c", species == 0 ? "e/#mu/#pi" : "K"));
  hB.SetMarkerStyle(20); hB.SetMarkerSize(0.6); hB.SetMinimum(0); hB.Draw("E");
  TLegend lg(0.40, 0.62, 0.89, 0.89); lg.SetBorderSize(0); lg.SetFillStyle(0);
  const int col[4] = {kBlue, kGreen + 2, kRed, kMagenta + 1};
  std::vector<TF1*> keep;
  for(int m = 0; m < 4; m++){
    // (1) solo el fondo
    std::function<double(double*, double*)> fb = [m](double* x, double* p) -> double { return F(m, x[0], p); };
    TF1* f = new TF1(Form("fb%d", m), fb, x0, x1, kNPar[m]); keep.push_back(f);
    Init(m, *f, hB.GetBinContent(1));
    const int st = hB.Fit(f, "LQNR0");
    const double chi2 = Pearson(hB, *f); const int ndf = nb - kNPar[m];
    f->SetLineColor(col[m]); f->SetNpx(500); f->Draw("same");
    lg.AddEntry(f, Form("%s: #chi^{2}/ndf %.0f/%d", kName[m], chi2, ndf), "l");
    // (2) espectro completo: forma de la señal del MC + esta función
    std::function<double(double*, double*)> fs = [m, &ts](double* x, double* p) -> double {
      const int b = std::min(nb, std::max(1, int((x[0] - x0)/bw) + 1));
      return p[0]*ts.GetBinContent(b) + F(m, x[0], p + 1); };
    TF1 g(Form("fs%d", m), fs, x0, x1, kNPar[m] + 1);
    for(int k = 0; k < kNPar[m]; k++){ double lo, hi; f->GetParLimits(k, lo, hi); g.SetParameter(k + 1, f->GetParameter(k)); if(lo < hi)g.SetParLimits(k + 1, lo, hi); }
    g.SetParameter(0, 0.1*hD.Integral()); g.SetParLimits(0, 0, hD.Integral());
    const int st2 = hD.Fit(&g, "LQNR0");
    printf("  %-40s solo fondo: chi2/ndf %5.1f/%d%s | con la señal: N_K = %6.1f +- %5.1f, frente al true %+5.1f %% (%+5.2f sigma), chi2/ndf %.1f/%d%s\n",
           kName[m], chi2, ndf, st ? " FALLO" : "", g.GetParameter(0), g.GetParError(0), 100*(g.GetParameter(0)/nK - 1), (g.GetParameter(0) - nK)/g.GetParError(0),
           Pearson(hD, g), nb - kNPar[m] - 1, st2 ? " FALLO" : "");
  }
  lg.Draw();
  // panel 2: residuos relativos del fondo
  cv.cd(2);
  TH1D frame("frame", "Background: (MC - fit)/fit;Momentum from range, #mu hypothesis (MeV/c);Relative residual", nb, x0, x1);
  frame.SetMinimum(-0.4); frame.SetMaximum(0.4); frame.Draw("AXIS");
  for(int m = 0; m < 4; m++){
    TH1D* r = new TH1D(Form("r%d", m), "", nb, x0, x1);
    for(int b = 1; b <= nb; b++){ const double mu = keep[m]->Eval(hB.GetBinCenter(b)); if(mu > 0){ r->SetBinContent(b, (hB.GetBinContent(b) - mu)/mu); r->SetBinError(b, std::sqrt(std::max(1., hB.GetBinContent(b)))/mu); } }
    r->SetLineColor(col[m]); r->SetMarkerColor(col[m]); r->SetMarkerStyle(20 + m); r->SetMarkerSize(0.5); r->Draw("E same");
  }
  cv.SaveAs(Form("%s.pdf", out)); cv.SaveAs(Form("%s.png", out));
  return 0;
}
