// Medida B (26-nov-2026): ajuste con el segundo pico de las K+ en reposo. El fondo tiene un pico de pi+ de K+ -> pi+ pi0
// en reposo (p = 205 MeV/c; en la hipótesis de mu, 160-180 MeV/c) que las funciones suaves no describen (parecía un
// hueco bajo el pico del mu). Modelo: forma del mu de K del MC (N_K) + forma del pi de K del MC (N_pi, libre) + fondo
// suave (exp(pol2) o subida saturada x exp). Misma muestra que measB_fit.C (en la beam slice, calidad provisional,
// cadena sin ventana, 150-400 MeV/c); desde el microtree directamente.
//
// sce = 1 (27-nov-2026, por defecto): calidad del haz beam_quality_sce de los xs_*d01b; sce = 0: la provisional.
//
//   root -l -b -q 'test/phase2/measB_fit2.C+("xs_v100to105d01b.root", 0, "figs/measB_fit2_pi")'
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TStyle.h"
#include <cmath>
#include <cstdio>
#include <functional>

int measB_fit2(const char* xs, int species, const char* out, int sce = 1){
  gStyle->SetOptStat(0);
  const double x0 = 150, x1 = 400, bw = 5; const int nb = 50;
  TFile fx(xs); TTree* a = (TTree*)fx.Get("ana");
  const TString q = sce ? "hasslice && beam_quality_sce==1"
                        : "hasslice && abs(seltrk_pos[0]+28.7)<3*4.17 && abs(seltrk_pos[1]-422.4)<3*3.91 && abs(seltrk_pos[2]-29.7)<3*0.59"
                          " && seltrk_dir[0]*beam_enddir[0]+seltrk_dir[1]*beam_enddir[1]+seltrk_dir[2]*beam_enddir[2]>0.95";
  const TString base = Form("beam_truespecies==%d && %s && cand_pop==0 && cand_level>=7 && cand_prange>=150 && cand_prange<400", species, q.Data());
  TH1D hD("hD", "", nb, x0, x1), hK("hK", "", nb, x0, x1), hP("hP", "", nb, x0, x1), hB("hB", "", nb, x0, x1);
  a->Project("hD", "cand_prange", base);
  a->Project("hK", "cand_prange", base + " && cand_truesignalmuon==1");
  a->Project("hP", "cand_prange", base + " && cand_truesignalmuon!=1 && cand_truepdg==211 && cand_trueparentpdg==321 && cand_trueproc==2");
  a->Project("hB", "cand_prange", base + " && cand_truesignalmuon!=1 && !(cand_truepdg==211 && cand_trueparentpdg==321 && cand_trueproc==2)");
  const double nK = hK.Integral(), nP = hP.Integral();
  printf("Haz %s: candidatos %.0f; mu de K %.0f, pi de K->pi pi0 %.0f, resto del fondo %.0f\n", species == 0 ? "e/mu/pi" : "K", hD.Integral(), nK, nP, hB.Integral());
  TH1D tK(hK), tP(hP); tK.Scale(1./tK.Integral()); tP.Scale(1./std::max(1., tP.Integral()));
  auto bin = [=](double x){ return std::min(nb, std::max(1, int((x - x0)/bw) + 1)); };
  for(int m = 0; m < 2; m++){
    for(int withPi = 0; withPi < 2; withPi++){
      std::function<double(double*, double*)> fn = [&, m, withPi](double* x, double* p) -> double {
        const int b = bin(x[0]);
        double bk;
        if(m == 0){ const double u = (x[0] - x0)/100; bk = std::exp(p[2] + p[3]*u + p[4]*u*u); }
        else{ const double d = x[0] - p[5]; bk = d > 0 ? std::exp(p[2])*(1 - std::exp(-d/p[3]))*std::exp(-d/p[4]) : 0; }
        return p[0]*tK.GetBinContent(b) + (withPi ? p[1]*tP.GetBinContent(b) : 0) + bk; };
      TF1 f(Form("f%d%d", m, withPi), fn, x0, x1, m == 0 ? 5 : 6);
      const double tot = hD.Integral(), l = std::log(hD.GetBinContent(1) + 1);
      f.SetParameter(0, 0.1*tot); f.SetParLimits(0, 0, tot);
      f.SetParameter(1, withPi ? 0.02*tot : 0); if(withPi)f.SetParLimits(1, 0, tot); else f.FixParameter(1, 0);
      if(m == 0){ f.SetParameter(2, l); f.SetParameter(3, -1); f.SetParameter(4, 0); }
      else{ f.SetParameter(2, l + 0.3); f.SetParameter(3, 20); f.SetParameter(4, 120); f.SetParameter(5, 120);
            f.SetParLimits(3, 1, 200); f.SetParLimits(4, 10, 1000); f.SetParLimits(5, 0, 149); }
      const int st = hD.Fit(&f, "LQNR0");
      double chi2 = 0; for(int b = 1; b <= nb; b++){ const double mu = f.Eval(hD.GetBinCenter(b)); if(mu > 0)chi2 += std::pow(hD.GetBinContent(b) - mu, 2)/mu; }
      printf("  fondo %-28s %-14s N_K = %6.1f +- %5.1f (true %.0f: %+5.1f %%, %+5.2f sigma)%s, chi2/ndf %.1f/%d%s\n", m == 0 ? "exp(pol2)" : "subida saturada x exp",
             withPi ? "+ pico K->pipi:" : "sin pico K->pipi:", f.GetParameter(0), f.GetParError(0), nK, 100*(f.GetParameter(0)/nK - 1), (f.GetParameter(0) - nK)/f.GetParError(0),
             withPi ? Form(", N_pi = %.0f +- %.0f (true %.0f)", f.GetParameter(1), f.GetParError(1), nP) : "", chi2, nb - (m == 0 ? 5 : 6) + (withPi ? 0 : 1), st ? " FALLO" : "");
    }
  }
  // figura: las tres componentes true apiladas (señal abajo) y el pico K->pipi
  TCanvas cv("cv", "", 800, 600);
  TH1D* d = (TH1D*)hD.Clone("d"); d->SetTitle(Form("%s beam, in the beam slice, beam quality;Momentum from range, #mu hypothesis (MeV/c);Candidates / 5 MeV/c", species == 0 ? "e/#mu/#pi" : "K"));
  d->SetMarkerStyle(20); d->SetMarkerSize(0.6); d->SetMinimum(0); d->Draw("E");
  TH1D* sK = (TH1D*)hK.Clone("sK"); sK->SetFillColor(kRed - 7); sK->SetLineColor(kRed - 7);
  TH1D* sP = (TH1D*)hP.Clone("sP"); sP->Add(sK); sP->SetFillColor(kOrange + 1); sP->SetLineColor(kOrange + 1);
  sP->Draw("hist same"); sK->Draw("hist same"); d->Draw("E same");
  TLatex tx; tx.SetNDC(); tx.SetTextSize(0.035);
  tx.SetTextColor(kRed - 7); tx.DrawLatex(0.55, 0.82, "K^{+} #rightarrow #mu^{+}#nu at rest");
  tx.SetTextColor(kOrange + 1); tx.DrawLatex(0.55, 0.77, "K^{+} #rightarrow #pi^{+}#pi^{0} at rest (#pi^{+})");
  cv.SaveAs(Form("%s.pdf", out)); cv.SaveAs(Form("%s.png", out));
  return 0;
}
