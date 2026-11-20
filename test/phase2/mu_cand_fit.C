// ¿Se puede sacar la señal ajustando el pico del p de rango? Candidatos a mu de K+ en reposo (árbol "cand"
// de mu_candidates_minitrees.C), por población, con todos los cortes de su cadena salvo la ventana.
//
//   root -l -b -q 'test/phase2/mu_cand_fit.C+("mu_cand_big19.root", "<directorio de salida>")'
//
// Modelo: gaussiana (señal) + exp(polinomio de grado 2) (fondo), ajuste de máxima verosimilitud de Poisson.
// Rango: en la beam slice 150-400 MeV/c; fuera 180-300 (el módulo solo guarda completos los candidatos en
// esa ventana). Se ajusta el MC tal cual (S + B escalado) y, con pseudoexperimentos de Poisson, se miden el
// sesgo y el error del número de señal con la estadística de esta muestra (x1) y con más (x6, x10).
// Guarda mu_cand_fit.pdf y un PNG por población (con broadSignal, señal ampliada: mucut::BroadSignal;
// mu_cand_fit_amplia.pdf).
#include "mu_cand_common.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TSystem.h"

namespace mufit {
  using namespace mucut;

  // n0 = número de señal (la gaussiana integra a n0 con el ancho de bin bw)
  TF1* Model(const char* name, double x0, double x1, double bw){
    TF1* f = new TF1(name, Form("[0]*%g*TMath::Gaus(x,[1],[2],true) + exp([3] + [4]*(x-%g)/100 + [5]*((x-%g)/100)^2)", bw, x0, x0), x0, x1);
    f->SetParNames("N_{señal}", "media", "sigma", "b0", "b1", "b2");
    return f;
  }

  struct FitRes { double n, en, mean, sigma; bool ok; };

  FitRes Fit(TH1D* h, double x0, double x1, double bw, bool quiet = true){
    TF1* f = Model(Form("f_%s", h->GetName()), x0, x1, bw);
    const double tot = h->Integral(h->FindBin(x0 + 1e-3), h->FindBin(x1 - 1e-3));
    f->SetParameters(0.1*tot, 228, 9, std::log(std::max(1., h->GetBinContent(h->FindBin(x0 + 1e-3)))), -1, 0);
    f->SetParLimits(0, 0, 10*tot + 10); f->SetParLimits(1, 210, 245); f->SetParLimits(2, 4, 20);
    const int st = h->Fit(f, quiet ? "LQNR0" : "LR0");
    FitRes r{f->GetParameter(0), f->GetParError(0), f->GetParameter(1), f->GetParameter(2), st == 0};
    if(quiet)delete f;
    return r;
  }
}

int mu_cand_fit(const char* in, const char* outdir, bool broadSignal = false, int nToys = 500, long nEvents = 73070, long nEventsNormal = 65977){
  using namespace mufit;
  gBroadSignal = broadSignal;
  if(broadSignal)std::cout << "Señal ampliada (mucut::BroadSignal)" << std::endl;
  const std::vector<C> cands = Load(in);
  if(cands.empty())return 1;
  const double scaleB = double(nEvents)/nEventsNormal;
  gSystem->mkdir(outdir, true);
  gStyle->SetOptStat(0);
  const std::string pdf = std::string(outdir) + (broadSignal ? "/mu_cand_fit_amplia.pdf" : "/mu_cand_fit.pdf");
  TCanvas cv("cvfit", "", 900, 650);
  cv.Print((pdf + "[").c_str());
  TRandom3 rnd(12345);

  struct Pop { const char* tag; const char* label; std::vector<Step> st; std::function<bool(const C&)> sel; double x0, x1; };
  const Pop pops[] = {{"bs", "en la beam slice", InSliceSteps(), [](const C& c){ return c.inSlice; }, 150, 400},
                      {"fuera", "fuera de la beam slice", OutSliceSteps(), [](const C& c){ return !c.inSlice; }, 180, 300}};
  const double bw = 5;
  for(const Pop& p : pops){
    int iw = -1; for(int i = 0; i < (int)p.st.size(); i++)if(p.st[i].name.rfind("p rango", 0) == 0)iw = i;
    const int nb = int((p.x1 - p.x0)/bw + 0.5);
    TH1D hS(Form("hS_%s", p.tag), "", nb, p.x0, p.x1), hB(Form("hB_%s", p.tag), "", nb, p.x0, p.x1);
    for(const C& c : cands){
      if(c.group < 0 || !p.sel(c) || !c.complete || !Pass(p.st, c, iw))continue;
      if(c.pRange < p.x0 || c.pRange >= p.x1)continue;
      if(c.group == 0)hS.Fill(c.pRange); else hB.Fill(c.pRange, scaleB);
    }
    TH1D hT(hS); hT.SetName(Form("hT_%s", p.tag)); hT.Add(&hB);
    // el "dato": S + B (B escalado, así que no es entero: para el ajuste de Poisson se redondea)
    TH1D hD(hT); hD.SetName(Form("hD_%s", p.tag));
    for(int b = 1; b <= nb; b++){ hD.SetBinContent(b, std::round(hT.GetBinContent(b))); hD.SetBinError(b, std::sqrt(hD.GetBinContent(b))); }
    const double Strue = hS.Integral();
    const double Swin = hS.Integral(hS.FindBin(205.001), hS.FindBin(239.999));

    // Ajuste del MC tal cual, con figura
    TF1* f = Model(Form("fMC_%s", p.tag), p.x0, p.x1, bw);
    const double tot = hD.Integral();
    f->SetParameters(0.1*tot, 228, 9, std::log(std::max(1., hD.GetBinContent(1))), -1, 0);
    f->SetParLimits(0, 0, 10*tot); f->SetParLimits(1, 210, 245); f->SetParLimits(2, 4, 20);
    hD.Fit(f, "LQR0");
    TF1* fb = new TF1(Form("fb_%s", p.tag), Form("exp([0] + [1]*(x-%g)/100 + [2]*((x-%g)/100)^2)", p.x0, p.x0), p.x0, p.x1);
    fb->SetParameters(f->GetParameter(3), f->GetParameter(4), f->GetParameter(5));
    std::cout << std::endl << "== " << p.label << " (rango " << p.x0 << "-" << p.x1 << " MeV/c, bins de " << bw << ") ==" << std::endl
              << Form("  señal true en el rango %.0f (en [205, 240]: %.0f); fondo %.0f", Strue, Swin, hB.Integral()) << std::endl
              << Form("  ajuste: N = %.0f ± %.0f, media %.1f, sigma %.1f MeV/c, chi2/ndf %.1f/%d", f->GetParameter(0), f->GetParError(0),
                      f->GetParameter(1), f->GetParameter(2), f->GetChisquare(), f->GetNDF()) << std::endl;
    cv.Clear();
    hD.SetTitle(Form("%s;p por rango (hip#acute{o}tesis #mu, MeV/c);candidatos / %g MeV/c", p.label, bw));
    hD.SetMarkerStyle(20); hD.SetLineColor(kBlack); hD.SetMinimum(0);
    hD.Draw("E");
    hB.SetLineColor(kAzure + 1); hB.SetLineWidth(2); hB.Draw("hist same");
    hS.SetLineColor(kRed + 1); hS.SetLineWidth(2); hS.Draw("hist same");
    f->SetLineColor(kBlack); f->SetLineWidth(2); f->Draw("same");
    fb->SetLineColor(kBlack); fb->SetLineStyle(2); fb->Draw("same");
    TLegend leg(0.55, 0.65, 0.89, 0.89); leg.SetBorderSize(0); leg.SetFillStyle(0);
    leg.AddEntry(&hD, "MC (S + B)", "pe"); leg.AddEntry(&hS, "se#tilde{n}al (truth)", "l"); leg.AddEntry(&hB, "fondo (truth)", "l");
    leg.AddEntry(f, "ajuste", "l"); leg.AddEntry(fb, "ajuste: fondo", "l");
    leg.Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.035);
    tx.DrawLatex(0.13, 0.85, Form("N_{ajuste} = %.0f #pm %.0f;  N_{true} = %.0f", f->GetParameter(0), f->GetParError(0), Strue));
    cv.Print(pdf.c_str()); cv.Print(Form("%s/fit_%s%s.png", outdir, p.tag, broadSignal ? "_amplia" : ""));

    // Pseudoexperimentos: Poisson sobre la expectativa (S + B), x factor
    for(double fac : {1., 6., 10.}){
      std::vector<double> n, pull; int fails = 0;
      for(int t = 0; t < nToys; t++){
        TH1D hToy(hT); hToy.SetName(Form("toy_%s_%g_%d", p.tag, fac, t));
        for(int b = 1; b <= nb; b++){ const double v = rnd.Poisson(fac*hT.GetBinContent(b)); hToy.SetBinContent(b, v); hToy.SetBinError(b, std::sqrt(v)); }
        const FitRes r = Fit(&hToy, p.x0, p.x1, bw);
        if(!r.ok || r.en <= 0){ fails++; continue; }
        n.push_back(r.n/fac); pull.push_back((r.n - fac*Strue)/r.en);
      }
      std::sort(n.begin(), n.end()); std::sort(pull.begin(), pull.end());
      auto q = [](const std::vector<double>& v, double x){ return v.empty() ? -999. : v[std::min<size_t>(v.size() - 1, size_t(x*v.size()))]; };
      const double med = q(n, 0.5), sig = (q(n, 0.84) - q(n, 0.16))/2;
      std::cout << Form("  toys x%-3g: N/factor mediana %.0f (true %.0f, sesgo %+.1f %%), error relativo %.1f %%; pull mediana %+.2f, anchura %.2f; fallos %d de %d",
                        fac, med, Strue, 100*(med/Strue - 1), 100*sig/med, q(pull, 0.5), (q(pull, 0.84) - q(pull, 0.16))/2, fails, nToys) << std::endl;
    }
  }
  cv.Print((pdf + "]").c_str());
  std::cout << std::endl << "Figuras en " << pdf << std::endl;
  return 0;
}
