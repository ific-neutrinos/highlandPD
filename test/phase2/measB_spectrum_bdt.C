// Nota técnica, medida B (26-nov-2026): como measB_spectrum.C (figura del espectro por especie del haz, en la beam
// slice, calidad del haz provisional) pero con la selección BDT (opt_bdt2.C, 10 variables) en lugar de la cadena de
// cortes. Solo la mitad de evaluación del BDT (half == 1), escalada x2. Categorías true como en la figura de cortes, y
// además el pi+ de K+ -> pi+ pi0 en reposo separado dentro del fondo (solo para verlo: no entra en la selección).
// Comprobación del fondo: observado en 200-245 MeV/c frente al interpolado con exp(pol3) ajustado fuera de esa región.
// Calidad del haz (27-nov-2026): sce = 1 (por defecto) usa beam_quality_sce de los xs nuevos (xs_*d01b, inicio del
// primario corregido por SCE y cortes de pdUtils.parameters.dat, como pdAnaUtils::GetBeamQuality); sce = 0, la provisional
// sin SCE (inicio a z ~ 30 cm, medias medidas en measB_beamq.C). Los xs_*d01b tienen los mismos eventos en el mismo orden
// que los xs_*d01 (comprobado), así que valen con opt_flat3 (entry).
//
//   root -l -b -q 'test/phase2/measB_spectrum_bdt.C+("opt_flat3_v100to105d01.root", "xs_v100to105d01b.root", "<pesos>", 0.19, "figs/measB_prange_in_bqsce_bdt")'
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TMVA/Reader.h"
#include <cmath>
#include <functional>
#include <algorithm>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

int measB_spectrum_bdt(const char* flat, const char* xs, const char* weights, double cut, const char* out,
                       const char* varList = "as_mean5,k_dedx,k_n,cnn_michel,michel_own,nmich10,cosbeam,dedx_stop,ar_charge10,intpeak",
                       int sce = 1){
  gStyle->SetOptStat(0);
  TFile fx(xs); TTree* ana = (TTree*)fx.Get("ana");
  Int_t sp, hs_, qsce = -1; Float_t pos[4], dir[3], bdir[3];
  ana->SetBranchStatus("*", 0);
  for(const char* b : {"beam_truespecies", "hasslice", "seltrk_pos", "seltrk_dir", "beam_enddir"})ana->SetBranchStatus(b, 1);
  if(sce){ ana->SetBranchStatus("beam_quality_sce", 1); ana->SetBranchAddress("beam_quality_sce", &qsce); }
  ana->SetBranchAddress("beam_truespecies", &sp); ana->SetBranchAddress("hasslice", &hs_); ana->SetBranchAddress("seltrk_pos", pos);
  ana->SetBranchAddress("seltrk_dir", dir); ana->SetBranchAddress("beam_enddir", bdir);
  std::vector<int> evSp(ana->GetEntries()); double nBeam[3] = {0};
  for(Long64_t e = 0; e < ana->GetEntries(); e++){
    ana->GetEntry(e); double c = 0; for(int i = 0; i < 3; i++)c += dir[i]*bdir[i];
    const bool q = sce ? (hs_ && qsce == 1)
                       : (hs_ && std::abs(pos[0] + 28.7) < 3*4.17 && std::abs(pos[1] - 422.4) < 3*3.91 && std::abs(pos[2] - 29.7) < 3*0.59 && c > 0.95);
    evSp[e] = q ? sp : -1; if(evSp[e] >= 0 && evSp[e] <= 2 && e % 2 == 1)nBeam[evSp[e]] += 2;   // mitad 1, x2
  }
  std::vector<std::string> vars; { std::string s = std::string(varList) + ","; for(size_t a = 0, b; (b = s.find(',', a)) != std::string::npos; a = b + 1)if(b > a)vars.push_back(s.substr(a, b - a)); }
  TFile ff(flat); TTree* t = (TTree*)ff.Get("flat");
  std::map<std::string, Float_t> val; for(auto& v : vars)val[v] = 0;
  for(const char* x : {"pop", "prange", "chi2mu", "faces", "pdg", "ppdg", "tproc"})val[x] = 0;
  Int_t cat; Long64_t entry;
  t->SetBranchStatus("*", 0);
  for(auto& kv : val){ t->SetBranchStatus(kv.first.c_str(), 1); t->SetBranchAddress(kv.first.c_str(), &kv.second); }
  t->SetBranchStatus("cat", 1); t->SetBranchAddress("cat", &cat); t->SetBranchStatus("entry", 1); t->SetBranchAddress("entry", &entry);
  TMVA::Reader rd("!Color:Silent"); std::vector<float> buf(vars.size());
  for(size_t i = 0; i < vars.size(); i++)rd.AddVariable(vars[i].c_str(), &buf[i]);
  rd.BookMVA("BDT", weights);
  // categorías: 0 A, 1 B no A, 2 C no B, 3 pi de K->pipi, 4 pi, 5 p, 6 mu, 7 e/gamma, 8 otro
  const char* name[9] = {"signal A", "signal B (not A)", "signal C (not B)", "#pi^{+} of K^{+}#rightarrow#pi^{+}#pi^{0}", "#pi^{#pm}", "p", "#mu (not signal)", "e, #gamma", "other"};
  const int col[9] = {kRed, kRed - 7, kRed - 10, kOrange + 1, kBlue, kGreen + 1, kCyan, kYellow, kGray + 1};
  TH1D* h[3][9];
  for(int s = 0; s < 3; s++)for(int k = 0; k < 9; k++){ h[s][k] = new TH1D(Form("h%d_%d", s, k), "", 50, 150, 400); h[s][k]->SetFillColor(col[k]); h[s][k]->SetLineColor(col[k]); }
  for(Long64_t i = 0; i < t->GetEntries(); i++){
    t->GetEntry(i);
    if((int)val["pop"] != 0 || entry % 2 != 1)continue;
    const double pr = val["prange"]; if(pr < 150 || pr >= 400)continue;
    if(val["chi2mu"] >= 20 || val["faces"] <= 10)continue;
    const int s = evSp[entry]; if(s < 0 || s > 2)continue;
    for(size_t k = 0; k < vars.size(); k++)buf[k] = val[vars[k]];
    if(rd.EvaluateMVA("BDT") <= cut)continue;
    int k;
    if(cat <= 2)k = cat;
    else if((int)val["pdg"] == 211 && (int)val["ppdg"] == 321 && (int)val["tproc"] == 2)k = 3;
    else k = cat + 1;
    h[s][k]->Fill(pr, 2.);
  }
  TCanvas cv("cv", "", 1500, 520); cv.Divide(3, 1);
  const char* spName[3] = {"e/#mu/#pi beam", "K beam", "p beam"};
  printf("BDT > %.2f (10 variables), en la beam slice, calidad del haz %s; mitad de evaluación x2\n", cut, sce ? "con SCE (beam_quality_sce)" : "provisional sin SCE");
  printf("  %-8s %9s | %6s %6s %6s %6s %6s %7s | ventana: B, total, pureza | fondo en 200-245: observado / interpolado\n", "haz", "partículas", "A", "B-A", "C-B", "Kpipi", "fondo", "total");
  for(int s = 0; s < 3; s++){
    cv.cd(s + 1); gPad->SetLeftMargin(0.14);
    THStack* hs = new THStack(Form("hs%d", s), Form("%s, in the beam slice, beam quality%s, BDT;Momentum from range, #mu hypothesis (MeV/c);Candidates / 5 MeV/c", spName[s], sce ? " (SCE)" : ""));
    for(int k = 0; k < 9; k++)hs->Add(h[s][k]);
    hs->Draw("hist");
    const double ymax = hs->GetMaximum();
    for(double x : {205., 240.}){ TLine* l = new TLine(x, 0, x, ymax); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw(); }
    double n[9], nw[9], tot = 0, totw = 0;
    for(int k = 0; k < 9; k++){ n[k] = h[s][k]->Integral(); nw[k] = h[s][k]->Integral(12, 18); tot += n[k]; totw += nw[k]; }
    // fondo (todo lo que no es mu de K, con el pi de K->pipi): ¿hueco bajo el pico?
    TH1D hb("hb", "", 50, 150, 400); for(int k = 3; k < 9; k++)hb.Add(h[s][k]);
    TH1D hf(hb); for(int b = 11; b <= 19; b++){ hf.SetBinContent(b, 0); hf.SetBinError(b, 1e9); }
    for(int b = 1; b <= 50; b++)if(b < 11 || b > 19)hf.SetBinError(b, std::sqrt(std::max(2., hf.GetBinContent(b))*2));   // pesos x2
    TF1 f(Form("f%d", s), "exp([0]+[1]*(x-150)/100+[2]*((x-150)/100)^2+[3]*((x-150)/100)^3)", 150, 400);
    f.SetParameters(std::log(hb.GetBinContent(1) + 1), -1, 0, 0); hf.Fit(&f, "QN0R");
    double obs = 0, pred = 0; for(int b = 11; b <= 19; b++){ obs += hb.GetBinContent(b); pred += f.Eval(hb.GetBinCenter(b)); }
    const double bkg = tot - n[0] - n[1] - n[2];
    printf("  %-8s %9.0f | %6.0f %6.0f %6.0f %6.0f %6.0f %7.0f | %4.0f %5.0f %5.1f %% | %5.0f / %6.1f (%+5.1f %%, %+4.1f sigma)\n", s == 0 ? "e/mu/pi" : s == 1 ? "K" : "p", nBeam[s],
           n[0], n[1], n[2], n[3], bkg, tot, nw[0] + nw[1], totw, totw > 0 ? 100*(nw[0] + nw[1])/totw : 0., obs, pred, 100*(obs/pred - 1), (obs - pred)/std::sqrt(2*pred));
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.04);
    tx.DrawLatex(0.50, 0.80, Form("%.0f beam particles", nBeam[s]));
    tx.DrawLatex(0.50, 0.75, Form("signal A+B: %.0f", n[0] + n[1]));
    if(s == 2){ TLegend* lg = new TLegend(0.50, 0.32, 0.89, 0.72); lg->SetBorderSize(0); lg->SetFillStyle(0); for(int k = 0; k < 9; k++)lg->AddEntry(h[s][k], name[k], "f"); lg->Draw(); }
  }
  cv.SaveAs(Form("%s.pdf", out)); cv.SaveAs(Form("%s.png", out));
  // ajuste (modelo de la nota: forma del mu de K + forma del pi de K->pipi, del MC, + exp(pol2)), sobre la mitad de
  // evaluación sin escalar (recuentos enteros); las formas salen de la misma muestra (paso 1, como measB_fit.C)
  for(int s = 0; s < 2; s++){
    TH1D hD("hD", "", 50, 150, 400), tK("tK", "", 50, 150, 400), tP("tP", "", 50, 150, 400);
    for(int k = 0; k < 9; k++)hD.Add(h[s][k], 0.5);
    for(int k = 0; k < 3; k++)tK.Add(h[s][k], 0.5);
    tP.Add(h[s][3], 0.5);
    const double nK = tK.Integral(), nP = tP.Integral();
    tK.Scale(1./nK); if(nP > 0)tP.Scale(1./nP);
    for(int withPi = 0; withPi < 2; withPi++){
      std::function<double(double*, double*)> fn = [&tK, &tP, withPi](double* x, double* p) -> double {
        const int b = std::min(50, std::max(1, int((x[0] - 150)/5) + 1)); const double u = (x[0] - 150)/100;
        return p[0]*tK.GetBinContent(b) + (withPi ? p[1]*tP.GetBinContent(b) : 0) + std::exp(p[2] + p[3]*u + p[4]*u*u); };
      TF1 f(Form("ff%d%d", s, withPi), fn, 150, 400, 5);
      const double tot = hD.Integral();
      f.SetParameters(0.2*tot, 0.02*tot, std::log(hD.GetBinContent(1) + 1), -1, 0); f.SetParLimits(0, 0, tot);
      if(withPi)f.SetParLimits(1, 0, tot); else f.FixParameter(1, 0);
      hD.Fit(&f, "LQNR0");
      printf("  ajuste %s, %s: N_K = %.1f +- %.1f (true %.0f: %+5.1f %%, %+5.2f sigma, error %.1f %%)%s\n", s == 0 ? "e/mu/pi" : "K", withPi ? "con pico K->pipi" : "sin pico K->pipi",
             f.GetParameter(0), f.GetParError(0), nK, 100*(f.GetParameter(0)/nK - 1), (f.GetParameter(0) - nK)/f.GetParError(0), 100*f.GetParError(0)/f.GetParameter(0),
             withPi ? Form(", N_pi %.0f +- %.0f (true %.0f)", f.GetParameter(1), f.GetParError(1), nP) : "");
    }
  }
  return 0;
}
