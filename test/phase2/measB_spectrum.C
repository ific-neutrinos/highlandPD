// Nota técnica, capítulo de la medida B (26-nov-2026): espectro de p de rango (hipótesis mu) de los candidatos que pasan
// la cadena de cortes de secondaryKaonProductionXS salvo la ventana, por especie del haz y población, apilado por
// categoría true con la señal abajo. Sobre el árbol plano de opt_flat.C (cadena: chain dentro, chain80 fuera) y el
// microtree de secondaryKaonProductionXS (especie true del haz por evento; en el MC no hay emulación de los Cherenkov).
//
//   root -l -b -q 'test/phase2/measB_spectrum.C+("opt_flat2_v100to105d01.root", "xs_v100to105d01.root", 0, "figs/measB_prange_in")'
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include <cstdio>
#include <vector>
#include <cmath>

// quality = 1: solo eventos que pasan la calidad del haz PROVISIONAL de measB_beamq.C (inicio del primario sin corregir por
// SCE a < 3 sigma de la media medida en el MC: x -28,7 +- 4,17, y 422,4 +- 3,91, z 29,7 +- 0,59 cm; cos con el BI > 0,95),
// también en el denominador. La definitiva será con las posiciones corregidas por SCE (pdAnaUtils::GetBeamQuality).
int measB_spectrum(const char* flat, const char* xs, int pop, const char* out, int quality = 0){
  gStyle->SetOptStat(0);
  TFile fx(xs); TTree* ana = (TTree*)fx.Get("ana");
  Int_t sp, hs_; Float_t pos[4], dir[3], bdir[3]; ana->SetBranchStatus("*", 0);
  for(const char* b : {"beam_truespecies", "hasslice", "seltrk_pos", "seltrk_dir", "beam_enddir"})ana->SetBranchStatus(b, 1);
  ana->SetBranchAddress("beam_truespecies", &sp); ana->SetBranchAddress("hasslice", &hs_); ana->SetBranchAddress("seltrk_pos", pos);
  ana->SetBranchAddress("seltrk_dir", dir); ana->SetBranchAddress("beam_enddir", bdir);
  std::vector<int> evSp(ana->GetEntries());
  double nBeam[3] = {0};
  for(Long64_t e = 0; e < ana->GetEntries(); e++){
    ana->GetEntry(e);
    double c = 0; for(int i = 0; i < 3; i++)c += dir[i]*bdir[i];
    const bool q = hs_ && std::abs(pos[0] + 28.7) < 3*4.17 && std::abs(pos[1] - 422.4) < 3*3.91 && std::abs(pos[2] - 29.7) < 3*0.59 && c > 0.95;
    evSp[e] = (quality && !q) ? -1 : sp; if(evSp[e] >= 0 && evSp[e] <= 2)nBeam[evSp[e]]++;
  }
  TFile ff(flat); TTree* t = (TTree*)ff.Get("flat");
  Float_t fpop, pr; Int_t cat, chain, chain80; Long64_t entry;
  t->SetBranchStatus("*", 0);
  for(const char* b : {"pop", "prange", "cat", "chain", "chain80", "entry"})t->SetBranchStatus(b, 1);
  t->SetBranchAddress("pop", &fpop); t->SetBranchAddress("prange", &pr); t->SetBranchAddress("cat", &cat);
  t->SetBranchAddress("chain", &chain); t->SetBranchAddress("chain80", &chain80); t->SetBranchAddress("entry", &entry);
  const char* spName[3] = {"e/#mu/#pi beam", "K beam", "p beam"};
  const char* name[8] = {"signal A", "signal B (not A)", "signal C (not B)", "#pi^{#pm}", "p", "#mu (not signal)", "e, #gamma", "other"};
  const int col[8] = {kRed, kRed - 7, kRed - 10, kBlue, kGreen + 1, kCyan, kYellow, kGray + 1};
  TH1D* h[3][8];
  for(int s = 0; s < 3; s++)for(int k = 0; k < 8; k++){
    h[s][k] = new TH1D(Form("h%d_%d", s, k), "", 50, 150, 400); h[s][k]->SetFillColor(col[k]); h[s][k]->SetLineColor(col[k]);
  }
  for(Long64_t i = 0; i < t->GetEntries(); i++){
    t->GetEntry(i);
    if((int)fpop != pop || !(pop == 0 ? chain : chain80) || pr < 150 || pr >= 400)continue;
    const int s = evSp[entry]; if(s < 0 || s > 2)continue;
    h[s][cat]->Fill(pr);
  }
  TCanvas cv("cv", "", 1500, 520); cv.Divide(3, 1);
  printf("Población %s, cadena de cortes sin la ventana, p de rango 150-400 MeV/c (v1_00-v1_05)%s\n", pop == 0 ? "dentro" : "fuera", quality ? ", con calidad del haz (provisional)" : "");
  printf("  %-8s %9s | %6s %6s %6s %6s %7s | en la ventana [205, 240]: B, total, pureza B\n", "haz", "partículas", "A", "B-A", "C-B", "fondo", "total");
  for(int s = 0; s < 3; s++){
    cv.cd(s + 1); gPad->SetLeftMargin(0.14);
    THStack* hs = new THStack(Form("hs%d", s), Form("%s, %s%s;Momentum from range, #mu hypothesis (MeV/c);Candidates / 5 MeV/c", spName[s],
                                                   pop == 0 ? "in the beam slice" : "out of the beam slice", quality ? ", beam quality" : ""));
    for(int k = 0; k < 8; k++)hs->Add(h[s][k]);   // señal abajo
    hs->Draw("hist");
    const double ymax = hs->GetMaximum();
    for(double x : {205., 240.}){ TLine* l = new TLine(x, 0, x, ymax); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw(); }
    double n[8], nw[8], tot = 0, totw = 0;
    for(int k = 0; k < 8; k++){ n[k] = h[s][k]->Integral(); nw[k] = h[s][k]->Integral(12, 18); tot += n[k]; totw += nw[k]; }
    const double bkg = tot - n[0] - n[1] - n[2];
    printf("  %-8s %9.0f | %6.0f %6.0f %6.0f %6.0f %7.0f | %4.0f %5.0f %5.1f %%\n", s == 0 ? "e/mu/pi" : s == 1 ? "K" : "p", nBeam[s], n[0], n[1], n[2], bkg, tot,
           nw[0] + nw[1], totw, totw > 0 ? 100*(nw[0] + nw[1])/totw : 0.);
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.04);
    tx.DrawLatex(0.50, 0.80, Form("%.0f beam particles", nBeam[s]));
    tx.DrawLatex(0.50, 0.75, Form("signal A+B: %.0f", n[0] + n[1]));
    if(s == 2){
      TLegend* lg = new TLegend(0.55, 0.35, 0.89, 0.72); lg->SetBorderSize(0); lg->SetFillStyle(0);
      for(int k = 0; k < 8; k++)lg->AddEntry(h[s][k], name[k], "f");
      lg->Draw();
    }
  }
  cv.SaveAs(Form("%s.pdf", out)); cv.SaveAs(Form("%s.png", out));
  return 0;
}
