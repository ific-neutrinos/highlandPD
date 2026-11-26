// Espectro de p de rango (hipótesis mu) tras la cadena de cortes y tras un BDT de opt_bdt*.C, apilado por categoría true
// y proyectado a los datos a 6 GeV (26-nov-2026). Solo la mitad de evaluación (half == 1), escalada x 2 x 1,196.
//
//   root -l -b -q 'test/phase2/opt_bdt_plot.C+("opt_flat2_v100to105d01.root", "opt_bdt_red/dataset/weights/bdt0_BDT.weights.xml", 0,
//                  "as_mean5,k_dedx,k_n,cnn_michel,michel_own,nmich10,cosbeam,dedx_stop", 0.20, "bdt_prange_in.pdf")'
//
// Izquierda: cadena de cortes sin la ventana (fuera, con el eje a 80 cm). Derecha: BDT > umbral con la preselección
// suelta de opt_bdt.C. Líneas: la ventana [205, 240] MeV/c.
#include "TFile.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TMVA/Reader.h"
#include <memory>
#include <string>
#include <vector>

int opt_bdt_plot(const char* in, const char* weights, int pop, const char* varList, double cut, const char* out){
  gStyle->SetOptStat(0);
  TFile* f = TFile::Open(in); TTree* t = (TTree*)f->Get("flat");
  t->SetAlias("gapmin", "min(min(abs(zstart-230),abs(zstart-462)),min(abs(zstop-230),abs(zstop-462)))");
  std::vector<std::string> vars;
  { std::string s = std::string(varList) + ","; for(size_t a = 0, b; (b = s.find(',', a)) != std::string::npos; a = b + 1)if(b > a)vars.push_back(s.substr(a, b - a)); }
  TMVA::Reader rd("!Color:Silent");
  std::vector<float> buf(vars.size());
  for(size_t i = 0; i < vars.size(); i++)rd.AddVariable(vars[i].c_str(), &buf[i]);
  rd.BookMVA("BDT", weights);
  std::vector<std::unique_ptr<TTreeFormula>> fv;
  for(auto& x : vars)fv.emplace_back(new TTreeFormula(x.c_str(), x.c_str(), t));
  const char* pre = pop == 0 ? "pop==0 && chi2mu<20 && faces>10" : "pop==1 && axis<150 && chi2mu<20";
  TTreeFormula fsel("sel", Form("pop==%d && half==1 && prange>=150 && prange<=400", pop), t), fpre("pre", pre, t),
               fch("ch", pop == 0 ? "chain" : "chain80", t), fcat("cat", "cat", t), fpr("pr", "prange", t);
  // categorías: A, B no A, C no B, pi, p, mu, e/gamma, otro (como technote_selection)
  const char* name[8] = {"signal A", "signal B (not A)", "signal C (not B)", "#pi^{#pm}", "p", "#mu (not signal)", "e, #gamma", "other"};
  const int col[8] = {kRed, kRed - 7, kRed - 10, kBlue, kGreen + 1, kCyan, kYellow, kGray + 1};
  TH1D* h[2][8];
  for(int s = 0; s < 2; s++)for(int k = 0; k < 8; k++){
    h[s][k] = new TH1D(Form("h%d_%d", s, k), "", 50, 150, 400);
    h[s][k]->SetFillColor(col[k]); h[s][k]->SetLineColor(col[k]);
  }
  for(Long64_t e = 0; e < t->GetEntries(); e++){
    t->LoadTree(e);
    if(fsel.EvalInstance() == 0)continue;
    const int k = (int)fcat.EvalInstance(); const double pr = fpr.EvalInstance();
    if(fch.EvalInstance() != 0)h[0][k]->Fill(pr);
    if(fpre.EvalInstance() != 0){
      for(size_t i = 0; i < vars.size(); i++)buf[i] = fv[i]->EvalInstance();
      if(rd.EvaluateMVA("BDT") > cut)h[1][k]->Fill(pr);
    }
  }
  const double scale = 2*1.196;
  TCanvas cv("cv", "", 1400, 600); cv.Divide(2, 1);
  const char* title[2] = {"cut chain (no window)", Form("BDT > %.2f (%zu variables)", cut, vars.size())};
  std::vector<THStack*> st;
  double ymax = 0;
  for(int s = 0; s < 2; s++){ double tot = 0; for(int k = 0; k < 8; k++){ h[s][k]->Scale(scale); } TH1D sum(*h[s][0]); for(int k = 1; k < 8; k++)sum.Add(h[s][k]); ymax = std::max(ymax, sum.GetMaximum()); (void)tot; }
  for(int s = 0; s < 2; s++){
    cv.cd(s + 1); gPad->SetLeftMargin(0.13);
    THStack* hs = new THStack(Form("hs%d", s), Form("%s, %s;Momentum from range, #mu hypothesis (MeV/c);Candidates / 5 MeV/c (projected to data)",
                                                   pop == 0 ? "in the beam slice" : "out of the beam slice", title[s]));
    for(int k = 0; k < 8; k++)hs->Add(h[s][k]);   // señal abajo (se ve su forma), fondo encima
    hs->SetMaximum(1.15*ymax);
    hs->Draw("hist");
    st.push_back(hs);
    for(double x : {205., 240.}){ TLine* l = new TLine(x, 0, x, 1.15*ymax); l->SetLineStyle(2); l->SetLineWidth(2); l->Draw(); }
    double sB = 0, nTot = 0, sBw = 0, nw = 0;
    for(int k = 0; k < 8; k++){ nTot += h[s][k]->Integral(); nw += h[s][k]->Integral(12, 18); }
    sB = h[s][0]->Integral() + h[s][1]->Integral(); sBw = h[s][0]->Integral(12, 18) + h[s][1]->Integral(12, 18);
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.035);
    tx.DrawLatex(0.45, 0.62, Form("150-400: B %.0f of %.0f (%.0f%%)", sB, nTot, 100*sB/nTot));
    tx.DrawLatex(0.45, 0.57, Form("A %.0f", h[s][0]->Integral()));
    tx.DrawLatex(0.45, 0.52, Form("window: B %.0f of %.0f (%.0f%%)", sBw, nw, 100*sBw/nw));
    if(s == 0){
      TLegend* lg = new TLegend(0.62, 0.68, 0.89, 0.89); lg->SetBorderSize(0); lg->SetFillStyle(0);
      for(int k = 0; k < 8; k++)lg->AddEntry(h[s][k], name[k], "f");
      lg->Draw();
    }
  }
  cv.SaveAs(out);
  cv.SaveAs(TString(out).ReplaceAll(".pdf", ".png"));
  return 0;
}
