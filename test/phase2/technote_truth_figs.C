// Figuras del capítulo 3 de la nota técnica (signal definition and truth studies), a partir del árbol kprod de
// kaon_acceptance_minitrees.C (23-nov-2026).
//
//   root -l -b -q 'test/phase2/technote_truth_figs.C+("kaon_acceptance_v100.root", "<technote>/3_truth/figs")'
//
// Muestra: K+ de la primera interacción inelástica de un pi+ del haz (level 0, A estricta) con el vértice en el
// volumen fiducial. Figuras (PDF):
//   kprod_T_fate.pdf      T_K de producción, apilado por destino de la K
//   kprod_theta_fate.pdf  theta_K, apilado por destino
//   acc_T.pdf             fracción que da mu+ nu en reposo con el mu contenido, frente a T_K
//   acc_theta.pdf         lo mismo frente a theta_K
//   acc_map.pdf           mapa en (T_K, theta_K), con el número de K producidas en cada celda
//   rmin.pdf              opción a: fracción mínima de momento conservada por el pi del haz (K de nivel >= 1)
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THStack.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include <iostream>

int technote_truth_figs(const char* in, const char* outdir){
  gStyle->SetOptStat(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
  gStyle->SetTitleFont(42, "XYZ"); gStyle->SetLabelFont(42, "XYZ"); gStyle->SetLegendFont(42);
  gSystem->mkdir(outdir, true);
  TFile f(in);
  TTree* t = (TTree*)f.Get("kprod");
  if(!t){ std::cerr << "No hay árbol kprod en " << in << std::endl; return 1; }
  Int_t beamPDG, level, fate, muContained, vtxInFV; Float_t kT, kTheta, rMin;
  t->SetBranchAddress("beamPDG", &beamPDG); t->SetBranchAddress("level", &level); t->SetBranchAddress("fate", &fate);
  t->SetBranchAddress("muContained", &muContained); t->SetBranchAddress("vtxInFV", &vtxInFV);
  t->SetBranchAddress("kT", &kT); t->SetBranchAddress("kTheta", &kTheta); t->SetBranchAddress("rMin", &rMin);

  // destinos (kaon_acceptance_minitrees_impl.C): 0 mu nu en reposo, 1 pi+ pi0 en reposo, 2 otro en reposo, 3 en vuelo,
  // 4 interacciona, 5 sale, 6 sin información
  const int NF = 6;
  const char* fateName[NF] = {"at rest #rightarrow #mu^{+}#nu", "at rest #rightarrow #pi^{+}#pi^{0}", "at rest #rightarrow other",
                              "decay in flight", "inelastic, no K^{+} left", "exits active volume"};
  const int fateCol[NF] = {kRed + 1, kOrange + 1, kYellow + 2, kAzure + 1, kGreen + 2, kGray + 1};
  const double Tbins[] = {0, 100, 200, 300, 400, 600, 800, 1200, 2000};
  const double thBins[] = {0, 15, 30, 45, 60, 90, 180};
  const int nT = 8, nTh = 6;
  TH1D* hT[NF]; TH1D* hTh[NF];
  for(int i = 0; i < NF; i++){
    hT[i] = new TH1D(Form("hT%d", i), ";T_{K} at production (MeV);K^{+} / 100 MeV", 20, 0, 2000);
    hTh[i] = new TH1D(Form("hTh%d", i), ";#theta_{K} with respect to the beam (deg);K^{+} / 5 deg", 36, 0, 180);
    hT[i]->SetFillColor(fateCol[i]); hT[i]->SetLineColor(fateCol[i]);
    hTh[i]->SetFillColor(fateCol[i]); hTh[i]->SetLineColor(fateCol[i]);
  }
  TEfficiency eT("eT", ";T_{K} at production (MeV);fraction giving a contained #mu^{+} from K^{+} #rightarrow #mu^{+}#nu at rest", nT, Tbins);
  TEfficiency eTh("eTh", ";#theta_{K} with respect to the beam (deg);fraction giving a contained #mu^{+} from K^{+} #rightarrow #mu^{+}#nu at rest", nTh, thBins);
  TH2D hProd("hProd", "", nT, Tbins, nTh, thBins), hGood("hGood", "", nT, Tbins, nTh, thBins);
  TH1D hR("hR", ";r_{min} (minimum momentum fraction kept by the beam #pi^{+});K^{+} / 0.02", 50, 0, 1);
  long nA = 0, nOptA = 0;
  for(Long64_t i = 0; i < t->GetEntries(); i++){
    t->GetEntry(i);
    if(beamPDG != 211 || !vtxInFV)continue;
    if(level >= 1){ hR.Fill(std::min(rMin, 0.9999f)); nOptA++; continue; }
    nA++;
    const bool good = fate == 0 && muContained;
    if(fate >= 0 && fate < NF){ hT[fate]->Fill(std::min(kT, 1999.f)); hTh[fate]->Fill(std::min(kTheta, 179.9f)); }
    eT.Fill(good, std::min(kT, 1999.f)); eTh.Fill(good, std::min(kTheta, 179.9f));
    hProd.Fill(std::min(kT, 1999.f), std::min(kTheta, 179.9f)); if(good)hGood.Fill(std::min(kT, 1999.f), std::min(kTheta, 179.9f));
  }
  std::cout << "K+ de la primera interacción (pi+, FV): " << nA << "; de la opción a (nivel >= 1): " << nOptA << std::endl;

  TCanvas c("c", "", 800, 600);
  c.SetLeftMargin(0.12); c.SetBottomMargin(0.12);
  auto stack = [&](TH1D** h, const char* name, const char* file){
    THStack s(name, h[0]->GetTitle());
    for(int i = NF - 1; i >= 0; i--)s.Add(h[i]);
    c.Clear(); s.Draw("hist"); s.GetXaxis()->SetTitle(h[0]->GetXaxis()->GetTitle()); s.GetYaxis()->SetTitle(h[0]->GetYaxis()->GetTitle());
    TLegend l(0.52, 0.55, 0.88, 0.88); l.SetBorderSize(0); l.SetFillStyle(0);
    for(int i = 0; i < NF; i++)l.AddEntry(h[i], fateName[i], "f");
    l.Draw();
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.035); tx.SetTextFont(42);
    tx.DrawLatex(0.14, 0.93, "MC 6 GeV, #pi^{+} beam: K^{+} from the first inelastic interaction (vertex in FV)");
    c.Print(Form("%s/%s", outdir, file));
  };
  stack(hT, "sT", "kprod_T_fate.pdf");
  stack(hTh, "sTh", "kprod_theta_fate.pdf");

  auto effplot = [&](TEfficiency& e, const char* file){
    c.Clear(); e.SetMarkerStyle(20); e.SetLineWidth(2); e.Draw("AP"); gPad->Update();
    e.GetPaintedGraph()->SetMinimum(0); e.GetPaintedGraph()->SetMaximum(0.8); gPad->Update();
    TLine bl; bl.SetLineStyle(2); bl.SetLineColor(kGray + 2);
    const double x0 = e.GetPaintedGraph()->GetXaxis()->GetXmin(), x1 = e.GetPaintedGraph()->GetXaxis()->GetXmax();
    bl.DrawLine(x0, 0.6356, x1, 0.6356);
    TLatex tx; tx.SetTextSize(0.03); tx.SetTextFont(42); tx.SetTextColor(kGray + 2);
    tx.DrawLatex(x0 + 0.03*(x1 - x0), 0.655, "BR(K^{+} #rightarrow #mu^{+}#nu) = 0.636");
    c.Print(Form("%s/%s", outdir, file));
  };
  effplot(eT, "acc_T.pdf");
  effplot(eTh, "acc_theta.pdf");

  // mapa: fracción en color y número de K producidas como texto
  TH2D hAcc(hGood); hAcc.SetName("hAcc"); hAcc.Divide(&hProd);
  hAcc.SetTitle(";T_{K} at production (MeV);#theta_{K} (deg)");
  c.Clear(); c.SetRightMargin(0.14); gStyle->SetPaintTextFormat(".2f");
  hAcc.SetMinimum(0); hAcc.SetMaximum(0.8); hAcc.Draw("COLZ TEXT");
  c.Print(Form("%s/acc_map.pdf", outdir));
  c.Clear(); hProd.SetTitle(";T_{K} at production (MeV);#theta_{K} (deg)"); gStyle->SetPaintTextFormat(".0f");
  hProd.Draw("COLZ TEXT"); c.Print(Form("%s/prod_map.pdf", outdir));
  c.SetRightMargin(0.05);

  c.Clear(); hR.SetLineWidth(2); hR.SetLineColor(kBlack); hR.Draw("hist");
  TLine l; l.SetLineStyle(2); l.SetLineColor(kRed + 1);
  for(double r : {0.9, 0.95})l.DrawLine(r, 0, r, hR.GetMaximum());
  c.Print(Form("%s/rmin.pdf", outdir));
  return 0;
}
