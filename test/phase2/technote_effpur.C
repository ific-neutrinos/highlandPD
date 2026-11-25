// Eficiencia y pureza tras cada corte de secondaryKaonProductionXS, por población (dentro y fuera de la beam slice) y
// por definición de la señal (A, B, C: capas acumuladas, A en B en C). Figura del capítulo de selección de la nota técnica.
//
//   root -l -b -q 'test/phase2/technote_effpur.C("xs_v100s_fix.root", "<technote>/4_selection/figs")'
//
// Eficiencia de la señal X en la población p tras el corte k: signal muons true de X (árbol truth) con un candidato
// emparejado que llega a cand_level >= k en p (accum_level[p] >= 1 + k), sobre todos los signal muons de X. Pureza:
// candidatos de p con cand_level >= k que son signal muons de X, sobre todos los candidatos de p con cand_level >= k
// (árbol ana). k = 1 es el corte de población. Imprime la tabla (con eficiencia x pureza) y guarda sel_effpur.pdf.
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include <cstdio>
#include <string>
#include <vector>

int technote_effpur(const char* in, const char* outdir){
  TFile* f = TFile::Open(in);
  if(!f || f->IsZombie())return -1;
  TTree* ana = (TTree*)f->Get("ana");
  TTree* tru = (TTree*)f->Get("truth");
  const std::vector<std::string> names[2] = {
    {"population", "faces", "vertex", "#chi^{2}_{#mu}", "dE/dx stop", "gaps", "Michel", "p range"},
    {"population", "not clear cosmic", "axis", "#chi^{2}_{#mu}", "dE/dx stop", "#chi^{2}_{p}", "dE/dx start", "gaps", "Michel", "p range"}};
  // capas acumuladas: A (KOrigin 0-1), B (0-3), C (cualquier signal muon)
  const char* layer[3] = {"A", "B", "C"};
  const char* korT[3] = {"(truemu_korigin==0||truemu_korigin==1)", "(truemu_korigin>=0&&truemu_korigin<=3)", "1"};
  const char* korC[3] = {"(cand_truekorigin==0||cand_truekorigin==1)", "(cand_truekorigin>=0&&cand_truekorigin<=3)", "1"};
  // paleta categórica de referencia (azul, naranja, aguamarina; validada para modo claro)
  const int col[3] = {TColor::GetColor("#2a78d6"), TColor::GetColor("#eb6834"), TColor::GetColor("#1baf7a")};
  const int mk[3] = {20, 21, 22};

  gStyle->SetOptStat(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
  TCanvas* cv = new TCanvas("cv_effpur", "", 1200, 900);
  cv->Divide(2, 2, 0.005, 0.005);
  std::vector<TObject*> keep;
  for(int p = 0; p < 2; p++){
    const int n = names[p].size();
    std::vector<double> eff[3], pur[3];
    printf("\n== %s ==\n  %-16s", p == 0 ? "Dentro de la beam slice" : "Fuera de la beam slice", "corte");
    for(int l = 0; l < 3; l++)printf("   eff %s   pur %s   e x p %s", layer[l], layer[l], layer[l]); printf("\n");
    for(int k = 1; k <= n; k++){
      const double nCand = ana->Draw("cand_level", Form("cand_pop==%d && cand_level>=%d", p, k), "goff");
      printf("  %-16s", names[p][k-1].c_str());
      for(int l = 0; l < 3; l++){
        const double nTot = tru->Draw("truemu_korigin", korT[l], "goff");
        const double nSel = tru->Draw("truemu_korigin", Form("%s && accum_level[%d]>=%d", korT[l], p, 1 + k), "goff");
        const double nSig = ana->Draw("cand_level", Form("cand_pop==%d && cand_level>=%d && cand_truesignalmuon==1 && %s", p, k, korC[l]), "goff");
        const double e = nTot > 0 ? nSel/nTot : 0, u = nCand > 0 ? nSig/nCand : 0;
        eff[l].push_back(e); pur[l].push_back(u);
        printf("   %6.4f   %6.4f   %7.5f", e, u, e*u);
      }
      printf("\n");
    }
    for(int q = 0; q < 2; q++){
      cv->cd(1 + 2*p + q);
      gPad->SetLogy(); gPad->SetBottomMargin(0.28); gPad->SetTopMargin(0.09); gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.03); gPad->SetGridy();
      TH1F* fr = new TH1F(Form("fr_%d_%d", p, q), "", n, 0.5, n + 0.5); keep.push_back(fr);
      for(int k = 1; k <= n; k++)fr->GetXaxis()->SetBinLabel(k, names[p][k-1].c_str());
      fr->GetXaxis()->SetLabelSize(0.055); fr->GetXaxis()->LabelsOption("v");
      fr->GetYaxis()->SetTitle(q == 0 ? "Efficiency" : "Purity"); fr->GetYaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetLabelSize(0.045);
      fr->SetMinimum(q == 0 ? 0.005 : 1e-4); fr->SetMaximum(1.2);
      fr->Draw("AXIS");
      TLegend* lg = q == 0 ? new TLegend(0.16, 0.31, 0.40, 0.50) : new TLegend(0.16, 0.66, 0.40, 0.85); lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.045); keep.push_back(lg);
      for(int l = 2; l >= 0; l--){
        TGraph* g = new TGraph(n); keep.push_back(g);
        for(int k = 0; k < n; k++)g->SetPoint(k, k + 1, q == 0 ? eff[l][k] : pur[l][k]);
        g->SetLineColor(col[l]); g->SetMarkerColor(col[l]); g->SetMarkerStyle(mk[l]); g->SetMarkerSize(1.1); g->SetLineWidth(2);
        g->Draw("LP same");
        lg->AddEntry(g, Form("signal %s", layer[l]), "lp");
      }
      lg->Draw();
      TLatex tt; tt.SetNDC(); tt.SetTextSize(0.045);
      tt.DrawLatex(0.13, 0.94, Form("#bf{DUNE:ProtoDUNE-SP} simulation, %s the beam slice", p == 0 ? "in" : "out of"));
      gPad->RedrawAxis();
    }
  }
  gSystem->mkdir(outdir, true);
  const TString pdf = Form("%s/sel_effpur.pdf", outdir);
  cv->Print(pdf);
  gSystem->Exec(Form("pdfcrop --margins 2 %s %s > /dev/null 2>&1", pdf.Data(), pdf.Data()));
  cv->Print(Form("%s/sel_effpur.png", outdir));
  return 0;
}
