// Medida B (26-nov-2026): ¿dónde aparece el hueco del fondo bajo el pico? Fondo true (no mu de K) de los candidatos en
// la beam slice, haz e/mu/pi, por nivel de la cadena de cortes (cand_level >= k): observado en 200-245 MeV/c frente al
// interpolado con exp(pol3) ajustado fuera de esa región (150-400). Por especie true del fondo.
//
//   root -l -b -q 'test/phase2/measB_dip.C+("opt_flat2_v100to105d01.root", "xs_v100to105d01.root")'
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include <cmath>
#include <cstdio>
#include <vector>

int measB_dip(const char* flat, const char* xs){
  TFile fx(xs); TTree* a = (TTree*)fx.Get("ana"); Int_t sp; a->SetBranchStatus("*", 0); a->SetBranchStatus("beam_truespecies", 1); a->SetBranchAddress("beam_truespecies", &sp);
  std::vector<char> es(a->GetEntries()); for(Long64_t e = 0; e < a->GetEntries(); e++){ a->GetEntry(e); es[e] = sp; }
  TFile ff(flat); TTree* t = (TTree*)ff.Get("flat"); Float_t pop, pr, level; Int_t cat; Long64_t entry;
  t->SetBranchStatus("*", 0); for(const char* b : {"pop", "prange", "cat", "level", "entry"})t->SetBranchStatus(b, 1);
  t->SetBranchAddress("pop", &pop); t->SetBranchAddress("prange", &pr); t->SetBranchAddress("cat", &cat); t->SetBranchAddress("level", &level); t->SetBranchAddress("entry", &entry);
  const char* lvName[8] = {"", "población", "caras", "vértice", "chi2_mu", "dE/dx parada", "huecos APA", "Michel"};
  const char* cn[4] = {"fondo", "pi", "p", "mu"};
  TH1D* h[8][4];
  for(int l = 1; l < 8; l++)for(int k = 0; k < 4; k++)h[l][k] = new TH1D(Form("h%d_%d", l, k), "", 50, 150, 400);
  for(Long64_t i = 0; i < t->GetEntries(); i++){
    t->GetEntry(i);
    if((int)pop != 0 || cat <= 2 || es[entry] != 0 || pr < 150 || pr >= 400)continue;
    const int k = cat == 3 ? 1 : cat == 4 ? 2 : cat == 5 ? 3 : -1;
    for(int l = 1; l < 8; l++){ if(level < l)break; h[l][0]->Fill(pr); if(k > 0)h[l][k]->Fill(pr); }
  }
  printf("Haz e/mu/pi, en la beam slice, fondo true con calidad o sin ella (todos los eventos): hueco en 200-245 MeV/c\n");
  printf("  %-14s", "nivel >="); for(int k = 0; k < 4; k++)printf(" | %-26s", cn[k]); printf("\n");
  for(int l = 1; l < 8; l++){
    printf("  %d %-12s", l, lvName[l]);
    for(int k = 0; k < 4; k++){
      TH1D hf(*h[l][k]); for(int b = 11; b <= 19; b++){ hf.SetBinContent(b, 0); hf.SetBinError(b, 1e9); }
      TF1 f(Form("f%d%d", l, k), "exp([0]+[1]*(x-150)/100+[2]*((x-150)/100)^2+[3]*((x-150)/100)^3)", 150, 400);
      f.SetParameters(std::log(h[l][k]->GetBinContent(1) + 1), -1, 0, 0);
      hf.Fit(&f, "QN0R");
      double obs = 0, pred = 0; for(int b = 11; b <= 19; b++){ obs += h[l][k]->GetBinContent(b); pred += f.Eval(h[l][k]->GetBinCenter(b)); }
      printf(" | %7.0f %+5.1f %% %+5.1f sig", obs, 100*(obs/pred - 1), (obs - pred)/std::sqrt(pred));
    }
    printf("\n");
  }
  return 0;
}
