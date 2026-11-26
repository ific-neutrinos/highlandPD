// Rendimientos de la selección con el BDT de opt_bdt.C frente a la cadena de cortes, proyectados a los datos a 6 GeV
// (26-nov-2026). Sobre la mitad de evaluación (half == 1) del árbol plano de opt_flat.C, escalada x 2 x 1,196.
//
//   root -l -b -q 'test/phase2/opt_bdt_yield.C+("opt_flat2_v100to105d01.root", "opt_bdt_nomom", 0.15, 0.18, "trmean,trmean_body,chi2mu,chi2mu_trimboth,chi2p")'
//
// Por población y selección (cadena de cortes sin la ventana; BDT > umbral): candidatos de cada categoría en la muestra
// del ajuste (p de rango 150-400 MeV/c, antes de la ventana) y en la ventana [205, 240]. Las variables del BDT se
// reconstruyen como en opt_bdt.C (mismo orden; exclude igual que al entrenar).
#include "TFile.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "TMVA/Reader.h"
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

namespace opty {
  const std::vector<std::string> kCommon = {"chi2mu", "chi2mu_trimboth", "chi2p", "dedx_stop", "dedx_start", "nmich10", "nmich5", "cnn_michel",
    "michel_own", "nmichel_own", "isshower", "cnn_trk", "cnn_em", "trmean", "trmean_body", "intpeak", "diry_stop", "cosbeam", "ar_charge10",
    "ar_n10", "as_charge10", "ndau_trk", "ndau_shw", "dau_nearstop", "dau_maxmichel", "k_n", "k_dedx", "faces", "gapmin", "as_n5", "as_n10", "as_mean5", "as_mean10"};
  const std::vector<std::string> kIn = {"vtx", "isprimary", "vc_dist", "vc_cover", "dstart_v", "dstop_v"};
  const std::vector<std::string> kOut = {"axis", "clear", "cosmiccrit", "slice_min"};
  const double kScale = 2*1.196;
  const char* kCat[5] = {"A", "B no A", "C no B", "fondo", "total"};
}

int opt_bdt_yield(const char* in, const char* dir, double cutIn, double cutOut, const char* exclude){
  using namespace opty;
  TFile* f = TFile::Open(in); TTree* t = (TTree*)f->Get("flat");
  t->SetAlias("gapmin", "min(min(abs(zstart-230),abs(zstart-462)),min(abs(zstop-230),abs(zstop-462)))");
  const std::string ex = std::string(",") + exclude + ",";
  printf("Proyección a los datos a 6 GeV (mitad de evaluación x %.3f). Muestra del ajuste: p de rango 150-400 MeV/c; ventana [205, 240].\n", kScale);
  for(int p = 0; p < 2; p++){
    std::vector<std::string> vars;
    for(auto& x : kCommon)if(ex.find("," + x + ",") == std::string::npos)vars.push_back(x);
    for(auto& x : (p == 0 ? kIn : kOut))if(ex.find("," + x + ",") == std::string::npos)vars.push_back(x);
    TMVA::Reader rd("!Color:Silent");
    std::vector<float> buf(vars.size());
    for(size_t i = 0; i < vars.size(); i++)rd.AddVariable(vars[i].c_str(), &buf[i]);
    rd.BookMVA("BDT", Form("%s/dataset/weights/bdt%d_BDT.weights.xml", dir, p));
    std::vector<std::unique_ptr<TTreeFormula>> fv;
    for(auto& x : vars)fv.emplace_back(new TTreeFormula(x.c_str(), x.c_str(), t));
    const char* pre = p == 0 ? "pop==0 && chi2mu<20 && faces>10 && half==1" : "pop==1 && axis<150 && chi2mu<20 && half==1";
    TTreeFormula fpre("pre", pre, t), fcat("cat", "cat", t), fpr("pr", "prange", t), fch("ch", p == 0 ? "chain" : "chain80", t),
                 fpop("pop", Form("pop==%d && half==1", p), t);
    double n[2][2][5] = {{{0}}};   // [selección 0 cadena, 1 BDT][0 ajuste, 1 ventana][categoría]
    for(Long64_t e = 0; e < t->GetEntries(); e++){
      t->LoadTree(e);
      if(fpop.EvalInstance() == 0)continue;
      const int cat = (int)fcat.EvalInstance(), k = cat <= 2 ? cat : 3;
      const double pr = fpr.EvalInstance();
      const bool win = pr >= 205 && pr <= 240;
      bool sel[2] = {fch.EvalInstance() != 0, false};
      if(fpre.EvalInstance() != 0){
        for(size_t i = 0; i < vars.size(); i++)buf[i] = fv[i]->EvalInstance();
        sel[1] = rd.EvaluateMVA("BDT") > (p == 0 ? cutIn : cutOut);
      }
      for(int s = 0; s < 2; s++){
        if(!sel[s])continue;
        n[s][0][k]++; n[s][0][4]++;
        if(win){ n[s][1][k]++; n[s][1][4]++; }
      }
    }
    printf("\n== %s (BDT > %.2f, %zu entradas) ==\n", p == 0 ? "dentro de la beam slice" : "fuera", p == 0 ? cutIn : cutOut, vars.size());
    printf("  %-22s %8s %8s %8s %8s %8s | %6s %6s %6s | pureza B, C\n", "", "A", "B no A", "C no B", "fondo", "total", "A", "B", "C");
    const char* sn[2] = {"cadena de cortes", "BDT"}, *wn[2] = {"ajuste (150-400)", "ventana"};
    for(int s = 0; s < 2; s++)for(int w = 0; w < 2; w++){
      const double* x = n[s][w];
      printf("  %-9s %-12s", sn[s], wn[w]);
      for(int k = 0; k < 5; k++)printf(" %8.0f", kScale*x[k]);
      printf(" | %6.0f %6.0f %6.0f | %5.1f%% %5.1f%%\n", kScale*x[0], kScale*(x[0] + x[1]), kScale*(x[0] + x[1] + x[2]),
             x[4] ? 100*(x[0] + x[1])/x[4] : 0., x[4] ? 100*(x[0] + x[1] + x[2])/x[4] : 0.);
    }
  }
  return 0;
}
