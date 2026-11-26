// Selección multivariable (BDT de TMVA) para los candidatos de secondaryKaonProductionXS, sobre el árbol plano de
// opt_flat.C (26-nov-2026). Un BDT por población; señal = signal muons de B (capas A y B no A), fondo = el resto (C no B
// incluido). p de rango NO es una entrada (es la variable del ajuste): la ventana [205, 240] se aplica después.
//
//   root -l -b -q 'test/phase2/opt_bdt.C+("opt_flat_v100to105d01.root", "opt_bdt")'
//
// Entrenamiento con la mitad 0 (half == 0), evaluación con la mitad 1. Figura de mérito como opt_cuts.C: error
// relativo de un recuento de B (y de A) en la ventana, con la estadística de los datos. El umbral del BDT se elige en
// la mitad 0 (sesgo de sobreentrenamiento) y en la mitad 1 (sesgo de elegir un parámetro: pequeño), y se dan las dos.
// Salida: <dir>/tmva_<pop>.root (TMVA) y pesos en <dir>/dataset/weights.
#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TTreeFormula.h"
#include <memory>
#include "TSystem.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

namespace optb2 {
  const double kHalf = 2*1.196;
  // entradas: comunes y por población (sin truth, sin longitud ni p de rango)
  const std::vector<std::string> kCommon = {"chi2mu", "chi2mu_trimboth", "chi2p", "dedx_stop", "dedx_start", "nmich10", "nmich5", "cnn_michel",
    "michel_own", "nmichel_own", "isshower", "cnn_trk", "cnn_em", "trmean", "trmean_body", "intpeak", "diry_stop", "cosbeam", "ar_charge10",
    "ar_n10", "as_charge10", "ndau_trk", "ndau_shw", "dau_nearstop", "dau_maxmichel", "k_n", "k_dedx", "faces", "gapmin", "as_n5", "as_n10", "as_mean5", "as_mean10"};
  const std::vector<std::string> kIn = {"vtx", "isprimary", "vc_dist", "vc_cover", "dstart_v", "dstop_v"};
  const std::vector<std::string> kOut = {"axis", "clear", "cosmiccrit", "slice_min"};
}

int opt_bdt2(const char* in, const char* dir = "opt_bdt", int popOnly = -1, const char* exclude = "", const char* bdtOpt = "!H:!V:NTrees=400:MinNodeSize=1%:MaxDepth=4:BoostType=AdaBoost:AdaBoostBeta=0.3:nCuts=40"){
  using namespace optb2;
  gSystem->mkdir(dir, true);
  TFile* f = TFile::Open(in);
  // TMVA no admite barras en el nombre del DataLoader: se trabaja dentro de dir (in tiene que ser una ruta absoluta)
  const TString here = gSystem->WorkingDirectory();
  gSystem->ChangeDirectory(dir); TTree* t = (TTree*)f->Get("flat");
  // gapmin no está en el árbol: se añade como alias
  t->SetAlias("gapmin", "min(min(abs(zstart-230),abs(zstart-462)),min(abs(zstop-230),abs(zstop-462)))");
  const char* popName[2] = {"dentro de la beam slice", "fuera"};
  for(int p = 0; p < 2; p++){
    if(popOnly >= 0 && p != popOnly)continue;
    std::vector<std::string> vars0 = kCommon, vars; for(auto& x : (p == 0 ? kIn : kOut))vars0.push_back(x);
    // exclude: variables separadas por comas que no entran (p. ej. las que dependen del momento y esculpen el fondo)
    // exclude que empieza por "+": lista de las únicas variables que entran (BDT reducido, en ese orden)
    if(exclude[0] == '+'){
      std::string inc = std::string(exclude + 1) + ",";
      for(size_t a = 0, b; (b = inc.find(',', a)) != std::string::npos; a = b + 1)if(b > a)vars.push_back(inc.substr(a, b - a));
    } else {
      const std::string ex = std::string(",") + exclude + ",";
      for(auto& x : vars0)if(ex.find("," + x + ",") == std::string::npos)vars.push_back(x);
    }
    // preselección suelta de entrenamiento (como opt_cuts); fuera, el fondo se submuestrea para que quepa
    const TCut pre = p == 0 ? "pop==0 && chi2mu<20 && faces>10" : "pop==1 && axis<150 && chi2mu<20";
    const TCut sig = "cat<=1", bkg = "cat>1";
    const TCut bkgSub = p == 0 ? TCut("1") : TCut("entry%3==0");
    TString outName = Form("tmva_%d.root", p);
    TFile* fo = TFile::Open(outName, "RECREATE");
    TMVA::Factory fac(Form("bdt%d", p), fo, "!V:!Silent:AnalysisType=Classification");
    TMVA::DataLoader dl("dataset");
    for(auto& x : vars)dl.AddVariable(x.c_str(), 'F');
    dl.AddTree(t, "Signal", 1., pre && sig && "half==0", TMVA::Types::kTraining);
    dl.AddTree(t, "Background", (p == 0 ? 1. : 3.), pre && bkg && bkgSub && "half==0", TMVA::Types::kTraining);
    dl.AddTree(t, "Signal", 1., pre && sig && "half==1", TMVA::Types::kTesting);
    dl.AddTree(t, "Background", (p == 0 ? 1. : 3.), pre && bkg && bkgSub && "half==1", TMVA::Types::kTesting);
    dl.PrepareTrainingAndTestTree("", "NormMode=None:!V");
    fac.BookMethod(&dl, TMVA::Types::kBDT, "BDT", bdtOpt);
    fac.TrainAllMethods(); fac.TestAllMethods(); fac.EvaluateAllMethods();
    fo->Close();

    // aplicación a todos los candidatos con la preselección (las dos mitades, sin submuestrear)
    TMVA::Reader rd("!Color:Silent");
    std::vector<float> buf(vars.size());
    for(size_t i = 0; i < vars.size(); i++)rd.AddVariable(vars[i].c_str(), &buf[i]);
    rd.BookMVA("BDT", Form("dataset/weights/bdt%d_BDT.weights.xml", p));
    std::vector<std::unique_ptr<TTreeFormula>> fv;
    for(auto& x : vars)fv.emplace_back(new TTreeFormula(x.c_str(), x.c_str(), t));
    TTreeFormula fpre("pre", pre.GetTitle(), t), fcat("cat", "cat", t), fhalf("half", "half", t), fpr("pr", "prange", t),
                 fch("ch", p == 0 ? "chain" : "chain80", t);
    struct Sc { float bdt; int cat, half; float prange; int chain; };
    std::vector<Sc> sc;
    for(Long64_t e = 0; e < t->GetEntries(); e++){
      t->LoadTree(e);
      if(fpre.EvalInstance() == 0)continue;
      for(size_t i = 0; i < vars.size(); i++)buf[i] = fv[i]->EvalInstance();
      sc.push_back({(float)rd.EvaluateMVA("BDT"), (int)fcat.EvalInstance(), (int)fhalf.EvalInstance(), (float)fpr.EvalInstance(), (int)fch.EvalInstance()});
    }
    auto count = [&](int h, double cut, bool useChain, double& sA, double& sB, double& n){
      sA = sB = n = 0;
      for(auto& c : sc){
        if(c.half != h || c.prange < 205 || c.prange > 240)continue;
        if(useChain ? !c.chain : c.bdt <= cut)continue;
        n++; if(c.cat <= 1)sB++; if(c.cat == 0)sA++;
      }
    };
    auto err = [](double s, double n){ return s > 0 ? 100.*std::sqrt(n)/s/std::sqrt(kHalf) : 1e9; };
    printf("\n== %s: BDT con %zu entradas ==\n", popName[p], vars.size());
    for(int which = 1; which >= 0; which--){
      double best[2] = {0, 0}, eBest[2] = {1e9, 1e9};
      for(int h = 0; h < 2; h++)for(double cut = -0.6; cut <= 0.6; cut += 0.01){
        double sA, sB, n; count(h, cut, false, sA, sB, n);
        const double e = err(which ? sB : sA, n);
        if(e < eBest[h]){ eBest[h] = e; best[h] = cut; }
      }
      double sA, sB, n;
      count(1, 0, true, sA, sB, n);
      printf("  error de %s, mitad 1: cortes actuales %.2f %% (S %.0f, N %.0f)\n", which ? "B" : "A", err(which ? sB : sA, n), which ? sB : sA, n);
      count(1, best[0], false, sA, sB, n);
      printf("    BDT > %.2f (elegido en la mitad 0, la de entrenamiento): %.2f %% (S %.0f, N %.0f)\n", best[0], err(which ? sB : sA, n), which ? sB : sA, n);
      count(1, best[1], false, sA, sB, n);
      printf("    BDT > %.2f (elegido en la mitad 1): %.2f %% (S %.0f, N %.0f)\n", best[1], err(which ? sB : sA, n), which ? sB : sA, n);
    }
    // espectro de p de rango del fondo con el BDT (mitad 1): ¿esculpe un pico?
    printf("  p de rango del fondo (mitad 1) en bins de 25 MeV/c de 150 a 400, cortes actuales / BDT > %.2f:\n   ", 0.);
    std::vector<int> hc(10, 0), hb(10, 0);
    double bestB = 0; { double e0 = 1e9; for(double cut = -0.6; cut <= 0.6; cut += 0.01){ double sA, sB, n; count(1, cut, false, sA, sB, n); if(err(sB, n) < e0){ e0 = err(sB, n); bestB = cut; } } }
    for(auto& c : sc){ if(c.half != 1 || c.cat <= 1)continue; const int b = (int)((c.prange - 150)/25); if(b < 0 || b > 9)continue; if(c.chain)hc[b]++; if(c.bdt > bestB)hb[b]++; }
    for(int b = 0; b < 10; b++)printf(" %d/%d", hc[b], hb[b]); printf("   (BDT > %.2f)\n", bestB);
  }
  gSystem->ChangeDirectory(here);
  return 0;
}
