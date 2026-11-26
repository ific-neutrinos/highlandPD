// Entrenamiento rápido de un BDT con una lista de variables, para la selección hacia delante de variables (26-nov-2026).
// Árbol plano preseleccionado por población (opt_fwd_slim.C). Entrena con la mitad 0 y evalúa con la mitad 1; la figura
// de mérito (error de un recuento de B y de A en la ventana [205, 240], con la estadística de los datos) se calcula
// sobre los árboles TrainTree/TestTree de TMVA, con p de rango y la categoría como espectadores (sin volver a aplicar
// el BDT). El umbral se elige en la mitad 0 (la de entrenamiento) y se evalúa en la 1.
//
//   root -l -b -q 'test/phase2/opt_fwd.C+("slim_0.root", 0, "as_mean5,k_dedx,...", "dir", "tag")'
//
// Añade una línea a <dir>/resultados_<pop>.txt: tag, número de variables, error de B, error de A, ROC, lista.
#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TSystem.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

int opt_fwd(const char* slim, int pop, const char* varList, const char* dir, const char* tag, int nTrees = 200){
  TFile* f = TFile::Open(slim); TTree* t = (TTree*)f->Get("flat");
  gSystem->mkdir(Form("%s/%s", dir, tag), true);
  const TString here = gSystem->WorkingDirectory();
  gSystem->ChangeDirectory(Form("%s/%s", dir, tag));
  std::vector<std::string> vars;
  { std::string s = std::string(varList) + ","; for(size_t a = 0, b; (b = s.find(',', a)) != std::string::npos; a = b + 1)if(b > a)vars.push_back(s.substr(a, b - a)); }
  TFile* fo = TFile::Open("tmva.root", "RECREATE");
  TMVA::Factory fac("fwd", fo, "!V:Silent:!DrawProgressBar:AnalysisType=Classification");
  TMVA::DataLoader dl("dataset");
  for(auto& x : vars)dl.AddVariable(x.c_str(), 'F');
  dl.AddSpectator("prange", 'F'); dl.AddSpectator("cat", 'I');
  const double wB = pop == 0 ? 1. : 3.;   // fuera el fondo está submuestreado 1 de 3 en el árbol preseleccionado
  dl.AddTree(t, "Signal", 1., "cat<=1 && half==0", TMVA::Types::kTraining);
  dl.AddTree(t, "Background", wB, "cat>1 && half==0", TMVA::Types::kTraining);
  dl.AddTree(t, "Signal", 1., "cat<=1 && half==1", TMVA::Types::kTesting);
  dl.AddTree(t, "Background", wB, "cat>1 && half==1", TMVA::Types::kTesting);
  dl.PrepareTrainingAndTestTree("", "NormMode=None:!V");
  fac.BookMethod(&dl, TMVA::Types::kBDT, "BDT", Form("!H:!V:NTrees=%d:MinNodeSize=1%%:MaxDepth=4:BoostType=AdaBoost:AdaBoostBeta=0.3:nCuts=40", nTrees));
  fac.TrainAllMethods(); fac.TestAllMethods(); fac.EvaluateAllMethods();
  const double roc = fac.GetROCIntegral(&dl, "BDT");
  fo->Close();
  // figura de mérito desde los árboles de TMVA
  TFile* fr = TFile::Open("tmva.root");
  TTree* tr[2] = {(TTree*)fr->Get("dataset/TrainTree"), (TTree*)fr->Get("dataset/TestTree")};
  const double kHalf = 2*1.196;
  struct E { float bdt, pr, w; int cat; };
  std::vector<E> ev[2];
  for(int h = 0; h < 2; h++){
    Float_t bdt, pr, w, cat;
    tr[h]->SetBranchAddress("BDT", &bdt); tr[h]->SetBranchAddress("prange", &pr); tr[h]->SetBranchAddress("weight", &w); tr[h]->SetBranchAddress("cat", &cat);
    for(Long64_t e = 0; e < tr[h]->GetEntries(); e++){ tr[h]->GetEntry(e); if(pr >= 205 && pr <= 240)ev[h].push_back({bdt, pr, w, (int)cat}); }
  }
  auto err = [&](int h, double cut, int which){
    double s = 0, n = 0;
    for(auto& x : ev[h]){ if(x.bdt <= cut)continue; n += x.w; if(which ? x.cat <= 1 : x.cat == 0)s += x.w; }
    return s > 0 ? 100.*std::sqrt(n)/s/std::sqrt(kHalf) : 1e9;
  };
  double res[2];
  for(int which = 1; which >= 0; which--){
    double best = 0, eb = 1e9;
    for(double c = -0.6; c <= 0.6; c += 0.01){ const double e = err(0, c, which); if(e < eb){ eb = e; best = c; } }
    res[which] = err(1, best, which);
  }
  gSystem->ChangeDirectory(here);
  FILE* out = fopen(Form("%s/resultados_%d.txt", dir, pop), "a");
  fprintf(out, "%s %zu %.3f %.3f %.4f %s\n", tag, vars.size(), res[1], res[0], roc, varList);
  fclose(out);
  printf("RESULTADO %s: %zu variables, error B %.2f %%, error A %.2f %%, ROC %.4f\n", tag, vars.size(), res[1], res[0], roc);
  gSystem->Exec(Form("rm -rf %s/%s", dir, tag));   // los pesos de cada prueba no hacen falta
  return 0;
}
