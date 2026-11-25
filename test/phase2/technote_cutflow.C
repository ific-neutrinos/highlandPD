// Tablas de cut flow del capítulo de selección de la nota técnica (filas LaTeX): para cada población y cada corte,
// candidatos con cand_level >= k por categoría (señal en capas A, B no A, C no B; fondo por especie true: pi, p, mu
// cósmico, e/gamma, otro), y la
// eficiencia (relativa al corte de población) y la pureza de la señal de la medida B, las dos en porcentaje.
//
//   root -l -b -q 'test/phase2/technote_cutflow.C("xs_v100s_fix.root")'
#include "TFile.h"
#include "TTree.h"
#include <cstdio>
#include <string>
#include <vector>

int technote_cutflow(const char* in){
  TFile* f = TFile::Open(in);
  if(!f || f->IsZombie())return -1;
  TTree* t = (TTree*)f->Get("ana");
  const std::vector<std::string> names[2] = {
    {"in the beam slice", "faces $>35$~cm", "vertex $>15$~cm", "$\\chi^2_\\mu<6$", "$dE/dx$ stop $>3$", "APA gaps", "Michel", "$p$ window"},
    {"out of the beam slice", "not \\code{IsClearCosmic}", "axis $<60$~cm", "$\\chi^2_\\mu<6$", "$dE/dx$ stop $>3$", "$\\chi^2_p>80$",
     "$dE/dx$ start $>1.8$", "APA gaps", "Michel", "$p$ window"}};
  const std::string sig = "cand_truesignalmuon==1";
  const std::vector<std::string> cat = {
    sig + " && (cand_truekorigin==0||cand_truekorigin==1)",
    sig + " && (cand_truekorigin==2||cand_truekorigin==3)",
    sig + " && !(cand_truekorigin>=0&&cand_truekorigin<=3)",
    "cand_truesignalmuon!=1 && abs(cand_truepdg)==211",
    "cand_truesignalmuon!=1 && abs(cand_truepdg)==2212",
    "cand_truesignalmuon!=1 && abs(cand_truepdg)==13 && cand_trueorigin==2",
    "cand_truesignalmuon!=1 && (abs(cand_truepdg)==11||cand_truepdg==22)",
    "cand_truesignalmuon!=1 && !(abs(cand_truepdg)==211||abs(cand_truepdg)==2212||(abs(cand_truepdg)==13&&cand_trueorigin==2)||abs(cand_truepdg)==11||cand_truepdg==22)"};
  auto num = [](double x){ // miles con \,
    char b[32]; snprintf(b, sizeof b, "%.0f", x); std::string s(b), o;
    for(size_t i = 0; i < s.size(); i++){ if(i > 0 && (s.size() - i) % 3 == 0)o += "\\,"; o += s[i]; }
    return o;
  };
  for(int p = 0; p < 2; p++){
    printf("\n%% %s\n", p == 0 ? "dentro de la beam slice" : "fuera");
    double sB0 = 0;
    for(size_t k = 1; k <= names[p].size(); k++){
      std::vector<double> n;
      for(const auto& c : cat)n.push_back(t->Draw("cand_level", Form("cand_pop==%d && cand_level>=%zu && %s", p, k, c.c_str()), "goff"));
      double tot = 0; for(double x : n)tot += x;
      const double sB = n[0] + n[1];
      if(k == 1)sB0 = sB;
      printf("%s", names[p][k-1].c_str());
      for(double x : n)printf(" & %s", num(x).c_str());
      printf(" & %.1f\\%%, %.1f\\%% \\\\\n", sB0 > 0 ? 100*sB/sB0 : 0., tot > 0 ? 100*sB/tot : 0.);
    }
  }
  return 0;
}
