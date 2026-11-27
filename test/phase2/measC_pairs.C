// AVISO (27-nov-2026): la etiqueta reconstruida de la K de la parte (2) sale mal (cand_k_dedx, cand_k_n leídos con
// TTreeFormula dan valores erróneos, igual que el BDT evaluado así); usar measC_compare.C, que lee con SetBranchAddress.
// Medida C (27-nov-2026): ¿cuántas parejas K -> mu tendríamos? Con la truth de C de secondaryKaonProductionXS
// (highlandPD 675d80c, en los xs_*d01b).
//
// (1) Truth (árbol truth, un registro por signal muon): ¿está reconstruida la K madre? (algún PFParticle cuya partícula
//     true es la K; el mejor con completitud y pureza > 0,5; traza o cascada; longitudes true y reconstruida).
// (2) Selección del mu (árbol ana, solo la mitad impar, x2: la par entrenó el BDT): BDT dentro de la beam slice (10
//     variables, > cut) y cadena de cortes dentro (chain) y fuera (eje a 80 cm, como opt_flat.C), en la ventana
//     [205, 240] MeV/c. Para los signal muons seleccionados: ¿el candidato a K al inicio es la madre true
//     (cand_k_truek)? Y una etiqueta reconstruida de la K: partícula con un extremo a < 10 cm del inicio del mu, traza,
//     con dE/dx en ese extremo > kDedx (pico de Bragg de la K que para): eficiencia en la señal y en el fondo.
// Proyección a los datos de 6 GeV: MC x 1/1,196 (xs_feasibility.C).
//
// La decisión del BDT se toma del árbol plano (opt_flat3, mismas entradas que el entrenamiento) y se cruza con el
// candidato del microtree por evento (entry) y p de rango; xs y flat deben ser de la misma muestra y en el mismo orden
// (vale xs_v100to105d01b con opt_flat3_v100to105d01: comprobado el 27-nov).
//
//   root -l -b -q 'test/phase2/measC_pairs.C+("xs_v100to105d01b.root", "opt_flat3_v100to105d01.root", "opt_bdt_final/dataset/weights/bdt0_BDT.weights.xml")'
#include "TFile.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "TMVA/Reader.h"
#include <cmath>
#include <cstdio>
#include <memory>
#include <set>
#include <utility>
#include <string>
#include <vector>

int measC_pairs(const char* xs, const char* flat, const char* weights, double cut = 0.19){
  const double kData = 1/1.196;
  TFile fx(xs);
  // ---------- (1) truth ----------
  TTree* tr = (TTree*)fx.Get("truth");
  const char* layer[3] = {"C", "B", "A"};
  const char* lsel[3] = {"truemu_korigin>=0", "truemu_korigin>=0 && truemu_korigin<=3", "truemu_korigin>=0 && truemu_korigin<=1"};
  printf("(1) Truth: signal muons por capa (C contiene a B, B a A) y su K madre reconstruida\n");
  printf("  %-2s %7s | K con algún PFParticle | mejor eff y pur > 0,5 | y es traza | longitud true de la K (mediana) | reco/true (mediana, traza buena)\n", "", "mu");
  for(int l = 0; l < 3; l++){
    const TString s = lsel[l];
    const double n = tr->GetEntries(s), n1 = tr->GetEntries(s + " && truek_nreco>0"), n2 = tr->GetEntries(s + " && truek_besteff>0.5 && truek_bestpur>0.5"),
                 n3 = tr->GetEntries(s + " && truek_besteff>0.5 && truek_bestpur>0.5 && truek_bestistrack==1");
    tr->Draw("truek_length", s + " && truek_length>=0", "goff"); std::vector<double> v(tr->GetV1(), tr->GetV1() + tr->GetSelectedRows()); std::sort(v.begin(), v.end());
    tr->Draw("truek_bestlength/truek_length", s + " && truek_besteff>0.5 && truek_bestpur>0.5 && truek_bestistrack==1 && truek_length>0", "goff");
    std::vector<double> r(tr->GetV1(), tr->GetV1() + tr->GetSelectedRows()); std::sort(r.begin(), r.end());
    printf("  %-2s %7.0f | %5.1f %% | %5.1f %% | %5.1f %% | %5.1f cm | %.2f\n", layer[l], n, 100*n1/n, 100*n2/n, 100*n3/n, v.empty() ? -1. : v[v.size()/2], r.empty() ? -1. : r[r.size()/2]);
  }
  // ---------- (2) selección ----------
  // decisión del BDT desde el árbol plano: (entry, p de rango) de los candidatos dentro, impares, que pasan
  std::set<std::pair<Long64_t, float>> bdtPass;
  { TFile ff(flat); TTree* tf = (TTree*)ff.Get("flat");
    const char* vs[10] = {"as_mean5", "k_dedx", "k_n", "cnn_michel", "michel_own", "nmich10", "cosbeam", "dedx_stop", "ar_charge10", "intpeak"};
    TMVA::Reader rf("!Color:Silent"); float bf[10]; Float_t vf[10], fpop, fpr, fchi2, ffaces; Long64_t fentry;
    tf->SetBranchStatus("*", 0);
    for(int i = 0; i < 10; i++){ rf.AddVariable(vs[i], &bf[i]); tf->SetBranchStatus(vs[i], 1); tf->SetBranchAddress(vs[i], &vf[i]); }
    for(const char* b : {"pop", "prange", "chi2mu", "faces", "entry"})tf->SetBranchStatus(b, 1);
    tf->SetBranchAddress("pop", &fpop); tf->SetBranchAddress("prange", &fpr); tf->SetBranchAddress("chi2mu", &fchi2); tf->SetBranchAddress("faces", &ffaces); tf->SetBranchAddress("entry", &fentry);
    rf.BookMVA("BDT", weights);
    for(Long64_t i = 0; i < tf->GetEntries(); i++){
      tf->GetEntry(i);
      if(fentry % 2 != 1 || (int)fpop != 0 || fpr < 205 || fpr > 240 || fchi2 >= 20 || ffaces <= 10)continue;
      for(int k = 0; k < 10; k++)bf[k] = vf[k];
      if(rf.EvaluateMVA("BDT") > cut)bdtPass.insert({fentry, fpr});
    }
    printf("\nBDT (árbol plano): %zu candidatos dentro, impares, en la ventana\n", bdtPass.size()); }
  TTree* t = (TTree*)fx.Get("ana");
  const std::vector<std::pair<std::string, std::string>> bdtVars = {
    {"as_mean5", "cand_as_n5>0 ? cand_as_charge5/cand_as_n5 : -1"}, {"k_dedx", "cand_k_dedx"}, {"k_n", "cand_k_n"}, {"cnn_michel", "cand_cnn_michel"},
    {"michel_own", "cand_michel_own"}, {"nmich10", "cand_ar_nmichel10"}, {"cosbeam", "cand_cosbeam"}, {"dedx_stop", "cand_dedx_stop"},
    {"ar_charge10", "cand_ar_charge10"}, {"intpeak", "cand_intpeak"}};
  const std::vector<std::pair<std::string, std::string>> other = {
    {"pop", "cand_pop"}, {"level", "cand_level"}, {"prange", "cand_prange"}, {"chi2mu", "cand_chi2mu"}, {"faces", "min(cand_faces_stop,cand_faces_start)"},
    {"sig", "cand_truesignalmuon"}, {"kor", "cand_truekorigin"}, {"ktrue", "cand_k_truek"}, {"kistrack", "cand_k_istrack"}, {"klen", "cand_k_length"},
    {"clear", "cand_clearcosmic"}, {"axis", "cand_axis_min"}, {"chi2p", "cand_chi2p"}, {"dedx_start", "cand_dedx_start"},
    {"zstart", "cand_start[][2]"}, {"zstop", "cand_stop[][2]"}, {"nmich", "cand_ar_nmichel10"}};
  std::vector<std::unique_ptr<TTreeFormula>> fb, fo;
  for(auto& v : bdtVars)fb.emplace_back(new TTreeFormula(v.first.c_str(), v.second.c_str(), t));
  for(auto& v : other)fo.emplace_back(new TTreeFormula(v.first.c_str(), v.second.c_str(), t));
  TTreeFormula nc("nc", "ncand", t);
  auto gap = [](double z){ return std::min(std::abs(z - 230.), std::abs(z - 462.)); };
  // selecciones: 0 BDT dentro, 1 cortes dentro, 2 cortes fuera (eje 80)
  const char* sn[3] = {"BDT dentro de la beam slice", "cortes dentro", "cortes fuera (eje 80 cm)"};
  const std::vector<double> kDedx = {0, 5, 8, 12};
  // [sel][capa 0 C,1 B,2 A][0 todos, 1 ktrue, 2.. etiqueta reco por umbral]; fondo [sel][umbral]
  double nS[3][3][2 + 4] = {{{0}}}, nB[3][1 + 4] = {{0}};
  for(Long64_t e = 1; e < t->GetEntries(); e += 2){
    t->LoadTree(e);
    for(auto& x : fb)x->GetNdata(); for(auto& x : fo)x->GetNdata();
    const int n = (int)nc.EvalInstance(0);
    for(int j = 0; j < n; j++){
      std::vector<double> o(other.size()); for(size_t i = 0; i < other.size(); i++)o[i] = fo[i]->EvalInstance(j);
      const int pop = (int)o[0], level = (int)o[1]; const double pr = o[2], chi2mu = o[3];
      if(pop < 0 || pop > 1 || level < 1 || chi2mu < 0 || pr < 205 || pr > 240)continue;
      bool pass[3] = {false, false, false};
      if(pop == 0){
        pass[1] = level >= 7;
        pass[0] = bdtPass.count({e, (float)pr}) > 0;
      } else {
        const bool michel = o[16] < 0 ? fb[3]->EvalInstance(j) > 0.02 : o[16] >= 2;   // como opt_flat.C
        pass[2] = !o[10] && o[11] < 80 && chi2mu < 6 && fb[7]->EvalInstance(j) > 3 && o[12] > 80 && o[13] > 1.8 && gap(o[14]) >= 5 && gap(o[15]) >= 5 && michel;
      }
      const int sig = (int)o[5], kor = (int)o[6], ktrue = (int)o[7], kistrack = (int)o[8];
      const double kd = fb[1]->EvalInstance(j), kn = fb[2]->EvalInstance(j);
      for(int s = 0; s < 3; s++){
        if(!pass[s])continue;
        bool tag[4]; for(int k = 0; k < 4; k++)tag[k] = kn > 0 && kistrack == 1 && kd > kDedx[k];
        if(sig == 1){
          for(int l = 0; l < 3; l++){
            const bool in = l == 0 ? true : l == 1 ? (kor >= 0 && kor <= 3) : (kor >= 0 && kor <= 1);
            if(!in)continue;
            nS[s][l][0] += 2; if(ktrue == 1)nS[s][l][1] += 2;
            for(int k = 0; k < 4; k++)if(tag[k])nS[s][l][2 + k] += 2;
          }
        } else { nB[s][0] += 2; for(int k = 0; k < 4; k++)if(tag[k])nB[s][1 + k] += 2; }
      }
    }
    if(e % 50001 == 0)fprintf(stderr, "[avance] %lld de %lld\n", e, t->GetEntries());
  }
  printf("\n(2) Signal muons seleccionados en la ventana (mitad impar x2; entre paréntesis, proyección a los datos de 6 GeV)\n");
  printf("    K madre = el candidato a K al inicio es la madre true; etiqueta = traza a < 10 cm del inicio con dE/dx en ese extremo > umbral\n");
  for(int s = 0; s < 3; s++){
    printf("  %s\n", sn[s]);
    for(int l = 0; l < 3; l++){
      const double* x = nS[s][l];
      printf("    %-2s mu %6.0f (%5.0f) | K madre al inicio %5.1f %% (%5.0f) | etiqueta dE/dx > 0/5/8/12: %5.1f %5.1f %5.1f %5.1f %%\n", layer[l], x[0], x[0]*kData,
             x[0] > 0 ? 100*x[1]/x[0] : 0., x[1]*kData, 100*x[2]/std::max(1., x[0]), 100*x[3]/std::max(1., x[0]), 100*x[4]/std::max(1., x[0]), 100*x[5]/std::max(1., x[0]));
    }
    printf("    fondo %6.0f | con etiqueta dE/dx > 0/5/8/12: %5.1f %5.1f %5.1f %5.1f %% | pureza (C) con etiqueta: %5.1f %5.1f %5.1f %5.1f %% (sin etiqueta %5.1f %%)\n",
           nB[s][0], 100*nB[s][1]/std::max(1., nB[s][0]), 100*nB[s][2]/std::max(1., nB[s][0]), 100*nB[s][3]/std::max(1., nB[s][0]), 100*nB[s][4]/std::max(1., nB[s][0]),
           100*nS[s][0][2]/std::max(1., nS[s][0][2] + nB[s][1]), 100*nS[s][0][3]/std::max(1., nS[s][0][3] + nB[s][2]),
           100*nS[s][0][4]/std::max(1., nS[s][0][4] + nB[s][3]), 100*nS[s][0][5]/std::max(1., nS[s][0][5] + nB[s][4]), 100*nS[s][0][0]/std::max(1., nS[s][0][0] + nB[s][0]));
  }
  return 0;
}
