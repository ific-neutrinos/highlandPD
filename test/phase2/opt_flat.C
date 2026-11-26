// Árbol plano de candidatos para optimizar la selección de secondaryKaonProductionXS (cortes y multivariable),
// 26-nov-2026. Un candidato por entrada: los que pasan el corte de población (cand_level >= 1), con calorimetría
// (cand_chi2mu >= 0) y p de rango en el rango del ajuste [150, 400] MeV/c. p de rango se guarda pero no es una
// variable de entrada: es la variable del ajuste.
//
//   root -l -b -q 'test/phase2/opt_flat.C+("xs_v100to105d01.root", "opt_flat_v100to105d01.root")'
//
// cat: 0 A, 1 B no A, 2 C no B, 3 pi, 4 p, 5 mu, 6 e/gamma, 7 otro. half: paridad del número de entrada del evento
// (mitades independientes para entrenar y evaluar). chain: pasa la cadena actual (cand_level de la cadena entera)
// salvo la ventana; chain80: con el eje a 80 cm emulado (fuera).
#include "TFile.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

int opt_flat(const char* in, const char* out){
  TFile* f = TFile::Open(in);
  if(!f || f->IsZombie())return -1;
  TTree* t = (TTree*)f->Get("ana");
  // variables (nombre en el árbol plano, expresión sobre el microtree)
  const std::vector<std::pair<std::string, std::string>> vars = {
    {"pop", "cand_pop"}, {"level", "cand_level"}, {"sig", "cand_truesignalmuon"}, {"kor", "cand_truekorigin"}, {"pdg", "cand_truepdg"},
    {"prange", "cand_prange"},
    {"faces", "min(cand_faces_stop,cand_faces_start)"}, {"vtx", "cand_vtx_min"}, {"axis", "cand_axis_min"}, {"clear", "cand_clearcosmic"},
    {"chi2mu", "cand_chi2mu"}, {"chi2mu_trimboth", "cand_chi2mu_trimboth"}, {"chi2p", "cand_chi2p"}, {"dedx_stop", "cand_dedx_stop"},
    {"dedx_start", "cand_dedx_start"}, {"zstart", "cand_start[][2]"}, {"zstop", "cand_stop[][2]"},
    {"nmich10", "cand_ar_nmichel10"}, {"nmich5", "cand_ar_nmichel5"}, {"cnn_michel", "cand_cnn_michel"}, {"michel_own", "cand_michel_own"},
    {"nmichel_own", "cand_nmichel_own"}, {"isshower", "cand_isshower"}, {"cnn_trk", "cand_cnn_trk"}, {"cnn_em", "cand_cnn_em"},
    {"trmean", "cand_trmean"}, {"trmean_body", "cand_trmean_body"}, {"intpeak", "cand_intpeak"}, {"diry_stop", "cand_diry_stop"},
    {"cosbeam", "cand_cosbeam"}, {"ar_charge10", "cand_ar_charge10"}, {"ar_n10", "cand_ar_n10"}, {"as_charge10", "cand_as_charge10"},
    {"ndau_trk", "cand_ndau_trk"}, {"ndau_shw", "cand_ndau_shw"}, {"dau_nearstop", "cand_dau_nearstop"}, {"dau_maxmichel", "cand_dau_maxmichel"},
    {"k_n", "cand_k_n"}, {"k_dedx", "cand_k_dedx"}, {"isprimary", "cand_isprimary"}, {"slice_min", "cand_slice_min"},
    {"vc_dist", "cand_vc_dist"}, {"vc_cover", "cand_vc_cover"}, {"dstart_v", "cand_dstart_v"}, {"dstop_v", "cand_dstop_v"},
    {"cosmiccrit", "cand_cosmiccrit"}, {"sh_lin", "cand_sh_lin"},
    // actividad alrededor del inicio del mu (idea del usuario, 26-nov: la K que para deja hits con más carga)
    {"as_n5", "cand_as_n5"}, {"as_n10", "cand_as_n10"}, {"as_charge5", "cand_as_charge5"},
    {"as_mean5", "cand_as_n5>0 ? cand_as_charge5/cand_as_n5 : -1"}, {"as_mean10", "cand_as_n10>0 ? cand_as_charge10/cand_as_n10 : -1"}};
  std::vector<std::unique_ptr<TTreeFormula>> fo;
  for(auto& v : vars){
    fo.emplace_back(new TTreeFormula(v.first.c_str(), v.second.c_str(), t));
    if(fo.back()->GetNdim() == 0){ printf("Expresión no válida: %s\n", v.second.c_str()); return -1; }
  }
  TFile fo_(out, "RECREATE");
  TTree o("flat", "candidatos");
  std::vector<float> val(vars.size());
  Int_t cat, half, chain, chain80; Long64_t entry;
  for(size_t i = 0; i < vars.size(); i++)o.Branch(vars[i].first.c_str(), &val[i]);
  o.Branch("cat", &cat); o.Branch("half", &half); o.Branch("chain", &chain); o.Branch("chain80", &chain80); o.Branch("entry", &entry);
  TTreeFormula nc("nc", "ncand", t);
  const int nCut[2] = {8, 10};
  auto gap = [](double z){ return std::min(std::abs(z - 230.), std::abs(z - 462.)); };
  for(Long64_t e = 0; e < t->GetEntries(); e++){
    t->LoadTree(e);
    for(auto& x : fo)x->GetNdata();
    const int n = (int)nc.EvalInstance(0);
    for(int j = 0; j < n; j++){
      for(size_t i = 0; i < vars.size(); i++)val[i] = fo[i]->EvalInstance(j);
      const int pop = (int)val[0], level = (int)val[1];
      if(pop < 0 || pop > 1 || level < 1)continue;
      const float prange = val[5], chi2mu = val[10];
      if(chi2mu < 0 || prange < 150 || prange > 400)continue;
      const int sig = (int)val[2], kor = (int)val[3], a = std::abs((int)val[4]);
      cat = sig == 1 ? ((kor == 0 || kor == 1) ? 0 : (kor == 2 || kor == 3) ? 1 : 2)
                     : (a == 211 ? 3 : a == 2212 ? 4 : a == 13 ? 5 : (a == 11 || a == 22) ? 6 : 7);
      half = e % 2; entry = e;
      chain = level >= nCut[pop] - 1;   // todos los cortes salvo la ventana
      // eje a 80 cm (fuera): la cadena sin el eje, emulada con las variables (como xs_nminus1)
      chain80 = chain;
      if(pop == 1){
        const float axis = val[8], clear = val[9], dstop = val[13], dstart = val[14], chi2p = val[12];
        const float nm = val[17], cnnm = val[19];
        const bool michel = nm < 0 ? cnnm > 0.02 : nm >= 2;
        const bool gaps = gap(val[15]) >= 5 && gap(val[16]) >= 5;
        chain80 = !clear && axis < 80 && chi2mu < 6 && dstop > 3 && chi2p > 80 && dstart > 1.8 && gaps && michel;
      }
      o.Fill();
    }
    if(e % 50000 == 0)fprintf(stderr, "[avance] %lld de %lld\n", e, t->GetEntries());
  }
  fo_.cd(); o.Write();
  printf("Candidatos: %lld en %s\n", o.GetEntries(), out);
  return 0;
}
