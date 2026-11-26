// Combinación de las dos selecciones (26-nov-2026): secondaryKaonProductionXS (mu solo) y las parejas K -> mu de la tesis
// (secondaryKaonPairXS + cortes de macros/tesis/xs/selection.C), sobre los mismos minitrees y con la misma señal.
//
//   root -l -b -q 'test/phase2/combine_pair_xs.C+("xs_v100d01.root", "pair_v100d01.root", 6)'
//
// k: factor de estadística para el error de un recuento (k = 6*73070/N_MC ~ datos a 6 GeV).
// Candidatos por evento (run, subrun, evt): los nuestros con la cadena hasta Michel (w0) o entera (w1), y las hijas de las
// parejas que pasan. Un candidato nuestro y una hija de pareja son el mismo objeto si sus extremos coinciden a < dMatch cm
// (en cualquier orientación) o, si los dos son signal muon, si es el mismo mu true (p true, |dp| < 0,5 MeV/c). Así la
// pureza de la unión no cuenta dos veces el mismo candidato.
// Selecciones: nuestra (w0, w1); nuestra con K en el inicio del mu (cand_k_n > 0; con traza; con dE/dx de K > 3); parejas;
// unión y solape de nuestra (w0, w1) con las parejas; parejas que no están en la nuestra. Por capa A y B: candidatos,
// signal muons de la capa (seleccionados entre los candidatos), eficiencia (sobre los signal muons true del árbol truth),
// pureza y error relativo de un recuento sqrt(k(S+B))/(kS).
#include "TFile.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include <cmath>
#include <cstdio>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace cmb {
  using Key = std::tuple<int, int, int>;
  bool InLayer(int kor, int l){ return l == 1 ? (kor >= 0 && kor <= 3) : (kor == 0 || kor == 1); }   // 0 A, 1 B
  struct Cand { float a[3], b[3]; int sig, kor; double pTrue; int lev, pop, kN, kTrk; float kdEdx; };
  struct PairC { float a[3], b[3]; int sig, kor; double pTrue; };
  double D(const float* x, const float* y){ return std::sqrt(std::pow(x[0]-y[0],2) + std::pow(x[1]-y[1],2) + std::pow(x[2]-y[2],2)); }
  bool Same(const Cand& c, const PairC& p, double dMatch){
    if(c.sig == 1 && p.sig == 1 && std::abs(c.pTrue - p.pTrue) < 0.5)return true;
    const double d1 = std::max(D(c.a, p.a), D(c.b, p.b)), d2 = std::max(D(c.a, p.b), D(c.b, p.a));
    return std::min(d1, d2) < dMatch;
  }
  struct Count { long cand = 0, sig[2] = {0, 0}; std::set<std::tuple<int,int,int,long>> mu[2]; };   // mu: signal muons distintos
  void Add(Count& n, const Key& k, int sig, int kor, double pTrue){
    n.cand++;
    if(sig != 1)return;
    for(int l = 0; l < 2; l++)if(InLayer(kor, l)){ n.sig[l]++; n.mu[l].insert({std::get<0>(k), std::get<1>(k), std::get<2>(k), std::lround(pTrue*10)}); }
  }
}

int combine_pair_xs(const char* xsFile, const char* pairFile, double k = 6, double dMatch = 3,
                    const char* pairCut = "candidates_dau_mom_muon>0.16 && candidates_dau_mom_muon<0.265 && candidates_distance_dau<5 && (candidates_truncated_dedx+candidates_dau_truncated_dedx)/2<4"){
  using namespace cmb;
  TFile* fx = TFile::Open(xsFile); TFile* fp = TFile::Open(pairFile);
  if(!fx || fx->IsZombie() || !fp || fp->IsZombie())return -1;
  TTree* tru = (TTree*)fx->Get("truth"); TTree* ax = (TTree*)fx->Get("ana"); TTree* ap = (TTree*)fp->Get("ana");
  if(!tru || !ax || !ap)return -1;
  const int nCut[2] = {8, 10};   // cortes de cada cadena contando el de población (cand_level de la cadena entera; technote_cutflow)

  long nTrue[2] = {0, 0};
  std::set<Key> evTruth;
  {
    Int_t run, sub, evt, kor;
    tru->SetBranchAddress("run", &run); tru->SetBranchAddress("subrun", &sub); tru->SetBranchAddress("evt", &evt); tru->SetBranchAddress("truemu_korigin", &kor);
    for(Long64_t e = 0; e < tru->GetEntries(); e++){ tru->GetEntry(e); for(int l = 0; l < 2; l++)if(InLayer(kor, l))nTrue[l]++; evTruth.insert({run, sub, evt}); }
  }

  // nuestros candidatos con la cadena hasta Michel (los de w1 son un subconjunto: lev mayor)
  std::map<Key, std::vector<Cand>> ours;
  {
    const int N = 100;
    Int_t run, sub, evt, nc, pop[N], lev[N], sig[N], kor[N], kn[N], ktrk[N]; Float_t st[N][3], sp[N][3], kd[N], pt[N];
    ax->SetBranchStatus("*", 0);
    for(const char* b : {"run", "subrun", "evt", "ncand", "cand_pop", "cand_level", "cand_truesignalmuon", "cand_truekorigin", "cand_start", "cand_stop",
                         "cand_k_n", "cand_k_istrack", "cand_k_dedx", "cand_truemom"})ax->SetBranchStatus(b, 1);
    ax->SetBranchAddress("run", &run); ax->SetBranchAddress("subrun", &sub); ax->SetBranchAddress("evt", &evt); ax->SetBranchAddress("ncand", &nc);
    ax->SetBranchAddress("cand_pop", pop); ax->SetBranchAddress("cand_level", lev); ax->SetBranchAddress("cand_truesignalmuon", sig);
    ax->SetBranchAddress("cand_truekorigin", kor); ax->SetBranchAddress("cand_start", st); ax->SetBranchAddress("cand_stop", sp);
    ax->SetBranchAddress("cand_k_n", kn); ax->SetBranchAddress("cand_k_istrack", ktrk); ax->SetBranchAddress("cand_k_dedx", kd); ax->SetBranchAddress("cand_truemom", pt);
    for(Long64_t e = 0; e < ax->GetEntries(); e++){
      ax->GetEntry(e);
      for(int j = 0; j < nc && j < N; j++){
        if(pop[j] < 0 || pop[j] > 1 || lev[j] < nCut[pop[j]] - 1)continue;
        Cand c; for(int i = 0; i < 3; i++){ c.a[i] = st[j][i]; c.b[i] = sp[j][i]; }
        c.sig = sig[j]; c.kor = kor[j]; c.pTrue = pt[j]; c.lev = lev[j]; c.pop = pop[j]; c.kN = kn[j]; c.kTrk = ktrk[j]; c.kdEdx = kd[j];
        ours[{run, sub, evt}].push_back(c);
      }
    }
  }
  // hijas de las parejas que pasan
  std::map<Key, std::vector<PairC>> pairs;
  long pairEvNotInXS = 0;
  {
    std::unique_ptr<TTreeFormula> fCut(new TTreeFormula("cut", pairCut, ap));
    if(fCut->GetNdim() == 0){ printf("Corte de parejas no válido\n"); return -1; }
    std::unique_ptr<TTreeFormula> fP(new TTreeFormula("p", "candidates_dau_truemom", ap));
    const int N = 200;
    Int_t run, sub, evt, acc[1], nc, sig[N], kor[N]; Float_t dp[N][4], de[N][4];
    ap->SetBranchAddress("run", &run); ap->SetBranchAddress("subrun", &sub); ap->SetBranchAddress("evt", &evt); ap->SetBranchAddress("accum_level", acc);
    ap->SetBranchAddress("ncandidates", &nc); ap->SetBranchAddress("candidatedau_truesignalmuon", sig); ap->SetBranchAddress("candidatedau_truekorigin", kor);
    ap->SetBranchAddress("candidates_dau_pos", dp); ap->SetBranchAddress("candidates_dau_endpos", de);
    for(Long64_t e = 0; e < ap->GetEntries(); e++){
      ap->LoadTree(e); ap->GetEntry(e);
      if(acc[0] < 3 || nc <= 0)continue;
      fCut->GetNdata(); fP->GetNdata();
      for(int j = 0; j < nc && j < N; j++){
        if(fCut->EvalInstance(j) == 0)continue;
        PairC p; for(int i = 0; i < 3; i++){ p.a[i] = dp[j][i]; p.b[i] = de[j][i]; }
        p.sig = sig[j]; p.kor = kor[j]; p.pTrue = 1000.*fP->EvalInstance(j);
        pairs[{run, sub, evt}].push_back(p);
      }
    }
    for(auto& kv : pairs)if(!evTruth.count(kv.first) && !ours.count(kv.first))pairEvNotInXS++;
  }

  // selecciones
  enum { kO0, kO1, kOK, kOKt, kOKd, kOK1, kP, kU0, kU1, kI0, kI1, kPnotO0, kN };
  const char* name[kN] = {"nuestra, hasta Michel", "nuestra, con ventana de p", "  ... w0 con K en el inicio", "  ... w0 con K traza",
                          "  ... w0 con K traza y dE/dx > 3", "  ... w1 con K en el inicio", "parejas (tesis)", "unión w0 + parejas", "unión w1 + parejas",
                          "solape w0 y parejas", "solape w1 y parejas", "parejas que no están en w0"};
  Count n[kN];
  std::set<Key> keys; for(auto& kv : ours)keys.insert(kv.first); for(auto& kv : pairs)keys.insert(kv.first);
  long nMatchBkg = 0;
  for(const Key& key : keys){
    static const std::vector<Cand> noC; static const std::vector<PairC> noP;
    const auto& vc = ours.count(key) ? ours.at(key) : noC;
    const auto& vp = pairs.count(key) ? pairs.at(key) : noP;
    std::vector<int> pairInW0(vp.size(), 0), pairInW1(vp.size(), 0);
    for(const Cand& c : vc){
      const bool w1 = c.lev >= nCut[c.pop];
      bool inP = false;
      for(size_t j = 0; j < vp.size(); j++)if(Same(c, vp[j], dMatch)){ inP = true; pairInW0[j] = 1; if(w1)pairInW1[j] = 1; if(c.sig != 1)nMatchBkg++; }
      Add(n[kO0], key, c.sig, c.kor, c.pTrue); Add(n[kU0], key, c.sig, c.kor, c.pTrue);
      if(c.kN > 0)Add(n[kOK], key, c.sig, c.kor, c.pTrue);
      if(c.kN > 0 && c.kTrk == 1)Add(n[kOKt], key, c.sig, c.kor, c.pTrue);
      if(c.kN > 0 && c.kTrk == 1 && c.kdEdx > 3)Add(n[kOKd], key, c.sig, c.kor, c.pTrue);
      if(inP)Add(n[kI0], key, c.sig, c.kor, c.pTrue);
      if(w1){ Add(n[kO1], key, c.sig, c.kor, c.pTrue); Add(n[kU1], key, c.sig, c.kor, c.pTrue); if(c.kN > 0)Add(n[kOK1], key, c.sig, c.kor, c.pTrue); if(inP)Add(n[kI1], key, c.sig, c.kor, c.pTrue); }
    }
    for(size_t j = 0; j < vp.size(); j++){
      const PairC& p = vp[j];
      Add(n[kP], key, p.sig, p.kor, p.pTrue);
      if(!pairInW0[j]){ Add(n[kU0], key, p.sig, p.kor, p.pTrue); Add(n[kPnotO0], key, p.sig, p.kor, p.pTrue); }
      if(!pairInW1[j])Add(n[kU1], key, p.sig, p.kor, p.pTrue);
    }
  }
  printf("Corte de parejas: %s\nEmparejamiento de candidatos: extremos a < %.1f cm o el mismo signal muon; k = %.3f\n", pairCut, dMatch, k);
  printf("Signal muons true: A %ld, B %ld. Candidatos de fondo emparejados entre las dos: %ld. Eventos con pareja fuera del microtree xs: %ld\n\n",
         nTrue[0], nTrue[1], nMatchBkg, pairEvNotInXS);
  printf("%-36s %7s |%6s %7s %7s %7s |%6s %7s %7s %7s\n", "selección", "cand.", "A", "eff", "pur", "error", "B", "eff", "pur", "error");
  for(int s = 0; s < kN; s++){
    printf("%-36s %7ld |", name[s], n[s].cand);
    for(int l = 0; l < 2; l++){
      const double S = n[s].sig[l], B = n[s].cand - S, eff = nTrue[l] ? 100.*n[s].mu[l].size()/nTrue[l] : 0., pur = n[s].cand ? 100.*S/n[s].cand : 0.;
      const double err = S > 0 ? 100.*std::sqrt(k*(S + B))/(k*S) : 0.;
      printf("%6.0f %6.2f%% %6.2f%% %6.1f%% |", S, eff, pur, err);
    }
    printf("\n");
  }
  // Categorías en lugar de cortes: la muestra se parte en subconjuntos excluyentes de distinta pureza y se combinan
  // (recuentos independientes, fondo conocido): 1/err^2 = sum_i k S_i^2/(S_i + B_i). Ganancia frente a la cadena entera sola.
  auto catErr = [&](const std::vector<std::pair<double, double>>& sc){   // (S, candidatos) por categoría
    double inv = 0; for(auto& x : sc)if(x.second > 0)inv += k*x.first*x.first/x.second;
    return inv > 0 ? 100./std::sqrt(inv) : 0.;
  };
  printf("\nCategorías (error relativo de A y B combinando subconjuntos excluyentes):\n");
  for(int l = 0; l < 2; l++){
    const double sO1 = n[kO1].sig[l], cO1 = n[kO1].cand, sK = n[kOK1].sig[l], cK = n[kOK1].cand, sI = n[kI1].sig[l], cI = n[kI1].cand;
    const double sPo = n[kU1].sig[l] - sO1, cPo = n[kU1].cand - cO1;   // parejas que no están en w1
    printf("  %s: w1 sola %.1f %% | w1 con K + sin K %.1f %% | w1 y parejas + w1 sin pareja %.1f %% | ... + parejas fuera de w1 %.1f %%\n",
           l == 0 ? "A" : "B", catErr({{sO1, cO1}}), catErr({{sK, cK}, {sO1 - sK, cO1 - cK}}), catErr({{sI, cI}, {sO1 - sI, cO1 - cI}}),
           catErr({{sI, cI}, {sO1 - sI, cO1 - cI}, {sPo, cPo}}));
  }
  return 0;
}
