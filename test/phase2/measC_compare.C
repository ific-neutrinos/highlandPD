// Medida C (27-nov-2026): parejas K -> mu con dos métodos, sobre los mismos eventos y con la misma truth.
//   (1) BDT + K: los mu que selecciona el BDT (dentro de la beam slice, ventana [205, 240] MeV/c, 10 variables) y, al
//       inicio del mu, el candidato a K (cand_k_*: partícula con un extremo a < 10 cm, la de mayor dE/dx). Pareja:
//       hay candidato a K; pareja con calorimetría: además es traza con dE/dx en ese extremo.
//   (2) Buscador de parejas de la tesis (secondaryKaonPairXS, forced daughters): eventos con accum_level >= 3 y, por
//       candidato a K, los cortes de pairCut (por defecto los de macros/tesis/xs/selection.C).
// Pareja correcta: (1) el mu es signal muon y el candidato a K es su madre true (cand_k_truek); (2) la hija es signal
// muon y el candidato a K es una K+ true. Solo eventos impares del microtree (la mitad par entrenó el BDT), x2.
// Eficiencia: signal muons true (árbol truth, eventos impares) con una pareja correcta, por capa. Solape por signal
// muon (run, subrun, evt, p true del mu a < 0,5 MeV/c). Proyección a los datos de 6 GeV: x 1/1,196.
//
//   root -l -b -q 'test/phase2/measC_compare.C+("xs_v100to105d01b.root", "opt_flat3_v100to105d01.root",
//       "opt_bdt_final/dataset/weights/bdt0_BDT.weights.xml", "pair_v100d01.root,pair_v101d01.root,...", "xs_v104d01b.root")'
// excludeXs: microtrees de los bloques sin parejas (pair_v104d01.root está vacío): sus eventos no cuentan en nada.
#include "TFile.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "TMVA/Reader.h"
#include <cmath>
#include <cstdio>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace mcc {
  using Key = std::tuple<int, int, int>;
  struct Mu { float p; int kor; int sel[4]; };   // sel: 0 BDT mu, 1 BDT pareja, 2 BDT pareja con calorimetría, 3 parejas tesis
  bool InLayer(int kor, int l){ return l == 2 ? true : l == 1 ? (kor >= 0 && kor <= 3) : (kor == 0 || kor == 1); }
  int Find(std::vector<Mu>& v, double p){ int b = -1; double d = 0.5; for(size_t i = 0; i < v.size(); i++)if(std::abs(v[i].p - p) < d){ d = std::abs(v[i].p - p); b = i; } return b; }
}

int measC_compare(const char* xs, const char* flat, const char* weights, const char* pairList, const char* excludeXs = "", double cut = 0.19,
                  const char* pairCut = "candidates_dau_mom_muon>0.16 && candidates_dau_mom_muon<0.265 && candidates_distance_dau<5 && (candidates_truncated_dedx+candidates_dau_truncated_dedx)/2<4"){
  using namespace mcc;
  const double kData = 2/1.196;   // mitad impar x2, y a datos
  TFile fx(xs); TTree* tru = (TTree*)fx.Get("truth"); TTree* ana = (TTree*)fx.Get("ana");
  // paridad de cada evento (por la entrada en el microtree)
  std::map<Key, int> odd;
  { Int_t run, sub, evt; ana->SetBranchStatus("*", 0); for(const char* b : {"run", "subrun", "evt"})ana->SetBranchStatus(b, 1);
    ana->SetBranchAddress("run", &run); ana->SetBranchAddress("subrun", &sub); ana->SetBranchAddress("evt", &evt);
    for(Long64_t e = 0; e < ana->GetEntries(); e++){ ana->GetEntry(e); odd[Key{run, sub, evt}] = e % 2; } }
  // bloques excluidos: sus eventos se marcan como no usables
  { std::string s = std::string(excludeXs) + ",";
    for(size_t a = 0, b; (b = s.find(',', a)) != std::string::npos; a = b + 1){
      if(b == a)continue;
      TFile fe(s.substr(a, b - a).c_str()); TTree* te = (TTree*)fe.Get("ana"); Int_t run, sub, evt; long n = 0;
      te->SetBranchStatus("*", 0); for(const char* c : {"run", "subrun", "evt"})te->SetBranchStatus(c, 1);
      te->SetBranchAddress("run", &run); te->SetBranchAddress("subrun", &sub); te->SetBranchAddress("evt", &evt);
      for(Long64_t e = 0; e < te->GetEntries(); e++){ te->GetEntry(e); auto it = odd.find(Key{run, sub, evt}); if(it != odd.end()){ it->second = -1; n++; } }
      printf("Excluidos %ld eventos de %s\n", n, s.substr(a, b - a).c_str());
    } }
  // signal muons true de los eventos impares
  std::map<Key, std::vector<Mu>> truth;
  { Int_t run, sub, evt, kor; Float_t p;
    tru->SetBranchAddress("run", &run); tru->SetBranchAddress("subrun", &sub); tru->SetBranchAddress("evt", &evt);
    tru->SetBranchAddress("truemu_korigin", &kor); tru->SetBranchAddress("truemu_mom", &p);
    for(Long64_t e = 0; e < tru->GetEntries(); e++){ tru->GetEntry(e); const Key k{run, sub, evt}; auto it = odd.find(k);
      if(it == odd.end() || it->second != 1)continue; truth[k].push_back({p, kor, {0, 0, 0, 0}}); } }
  long nTrue[3] = {0}; for(auto& kv : truth)for(auto& m : kv.second)for(int l = 0; l < 3; l++)if(InLayer(m.kor, l))nTrue[l]++;
  // decisión del BDT desde el árbol plano (como measC_pairs.C)
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
    } }
  // (1) BDT + K
  double b1[3] = {0}, b1ok[3] = {0};   // [0 mu, 1 pareja, 2 pareja con calorimetría]: candidatos y correctos
  { const int N = 1000;
    if(ana->GetMaximum("ncand") > N){ printf("ERROR: ncand > %d\n", N); return 1; }
    Int_t run, sub, evt, nc, pop[N], sig[N], kn[N], kist[N], ktrue[N]; Float_t pr[N], kd[N], tp[N];
    ana->SetBranchStatus("*", 0);
    for(const char* b : {"run", "subrun", "evt", "ncand", "cand_pop", "cand_prange", "cand_truesignalmuon", "cand_truemom", "cand_k_n", "cand_k_istrack", "cand_k_dedx", "cand_k_truek"})ana->SetBranchStatus(b, 1);
    ana->SetBranchAddress("run", &run); ana->SetBranchAddress("subrun", &sub); ana->SetBranchAddress("evt", &evt); ana->SetBranchAddress("ncand", &nc);
    ana->SetBranchAddress("cand_pop", pop); ana->SetBranchAddress("cand_prange", pr); ana->SetBranchAddress("cand_truesignalmuon", sig); ana->SetBranchAddress("cand_truemom", tp);
    ana->SetBranchAddress("cand_k_n", kn); ana->SetBranchAddress("cand_k_istrack", kist); ana->SetBranchAddress("cand_k_dedx", kd); ana->SetBranchAddress("cand_k_truek", ktrue);
    for(Long64_t e = 1; e < ana->GetEntries(); e += 2){
      ana->GetEntry(e);
      if(odd[Key{run, sub, evt}] != 1)continue;
      for(int j = 0; j < nc; j++){
        if(pop[j] != 0 || !bdtPass.count({e, pr[j]}))continue;
        const bool lv[3] = {true, kn[j] > 0, kn[j] > 0 && kist[j] == 1 && kd[j] > 0};
        const bool isSig = sig[j] == 1;
        Mu* m = nullptr;
        if(isSig){ auto it = truth.find(Key{run, sub, evt}); if(it != truth.end()){ const int i = Find(it->second, tp[j]); if(i >= 0)m = &it->second[i]; } }
        for(int w = 0; w < 3; w++){
          if(!lv[w])continue;
          b1[w]++;
          const bool ok = w == 0 ? isSig : isSig && ktrue[j] == 1;
          if(ok){ b1ok[w]++; if(m)m->sel[w] = 1; }
        }
      }
    } }
  // (2) parejas de la tesis
  double p2 = 0, p2ok = 0, p2f = 0, p2fok = 0; long p2noTruth = 0;
  { std::string s = std::string(pairList) + ",";
    for(size_t a = 0, b; (b = s.find(',', a)) != std::string::npos; a = b + 1){
      if(b == a)continue;
      TFile fp(s.substr(a, b - a).c_str()); TTree* ap = fp.IsZombie() ? nullptr : (TTree*)fp.Get("ana");
      if(!ap || !ap->GetBranch("candidates_dau_forced")){ printf("AVISO: %s sin árbol de parejas: se salta\n", s.substr(a, b - a).c_str()); continue; }
      TTreeFormula fCut("cut", pairCut, ap), fP("p", "candidates_dau_truemom", ap), fSig("s", "candidatedau_truesignalmuon", ap),
                   fK("k", "candidates_truepdg", ap), fF("f", "candidates_dau_forced", ap);
      Int_t run, sub, evt, acc[1], nc;
      ap->SetBranchAddress("run", &run); ap->SetBranchAddress("subrun", &sub); ap->SetBranchAddress("evt", &evt);
      ap->SetBranchAddress("accum_level", acc); ap->SetBranchAddress("ncandidates", &nc);
      for(Long64_t e = 0; e < ap->GetEntries(); e++){
        ap->LoadTree(e); ap->GetEntry(e);
        if(acc[0] < 3 || nc <= 0)continue;
        const Key k{run, sub, evt}; auto io = odd.find(k); if(io == odd.end() || io->second != 1)continue;
        fCut.GetNdata(); fP.GetNdata(); fSig.GetNdata(); fK.GetNdata(); fF.GetNdata();
        for(int j = 0; j < nc; j++){
          if(fCut.EvalInstance(j) == 0)continue;
          const bool forced = fF.EvalInstance(j) == 1;
          p2++; if(forced)p2f++;
          if((int)fSig.EvalInstance(j) != 1 || (int)fK.EvalInstance(j) != 321)continue;
          p2ok++; if(forced)p2fok++;
          auto it = truth.find(k); const int i = it == truth.end() ? -1 : Find(it->second, 1000.*fP.EvalInstance(j));
          if(i < 0){ p2noTruth++; continue; }
          it->second[i].sel[3] = 1;
        }
      }
    } }
  // resultados
  long eff[5][3] = {{0}}, both[3] = {0}, onlyB[3] = {0}, onlyP[3] = {0};
  for(auto& kv : truth)for(auto& m : kv.second)for(int l = 0; l < 3; l++){
    if(!InLayer(m.kor, l))continue;
    for(int w = 0; w < 4; w++)if(m.sel[w])eff[w][l]++;
    if(m.sel[1] && m.sel[3])both[l]++; else if(m.sel[1])onlyB[l]++; else if(m.sel[3])onlyP[l]++;
  }
  printf("Signal muons true en los eventos impares: A %ld, B %ld, C %ld\n", nTrue[0], nTrue[1], nTrue[2]);
  printf("Corte de parejas (tesis): %s\n\n", pairCut);
  printf("%-38s %8s %8s %7s | eficiencia A / B / C (pareja correcta)   | en datos de 6 GeV: total, correctas\n", "método", "cand.", "correct.", "pureza");
  auto row = [&](const char* n, double c, double ok, int w){
    printf("%-38s %8.0f %8.0f %6.1f%% | %5.2f%% %5.2f%% %5.2f%%                     | %6.0f %6.0f\n", n, c, ok, c > 0 ? 100*ok/c : 0.,
           100.*eff[w][0]/nTrue[0], 100.*eff[w][1]/nTrue[1], 100.*eff[w][2]/nTrue[2], c*kData, ok*kData); };
  row("BDT: mu (sin K)", b1[0], b1ok[0], 0);
  row("BDT + K al inicio", b1[1], b1ok[1], 1);
  row("BDT + K traza con dE/dx", b1[2], b1ok[2], 2);
  row("parejas de la tesis (forced daughters)", p2, p2ok, 3);
  printf("  de ellas con forced daughter: %.0f (correctas %.0f); parejas correctas sin su mu en el árbol truth: %ld\n", p2f, p2fok, p2noTruth);
  printf("\nSolape de parejas correctas por signal muon (BDT + K al inicio frente a parejas de la tesis):\n");
  for(int l = 0; l < 3; l++)printf("  %s: las dos %ld, solo BDT %ld, solo tesis %ld; unión %ld (%.2f %% de los true; datos %.0f)\n", l == 0 ? "A" : l == 1 ? "B" : "C",
                                   both[l], onlyB[l], onlyP[l], both[l] + onlyB[l] + onlyP[l], 100.*(both[l] + onlyB[l] + onlyP[l])/nTrue[l], (both[l] + onlyB[l] + onlyP[l])*kData);
  return 0;
}
