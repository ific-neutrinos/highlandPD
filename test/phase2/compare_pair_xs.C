// Comparación de la selección de parejas K -> mu de la tesis (secondaryKaonPairXS + cortes de macro) con
// secondaryKaonProductionXS, con la misma definición de señal (signal muon, capas A en B en C) y sobre los mismos
// minitrees. Denominador común: los signal muons del árbol truth de secondaryKaonProductionXS.
//
//   root -l -b -q 'test/phase2/compare_pair_xs.C+("xs_v100d01.root", "pair_v100d01.root")'
//   root -l -b -q 'test/phase2/compare_pair_xs.C+("xs.root", "pair.root", "<cortes de la tesis sobre los candidatos>")'
//
// Selección de parejas: eventos con accum_level[0] >= 3 (traza del haz, cualquier especie, hay candidato a K; forced
// daughter buscada) y, por candidato a K, los cortes de pairCut (por defecto los de macros/tesis/xs/selection.C:
// p de rango del mu en [0.16, 0.265] GeV/c, distancia K-hija < 5 cm, dE/dx truncado medio de K e hija < 4 MeV/cm).
// Nuestra selección, en dos niveles: hasta Michel (nivel del ajuste, sin ventana de p) y con la ventana de p.
//
// Eficiencia: signal muons true de la capa seleccionados (en la de parejas, hija de un candidato que pasa) sobre todos.
// Pureza: candidatos seleccionados cuya traza (la hija, en la de parejas) es signal muon de la capa, sobre todos.
// Un signal muon se identifica por (run, subrun, evt, p true del mu): |dp| < 0.5 MeV/c. Imprime también el solape
// entre las dos selecciones y el reparto por especie true del haz.
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

namespace cpx {
  using Key = std::tuple<int, int, int>;   // run, subrun, evt
  struct Mu { double p; int kor; int beam; };
  // capa: 0 A (KOrigin 0-1), 1 B (0-3), 2 C (cualquiera); acumuladas
  bool InLayer(int kor, int l){ return l == 2 ? true : l == 1 ? (kor >= 0 && kor <= 3) : (kor == 0 || kor == 1); }
  const char* kLayer[3] = {"A", "B", "C"};
  int Species(int pdg){ const int a = std::abs(pdg); return a == 211 ? 0 : a == 2212 ? 1 : a == 321 ? 2 : (a == 13 || a == 11) ? 3 : 4; }
  const char* kSpecies[5] = {"pi", "p", "K", "mu/e", "otra"};
  // índice del signal muon true del evento con p (MeV/c) a menos de 0.5; -1 si no hay
  int Find(const std::vector<Mu>& v, double p){
    int best = -1; double d = 0.5;
    for(size_t i = 0; i < v.size(); i++)if(std::abs(v[i].p - p) < d){ d = std::abs(v[i].p - p); best = i; }
    return best;
  }
}

int compare_pair_xs(const char* xsFile, const char* pairFile,
                    const char* pairCut = "candidates_dau_mom_muon>0.16 && candidates_dau_mom_muon<0.265 && candidates_distance_dau<5 && (candidates_truncated_dedx+candidates_dau_truncated_dedx)/2<4"){
  using namespace cpx;
  TFile* fx = TFile::Open(xsFile); TFile* fp = TFile::Open(pairFile);
  if(!fx || fx->IsZombie() || !fp || fp->IsZombie())return -1;
  TTree* tru = (TTree*)fx->Get("truth"); TTree* ax = (TTree*)fx->Get("ana"); TTree* ap = (TTree*)fp->Get("ana");
  if(!tru || !ax || !ap)return -1;
  // niveles de nuestra cadena (secondaryKaonProductionXSSelection): dentro 7 cortes y fuera 9 más el de población; el último es la ventana de p
  const int nCut[2] = {8, 10};   // cortes de cada cadena contando el de población (cand_level de la cadena entera; technote_cutflow)

  // ---- signal muons true (denominador) y los que selecciona nuestra cadena
  std::map<Key, std::vector<Mu>> truth;
  std::map<Key, std::vector<int>> oursSel[2];   // [0 hasta Michel, 1 con ventana de p]: 1 si seleccionado
  {
    Int_t run, sub, evt, kor, acc[2], bpdg; Float_t p;
    tru->SetBranchAddress("run", &run); tru->SetBranchAddress("subrun", &sub); tru->SetBranchAddress("evt", &evt);
    tru->SetBranchAddress("truemu_korigin", &kor); tru->SetBranchAddress("truemu_mom", &p); tru->SetBranchAddress("accum_level", acc);
    tru->SetBranchAddress("beam_truepdg", &bpdg);
    for(Long64_t e = 0; e < tru->GetEntries(); e++){
      tru->GetEntry(e);
      const Key k{run, sub, evt};
      truth[k].push_back({p, kor, Species(bpdg)});
      for(int w = 0; w < 2; w++){
        const int l0 = nCut[0] - 1 + w, l1 = nCut[1] - 1 + w;   // cand_level; accum_level = 1 + cand_level
        oursSel[w][k].push_back(acc[0] >= 1 + l0 || acc[1] >= 1 + l1);
      }
    }
  }
  long nTot[3] = {0, 0, 0}, nTotS[5][3] = {{0}};
  for(auto& kv : truth)for(auto& m : kv.second)for(int l = 0; l < 3; l++)if(InLayer(m.kor, l)){ nTot[l]++; nTotS[m.beam][l]++; }

  // ---- nuestra selección: candidatos (pureza)
  long oCand[2] = {0, 0}, oSig[2][3] = {{0}}, oCandS[2][5] = {{0}}, oSigBS[2][5] = {{0}};
  {
    const int N = 100;
    Int_t nc, pop[N], lev[N], sig[N], kor[N], bpdg;
    ax->SetBranchStatus("*", 0);
    for(const char* b : {"ncand", "cand_pop", "cand_level", "cand_truesignalmuon", "cand_truekorigin", "beam_truepdg"})ax->SetBranchStatus(b, 1);
    ax->SetBranchAddress("ncand", &nc); ax->SetBranchAddress("cand_pop", pop); ax->SetBranchAddress("cand_level", lev);
    ax->SetBranchAddress("cand_truesignalmuon", sig); ax->SetBranchAddress("cand_truekorigin", kor); ax->SetBranchAddress("beam_truepdg", &bpdg);
    for(Long64_t e = 0; e < ax->GetEntries(); e++){
      ax->GetEntry(e);
      for(int j = 0; j < nc && j < N; j++){
        if(pop[j] < 0 || pop[j] > 1)continue;
        for(int w = 0; w < 2; w++){
          if(lev[j] < nCut[pop[j]] - 1 + w)continue;
          oCand[w]++; oCandS[w][Species(bpdg)]++;
          if(sig[j] == 1){ for(int l = 0; l < 3; l++)if(InLayer(kor[j], l))oSig[w][l]++; if(InLayer(kor[j], 1))oSigBS[w][Species(bpdg)]++; }
        }
      }
    }
  }

  // ---- selección de parejas
  long pCand = 0, pSig[3] = {0, 0, 0}, pCandS[5] = {0}, pSigBS[5] = {0}, pNoTruth = 0;
  std::map<Key, std::vector<int>> pairSel;   // por signal muon true: 1 si alguna pareja lo selecciona
  for(auto& kv : truth)pairSel[kv.first].assign(kv.second.size(), 0);
  {
    std::unique_ptr<TTreeFormula> fCut(new TTreeFormula("cut", pairCut, ap));
    std::unique_ptr<TTreeFormula> fP(new TTreeFormula("p", "candidates_dau_truemom", ap));
    std::unique_ptr<TTreeFormula> fSig(new TTreeFormula("s", "candidatedau_truesignalmuon", ap));
    std::unique_ptr<TTreeFormula> fKor(new TTreeFormula("k", "candidatedau_truekorigin", ap));
    if(fCut->GetNdim() == 0){ printf("Corte de parejas no válido: %s\n", pairCut); return -1; }
    Int_t run, sub, evt, acc[1], bpdg, nc;
    ap->SetBranchAddress("run", &run); ap->SetBranchAddress("subrun", &sub); ap->SetBranchAddress("evt", &evt);
    ap->SetBranchAddress("accum_level", acc); ap->SetBranchAddress("beam_truepdg", &bpdg); ap->SetBranchAddress("ncandidates", &nc);
    for(Long64_t e = 0; e < ap->GetEntries(); e++){
      ap->LoadTree(e); ap->GetEntry(e);
      if(acc[0] < 3 || nc <= 0)continue;
      fCut->GetNdata(); fP->GetNdata(); fSig->GetNdata(); fKor->GetNdata();
      const Key k{run, sub, evt};
      for(int j = 0; j < nc; j++){
        if(fCut->EvalInstance(j) == 0)continue;
        pCand++; pCandS[Species(bpdg)]++;
        if((int)fSig->EvalInstance(j) != 1)continue;
        const int kor = (int)fKor->EvalInstance(j);
        for(int l = 0; l < 3; l++)if(InLayer(kor, l))pSig[l]++;
        if(InLayer(kor, 1))pSigBS[Species(bpdg)]++;
        const double pMeV = 1000.*fP->EvalInstance(j);   // kaonTree: GeV/c
        auto it = truth.find(k); const int i = it == truth.end() ? -1 : Find(it->second, pMeV);
        if(i < 0){ pNoTruth++; continue; }
        pairSel[k][i] = 1;
      }
    }
  }

  // ---- eficiencias y solape
  long oEff[2][3] = {{0}}, pEff[3] = {0}, both[3] = {0}, onlyO[3] = {0}, onlyP[3] = {0}, oEffS[2][5] = {{0}}, pEffS[5] = {0};
  for(auto& kv : truth){
    for(size_t i = 0; i < kv.second.size(); i++){
      const Mu& m = kv.second[i];
      const bool p = pairSel[kv.first][i], o = oursSel[0][kv.first][i];
      for(int l = 0; l < 3; l++){
        if(!InLayer(m.kor, l))continue;
        for(int w = 0; w < 2; w++)if(oursSel[w][kv.first][i])oEff[w][l]++;
        if(p)pEff[l]++;
        if(p && o)both[l]++; else if(o)onlyO[l]++; else if(p)onlyP[l]++;
      }
      if(InLayer(m.kor, 1)){ for(int w = 0; w < 2; w++)if(oursSel[w][kv.first][i])oEffS[w][m.beam]++; if(p)pEffS[m.beam]++; }
    }
  }
  auto pc = [](double a, double b){ return b > 0 ? 100.*a/b : 0.; };
  printf("Corte de parejas: %s\n", pairCut);
  printf("Signal muons true: A %ld, B %ld, C %ld\n\n", nTot[0], nTot[1], nTot[2]);
  printf("%-30s %8s", "selección", "cand.");
  for(int l = 0; l < 3; l++)printf("   eff %s   pur %s  e x p %s", kLayer[l], kLayer[l], kLayer[l]); printf("\n");
  const char* oName[2] = {"nuestra, hasta Michel", "nuestra, con ventana de p"};
  for(int w = 0; w < 2; w++){
    printf("%-30s %8ld", oName[w], oCand[w]);
    for(int l = 0; l < 3; l++)printf("  %6.2f%%  %6.2f%%  %6.3f%%", pc(oEff[w][l], nTot[l]), pc(oSig[w][l], oCand[w]), pc(oEff[w][l], nTot[l])*pc(oSig[w][l], oCand[w])/100.);
    printf("\n");
  }
  printf("%-30s %8ld", "parejas (tesis)", pCand);
  for(int l = 0; l < 3; l++)printf("  %6.2f%%  %6.2f%%  %6.3f%%", pc(pEff[l], nTot[l]), pc(pSig[l], pCand), pc(pEff[l], nTot[l])*pc(pSig[l], pCand)/100.);
  printf("\n  (parejas con hija signal muon sin su mu en el árbol truth: %ld)\n", pNoTruth);
  printf("\nSolape de signal muons seleccionados (nuestra hasta Michel frente a parejas):\n");
  for(int l = 0; l < 3; l++)printf("  %s: los dos %ld, solo nuestra %ld, solo parejas %ld (unión %ld = %.2f%% de los true)\n", kLayer[l], both[l], onlyO[l], onlyP[l],
                                   both[l] + onlyO[l] + onlyP[l], pc(both[l] + onlyO[l] + onlyP[l], nTot[l]));
  printf("\nPor especie true del haz (señal B): true / eff nuestra hasta Michel, con p / eff parejas | pureza B nuestra hasta Michel, con p / parejas\n");
  for(int s = 0; s < 5; s++){
    if(!nTotS[s][1] && !pCandS[s] && !oCandS[0][s])continue;
    printf("  %-5s %6ld / %6.2f%% %6.2f%% / %6.2f%% | %6.2f%% %6.2f%% / %6.2f%%\n", kSpecies[s], nTotS[s][1], pc(oEffS[0][s], nTotS[s][1]), pc(oEffS[1][s], nTotS[s][1]),
           pc(pEffS[s], nTotS[s][1]), pc(oSigBS[0][s], oCandS[0][s]), pc(oSigBS[1][s], oCandS[1][s]), pc(pSigBS[s], pCandS[s]));
  }
  return 0;
}
