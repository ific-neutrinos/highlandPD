// Composición de los signal muons (t_sigMu: mu+ de K+ -> mu+ nu de una K+ en reposo) por especie de la partícula true
// del haz que dispara, en el truth skim. Rehace con la definición de señal actual lo que daba kaon_origin_minitrees
// (ventana de p del mu, definición anterior): nota técnica, cap. de truth, composición (items 2 y 3).
//
//   root -l -b -q 'test/phase2/skim_composition.C+("truth_skim_v100s.root")'
//
// Por especie del haz (t_pdg de beamIdx): signal muons por KOrigin (0 1a interacción, 1 con K+ -> K+, 2 secundaria,
// 3 K del haz, 4 cadena cortada, 6 otra raíz, 5 madre no guardada). En las secundarias (2), PDG de la madre de la
// primera K de la cadena (subiendo las K+ -> K+) y su momento al final (la interacción que crea la K, si la K nace en
// su final). En A (0-1), momento final de la partícula del haz.
#include "skim_common.h"
#include <algorithm>
#include <map>
#include <string>

namespace skco {
  const char* Name(int pdg){
    switch(pdg){
      case 211: return "pi+"; case -211: return "pi-"; case 321: return "K+"; case -321: return "K-"; case 2212: return "p";
      case 2112: return "n"; case 130: return "K0L"; case 310: return "K0S"; case -11: return "e+"; case 11: return "e-";
      case -13: return "mu+"; case 13: return "mu-"; default: return "otro";
    }
  }
  double Median(std::vector<double> v){ if(v.empty())return -1; std::sort(v.begin(), v.end()); return v[v.size()/2]; }
}

void skim_composition(const char* input){
  using namespace skim;
  SkimReader r(input);
  std::map<int, std::map<int, long>> n;            // [PDG del haz][KOrigin]
  std::map<int, std::map<int, long>> secMother;    // [PDG del haz][PDG de la madre de la primera K]
  std::map<int, std::vector<double>> pBeamA;       // [PDG del haz]: momento final del haz en A
  std::map<int, std::map<int, std::vector<double>>> pMother;   // [PDG del haz][PDG madre]: momento final de la madre
  long nSig = 0, nNoBeam = 0, secMid = 0, secEnd = 0;
  for(Long64_t e = 0; r.Next(e); e++){
    const Event& ev = r.ev;
    for(size_t i = 0; i < ev.t_id->size(); i++){
      if((*ev.t_sigMu)[i] != 1)continue;
      nSig++;
      if(ev.beamIdx < 0){ nNoBeam++; continue; }
      const int bpdg = (*ev.t_pdg)[ev.beamIdx];
      const int ko = (*ev.t_kOrigin)[i];
      n[bpdg][ko]++;
      if(ko == 0 || ko == 1)pBeamA[bpdg].push_back((*ev.t_pEnd)[ev.beamIdx]);
      if(ko != 2)continue;
      // primera K de la cadena: desde la K madre del mu, se suben las K+ -> K+
      int k = r.Parent(i);
      while(k >= 0){ const int m = r.Parent(k); if(m < 0 || (*ev.t_pdg)[m] != 321)break; k = m; }
      const int m = k >= 0 ? r.Parent(k) : -1;
      if(m < 0)continue;
      const int mpdg = (*ev.t_pdg)[m];
      secMother[bpdg][mpdg]++;
      const double d = Dist((*ev.t_x)[k], (*ev.t_y)[k], (*ev.t_z)[k], (*ev.t_xe)[m], (*ev.t_ye)[m], (*ev.t_ze)[m]);
      if(d < 1.){ secEnd++; pMother[bpdg][mpdg].push_back((*ev.t_pEnd)[m]); } else secMid++;
    }
  }
  printf("Signal muons: %ld (sin partícula true del haz: %ld)\n", nSig, nNoBeam);
  printf("Por especie del haz que dispara; KOrigin: 0 1a interacción, 1 con K+ -> K+, 2 secundaria, 3 K del haz, 4 cadena cortada, 5 madre no guardada, 6 otra raíz\n");
  std::vector<std::pair<long, int>> order;
  for(auto& b : n){ long t = 0; for(auto& o : b.second)t += o.second; order.push_back({t, b.first}); }
  std::sort(order.rbegin(), order.rend());
  for(auto& ob : order){
    const int b = ob.second; auto& nb = n[b];
    const long chain = nb[0] + nb[1] + nb[2] + nb[3];
    printf("\n== Haz %s (%d): %ld signal muons; con la cadena hasta el haz que dispara (0-3) %ld\n", skco::Name(b), b, ob.first, chain);
    for(int o : {0, 1, 2, 3, 4, 5, 6})if(nb.count(o)){
      printf("    KOrigin %d: %6ld (%5.1f %% del total", o, nb[o], 100.*nb[o]/ob.first);
      if(o <= 3 && chain)printf(", %5.1f %% de 0-3", 100.*nb[o]/chain);
      printf(")\n");
    }
    if(chain)printf("    A (0-1) %ld (%.1f %%), secundaria %ld (%.1f %%), K del haz %ld (%.1f %%); p final del haz en A: mediana %.2f GeV/c\n",
                    nb[0] + nb[1], 100.*(nb[0] + nb[1])/chain, nb[2], 100.*nb[2]/chain, nb[3], 100.*nb[3]/chain, skco::Median(pBeamA[b]));
    if(nb[2]){
      std::vector<std::pair<long, int>> mo; for(auto& kv : secMother[b])mo.push_back({kv.second, kv.first});
      std::sort(mo.rbegin(), mo.rend());
      printf("    madre de la primera K (secundarias), y su momento final (GeV/c: < 0,5 / 0,5-1 / 1-2 / 2-4 / > 4):\n");
      for(auto& x : mo){
        int bin[5] = {0, 0, 0, 0, 0};
        for(double p : pMother[b][x.second])bin[p < 0.5 ? 0 : p < 1 ? 1 : p < 2 ? 2 : p < 4 ? 3 : 4]++;
        printf("      %-4s %5ld   %4d %4d %4d %4d %4d\n", skco::Name(x.second), x.first, bin[0], bin[1], bin[2], bin[3], bin[4]);
      }
    }
  }
  printf("\nSecundarias: la K nace en el final de la madre %ld, a mitad de la madre %ld (sin momento)\n", secEnd, secMid);
}
