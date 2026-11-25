// Compara la truth de dos truth skims de los mismos eventos (p. ej. 23480 frente a 23527: mismo MC, módulo distinto,
// la truth debería ser idéntica salvo el tiempo, que el 23480 no tenía). Por evento (run, subrun, event): número de
// partículas true, PDG de la partícula del haz, signal muons y su KOrigin.
//
//   root -l -b -q 'test/phase2/skim_compare.C+("truth_skim_v100s.root", "truth_skim_v100.root")'
#include "skim_common.h"
#include <map>
#include <tuple>

void skim_compare(const char* refFile, const char* newFile){
  using namespace skim;
  struct Sum { int nTrue = 0, beamPDG = 0, nSig = 0; long korSum = 0; };
  auto summarize = [](SkimReader& r, std::map<std::tuple<int,int,int>, Sum>& m){
    for(Long64_t e = 0; r.Next(e); e++){
      const Event& v = r.ev; Sum s;
      s.nTrue = v.t_id->size(); s.beamPDG = v.beamIdx >= 0 ? (*v.t_pdg)[v.beamIdx] : 0;
      for(size_t i = 0; i < v.t_id->size(); i++)if((*v.t_sigMu)[i] == 1){ s.nSig++; s.korSum += 10 + (*v.t_kOrigin)[i]; }
      m[std::make_tuple(v.run, v.subrun, v.event)] = s;
    }
  };
  std::map<std::tuple<int,int,int>, Sum> a, b;
  { SkimReader r(refFile); summarize(r, a); }
  { SkimReader r(newFile); summarize(r, b); }
  long common = 0, onlyA = 0, onlyB = 0, dTrue = 0, dBeam = 0, dSig = 0, sigA = 0, sigB = 0;
  for(auto& kv : a){
    auto it = b.find(kv.first);
    if(it == b.end()){ onlyA++; continue; }
    common++; sigA += kv.second.nSig; sigB += it->second.nSig;
    if(kv.second.nTrue != it->second.nTrue)dTrue++;
    if(kv.second.beamPDG != it->second.beamPDG)dBeam++;
    if(kv.second.nSig != it->second.nSig || kv.second.korSum != it->second.korSum)dSig++;
  }
  for(auto& kv : b)if(!a.count(kv.first))onlyB++;
  printf("Eventos: referencia %zu, nueva %zu, comunes %ld (solo en la referencia %ld, solo en la nueva %ld)\n", a.size(), b.size(), common, onlyA, onlyB);
  printf("En los comunes, eventos distintos: número de partículas true %ld, PDG del haz %ld, signal muons o su KOrigin %ld\n", dTrue, dBeam, dSig);
  printf("Signal muons en los comunes: referencia %ld, nueva %ld\n", sigA, sigB);
}
