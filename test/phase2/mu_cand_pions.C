// Qué son los pi+- que quedan tras las cadenas de cortes de los candidatos a mu de K+ en reposo
// (InSliceSteps y OutSliceSteps de mu_cand_common.h), comparados con la señal.
//
//   root -l -b -q 'test/phase2/mu_cand_pions.C+("mu_cand_big19.root")'
//
// Por población: carga, PDG de la madre, proceso de creación, si descienden del haz; pureza y completitud
// del emparejamiento; longitud reco / distancia true inicio-fin (¿la traza es el pion entero?); distancia
// true inicio-fin. Sin escalar (fondo de eventos normales).
#include "mu_cand_common.h"
#include <map>

namespace mupi {
  using namespace mucut;

  const char* ProcName(int p){
    static const char* n[] = {"Unknown", "primary", "Decay", "kaonplusInelastic", "kaonminusInelastic", "neutronInelastic", "hadElastic",
                              "nCapture", "protonInelastic", "piplusInelastic", "piminusInelastic", "hBertiniCaptureAtRest", "CoulombScat",
                              "kaon0LInelastic", "electronNuclear", "muMinusCaptureAtRest", "dInelastic", "kaon0SInelastic", "positronNuclear",
                              "lambdaInelastic", "tInelastic"};
    return p >= 0 && p <= 20 ? n[p] : "?";
  }

  void Counts(const std::vector<const C*>& v, const std::function<std::string(const C&)>& key, const char* title){
    std::map<std::string, int> m;
    for(const C* c : v)m[key(*c)]++;
    std::vector<std::pair<int, std::string>> s; for(auto& kv : m)s.push_back({kv.second, kv.first});
    std::sort(s.rbegin(), s.rend());
    std::cout << "  " << title << ":";
    for(size_t i = 0; i < std::min<size_t>(s.size(), 6); i++)std::cout << Form("  %s %.0f %%", s[i].second.c_str(), 100.*s[i].first/v.size());
    std::cout << std::endl;
  }

  void Quant(const std::vector<const C*>& v, const std::function<double(const C&)>& f, const char* title){
    std::vector<double> x;
    for(const C* c : v){ const double y = f(*c); if(y > -900)x.push_back(y); }
    if(x.empty())return;
    std::sort(x.begin(), x.end());
    auto q = [&](double p){ return x[std::min<size_t>(x.size() - 1, size_t(p*x.size()))]; };
    std::cout << Form("  %-38s %7.3g [%7.3g, %7.3g]", title, q(0.5), q(0.16), q(0.84)) << std::endl;
  }

  void Describe(const std::vector<const C*>& v, const char* label){
    std::cout << std::endl << "-- " << label << ": " << v.size() << " candidatos" << std::endl;
    if(v.empty())return;
    Counts(v, [](const C& c){ return std::string(Form("%d", c.tPDG)); }, "PDG");
    Counts(v, [](const C& c){ return std::string(Form("%d", c.tParentPDG)); }, "PDG madre");
    Counts(v, [](const C& c){ return std::string(ProcName(c.tProcStart)); }, "proceso de creación");
    Counts(v, [](const C& c){ return std::string(c.tFromBeam ? "del haz" : (c.tOrigin == 2 ? "cósmico" : "otro")); }, "origen");
    Quant(v, [](const C& c){ return c.tPur; }, "pureza (mediana [16, 84 %])");
    Quant(v, [](const C& c){ return c.tEff; }, "completitud");
    Quant(v, [](const C& c){ return c.tLength > 0 ? c.length/c.tLength : -999; }, "longitud reco / true inicio-fin");
    Quant(v, [](const C& c){ return c.tLength; }, "distancia true inicio-fin (cm)");
    Quant(v, [](const C& c){ return c.length; }, "longitud reco (cm)");
    int whole = 0;
    for(const C* c : v)if(c->tEff > 0.8 && c->tLength > 0 && std::abs(c->length/c->tLength - 1) < 0.15)whole++;
    std::cout << Form("  traza = partícula entera (completitud > 0,8 y longitud a < 15 %%): %.0f %%", 100.*whole/v.size()) << std::endl;
  }
}

int mu_cand_pions(const char* in){
  using namespace mupi;
  const std::vector<C> cands = Load(in);
  if(cands.empty())return 1;
  struct Pop { const char* label; std::vector<Step> st; std::function<bool(const C&)> sel; };
  const Pop pops[] = {{"en la beam slice", InSliceSteps(), [](const C& c){ return c.inSlice; }},
                      {"fuera de la beam slice", OutSliceSteps(), [](const C& c){ return !c.inSlice; }}};
  for(const Pop& p : pops){
    std::vector<const C*> sig, pi, piPlus, piMinus;
    for(const C& c : cands){
      if(!p.sel(c) || !InWindow(c) || !Pass(p.st, c))continue;
      if(c.group == 0)sig.push_back(&c);
      if(c.group == 3){ pi.push_back(&c); (c.tPDG > 0 ? piPlus : piMinus).push_back(&c); }
    }
    std::cout << std::endl << "==== " << p.label << " (todos los cortes) ====" << std::endl;
    Describe(sig, "señal");
    Describe(pi, "pi+-");
    Describe(piPlus, "pi+");
    Describe(piMinus, "pi-");
  }
  return 0;
}
