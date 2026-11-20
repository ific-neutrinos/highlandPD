// La señal (mu+ de K+ en reposo) en los candidatos: por qué se pierde en el corte de chi2 mu y escala y
// resolución del rango. Árbol "cand" de mu_candidates_minitrees.C (ver mu_cand_common.h).
//
//   root -l -b -q 'test/phase2/mu_cand_signal.C+("mu_cand_big19.root")'
//
// 1. Candidatos de señal con los cortes del módulo, por clases de chi2 mu (extremo de parada = mayor
//    dE/dx): pureza, completitud, longitud reco / true, si la otra orientación da mejor chi2 mu.
// 2. Rango: p de rango, longitud reco - true (true = distancia inicio-fin del mu; el mu de K+ en reposo
//    tiene 235,5 MeV/c y un rango CSDA de ~54 cm), por población y por completitud.
#include "mu_cand_common.h"

namespace musig {
  using namespace mucut;

  double Q(std::vector<double> v, double p){
    if(v.empty())return -999;
    std::sort(v.begin(), v.end());
    return v[std::min<size_t>(v.size() - 1, size_t(p*v.size()))];
  }
  std::string QS(const std::vector<double>& v, int prec = 3){
    return Form("%.*g [%.*g, %.*g]", prec, Q(v, 0.5), prec, Q(v, 0.16), prec, Q(v, 0.84));
  }

  void Class(const std::vector<const C*>& v, const char* label){
    std::vector<double> pur, eff, ratio, dStop, dOther, nh, pr;
    int whole = 0, otherBetter = 0, otherGood = 0, inS = 0, noChi = 0;
    for(const C* c : v){
      pur.push_back(c->tPur); eff.push_back(c->tEff);
      if(c->tLength > 0)ratio.push_back(c->length/c->tLength);
      if(c->dEdxStop > 0)dStop.push_back(c->dEdxStop);
      if(c->dEdxStart > 0)dOther.push_back(c->dEdxStart);
      nh.push_back(c->nHits2); pr.push_back(c->pRange);
      if(c->tEff > 0.8 && c->tLength > 0 && std::abs(c->length/c->tLength - 1) < 0.15)whole++;
      if(c->chi2muOther >= 0 && (c->chi2muStop < 0 || c->chi2muOther < c->chi2muStop))otherBetter++;
      if(c->chi2muOther >= 0 && c->chi2muOther < 6)otherGood++;
      if(c->inSlice)inS++;
      if(c->chi2muStop < 0)noChi++;
    }
    const double n = v.size();
    std::cout << Form("-- %s: %.0f candidatos (en la beam slice %.0f %%; sin chi2 %.0f %%)", label, n, 100*inS/n, 100*noChi/n) << std::endl
              << "   pureza                     " << QS(pur) << std::endl
              << "   completitud                " << QS(eff) << std::endl
              << "   longitud reco / true       " << QS(ratio) << std::endl
              << "   p de rango (MeV/c)         " << QS(pr) << std::endl
              << "   dE/dx parada / inicio      " << QS(dStop) << "  /  " << QS(dOther) << std::endl
              << "   hits plano 2               " << QS(nh) << std::endl
              << Form("   traza = mu entero: %.0f %%; la otra orientación da mejor chi2 mu: %.0f %% (y < 6: %.0f %%)",
                      100*whole/n, 100*otherBetter/n, 100*otherGood/n) << std::endl;
  }
}

int mu_cand_signal(const char* in){
  using namespace musig;
  const std::vector<C> cands = Load(in);
  if(cands.empty())return 1;

  // 1. Pérdidas en chi2 mu
  std::vector<const C*> good, mid, bad, all;
  for(const C& c : cands){
    if(c.group != 0 || !InWindow(c))continue;
    all.push_back(&c);
    if(c.chi2muStop >= 0 && c.chi2muStop < 6)good.push_back(&c);
    else if(c.chi2muStop >= 6 && c.chi2muStop < 15)mid.push_back(&c);
    else bad.push_back(&c);
  }
  std::cout << std::endl << "== 1. Señal con los cortes del módulo, por chi2 mu ==" << std::endl;
  Class(all, "todos");
  Class(good, "chi2 mu < 6");
  Class(mid, "6 <= chi2 mu < 15");
  Class(bad, "chi2 mu >= 15 o sin chi2");
  // ¿y si se eligiera la orientación por el mejor chi2 mu?
  {
    int n6 = 0, n6best = 0;
    for(const C* c : all){
      if(c->chi2muStop >= 0 && c->chi2muStop < 6)n6++;
      const double best = std::min(c->chi2muStop < 0 ? 1e9 : c->chi2muStop, c->chi2muOther < 0 ? 1e9 : c->chi2muOther);
      if(best < 6)n6best++;
    }
    std::cout << Form("Con chi2 mu < 6: extremo de mayor dE/dx %d; mejor de las dos orientaciones %d (de %zu)", n6, n6best, all.size()) << std::endl;
    // lo mismo en el fondo (eventos normales), para ver el coste
    long b6 = 0, b6best = 0;
    for(const C& c : cands){
      if(c.group <= 0 || !InWindow(c))continue;
      if(c.chi2muStop >= 0 && c.chi2muStop < 6)b6++;
      const double best = std::min(c.chi2muStop < 0 ? 1e9 : c.chi2muStop, c.chi2muOther < 0 ? 1e9 : c.chi2muOther);
      if(best < 6)b6best++;
    }
    std::cout << Form("  fondo (sin escalar): %ld y %ld", b6, b6best) << std::endl;
  }

  // 2. Escala y resolución del rango (señal, p de rango en [100, 600], completos)
  std::cout << std::endl << "== 2. Rango de la señal (candidatos completos, p de rango en [100, 600]) ==" << std::endl;
  struct G { const char* label; std::function<bool(const C&)> sel; };
  const G groups[] = {
    {"todos", [](const C&){ return true; }},
    {"en la beam slice", [](const C& c){ return c.inSlice; }},
    {"fuera", [](const C& c){ return !c.inSlice; }},
    {"completitud > 0,9", [](const C& c){ return c.tEff > 0.9; }},
    {"completitud 0,7-0,9", [](const C& c){ return c.tEff > 0.7 && c.tEff <= 0.9; }},
    {"completitud < 0,7", [](const C& c){ return c.tEff <= 0.7; }},
    {"compl. > 0,9 y pureza > 0,9", [](const C& c){ return c.tEff > 0.9 && c.tPur > 0.9; }},
  };
  std::cout << Form("  %-30s %6s  %-26s %-26s %-26s", "", "N", "p de rango (MeV/c)", "reco - true (cm)", "true inicio-fin (cm)") << std::endl;
  for(const G& g : groups){
    std::vector<double> pr, d, tl;
    for(const C& c : cands){
      if(c.group != 0 || !c.complete || !g.sel(c))continue;
      pr.push_back(c.pRange);
      if(c.tLength > 0){ d.push_back(c.length - c.tLength); tl.push_back(c.tLength); }
    }
    std::cout << Form("  %-30s %6zu  %-26s %-26s %-26s", g.label, pr.size(), QS(pr).c_str(), QS(d).c_str(), QS(tl).c_str()) << std::endl;
  }
  return 0;
}
