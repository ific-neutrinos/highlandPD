// Cortes en cadena sobre el árbol "cand" de mu_candidates_minitrees.C (candidatos a mu de K+ en reposo).
//
//   root -l -b -q 'test/phase2/mu_cand_cutflow.C+("mu_cand_big19.root", "cutflow.root")'
//
// Señal: tSignal con la K en reposo (tKMomEnd < 50 MeV/c), en todos los eventos. Fondo: todo lo que no es
// señal, solo en eventos normales (sin volcado), escalado al número total de eventos. "Extremo de parada"
// = el de mayor dE/dx (stopAtB), como en mu_candidates_minitrees.C.
//
// Una sola pasada por el árbol (mu_cand_common.h): se guardan en memoria los candidatos que pasan los
// cortes del módulo sin la ventana y con p de rango en [100, 600]; todo lo demás se calcula sobre ellos.
// Los cortes están en mucut::DefaultSteps().
//
// Imprime: corte a corte (acumulado) por grupo, S y B por evento, S/B y S/sqrt(S+B); N-1 de cada corte;
// cuantiles de las variables tras todos los cortes; composición del fondo que queda; espectro de p de
// rango con los cortes salvo la ventana (solo candidatos completos: los de fuera de la beam slice fuera de
// la ventana van en resumen, sin calorimetría). Guarda los histogramas.
#include "mu_cand_common.h"
#include "TH1D.h"
#include <map>

namespace mucut {

  struct Tally { double n[6] = {0}; };

  void Header(){
    std::cout << Form("  %-28s %7s", "", kGroupName[0]);
    for(int g = 1; g < 6; g++)std::cout << Form(" %10s", kGroupName[g]);
    std::cout << Form("  %8s %8s  %8s %8s", "S/ev", "B/ev", "S/B", "S/√(S+B)") << std::endl;
  }

  void Row(const std::string& name, const Tally& t, double scaleB, double nEv){
    double B = 0; for(int g = 1; g < 6; g++)B += t.n[g];
    const double S = t.n[0], Bs = B*scaleB;
    std::cout << Form("  %-28s %7.0f", name.c_str(), S);
    for(int g = 1; g < 6; g++)std::cout << Form(" %10.0f", t.n[g]);
    std::cout << Form("  %8.4f %8.4f  %8.4f %8.2f", S/nEv, Bs/nEv, Bs > 0 ? S/Bs : 0., S + Bs > 0 ? S/std::sqrt(S + Bs) : 0.) << std::endl;
  }

  double Q(std::vector<double>& v, double q){ return v[std::min<size_t>(v.size() - 1, size_t(q*v.size()))]; }
}

int mu_cand_cutflow(const char* in, const char* out = "mu_cand_cutflow.root", bool orientByChi2 = false, long nEvents = 73070, long nEventsNormal = 65977){
  using namespace mucut;
  const std::vector<C> cands = Load(in, orientByChi2);
  if(orientByChi2)std::cout << "Extremo de parada: orientación con mejor chi2 mu" << std::endl;
  if(cands.empty())return 1;

  const double scaleB = double(nEvents)/nEventsNormal;
  const std::vector<Step> steps = DefaultSteps();
  const int nS = steps.size();
  std::cout << "Eventos " << nEvents << " (normales " << nEventsNormal << "); fondo escalado x" << scaleB << std::endl;

  // Cortes acumulados y N-1, en la ventana del módulo: cadena común antigua y una por población
  auto Tables = [&](const std::vector<Step>& st, const std::function<bool(const C&)>& pop, const std::string& label){
    const int n = st.size();
    std::vector<Tally> cum(n + 1), nm1(n + 1);
    for(const C& c : cands){
      if(c.group < 0 || !pop(c))continue;
      // la ventana estrecha ya está dentro de la del módulo; sin ella, la base es la ventana del módulo
      if(!InWindow(c))continue;
      cum[0].n[c.group]++;
      for(int i = 0; i < n; i++){ if(!st[i].pass(c))break; cum[i + 1].n[c.group]++; }
      for(int i = 0; i < n; i++)if(Pass(st, c, i))nm1[i].n[c.group]++;
      if(Pass(st, c))nm1[n].n[c.group]++;
    }
    std::cout << std::endl << "== " << label << " ==" << std::endl << "Cortes acumulados (candidatos):" << std::endl;
    Header();
    Row("cortes del módulo", cum[0], scaleB, nEvents);
    for(int i = 0; i < n; i++)Row("+ " + st[i].name, cum[i + 1], scaleB, nEvents);
    std::cout << std::endl << "N-1 (todos los cortes menos uno):" << std::endl;
    Header();
    Row("todos", nm1[n], scaleB, nEvents);
    for(int i = 0; i < n; i++)Row("sin " + st[i].name, nm1[i], scaleB, nEvents);
  };
  Tables(steps, [](const C&){ return true; }, "Cadena común (DefaultSteps), todos");
  Tables(InSliceSteps(), [](const C& c){ return c.inSlice; }, "En la beam slice (InSliceSteps)");
  Tables(OutSliceSteps(), [](const C& c){ return !c.inSlice; }, "Fuera de la beam slice (OutSliceSteps)");

  // Variables tras todos los cortes
  struct V { const char* name; std::function<double(const C&)> get; };
  const std::vector<V> vars = {
    {"dEdxStop",   [](const C& c){ return c.dEdxStop; }},   {"dEdxStart", [](const C& c){ return c.dEdxStart; }},
    {"chi2muStop", [](const C& c){ return c.chi2muStop; }}, {"chi2pStop", [](const C& c){ return c.chi2pStop; }},
    {"michelStop", [](const C& c){ return c.michelStop; }}, {"axisMin",   [](const C& c){ return c.axisMin; }},
    {"vtxMin",     [](const C& c){ return c.vtxMin; }},     {"cnnEM",     [](const C& c){ return c.cnnEM; }},
    {"cnnTrk",     [](const C& c){ return c.cnnTrk; }},     {"cnnMichel", [](const C& c){ return c.cnnMichel; }},
    {"actStart",   [](const C& c){ return c.actStart; }},   {"intPeak",   [](const C& c){ return c.intPeak; }},
    {"dirYStop",   [](const C& c){ return c.dirYStop; }},   {"length",    [](const C& c){ return c.length; }},
    {"nHits2",     [](const C& c){ return c.nHits2; }},     {"inSlice",   [](const C& c){ return c.inSlice; }},
    {"hasSlice",   [](const C& c){ return c.hasSlice; }}};
  std::cout << std::endl << "Variables tras todos los cortes (mediana [16 %, 84 %]; vtxMin -1 = sin vértice):" << std::endl;
  std::cout << Form("  %-12s", "");
  for(int g = 0; g < 6; g++)std::cout << Form(" %24s", kGroupName[g]);
  std::cout << std::endl;
  for(const V& v : vars){
    std::cout << Form("  %-12s", v.name);
    for(int g = 0; g < 6; g++){
      std::vector<double> val;
      for(const C& c : cands)if(c.group == g && InWindow(c) && Pass(steps, c)){ const double y = v.get(c); if(y > -900)val.push_back(y); }
      if(val.empty()){ std::cout << Form(" %24s", "-"); continue; }
      std::sort(val.begin(), val.end());
      std::cout << Form(" %7.3g [%7.3g,%7.3g]", Q(val, 0.5), Q(val, 0.16), Q(val, 0.84));
    }
    std::cout << std::endl;
  }

  // Composición del fondo que queda
  std::cout << std::endl << "Fondo tras todos los cortes, por PDG true y madre (sin escalar):" << std::endl;
  {
    std::map<std::string, int> m; int n = 0;
    for(const C& c : cands)if(c.group > 0 && InWindow(c) && Pass(steps, c)){
      m[Form("PDG %6d madre %6d %-8s %s", c.tPDG, c.tParentPDG, c.tOrigin == 2 ? "cósmico" : "haz/otro", c.inSlice ? "en BS" : "fuera")]++; n++;
    }
    std::vector<std::pair<int, std::string>> v; for(auto& kv : m)v.push_back({kv.second, kv.first});
    std::sort(v.rbegin(), v.rend());
    for(size_t i = 0; i < std::min<size_t>(v.size(), 20); i++)std::cout << Form("  %-50s %6d (%5.1f %%)", v[i].second.c_str(), v[i].first, 100.*v[i].first/n) << std::endl;
  }

  // Espectro de p de rango: completos, todos los cortes salvo la ventana
  int iWindow = -1;
  for(int i = 0; i < nS; i++)if(steps[i].name.rfind("p rango", 0) == 0)iWindow = i;
  TFile fo(out, "RECREATE");
  TH1D hS("hS", "señal;p de rango (MeV/c)", 25, 100, 600), hB("hB", "fondo (escalado);p de rango (MeV/c)", 25, 100, 600);
  TH1D hBin("hBin", "fondo en la beam slice (escalado)", 25, 100, 600), hBout("hBout", "fondo fuera (escalado)", 25, 100, 600);
  TH1D hF("hF", "señal con la K en vuelo", 25, 100, 600);
  for(const C& c : cands){
    if(!c.complete || !Pass(steps, c, iWindow))continue;
    if(c.group == 0)hS.Fill(c.pRange);
    else if(c.group == -1)hF.Fill(c.pRange);
    else { hB.Fill(c.pRange, scaleB); (c.inSlice ? hBin : hBout).Fill(c.pRange, scaleB); }
  }
  std::cout << std::endl << "Espectro de p de rango (completos, todos los cortes salvo la ventana; fondo escalado), por 20 MeV/c:" << std::endl;
  std::cout << Form("  %-10s %7s %9s %9s %11s %10s", "p (MeV/c)", "señal", "fondo", "fondo BS", "fondo fuera", "K en vuelo") << std::endl;
  for(int b = 1; b <= hS.GetNbinsX(); b++)
    std::cout << Form("  %4.0f-%4.0f  %7.0f %9.0f %9.0f %11.0f %10.0f", hS.GetBinLowEdge(b), hS.GetBinLowEdge(b + 1), hS.GetBinContent(b),
                      hB.GetBinContent(b), hBin.GetBinContent(b), hBout.GetBinContent(b), hF.GetBinContent(b)) << std::endl;
  fo.Write();
  std::cout << std::endl << "Histogramas en " << out << std::endl;
  return 0;
}
