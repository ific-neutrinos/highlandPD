// Cortes en cadena sobre el árbol "cand" de mu_candidates_minitrees.C (candidatos a mu de K+ en reposo).
//
//   root -l -b -q 'test/phase2/mu_cand_cutflow.C+("mu_cand_big19.root", "cutflow.root")'
//
// Señal: tSignal con la K en reposo (tKMomEnd < 50 MeV/c), en todos los eventos. Fondo: todo lo que no es
// señal, solo en eventos normales (sin volcado), escalado al número total de eventos. "Extremo de parada"
// = el de mayor dE/dx (stopAtB), como en mu_candidates_minitrees.C.
//
// Una sola pasada por el árbol: se guardan en memoria los candidatos que pasan los cortes del módulo sin
// la ventana (o en la ventana) y con p de rango en [100, 600]; todo lo demás se calcula sobre ellos.
//
// Imprime: corte a corte (acumulado) por grupo, S y B por evento, S/B y S/sqrt(S+B); N-1 de cada corte;
// cuantiles de las variables tras todos los cortes; composición del fondo que queda; espectro de p de
// rango con los cortes salvo la ventana (solo candidatos completos: los de fuera de la beam slice fuera de
// la ventana van en resumen, sin calorimetría). Guarda los histogramas.
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace mucut {

  // Lo que se lee del árbol (mismos tipos que mu_candidates_minitrees_impl.C)
  struct In {
    Int_t fullDump, hasSlice, inSlice, complete, stopAtB, clearCosmic, actPartsA, actPartsB, nHits2;
    Int_t tPDG, tParentPDG, tOrigin, tSignal;
    Float_t pRange, length, cnnTrk, cnnEM, cnnMichel, facesA, facesB, axisA, axisB, vtxA, vtxB, dEdxA, dEdxB;
    Float_t chi2muF, chi2muR, chi2pF, chi2pR, michelOwnA, michelOwnB, intPeak, dirYStop, tKMomEnd;
  };

  // Candidato en memoria, con las variables ya referidas al extremo de parada
  struct C {
    int group;          // 0 señal, 1 mu cósmico, 2 p, 3 pi, 4 EM, 5 resto, -1 señal en vuelo, -2 fondo de eventos con volcado
    bool inSlice, hasSlice, complete, clearCosmic;
    int tPDG, tParentPDG, tOrigin, actStart, nHits2;
    double pRange, length, axisMin, vtxMin, facesMax, dEdxStop, dEdxStart, chi2muStop, chi2pStop, michelStop;
    double cnnTrk, cnnEM, cnnMichel, intPeak, dirYStop;
  };

  const char* kGroupName[6] = {"señal", "mu cósmico", "protón", "pi+-", "EM", "resto"};

  int Group(const In& x){
    if(x.tSignal)return x.tKMomEnd < 50 ? 0 : -1;
    if(x.fullDump)return -2;
    if(std::abs(x.tPDG) == 13 && x.tOrigin == 2)return 1;
    if(std::abs(x.tPDG) == 2212)return 2;
    if(std::abs(x.tPDG) == 211)return 3;
    if(std::abs(x.tPDG) == 11 || x.tPDG == 22)return 4;
    return 5;
  }

  bool InWindow(const C& c){ return c.pRange >= 180 && c.pRange <= 300; }
  bool ModuleNoWindow(const C& c){ return c.facesMax >= 10 && c.axisMin <= 150; }

  struct Step { std::string name; std::function<bool(const C&)> pass; };

  std::vector<Step> DefaultSteps(){
    return {
      {"no IsClearCosmic",         [](const C& c){ return !c.clearCosmic; }},
      {"dist. al eje < 60 cm",     [](const C& c){ return c.axisMin < 60; }},
      {"dE/dx parada > 3",         [](const C& c){ return c.dEdxStop > 3; }},
      {"chi2 mu < 12",             [](const C& c){ return c.chi2muStop >= 0 && c.chi2muStop < 12; }},
      {"chi2 p > 80",              [](const C& c){ return c.chi2pStop > 80; }},
      {"dist. al vértice > 15 cm", [](const C& c){ return c.vtxMin < 0 || c.vtxMin > 15; }},
      {"cnnEM < 0.5",              [](const C& c){ return c.cnnEM < 0.5; }},
    };
  }

  // Pasa los pasos [0, upto] salvo skip
  bool Pass(const std::vector<Step>& s, const C& c, int skip = -1, int upto = 1000){
    for(int i = 0; i < (int)s.size() && i <= upto; i++)if(i != skip && !s[i].pass(c))return false;
    return true;
  }

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

int mu_cand_cutflow(const char* in, const char* out = "mu_cand_cutflow.root", long nEvents = 73070, long nEventsNormal = 65977){
  using namespace mucut;
  TFile f(in);
  TTree* t = (TTree*)f.Get("cand");
  if(!t){ std::cerr << "No hay árbol cand en " << in << std::endl; return 1; }

  // Lectura: solo las ramas que hacen falta
  In x;
  t->SetBranchStatus("*", 0);
#define MC_BR(v) t->SetBranchStatus(#v, 1); t->SetBranchAddress(#v, &x.v);
  MC_BR(fullDump) MC_BR(hasSlice) MC_BR(inSlice) MC_BR(complete) MC_BR(stopAtB) MC_BR(clearCosmic) MC_BR(actPartsA) MC_BR(actPartsB) MC_BR(nHits2)
  MC_BR(tPDG) MC_BR(tParentPDG) MC_BR(tOrigin) MC_BR(tSignal)
  MC_BR(pRange) MC_BR(length) MC_BR(cnnTrk) MC_BR(cnnEM) MC_BR(cnnMichel) MC_BR(facesA) MC_BR(facesB) MC_BR(axisA) MC_BR(axisB) MC_BR(vtxA) MC_BR(vtxB)
  MC_BR(dEdxA) MC_BR(dEdxB) MC_BR(chi2muF) MC_BR(chi2muR) MC_BR(chi2pF) MC_BR(chi2pR) MC_BR(michelOwnA) MC_BR(michelOwnB) MC_BR(intPeak) MC_BR(dirYStop) MC_BR(tKMomEnd)
#undef MC_BR

  std::vector<C> cands;
  const Long64_t nEntries = t->GetEntries();
  for(Long64_t i = 0; i < nEntries; i++){
    t->GetEntry(i);
    if(x.pRange < 100 || x.pRange > 600)continue;
    C c;
    c.group = Group(x);
    if(c.group == -2)continue;   // fondo de eventos con volcado: no se usa
    const bool b = x.stopAtB;
    c.inSlice = x.inSlice; c.hasSlice = x.hasSlice; c.complete = x.complete; c.clearCosmic = x.clearCosmic;
    c.tPDG = x.tPDG; c.tParentPDG = x.tParentPDG; c.tOrigin = x.tOrigin; c.nHits2 = x.nHits2;
    c.actStart = b ? x.actPartsA : x.actPartsB;
    c.pRange = x.pRange; c.length = x.length;
    c.axisMin = std::min(x.axisA, x.axisB); c.vtxMin = std::min(x.vtxA, x.vtxB); c.facesMax = std::max(x.facesA, x.facesB);
    c.dEdxStop = b ? x.dEdxB : x.dEdxA; c.dEdxStart = b ? x.dEdxA : x.dEdxB;
    c.chi2muStop = b ? x.chi2muF : x.chi2muR; c.chi2pStop = b ? x.chi2pF : x.chi2pR;
    c.michelStop = b ? x.michelOwnB : x.michelOwnA;
    c.cnnTrk = x.cnnTrk; c.cnnEM = x.cnnEM; c.cnnMichel = x.cnnMichel; c.intPeak = x.intPeak; c.dirYStop = x.dirYStop;
    if(!ModuleNoWindow(c))continue;
    cands.push_back(c);
  }
  std::cout << "Leídas " << nEntries << " entradas; en memoria " << cands.size() << " candidatos" << std::endl;

  const double scaleB = double(nEvents)/nEventsNormal;
  const std::vector<Step> steps = DefaultSteps();
  const int nS = steps.size();
  std::cout << "Eventos " << nEvents << " (normales " << nEventsNormal << "); fondo escalado x" << scaleB << std::endl;

  // Cortes acumulados y N-1, en la ventana del módulo
  std::vector<Tally> cum(nS + 1), nm1(nS + 1);
  for(const C& c : cands){
    if(c.group < 0 || !InWindow(c))continue;
    cum[0].n[c.group]++;
    for(int i = 0; i < nS; i++){ if(!steps[i].pass(c))break; cum[i + 1].n[c.group]++; }
    for(int i = 0; i < nS; i++)if(Pass(steps, c, i))nm1[i].n[c.group]++;
    if(Pass(steps, c))nm1[nS].n[c.group]++;
  }
  std::cout << std::endl << "Cortes acumulados (candidatos):" << std::endl;
  Header();
  Row("cortes del módulo", cum[0], scaleB, nEvents);
  for(int i = 0; i < nS; i++)Row("+ " + steps[i].name, cum[i + 1], scaleB, nEvents);
  std::cout << std::endl << "N-1 (todos los cortes menos uno):" << std::endl;
  Header();
  Row("todos", nm1[nS], scaleB, nEvents);
  for(int i = 0; i < nS; i++)Row("sin " + steps[i].name, nm1[i], scaleB, nEvents);

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
  TFile fo(out, "RECREATE");
  TH1D hS("hS", "señal;p de rango (MeV/c)", 25, 100, 600), hB("hB", "fondo (escalado);p de rango (MeV/c)", 25, 100, 600);
  TH1D hBin("hBin", "fondo en la beam slice (escalado)", 25, 100, 600), hBout("hBout", "fondo fuera (escalado)", 25, 100, 600);
  TH1D hF("hF", "señal con la K en vuelo", 25, 100, 600);
  for(const C& c : cands){
    if(!c.complete || !Pass(steps, c))continue;
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
