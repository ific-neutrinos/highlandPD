// Común a mu_cand_cutflow.C y mu_cand_plots.C: lectura del árbol "cand" de mu_candidates_minitrees.C,
// categorías y cortes.
//
// Señal: tSignal con la K en reposo (tKMomEnd < 50 MeV/c), en todos los eventos. Fondo: todo lo que no es
// señal, solo en eventos normales (sin volcado), escalado al número total de eventos. "Extremo de parada"
// = el de mayor dE/dx (stopAtB), como en mu_candidates_minitrees.C.
#ifndef MU_CAND_COMMON_H
#define MU_CAND_COMMON_H
#include "TFile.h"
#include "TTree.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
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
    int group;          // 0 señal, 1 mu cósmico, 2 p, 3 pi, 4 EM, 5 resto, -1 señal en vuelo
    bool inSlice, hasSlice, complete, clearCosmic;
    int tPDG, tParentPDG, tOrigin, actStart, nHits2;
    double pRange, length, axisMin, vtxMin, facesMax, dEdxStop, dEdxStart, chi2muStop, chi2pStop, michelStop;
    double cnnTrk, cnnEM, cnnMichel, intPeak, dirYStop;
  };

  constexpr int kNGroups = 6;
  const char* kGroupName[kNGroups] = {"señal", "mu cósmico", "protón", "pi+-", "EM", "resto"};

  // -2: fondo de eventos con volcado (no se usa)
  int Group(const In& x){
    if(x.tSignal)return x.tKMomEnd < 50 ? 0 : -1;
    if(x.fullDump)return -2;
    if(std::abs(x.tPDG) == 13 && x.tOrigin == 2)return 1;
    if(std::abs(x.tPDG) == 2212)return 2;
    if(std::abs(x.tPDG) == 211)return 3;
    if(std::abs(x.tPDG) == 11 || x.tPDG == 22)return 4;
    return 5;
  }

  bool InWindow(const C& c){ return c.pRange >= 180 && c.pRange <= 300; }   // ventana del módulo
  bool ModuleNoWindow(const C& c){ return c.facesMax >= 10 && c.axisMin <= 150; }

  // Lee los candidatos con p de rango en [100, 600] que pasan los cortes del módulo salvo la ventana
  std::vector<C> Load(const char* in){
    std::vector<C> cands;
    TFile f(in);
    TTree* t = (TTree*)f.Get("cand");
    if(!t){ std::cerr << "No hay árbol cand en " << in << std::endl; return cands; }
    In x;
    t->SetBranchStatus("*", 0);
#define MC_BR(v) t->SetBranchStatus(#v, 1); t->SetBranchAddress(#v, &x.v);
    MC_BR(fullDump) MC_BR(hasSlice) MC_BR(inSlice) MC_BR(complete) MC_BR(stopAtB) MC_BR(clearCosmic) MC_BR(actPartsA) MC_BR(actPartsB) MC_BR(nHits2)
    MC_BR(tPDG) MC_BR(tParentPDG) MC_BR(tOrigin) MC_BR(tSignal)
    MC_BR(pRange) MC_BR(length) MC_BR(cnnTrk) MC_BR(cnnEM) MC_BR(cnnMichel) MC_BR(facesA) MC_BR(facesB) MC_BR(axisA) MC_BR(axisB) MC_BR(vtxA) MC_BR(vtxB)
    MC_BR(dEdxA) MC_BR(dEdxB) MC_BR(chi2muF) MC_BR(chi2muR) MC_BR(chi2pF) MC_BR(chi2pR) MC_BR(michelOwnA) MC_BR(michelOwnB) MC_BR(intPeak) MC_BR(dirYStop) MC_BR(tKMomEnd)
#undef MC_BR
    const Long64_t nEntries = t->GetEntries();
    for(Long64_t i = 0; i < nEntries; i++){
      t->GetEntry(i);
      if(x.pRange < 100 || x.pRange > 600)continue;
      C c;
      c.group = Group(x);
      if(c.group == -2)continue;
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
    return cands;
  }

  // Un corte: nombre, condición, y para dibujarlo la variable, el binning y dónde se corta
  struct Step {
    std::string name;
    std::function<bool(const C&)> pass;
    std::function<double(const C&)> var;
    int nbins; double xmin, xmax;
    std::string xtitle;
    std::vector<double> lines;   // posiciones del corte
    bool logy;
  };

  std::vector<Step> DefaultSteps(){
    return {
      {"no IsClearCosmic", [](const C& c){ return !c.clearCosmic; }, [](const C& c){ return double(c.clearCosmic); },
       2, -0.5, 1.5, "IsClearCosmic", {0.5}, false},
      {"dist. al eje < 60 cm", [](const C& c){ return c.axisMin < 60; }, [](const C& c){ return c.axisMin; },
       50, 0, 150, "distancia mínima al eje del haz (cm)", {60}, false},
      {"dE/dx parada > 3", [](const C& c){ return c.dEdxStop > 3; }, [](const C& c){ return c.dEdxStop; },
       60, 0, 20, "dE/dx medio en los 5 cm del extremo de parada (MeV/cm)", {3}, false},
      {"chi2 mu < 12", [](const C& c){ return c.chi2muStop >= 0 && c.chi2muStop < 12; }, [](const C& c){ return c.chi2muStop; },
       60, 0, 60, "#chi^{2}_{#mu}/ndf (parada en el extremo de mayor dE/dx)", {12}, false},
      {"chi2 p > 80", [](const C& c){ return c.chi2pStop > 80; }, [](const C& c){ return c.chi2pStop; },
       60, 0, 400, "#chi^{2}_{p}/ndf (misma orientación)", {80}, false},
      {"dist. al vértice > 15 cm", [](const C& c){ return c.vtxMin < 0 || c.vtxMin > 15; }, [](const C& c){ return c.vtxMin; },
       52, -4, 300, "distancia mínima al vértice del haz (cm; < 0: sin vértice)", {15}, false},
      {"cnnEM < 0.5", [](const C& c){ return c.cnnEM < 0.5; }, [](const C& c){ return c.cnnEM; },
       50, 0, 1, "puntuación CNN EM de la partícula", {0.5}, false},
      {"p rango en [200, 245]", [](const C& c){ return c.pRange >= 200 && c.pRange <= 245; }, [](const C& c){ return c.pRange; },
       60, 100, 400, "p por rango (hipótesis #mu, MeV/c)", {200, 245}, false},
      {"cnnMichel > 0.02", [](const C& c){ return c.cnnMichel > 0.02; }, [](const C& c){ return c.cnnMichel; },
       50, 0, 0.5, "puntuación CNN Michel de la partícula", {0.02}, true},
    };
  }

  // Pasa los pasos [0, upto] salvo skip
  bool Pass(const std::vector<Step>& s, const C& c, int skip = -1, int upto = 1000){
    for(int i = 0; i < (int)s.size() && i <= upto; i++)if(i != skip && !s[i].pass(c))return false;
    return true;
  }
}
#endif
