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
    Int_t fullDump, hasSlice, inSlice, complete, stopAtB, clearCosmic, actPartsA, actPartsB, nHits2, nDaughters;
    Int_t tPDG, tParentPDG, tOrigin, tSignal, tProcStart, tFromBeam;
    Float_t tLength, tPur, tEff;
    Float_t pRange, length, cnnTrk, cnnEM, cnnMichel, facesA, facesB, axisA, axisB, vtxA, vtxB, dEdxA, dEdxB;
    Float_t chi2muF, chi2muR, chi2pF, chi2pR, michelOwnA, michelOwnB, intPeak, dirYStop, tKMomEnd, actMichelA, actMichelB;
  };

  // Candidato en memoria, con las variables ya referidas al extremo de parada
  struct C {
    int group;          // 0 señal, 1 mu cósmico, 2 p, 3 pi, 4 EM, 5 resto, -1 señal en vuelo
    bool inSlice, hasSlice, complete, clearCosmic;
    int tProcStart, tFromBeam; double tLength, tPur, tEff;   // truth: proceso de creación, longitud true inicio-fin, pureza, completitud
    int tPDG, tParentPDG, tOrigin, actStart, actStop, nHits2, nDaughters;
    double pRange, length, axisMin, vtxMin, facesMax, dEdxStop, dEdxStart, chi2muStop, chi2pStop, michelStop;
    double cnnTrk, cnnEM, cnnMichel, intPeak, dirYStop;
    double facesStop, facesStart, actMichelStop;
    double chi2muOther, dEdxA, dEdxB;   // chi2 mu con la orientación contraria; dE/dx en PositionStart y PositionEnd de Pandora
    bool stopAtB;   // distancia de cada extremo a las caras; CNN Michel de las partículas a < 10 cm de la parada
  };

  constexpr int kNGroups = 6;
  const char* kGroupName[kNGroups] = {"señal", "mu cósmico", "protón", "pi+-", "EM", "resto"};

  // Señal ampliada: además de los mu+ de K+ de la cadena del haz (tSignal), cualquier mu+ no cósmico creado
  // en una desintegración cuya madre es un K+ o no está guardada (K de interacciones secundarias cuya
  // cadena no se guarda en la truth). Para estos no se sabe si la K estaba en reposo.
  bool gBroadSignal = false;
  bool BroadSignal(const In& x){
    return x.tPDG == -13 && x.tOrigin != 2 && x.tProcStart == 2 && (x.tParentPDG == 321 || x.tParentPDG == 0);
  }

  // -2: fondo de eventos con volcado (no se usa)
  int Group(const In& x){
    if(x.tSignal)return x.tKMomEnd < 50 ? 0 : -1;
    if(gBroadSignal && BroadSignal(x))return 0;
    if(x.fullDump)return -2;
    if(std::abs(x.tPDG) == 13 && x.tOrigin == 2)return 1;
    if(std::abs(x.tPDG) == 2212)return 2;
    if(std::abs(x.tPDG) == 211)return 3;
    if(std::abs(x.tPDG) == 11 || x.tPDG == 22)return 4;
    return 5;
  }

  bool InWindow(const C& c){ return c.pRange >= 180 && c.pRange <= 300; }   // ventana del módulo
  bool ModuleNoWindow(const C& c){ return c.facesMax >= 10 && c.axisMin <= 150; }

  // Lee los candidatos con p de rango en [100, 600] que pasan los cortes del módulo salvo la ventana.
  // orientByChi2: extremo de parada = el de la orientación con mejor chi2 mu (por defecto, el de mayor dE/dx)
  std::vector<C> Load(const char* in, bool orientByChi2 = false){
    std::vector<C> cands;
    TFile f(in);
    TTree* t = (TTree*)f.Get("cand");
    if(!t){ std::cerr << "No hay árbol cand en " << in << std::endl; return cands; }
    In x;
    t->SetBranchStatus("*", 0);
#define MC_BR(v) t->SetBranchStatus(#v, 1); t->SetBranchAddress(#v, &x.v);
    MC_BR(fullDump) MC_BR(hasSlice) MC_BR(inSlice) MC_BR(complete) MC_BR(stopAtB) MC_BR(clearCosmic) MC_BR(actPartsA) MC_BR(actPartsB) MC_BR(nHits2) MC_BR(nDaughters)
    MC_BR(tPDG) MC_BR(tParentPDG) MC_BR(tOrigin) MC_BR(tSignal) MC_BR(tProcStart) MC_BR(tFromBeam) MC_BR(tLength) MC_BR(tPur) MC_BR(tEff)
    MC_BR(pRange) MC_BR(length) MC_BR(cnnTrk) MC_BR(cnnEM) MC_BR(cnnMichel) MC_BR(facesA) MC_BR(facesB) MC_BR(axisA) MC_BR(axisB) MC_BR(vtxA) MC_BR(vtxB)
    MC_BR(dEdxA) MC_BR(dEdxB) MC_BR(chi2muF) MC_BR(chi2muR) MC_BR(chi2pF) MC_BR(chi2pR) MC_BR(michelOwnA) MC_BR(michelOwnB) MC_BR(intPeak) MC_BR(dirYStop) MC_BR(tKMomEnd) MC_BR(actMichelA) MC_BR(actMichelB)
#undef MC_BR
    const Long64_t nEntries = t->GetEntries();
    for(Long64_t i = 0; i < nEntries; i++){
      t->GetEntry(i);
      if(x.pRange < 100 || x.pRange > 600)continue;
      C c;
      c.group = Group(x);
      if(c.group == -2)continue;
      bool b = x.stopAtB;
      if(orientByChi2){
        const double cB = x.chi2muF < 0 ? 1e9 : x.chi2muF, cA = x.chi2muR < 0 ? 1e9 : x.chi2muR;   // F: parada en B
        if(cA < 1e9 || cB < 1e9)b = cB <= cA;
      }
      c.inSlice = x.inSlice; c.hasSlice = x.hasSlice; c.complete = x.complete; c.clearCosmic = x.clearCosmic;
      c.tProcStart = x.tProcStart; c.tFromBeam = x.tFromBeam; c.tLength = x.tLength; c.tPur = x.tPur; c.tEff = x.tEff;
      c.tPDG = x.tPDG; c.tParentPDG = x.tParentPDG; c.tOrigin = x.tOrigin; c.nHits2 = x.nHits2;
      c.actStart = b ? x.actPartsA : x.actPartsB; c.actStop = b ? x.actPartsB : x.actPartsA; c.nDaughters = x.nDaughters;
      c.chi2muOther = b ? x.chi2muR : x.chi2muF; c.dEdxA = x.dEdxA; c.dEdxB = x.dEdxB; c.stopAtB = b;
      c.facesStop = b ? x.facesB : x.facesA; c.facesStart = b ? x.facesA : x.facesB; c.actMichelStop = b ? x.actMichelB : x.actMichelA;
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

  // Cortes por población (la beam slice y lo de fuera los reconstruyen pasadas distintas de Pandora y tienen
  // fondos distintos: dentro la traza del haz y sus pi y p; fuera cósmicos y EM). Cada uno se define una vez.
  Step MakeStep(const std::string& key){
    if(key == "noclear")   return {"no IsClearCosmic", [](const C& c){ return !c.clearCosmic; }, [](const C& c){ return double(c.clearCosmic); },
                                   2, -0.5, 1.5, "IsClearCosmic", {0.5}, false};
    if(key == "caras20")   return {"extremos a > 20 cm de las caras", [](const C& c){ return std::min(c.facesStop, c.facesStart) > 20; },
                                   [](const C& c){ return std::min(c.facesStop, c.facesStart); }, 50, 0, 150, "distancia mínima de los dos extremos a las caras (cm)", {20}, false};
    if(key == "caras35")   return {"extremos a > 35 cm de las caras", [](const C& c){ return std::min(c.facesStop, c.facesStart) > 35; },
                                   [](const C& c){ return std::min(c.facesStop, c.facesStart); }, 50, 0, 150, "distancia mínima de los dos extremos a las caras (cm)", {35}, false};
    if(key == "eje<60")    return {"dist. al eje < 60 cm", [](const C& c){ return c.axisMin < 60; }, [](const C& c){ return c.axisMin; },
                                   50, 0, 150, "distancia mínima al eje del haz (cm)", {60}, false};
    if(key == "eje>10")    return {"dist. al eje > 10 cm", [](const C& c){ return c.axisMin > 10; }, [](const C& c){ return c.axisMin; },
                                   50, 0, 150, "distancia mínima al eje del haz (cm)", {10}, false};
    if(key == "vtx>15")    return {"dist. al vértice > 15 cm", [](const C& c){ return c.vtxMin < 0 || c.vtxMin > 15; }, [](const C& c){ return c.vtxMin; },
                                   52, -4, 300, "distancia mínima al vértice del haz (cm; < 0: sin vértice)", {15}, false};
    if(key == "chi2mu<5")  return {"chi2 mu < 5", [](const C& c){ return c.chi2muStop >= 0 && c.chi2muStop < 5; }, [](const C& c){ return c.chi2muStop; },
                                   60, 0, 60, "#chi^{2}_{#mu}/ndf (parada en el extremo de mayor dE/dx)", {5}, false};
    if(key == "chi2mu<6")  return {"chi2 mu < 6", [](const C& c){ return c.chi2muStop >= 0 && c.chi2muStop < 6; }, [](const C& c){ return c.chi2muStop; },
                                   60, 0, 60, "#chi^{2}_{#mu}/ndf (parada en el extremo de mayor dE/dx)", {6}, false};
    if(key == "dedxStop>3")return {"dE/dx parada > 3", [](const C& c){ return c.dEdxStop > 3; }, [](const C& c){ return c.dEdxStop; },
                                   60, 0, 20, "dE/dx medio en los 5 cm del extremo de parada (MeV/cm)", {3}, false};
    if(key == "dedxStart>1.8")return {"dE/dx inicio > 1.8", [](const C& c){ return c.dEdxStart > 1.8; }, [](const C& c){ return c.dEdxStart; },
                                   50, 0, 10, "dE/dx medio en los 5 cm del extremo de inicio (MeV/cm)", {1.8}, false};
    if(key == "chi2p>80")  return {"chi2 p > 80", [](const C& c){ return c.chi2pStop > 80; }, [](const C& c){ return c.chi2pStop; },
                                   60, 0, 400, "#chi^{2}_{p}/ndf (misma orientación)", {80}, false};
    if(key == "cnnEM")     return {"cnnEM < 0.5", [](const C& c){ return c.cnnEM < 0.5; }, [](const C& c){ return c.cnnEM; },
                                   50, 0, 1, "puntuación CNN EM de la partícula", {0.5}, false};
    if(key == "ventana")   return {"p rango en [200, 245]", [](const C& c){ return c.pRange >= 200 && c.pRange <= 245; }, [](const C& c){ return c.pRange; },
                                   60, 100, 400, "p por rango (hipótesis #mu, MeV/c)", {200, 245}, false};
    if(key == "ventana205")return {"p rango en [205, 240]", [](const C& c){ return c.pRange >= 205 && c.pRange <= 240; }, [](const C& c){ return c.pRange; },
                                   60, 100, 400, "p por rango (hipótesis #mu, MeV/c)", {205, 240}, false};
    if(key == "michel")    return {"cnnMichel > 0.02", [](const C& c){ return c.cnnMichel > 0.02; }, [](const C& c){ return c.cnnMichel; },
                                   50, 0, 0.5, "puntuación CNN Michel de la partícula", {0.02}, true};
    std::cerr << "MakeStep: corte desconocido " << key << std::endl;
    exit(1);
  }

  std::vector<Step> Steps(const std::vector<std::string>& keys){
    std::vector<Step> s;
    for(const std::string& k : keys)s.push_back(MakeStep(k));
    return s;
  }
  std::vector<Step> InSliceSteps(){
    // 21-nov: tras el N-1 de la cadena larga (fuera eje > 10, chi2 p, dE/dx inicio y cnnEM) y el escaneo de
    // mu_cand_scan.C (chi2 mu, Michel y ventana; superficie plana alrededor de estos valores)
    return Steps({"caras35", "vtx>15", "chi2mu<6", "dedxStop>3", "ventana205", "michel"});
  }
  std::vector<Step> OutSliceSteps(){
    // 21-nov: fuera caras > 20 y cnnEM (N-1); chi2 mu, Michel y ventana de mu_cand_scan.C
    return Steps({"noclear", "eje<60", "chi2mu<6", "dedxStop>3", "chi2p>80", "dedxStart>1.8", "ventana205", "michel"});
  }

  // Pasa los pasos [0, upto] salvo skip
  bool Pass(const std::vector<Step>& s, const C& c, int skip = -1, int upto = 1000){
    for(int i = 0; i < (int)s.size() && i <= upto; i++)if(i != skip && !s[i].pass(c))return false;
    return true;
  }
}
#endif
