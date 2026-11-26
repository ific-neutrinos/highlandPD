// Viabilidad de la medida A (producción de K+ en la primera interacción del haz) con la salida de
// secondaryKaonProductionXS (árbol ana), 23-nov-2026.
//
//   root -l -b -q 'test/phase2/xs_feasibility.C+("xs_v100.root")'
//
// Para los candidatos que pasan la cadena entera de su población (cand_level máximo), por clase:
//   - mu de K+ en reposo por el origen de su K (A, secundaria, K+ del haz, pileup, cadena cortada), de todos los eventos;
//   - fondo que no es de K (pi, p, mu cósmico, EM, resto), solo de eventos normales (fulldump == 0), escalado.
// Conexión con el vértice (un PFParticle une el vértice con un extremo del mu, o un extremo del mu está a < 10 cm del
// vértice) con dos vértices: el de Pandora (final del primario) y el perfecto (MC: último hit del haz true, en el
// sistema de la reco). Con y sin exigir que el haz interaccione a z true > 30 cm (beam_trueendpos o tbeamend).
// Proyección: error relativo de un recuento de la señal A con k veces la estadística del MC (k = 6 ~ datos a 6 GeV con
// un bloque de ~73.000 eventos de MC; con más MC, k = 6*73070/N_MC: segundo argumento, 26-nov-2026).
//   root -l -b -q 'test/phase2/xs_feasibility.C+("xs_v100to105d01.root", 1.196)'   // v1_00-v1_05, 366.460 eventos
// Lee los nombres antiguos (cand_tpdg, tbeamend) y los nuevos (cand_truepdg, beam_trueendpos).
#include "TFile.h"
#include "TTree.h"
#include <cmath>
#include <iostream>
#include <string>

namespace feas {
  const int NMAX = 300;
  enum Cls { kA = 0, kSec, kBeamK, kPileup, kCut, kPi, kP, kMuCos, kEM, kRest, kNCls };
  const char* kName[kNCls] = {"mu de K: A", "mu de K: secundaria", "mu de K: K+ del haz", "mu de K: pileup", "mu de K: cadena cortada",
                              "pi", "p", "mu cósmico", "EM", "resto"};
  // nombre nuevo si está, si no el antiguo
  const char* Br(TTree* t, const char* nuevo, const char* viejo){ return t->GetBranch(nuevo) ? nuevo : viejo; }
}

int xs_feasibility(const char* in, double kData = 6){
  using namespace feas;
  TFile f(in);
  TTree* t = (TTree*)f.Get("ana");
  if(!t){ std::cerr << "No hay árbol ana" << std::endl; return 1; }
  Int_t ncand, fulldump, pop[NMAX], level[NMAX], tpdg[NMAX], tproc[NMAX], torig[NMAX], tpar[NMAX], tkor[NMAX];
  Float_t tmom[NMAX], lkSdV[NMAX], lkSdE[NMAX], lkPdV[NMAX], lkPdE[NMAX], dS[NMAX], dP[NMAX];
  Float_t lktSdV[NMAX], lktSdE[NMAX], lktPdV[NMAX], lktPdE[NMAX], dSt[NMAX], dPt[NMAX], beamEnd[4] = {-999, -999, -999, -999};
  t->SetBranchStatus("*", 0);
  auto B = [&](const char* n, void* a){ t->SetBranchStatus(n, 1); t->SetBranchAddress(n, a); };
  B("ncand", &ncand); B("fulldump", &fulldump); B("cand_pop", pop); B("cand_level", level);
  B(Br(t, "cand_truepdg", "cand_tpdg"), tpdg); B(Br(t, "cand_trueproc", "cand_tprocstart"), tproc); B(Br(t, "cand_trueorigin", "cand_torigin"), torig);
  B(Br(t, "cand_trueparentpdg", "cand_tparentpdg"), tpar); B(Br(t, "cand_truekorigin", "cand_tkorigin"), tkor); B(Br(t, "cand_truemom", "cand_tmom"), tmom);
  B("cand_lk_s_dv", lkSdV); B("cand_lk_s_de", lkSdE); B("cand_lk_p_dv", lkPdV); B("cand_lk_p_de", lkPdE); B("cand_dstart_v", dS); B("cand_dstop_v", dP);
  B("cand_lkt_s_dv", lktSdV); B("cand_lkt_s_de", lktSdE); B("cand_lkt_p_dv", lktPdV); B("cand_lkt_p_de", lktPdE); B("cand_dstart_tv", dSt); B("cand_dstop_tv", dPt);
  const bool newBeam = t->GetBranch("beam_trueendpos");
  B(newBeam ? "beam_trueendpos" : "tbeamend", beamEnd);   // beam_trueendpos tiene 4 componentes, tbeamend 3: la z es [2] en los dos

  const Long64_t n = t->GetEntries();
  Long64_t nNormal = 0;
  for(Long64_t i = 0; i < n; i++){ t->GetEntry(i); if(!fulldump)nNormal++; }
  const double sc = nNormal ? double(n)/nNormal : 1.;
  const int fullLevel[2] = {8, 10};   // cadena entera: población + 7 cortes dentro, + 9 fuera
  auto conn = [](double sdv, double sde, double pdv, double pde, double ds, double dp){
    if(ds < 0)return false;
    return (sdv >= 0 && sdv < 10 && sde < 10) || (pdv >= 0 && pdv < 10 && pde < 10) || ds < 10 || dp < 10;
  };

  // [z: 0 todas, 1 z > 30][pop][clase]: total, conectados (Pandora), conectados (perfecto)
  double N[2][2][kNCls] = {{{0}}}, CP[2][2][kNCls] = {{{0}}}, CT[2][2][kNCls] = {{{0}}};
  for(Long64_t i = 0; i < n; i++){
    t->GetEntry(i);
    const bool zok = beamEnd[2] > 30;
    for(int k = 0; k < ncand && k < NMAX; k++){
      const int p = pop[k];
      if(level[k] < fullLevel[p])continue;
      const bool muK = tpdg[k] == -13 && tproc[k] == 2 && torig[k] != 2 && tmom[k] > 230 && tmom[k] < 241 && (tpar[k] == 321 || tpar[k] == 0);
      int c;
      if(muK)c = (tkor[k] == 0 || tkor[k] == 1) ? kA : tkor[k] == 2 ? kSec : tkor[k] == 3 ? kBeamK : tkor[k] == 6 ? kPileup : kCut;
      else {
        if(fulldump)continue;
        const int a = std::abs(tpdg[k]);
        c = a == 211 ? kPi : a == 2212 ? kP : (a == 13 && torig[k] == 2) ? kMuCos : (a == 11 || tpdg[k] == 22) ? kEM : kRest;
      }
      const double w = c >= kPi ? sc : 1.;
      const bool cp = conn(lkSdV[k], lkSdE[k], lkPdV[k], lkPdE[k], dS[k], dP[k]);
      const bool ct = conn(lktSdV[k], lktSdE[k], lktPdV[k], lktPdE[k], dSt[k], dPt[k]);
      for(int z = 0; z < 2; z++){
        if(z == 1 && !zok)continue;
        N[z][p][c] += w; if(cp)CP[z][p][c] += w; if(ct)CT[z][p][c] += w;
      }
    }
  }

  std::cout << "Eventos " << n << " (normales " << nNormal << ", fondo escalado x" << sc << "); "
            << (newBeam ? "beam_trueendpos" : "tbeamend") << " para la z de la interacción" << std::endl;
  const char* popName[2] = {"en la beam slice", "fuera de la beam slice"};
  for(int z = 0; z < 2; z++)for(int p = 0; p < 2; p++){
    std::cout << std::endl << "== " << popName[p] << ", cadena entera" << (z ? ", interacción del haz a z > 30 cm" : "") << " ==" << std::endl
              << Form("  %-26s %10s %22s %22s", "clase", "candidatos", "conectados (Pandora)", "conectados (perfecto)") << std::endl;
    double tot = 0, totP = 0, totT = 0;
    for(int c = 0; c < kNCls; c++){
      tot += N[z][p][c]; totP += CP[z][p][c]; totT += CT[z][p][c];
      std::cout << Form("  %-26s %10.1f %13.1f (%3.0f %%) %13.1f (%3.0f %%)", kName[c], N[z][p][c], CP[z][p][c], N[z][p][c] ? 100*CP[z][p][c]/N[z][p][c] : 0.,
                        CT[z][p][c], N[z][p][c] ? 100*CT[z][p][c]/N[z][p][c] : 0.) << std::endl;
    }
    const double S = N[z][p][kA], SP = CP[z][p][kA], ST = CT[z][p][kA];
    auto err = [](double s, double b, double k){ return s > 0 ? std::sqrt(k*(s + b))/(k*s) : 0.; };
    std::cout << Form("  pureza para A: sin conexión %.3f, con la de Pandora %.3f, con la perfecta %.3f", tot ? S/tot : 0., totP ? SP/totP : 0., totT ? ST/totT : 0.) << std::endl
              << Form("  error relativo de un recuento de A con %.2f veces esta estadística: sin conexión %.1f %%, Pandora %.1f %%, perfecta %.1f %%", kData,
                      100*err(S, tot - S, kData), 100*err(SP, totP - SP, kData), 100*err(ST, totT - ST, kData)) << std::endl;
  }
  return 0;
}
