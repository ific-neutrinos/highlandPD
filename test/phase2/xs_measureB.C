// Medida B (blanco grueso): mu+ de K+ en reposo por partícula del haz de cada especie, con la salida de
// secondaryKaonProductionXS (árboles ana y truth), 23-nov-2026.
//
//   root -l -b -q 'test/phase2/xs_measureB.C+("xs_v100.root")'
//
// Especie (MC): la true del haz agrupada como la separan los Cherenkov a 6 GeV (e/mu/pi, K, p); de beam_truespecies o,
// en ficheros anteriores, de beam_truepdg. Señal B: mu+ de K+ en reposo (mu+ de desintegración no cósmico, p true en
// [230, 241] MeV/c, madre K+) cuya K viene de la cadena del haz que dispara (origen 0-3: primera interacción,
// secundarias, K+ del haz). El pileup (otra partícula del haz, origen 6) es fondo, como lo que no es de K. Señal de
// todos los eventos; fondo de los eventos normales, escalado.
// Por especie y población, con la cadena entera: partículas del haz, señal y fondo seleccionados, rendimiento por
// partícula, eficiencia (frente a los mu true del árbol truth de esa especie), pureza y error relativo de un recuento
// con k veces la estadística del MC.
#include "TFile.h"
#include "TTree.h"
#include <cmath>
#include <iostream>
#include <map>

namespace measB {
  const int NMAX = 300;
  const char* kSpName[3] = {"e/mu/pi", "K", "p"};
  int SpeciesFromPDG(int pdg){ const int a = std::abs(pdg); return (a == 11 || a == 13 || a == 211) ? 0 : a == 321 ? 1 : a == 2212 ? 2 : -1; }
  const char* Br(TTree* t, const char* nuevo, const char* viejo){ return t->GetBranch(nuevo) ? nuevo : viejo; }
  enum Cls { kSig = 0, kPileup, kOtherK, kPi, kP, kMuCos, kEM, kRest, kNCls };
  const char* kClsName[kNCls] = {"señal B (mu de K, haz que dispara)", "mu de K: pileup", "mu de K: otro", "pi", "p", "mu cósmico", "EM", "resto"};
}

int xs_measureB(const char* in, double k = 6){
  using namespace measB;
  TFile f(in);
  TTree* t = (TTree*)f.Get("ana");
  TTree* tt = (TTree*)f.Get("truth");
  if(!t){ std::cerr << "No hay árbol ana" << std::endl; return 1; }
  Int_t ncand, fulldump, beamTruePDG = 0, beamTrueSp = -1, pop[NMAX], level[NMAX], tpdg[NMAX], tproc[NMAX], torig[NMAX], tpar[NMAX], tkor[NMAX];
  Float_t tmom[NMAX];
  t->SetBranchStatus("*", 0);
  auto B = [&](TTree* tr, const char* n, void* a){ tr->SetBranchStatus(n, 1); tr->SetBranchAddress(n, a); };
  B(t, "ncand", &ncand); B(t, "fulldump", &fulldump); B(t, "cand_pop", pop); B(t, "cand_level", level);
  const bool hasSp = t->GetBranch("beam_truespecies");
  if(hasSp)B(t, "beam_truespecies", &beamTrueSp); else B(t, "beam_truepdg", &beamTruePDG);
  B(t, Br(t, "cand_truepdg", "cand_tpdg"), tpdg); B(t, Br(t, "cand_trueproc", "cand_tprocstart"), tproc); B(t, Br(t, "cand_trueorigin", "cand_torigin"), torig);
  B(t, Br(t, "cand_trueparentpdg", "cand_tparentpdg"), tpar); B(t, Br(t, "cand_truekorigin", "cand_tkorigin"), tkor); B(t, Br(t, "cand_truemom", "cand_tmom"), tmom);

  const Long64_t n = t->GetEntries();
  Long64_t nNormal = 0;
  for(Long64_t i = 0; i < n; i++){ t->GetEntry(i); if(!fulldump)nNormal++; }
  const double sc = nNormal ? double(n)/nNormal : 1.;
  const int fullLevel[2] = {8, 10};

  double nBeam[3] = {0}, sel[3][2][kNCls] = {{{0}}};
  for(Long64_t i = 0; i < n; i++){
    t->GetEntry(i);
    const int sp = hasSp ? beamTrueSp : SpeciesFromPDG(beamTruePDG);
    if(sp < 0 || sp > 2)continue;
    nBeam[sp]++;
    for(int c = 0; c < ncand && c < NMAX; c++){
      const int p = pop[c];
      if(level[c] < fullLevel[p])continue;
      const bool muK = tpdg[c] == -13 && tproc[c] == 2 && torig[c] != 2 && tmom[c] > 230 && tmom[c] < 241 && (tpar[c] == 321 || tpar[c] == 0);
      int cls;
      if(muK)cls = (tkor[c] >= 0 && tkor[c] <= 3) ? kSig : tkor[c] == 6 ? kPileup : kOtherK;
      else {
        if(fulldump)continue;
        const int a = std::abs(tpdg[c]);
        cls = a == 211 ? kPi : a == 2212 ? kP : (a == 13 && torig[c] == 2) ? kMuCos : (a == 11 || tpdg[c] == 22) ? kEM : kRest;
      }
      sel[sp][p][cls] += (cls >= kPi) ? sc : 1.;
    }
  }

  // mu true de la señal B por especie (árbol truth: mu+ de K+ en reposo con origen 0-3)
  double nTrue[3] = {0};
  if(tt){
    Int_t kor = -999, btp = 0, bts = -1;
    tt->SetBranchStatus("*", 0);
    B(tt, "truemu_korigin", &kor);
    const bool ttSp = tt->GetBranch("beam_truespecies");
    if(ttSp)B(tt, "beam_truespecies", &bts); else if(tt->GetBranch("beam_truepdg"))B(tt, "beam_truepdg", &btp); else if(tt->GetBranch("truemu_beampdg"))B(tt, "truemu_beampdg", &btp);
    for(Long64_t i = 0; i < tt->GetEntries(); i++){
      tt->GetEntry(i);
      const int sp = ttSp ? bts : SpeciesFromPDG(btp);
      if(sp >= 0 && sp <= 2 && kor >= 0 && kor <= 3)nTrue[sp]++;
    }
  }

  std::cout << "Eventos " << n << " (normales " << nNormal << ", fondo escalado x" << sc << "); error proyectado con " << k << " veces la estadística del MC" << std::endl;
  const char* popName[2] = {"en la beam slice", "fuera"};
  for(int sp = 0; sp < 3; sp++){
    std::cout << std::endl << Form("== Haz %s: %.0f partículas; mu de señal B true (K en reposo, todo el detector): %.0f (%.4f por partícula) ==",
                                   kSpName[sp], nBeam[sp], nTrue[sp], nBeam[sp] ? nTrue[sp]/nBeam[sp] : 0.) << std::endl;
    double S2 = 0, T2 = 0;
    for(int p = 0; p < 2; p++){
      double tot = 0; for(int c = 0; c < kNCls; c++)tot += sel[sp][p][c];
      const double S = sel[sp][p][kSig];
      S2 += S; T2 += tot;
      std::cout << Form("  %-17s seleccionados %7.1f: señal %5.0f (%.5f por partícula; eficiencia %.3f), pureza %.3f, error %5.1f %%",
                        popName[p], tot, S, nBeam[sp] ? S/nBeam[sp] : 0., nTrue[sp] ? S/nTrue[sp] : 0., tot ? S/tot : 0., S > 0 ? 100*std::sqrt(k*tot)/(k*S) : 0.) << std::endl;
      std::cout << "      fondo:";
      for(int c = 1; c < kNCls; c++)if(sel[sp][p][c] > 0)std::cout << Form("  %s %.1f", kClsName[c], sel[sp][p][c]);
      std::cout << std::endl;
    }
    std::cout << Form("  %-17s seleccionados %7.1f: señal %5.0f (eficiencia %.3f), pureza %.3f, error %5.1f %%", "las dos", T2, S2,
                      nTrue[sp] ? S2/nTrue[sp] : 0., T2 ? S2/T2 : 0., S2 > 0 ? 100*std::sqrt(k*T2)/(k*S2) : 0.) << std::endl;
  }
  return 0;
}
