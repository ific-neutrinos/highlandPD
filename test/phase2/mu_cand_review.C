// Revisión de la selección de candidatos a mu de K+ en reposo (23-nov-2026), sobre el árbol cand de
// mu_candidates_minitrees.C (ficheros desde el 20-nov):
//   A. corte vtx>15 (beam slice) por origen de la K; B. corte de Michel: CNN de la partícula frente a hits Michel
//   alrededor de la parada; C. candidatos tipo shower; D. orientación (dE/dx, chi2, cualquiera); E. K y mu
//   fundidos en una traza; F. qué hay bajo el pico; G. resolución en longitud y en qué extremo se pierde;
//   H. candidatos repetidos; K. pi+ de K+ -> pi+ pi0.
//
//   root -l -b -q 'test/phase2/mu_cand_review.C+("mu_cand_big19v2.root")'
#include "mu_cand_common.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TMath.h"
#include <map>
#include <set>
using namespace mucut;

const double kSc = 73070./65977.;   // fondo de eventos normales escalado a todos
struct Tally { double S = 0, B = 0; void add(const C& c){ if(c.group == 0)S++; else if(c.group > 0)B += kSc; } };
std::string Fmt(const Tally& t){
  return Form("S %5.0f  B %7.1f  S/B %.3f  S/sqrt(S+B) %.2f", t.S, t.B, t.B > 0 ? t.S/t.B : 0., (t.S + t.B) > 0 ? t.S/std::sqrt(t.S + t.B) : 0.);
}
int OClass(const C& c){ if(c.tKOrigin == 0 || c.tKOrigin == 1)return 0; if(c.tKOrigin == 2)return 1; if(c.tKOrigin == 3)return 2; return 3; }
const char* kOName[4] = {"A (1a interacción)", "secundaria", "K+ del haz", "otro"};

void Quant(TTree* t, const char* var, const char* sel, const char* label){
  const Long64_t n = t->Draw(var, sel, "goff");
  if(n <= 0){ printf("    %-44s n 0\n", label); return; }
  std::vector<double> v(t->GetV1(), t->GetV1() + n);
  std::sort(v.begin(), v.end());
  auto q = [&](double f){ return v[std::min<size_t>(v.size() - 1, size_t(f*v.size()))]; };
  printf("    %-44s n %5lld  mediana %6.2f  [16 %%, 84 %%] = [%6.2f, %6.2f]  (84-16)/2 = %.2f\n", label, n, q(0.5), q(0.16), q(0.84), 0.5*(q(0.84) - q(0.16)));
}

void ReviewPart1(const char* in){
  std::vector<C> cands = Load(in);
  const std::vector<Step> stIn = InSliceSteps(), stOut = OutSliceSteps();
  const int iVtx = 1, iMichelIn = int(stIn.size()) - 1, iMichelOut = int(stOut.size()) - 1;

  // ---------------- A ----------------
  printf("\n== A. Corte vtx>15 (dentro de la beam slice), señal por origen de la K ==\n");
  for(int stage = 0; stage < 2; stage++){
    double tot[4] = {0}, fail[4] = {0}; Tally all, rem;
    for(const C& c : cands){
      if(!c.inSlice)continue;
      const bool base = stage == 0 ? (stIn[0].pass(c) && c.pRange >= 205 && c.pRange <= 240) : Pass(stIn, c, iVtx);
      if(!base)continue;
      const bool f = !stIn[iVtx].pass(c);
      if(c.group == 0){ const int o = OClass(c); tot[o]++; if(f)fail[o]++; }
      all.add(c); if(f)rem.add(c);
    }
    printf("  %s\n", stage == 0 ? "base: caras35 + ventana [205, 240]" : "base: toda la cadena salvo vtx>15 (N-1)");
    for(int o = 0; o < 4; o++)printf("    señal %-20s %4.0f, la quita vtx>15: %4.0f (%3.0f %%)\n", kOName[o], tot[o], fail[o], tot[o] ? 100*fail[o]/tot[o] : 0.);
    printf("    fondo %7.1f, lo quita: %7.1f (%3.0f %%)\n", all.B, rem.B, all.B ? 100*rem.B/all.B : 0.);
  }

  // ---------------- B ----------------
  printf("\n== B. Michel: corte actual (CNN Michel media de la partícula > 0,02) frente a hits Michel alrededor de la parada ==\n");
  struct Opt { const char* n; std::function<bool(const C&)> f; };
  const std::vector<Opt> opts = {
    {"sin corte de Michel",          [](const C&){ return true; }},
    {"actual: cnnMichel > 0,02",     [](const C& c){ return c.cnnMichel > 0.02; }},
    {"arNMichel5 >= 1",              [](const C& c){ return c.arNMichel[0] >= 1; }},
    {"arNMichel5 >= 2",              [](const C& c){ return c.arNMichel[0] >= 2; }},
    {"arNMichel10 >= 2",             [](const C& c){ return c.arNMichel[1] >= 2; }},
    {"arNMichel10 >= 3",             [](const C& c){ return c.arNMichel[1] >= 3; }},
    {"arNMichel10 >= 5",             [](const C& c){ return c.arNMichel[1] >= 5; }},
    {"arNMichel20 >= 3",             [](const C& c){ return c.arNMichel[2] >= 3; }},
    {"arMichel10 > 0,1",             [](const C& c){ return c.arMichel[1] > 0.1; }},
    {"actual O arNMichel10 >= 3",    [](const C& c){ return c.cnnMichel > 0.02 || c.arNMichel[1] >= 3; }},
    {"actual Y arNMichel10 >= 2",    [](const C& c){ return c.cnnMichel > 0.02 && c.arNMichel[1] >= 2; }},
  };
  for(int pop = 1; pop >= 0; pop--){
    const std::vector<Step>& st = pop ? stIn : stOut;
    printf("  %s (resto de la cadena aplicado)\n", pop ? "dentro de la beam slice" : "fuera");
    for(const Opt& o : opts){
      Tally t;
      for(const C& c : cands)if(c.inSlice == (bool)pop && Pass(st, c, pop ? iMichelIn : iMichelOut) && o.f(c))t.add(c);
      printf("    %-28s %s\n", o.n, Fmt(t).c_str());
    }
  }

  // ---------------- C ----------------
  printf("\n== C. Candidatos tipo shower (Pandora no hace traza; longitud de la traza forzada), hoy excluidos ==\n");
  TFile f(in);
  TTree* t = (TTree*)f.Get("cand");
  const char* mod = "max(facesA,facesB)>=10&&min(axisA,axisB)<=150";
  const char* win = "pRange>=205&&pRange<=240";
  const char* sig = "tSignal&&tKMomEnd<50";
  const char* bkg = "!tSignal&&!fullDump";
  for(int pop = 1; pop >= 0; pop--){
    for(int shw = 0; shw < 2; shw++){
      const double S = t->GetEntries(Form("%s&&%s&&%s&&inSlice==%d&&isShower==%d", mod, win, sig, pop, shw));
      const double B = kSc*t->GetEntries(Form("%s&&%s&&%s&&inSlice==%d&&isShower==%d", mod, win, bkg, pop, shw));
      printf("  %-8s %-7s módulo + ventana: S %5.0f  B %8.0f  S/B %.4f\n", pop ? "dentro" : "fuera", shw ? "shower" : "traza", S, B, B > 0 ? S/B : 0.);
    }
  }
  printf("  p de rango (traza forzada) de los mu de señal tipo shower, cortes del módulo:\n");
  Quant(t, "pRange", Form("%s&&%s&&isShower==1&&pRange>100&&pRange<600", sig, mod), "señal shower: p de rango (MeV/c)");
  Quant(t, "pRange", Form("%s&&%s&&isShower==0&&pRange>100&&pRange<600", sig, mod), "señal traza: p de rango (MeV/c)");
  // selección mínima para showers: ventana, geometría y Michel de la partícula, chi2 de la traza forzada
  for(int pop = 1; pop >= 0; pop--){
    const char* geo = pop ? "min(facesA,facesB)>35&&(min(vtxA,vtxB)<0||min(vtxA,vtxB)>15)" : "!clearCosmic&&min(axisA,axisB)<60";
    const char* sel = Form("%s&&%s&&isShower==1&&inSlice==%d&&%s&&cnnMichel>0.02&&ndf>0&&chi2mu/ndf<6", mod, win, pop, geo);
    const double S = t->GetEntries(Form("%s&&%s", sel, sig)), B = kSc*t->GetEntries(Form("%s&&%s", sel, bkg));
    printf("  showers %-7s con geometría de la cadena + cnnMichel > 0,02 + chi2 mu (traza forzada) < 6: S %4.0f  B %6.1f\n", pop ? "dentro" : "fuera", S, B);
  }

  // ---------------- D ----------------
  printf("\n== D. Orientación (extremo de parada): dE/dx (por defecto) frente a chi2 mu; señal con cortes del módulo y ventana ==\n");
  std::vector<C> candsChi2 = Load(in, true);
  for(int which = 0; which < 2; which++){
    const std::vector<C>& cc = which ? candsChi2 : cands;
    for(int pop = 1; pop >= 0; pop--){
      double n[2] = {0, 0}, ok[2] = {0, 0};
      for(const C& c : cc){
        if(c.group != 0 || c.inSlice != (bool)pop || c.pRange < 205 || c.pRange > 240)continue;
        if(c.tStartZ < -900 || c.tEndZ < -900 || std::abs(c.tEndZ - c.tStartZ) < 10)continue;
        const int k = OClass(c) == 0 ? 0 : 1;
        n[k]++; if(std::abs(c.stopZ - c.tEndZ) < std::abs(c.stopZ - c.tStartZ))ok[k]++;
      }
      printf("  %-6s %-7s: bien orientados  A %3.0f de %3.0f (%3.0f %%)   otros orígenes %3.0f de %3.0f (%3.0f %%)\n",
             which ? "chi2" : "dE/dx", pop ? "dentro" : "fuera", ok[0], n[0], n[0] ? 100*ok[0]/n[0] : 0., ok[1], n[1], n[1] ? 100*ok[1]/n[1] : 0.);
    }
  }
  // efecto en la cadena: señal que pasa con cada orientación
  for(int which = 0; which < 2; which++){
    const std::vector<C>& cc = which ? candsChi2 : cands;
    Tally tin, tout;
    for(const C& c : cc){ if(c.inSlice){ if(Pass(stIn, c))tin.add(c); } else if(Pass(stOut, c))tout.add(c); }
    printf("  cadena corta con orientación por %-5s dentro: %s\n", which ? "chi2" : "dE/dx", Fmt(tin).c_str());
    printf("  %38s fuera:  %s\n", "", Fmt(tout).c_str());
  }

  // ---------------- E ----------------
  printf("\n== E. K y mu fundidos en una traza (p de rango por encima de la ventana) ==\n");
  for(int pop = 1; pop >= 0; pop--){
    const double all = t->GetEntries(Form("%s&&%s&&isShower==0&&inSlice==%d&&pRange>100&&pRange<600", sig, mod, pop));
    const double hi = t->GetEntries(Form("%s&&%s&&isShower==0&&inSlice==%d&&pRange>260&&pRange<600", sig, mod, pop));
    const char* pk = "intPeak>6&&((intPeakRR>40&&intPeakRR<65)||(length-intPeakRR>40&&length-intPeakRR<65))";
    const double S = t->GetEntries(Form("%s&&%s&&isShower==0&&inSlice==%d&&pRange>260&&pRange<600&&%s", sig, mod, pop, pk));
    const double B = kSc*t->GetEntries(Form("%s&&%s&&isShower==0&&inSlice==%d&&pRange>260&&pRange<600&&%s", bkg, mod, pop, pk));
    printf("  %-7s señal (traza, módulo) %4.0f; con p > 260: %4.0f (%2.0f %%); de ellas con pico de Bragg interior a 40-65 cm de un extremo: S %4.0f (B %7.0f)\n",
           pop ? "dentro" : "fuera", all, hi, all ? 100*hi/all : 0., S, B);
  }

  // ---------------- F ----------------
  printf("\n== F. mu+ de K+ que no son señal (pileup o cadena cortada en este fichero) ==\n");
  for(int pop = 1; pop >= 0; pop--){
    const std::vector<Step>& st = pop ? stIn : stOut;
    const int iWin = pop ? 4 : 6;
    double sAll = 0, kmu = 0, kmuWin = 0, kmuNoWin = 0; Tally tt;
    for(const C& c : cands){
      if(c.inSlice != (bool)pop)continue;
      const bool isKmu = c.group > 0 && c.tPDG == -13 && c.tParentPDG == 321;
      if(Pass(st, c)){ tt.add(c); if(c.group == 0)sAll++; if(isKmu)kmu += kSc; }
      if(Pass(st, c, iWin) && isKmu){ kmuNoWin += kSc; if(c.pRange >= 205 && c.pRange <= 240)kmuWin += kSc; }
    }
    printf("  %-7s tras la cadena: señal %3.0f, fondo %6.1f, de él mu+ de K+ no señal %5.1f (%2.0f %% del fondo); sin la ventana %5.1f, en la ventana %5.1f\n",
           pop ? "dentro" : "fuera", sAll, tt.B, kmu, tt.B ? 100*kmu/tt.B : 0., kmuNoWin, kmuWin);
  }

  // ---------------- G ----------------
  printf("\n== G. Resolución en longitud (reco - true; true = distancia inicio-fin), señal traza con cortes del módulo ==\n");
  const char* gs = Form("%s&&%s&&isShower==0&&pRange>150&&pRange<350&&tLength>0", sig, mod);
  Quant(t, "length-tLength", gs, "todas");
  Quant(t, "length-tLength", Form("%s&&inSlice", gs), "dentro");
  Quant(t, "length-tLength", Form("%s&&!inSlice", gs), "fuera");
  Quant(t, "length-tLength", Form("%s&&abs(stopX)<120", gs), "parada |x| < 120 (cerca del cátodo)");
  Quant(t, "length-tLength", Form("%s&&abs(stopX)>=120&&abs(stopX)<240", gs), "parada 120 < |x| < 240");
  Quant(t, "length-tLength", Form("%s&&abs(stopX)>=240", gs), "parada |x| > 240 (cerca de las APA)");
  Quant(t, "length-tLength", Form("%s&&kN>0", gs), "con algo que acaba en el inicio (kN > 0)");
  Quant(t, "length-tLength", Form("%s&&kN==0", gs), "sin nada en el inicio (kN = 0)");
  Quant(t, "length-tLength", Form("%s&&arNMichel[1]>=3", gs), "con Michel alrededor de la parada");
  Quant(t, "length-tLength", Form("%s&&cnnMichel>0.1", gs), "Michel dentro de la traza (cnnMichel > 0,1)");

  // ---------------- H ----------------
  printf("\n== H. Señal seleccionada repetida (dos candidatos del mismo mu true), tras la cadena ==\n");
  for(int pop = 1; pop >= 0; pop--){
    std::map<std::pair<int, float>, int> m;
    for(const C& c : cands)if(c.group == 0 && c.inSlice == (bool)pop && Pass(pop ? stIn : stOut, c))m[{c.entry, (float)c.tMom}]++;
    int dup = 0; for(auto& kv : m)if(kv.second > 1)dup += kv.second - 1;
    printf("  %-7s mu distintos %zu, candidatos de más %d\n", pop ? "dentro" : "fuera", m.size(), dup);
  }
}

void ReviewPart2(const char* in){
  std::vector<C> cands = Load(in);
  const std::vector<Step> stIn = InSliceSteps(), stOut = OutSliceSteps();

  // D2: regla sin orientación: pasa si alguna de las dos orientaciones cumple chi2 mu < 6 y dE/dx de parada > 3
  printf("\n== D2. Cadena con chi2 mu y dE/dx de parada en cualquiera de las dos orientaciones ==\n");
  // índices: dentro {caras35, vtx>15, chi2mu<6, dedxStop>3, ventana205, michel}; fuera {noclear, eje<60, chi2mu<6, dedxStop>3, chi2p>80, dedxStart>1.8, ventana205, michel}
  for(int pop = 1; pop >= 0; pop--){
    const std::vector<Step>& st = pop ? stIn : stOut;
    Tally now, either, misS;
    for(const C& c : cands){
      if(c.inSlice != (bool)pop)continue;
      bool rest = true;
      for(int i = 0; i < (int)st.size(); i++){
        if(i == 2 || i == 3)continue;
        if(!pop && (i == 4 || i == 5))continue;   // chi2p y dE/dx de inicio dependen de la orientación: fuera, se quitan en las dos reglas
        if(!st[i].pass(c)){ rest = false; break; }
      }
      if(!rest)continue;
      const bool fwd = c.chi2muStop >= 0 && c.chi2muStop < 6 && c.dEdxStop > 3;
      const bool rev = c.chi2muOther >= 0 && c.chi2muOther < 6 && c.dEdxStart > 3;
      if(fwd)now.add(c);
      if(fwd || rev)either.add(c);
      if(!fwd && rev)misS.add(c);
    }
    printf("  %-7s (sin chi2p ni dE/dx de inicio fuera) solo la orientación de dE/dx: %s\n", pop ? "dentro" : "fuera", Fmt(now).c_str());
    printf("  %-7s cualquiera de las dos:                                %s\n", "", Fmt(either).c_str());
    printf("  %-7s lo que solo entra al revés:                            %s\n", "", Fmt(misS).c_str());
  }

  // F2: el bulto del ajuste: mu+ de desintegración no cósmicos con madre no guardada, y la señal en vuelo
  printf("\n== F2. Qué hay bajo el pico tras la cadena salvo la ventana (fondo escalado; candidatos en [205, 240]) ==\n");
  for(int pop = 1; pop >= 0; pop--){
    const std::vector<Step>& st = pop ? stIn : stOut;
    const int iWin = pop ? 4 : 6;
    std::map<std::string, double> m; double S = 0, inFlight = 0;
    for(const C& c : cands){
      if(c.inSlice != (bool)pop || !Pass(st, c, iWin) || c.pRange < 205 || c.pRange > 240)continue;
      if(c.group == 0){ S++; continue; }
      if(c.group == -1){ inFlight++; continue; }
      std::string k;
      if(c.tPDG == -13 && c.tProcStart == 2 && c.tOrigin != 2)k = c.tParentPDG == 321 ? "mu+ de desintegración, madre K+ (no señal)" : (c.tParentPDG == 0 ? "mu+ de desintegración, madre no guardada" : Form("mu+ de desintegración, madre %d", c.tParentPDG));
      else if(std::abs(c.tPDG) == 211)k = "pi";
      else if(std::abs(c.tPDG) == 2212)k = "p";
      else if(std::abs(c.tPDG) == 13 && c.tOrigin == 2)k = "mu cósmico";
      else if(std::abs(c.tPDG) == 13)k = Form("mu%s no cósmico, otro (proc %d, madre %d)", c.tPDG < 0 ? "+" : "-", c.tProcStart, c.tParentPDG);
      else k = "resto";
      m[k] += kSc;
    }
    printf("  %-7s señal %3.0f, señal con K en vuelo (fuera de S y B) %3.0f\n", pop ? "dentro" : "fuera", S, inFlight);
    std::vector<std::pair<double, std::string>> v; for(auto& kv : m)v.push_back({kv.second, kv.first});
    std::sort(v.rbegin(), v.rend());
    for(auto& x : v)printf("      %-60s %6.1f\n", x.second.c_str(), x.first);
  }

  // G2: dónde se pierde longitud (proyectando en z: fracción del recorrido true)
  printf("\n== G2. Dónde se pierde la longitud: desplazamiento de cada extremo reco respecto al true a lo largo del mu (cm, + = hacia dentro) ==\n");
  TFile f(in); TTree* t = (TTree*)f.Get("cand");
  const char* base = "tSignal&&tKMomEnd<50&&max(facesA,facesB)>=10&&min(axisA,axisB)<=150&&isShower==0&&pRange>150&&pRange<350&&tLength>0&&abs(tEndZ-tStartZ)>20";
  // bien orientados (la parada reco más cerca del final true en z)
  const char* ok = "abs(stopZ-tEndZ)<abs(stopZ-tStartZ)";
  // a lo largo del mu: s = (z - tStartZ)/(tEndZ - tStartZ) * tLength (aprox. recta)
  Quant(t, "(startZ-tStartZ)/(tEndZ-tStartZ)*tLength", Form("%s&&%s", base, ok), "inicio reco - inicio true (cm; + = empieza más tarde)");
  Quant(t, "(tEndZ-stopZ)/(tEndZ-tStartZ)*tLength", Form("%s&&%s", base, ok), "final true - parada reco (cm; + = acaba antes)");
  Quant(t, "(startZ-tStartZ)/(tEndZ-tStartZ)*tLength", Form("%s&&%s&&kN>0", base, ok), "  inicio, con algo que acaba en el inicio");
  Quant(t, "(startZ-tStartZ)/(tEndZ-tStartZ)*tLength", Form("%s&&%s&&kN==0", base, ok), "  inicio, sin nada en el inicio");
  Quant(t, "(tEndZ-stopZ)/(tEndZ-tStartZ)*tLength", Form("%s&&%s&&arNMichel[1]>=3", base, ok), "  parada, con Michel alrededor");

  // C2: señal tipo shower por población, cortes del módulo y cualquier p
  printf("\n== C2. Señal (mu de K en reposo) con cortes del módulo y p de rango 100-600: traza frente a shower ==\n");
  for(int pop = 1; pop >= 0; pop--){
    const double tr = t->GetEntries(Form("tSignal&&tKMomEnd<50&&max(facesA,facesB)>=10&&min(axisA,axisB)<=150&&inSlice==%d&&isShower==0", pop));
    const double sh = t->GetEntries(Form("tSignal&&tKMomEnd<50&&max(facesA,facesB)>=10&&min(axisA,axisB)<=150&&inSlice==%d&&isShower==1", pop));
    printf("  %-7s traza %4.0f, shower %4.0f (%2.0f %% del total)\n", pop ? "dentro" : "fuera", tr, sh, 100*sh/(tr + sh));
  }
  Quant(t, "length/tLength", "tSignal&&tKMomEnd<50&&max(facesA,facesB)>=10&&min(axisA,axisB)<=150&&isShower==1&&tLength>0", "showers de señal: longitud forzada / true");
  Quant(t, "length/tLength", "tSignal&&tKMomEnd<50&&max(facesA,facesB)>=10&&min(axisA,axisB)<=150&&isShower==0&&tLength>0", "trazas de señal: longitud / true");
}

void ReviewKPiPi(const char* in){
  printf("\n== K. pi+ de K+ -> pi+ pi0 entre los candidatos ==\n");
  TFile f(in); TTree* t=(TTree*)f.Get("cand");
  const char* mod="max(facesA,facesB)>=10&&min(axisA,axisB)<=150&&isShower==0";
  const char* kpi="tPDG==211&&tParentPDG==321&&tProcStart==2";
  printf("pi+ de desintegración de K+ (trazas, cortes del módulo): %lld (dentro %lld, fuera %lld); mu de K en reposo (señal): %lld\n",
    t->GetEntries(Form("%s&&%s",mod,kpi)), t->GetEntries(Form("%s&&%s&&inSlice",mod,kpi)), t->GetEntries(Form("%s&&%s&&!inSlice",mod,kpi)), t->GetEntries(Form("%s&&tSignal&&tKMomEnd<50",mod)));
  // espectro de p de rango (hipótesis mu), por 10 MeV/c
  TH1D h("h","",30,100,400), hs("hs","",30,100,400);
  t->Draw("pRange>>h", Form("%s&&%s",mod,kpi), "goff");
  t->Draw("pRange>>hs", Form("%s&&tSignal&&tKMomEnd<50",mod), "goff");
  printf("p de rango (MeV/c), por 10:   pi+ de K+   |   mu de K en reposo\n");
  for(int i=1;i<=30;i++) printf("  %3.0f-%3.0f  %5.0f   %5.0f\n", h.GetBinLowEdge(i), h.GetBinLowEdge(i)+10, h.GetBinContent(i), hs.GetBinContent(i));
  // momento true de esos pi (¿205 MeV/c, de K en reposo?)
  Long64_t n=t->Draw("tMom", Form("%s&&%s",mod,kpi), "goff"); std::vector<double> v(t->GetV1(), t->GetV1()+n); std::sort(v.begin(),v.end());
  long at=0; for(double x: v) if(x>195&&x<215) at++;
  printf("p true de los pi+: mediana %.0f MeV/c; en [195, 215] (K en reposo): %ld de %lld\n", v[v.size()/2], at, n);
  // fracción con Michel alrededor de la parada (pi+ -> mu+ -> e+)
  printf("con arNMichel10 >= 2: pi+ de K %lld de %lld; señal %lld de %lld\n", t->GetEntries(Form("%s&&%s&&arNMichel[1]>=2",mod,kpi)), t->GetEntries(Form("%s&&%s",mod,kpi)),
    t->GetEntries(Form("%s&&tSignal&&tKMomEnd<50&&arNMichel[1]>=2",mod)), t->GetEntries(Form("%s&&tSignal&&tKMomEnd<50",mod)));
}

int mu_cand_review(const char* in){
  ReviewPart1(in);
  ReviewPart2(in);
  ReviewKPiPi(in);
  return 0;
}
