// Variables nuevas del árbol "cand" de mu_candidates_minitrees.C (20-nov: media truncada, actividad alrededor
// de la parada y del inicio, lo que acaba en el inicio (la K), hijas de Pandora, z de los extremos): ¿separan
// la señal (mu+ de K+ en reposo) del fondo, y lo que añaden a la cadena corta?
//
//   root -l -b -q 'test/phase2/mu_cand_newvars.C+("mu_cand_big19v2.root", "<directorio de salida>")'
//
// Para cada población (en la beam slice y fuera: pasadas distintas de Pandora) y dos etapas:
//   - "suelta": cortes del módulo y p de rango en [205, 240] MeV/c;
//   - "cadena": la cadena corta de la población (InSliceSteps / OutSliceSteps de mu_cand_common.h): lo que
//     queda de fondo, para ver si la variable añade algo a los cortes de ahora.
// Separación: AUC de señal frente a fondo (fondo con los pesos de los eventos normales; la escala no cambia el
// AUC), dada como max(AUC, 1 - AUC) con el sentido ("señal >" o "señal <"); los centinelas (< -900, -1 de los
// que no tienen hits) cuentan como el valor más bajo. 0,5 = no separa.
// Además: la señal partida por el origen de la K (tKOrigin; la medida A es 0 y 1) y la AUC de la señal de la
// medida A frente a la de otros orígenes; y la z de los extremos en el pico de ~570 MeV/c (diverters).
// Figuras: <dir>/mu_cand_newvars.pdf, una página por variable (formas normalizadas, arriba en la beam slice,
// abajo fuera; izquierda suelta, derecha cadena).
#include "mu_cand_common.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TSystem.h"
#include <map>
#include <memory>

namespace newvars {
  using mucut::C;

  struct Var {
    std::string name, title;
    std::function<double(const C&)> f;
    int nbins; double xmin, xmax;
  };

  std::vector<Var> Vars(){
    std::vector<Var> v = {
      {"trMean",      "media truncada del dE/dx (MeV/cm)",                         [](const C& c){ return c.trMean; },      60, 0, 6},
      {"trMeanBody",  "media truncada sin 5 cm de cada extremo (MeV/cm)",          [](const C& c){ return c.trMeanBody; },  60, 0, 6},
      {"asN5",        "hits a < 5 cm del inicio",                                  [](const C& c){ return c.asN5; },        40, 0, 40},
      {"asN10",       "hits a < 10 cm del inicio",                                 [](const C& c){ return c.asN10; },       50, 0, 100},
      {"asCharge5",   "carga a < 5 cm del inicio",                                 [](const C& c){ return c.asCharge5; },   50, 0, 20000},
      {"asCharge10",  "carga a < 10 cm del inicio",                                [](const C& c){ return c.asCharge10; },  50, 0, 40000},
      {"kN",          "PFParticles que acaban en el inicio (< 10 cm)",             [](const C& c){ return c.kN; },          10, 0, 10},
      {"kIsTrack",    "lo que acaba en el inicio es traza",                        [](const C& c){ return c.kIsTrack; },    3, 0, 3},
      {"kDEdx",       "dE/dx del que acaba en el inicio, en ese extremo (MeV/cm)", [](const C& c){ return c.kDEdx; },       50, 0, 25},
      {"kLength",     "longitud del que acaba en el inicio (cm)",                  [](const C& c){ return c.kLength; },     50, 0, 200},
      {"kPRange",     "p de rango (K) del que acaba en el inicio (MeV/c)",         [](const C& c){ return c.kPRange; },     50, 0, 1000},
      {"nDauTrk",     "hijas traza (Pandora)",                                     [](const C& c){ return c.nDauTrk; },     8, 0, 8},
      {"nDauShw",     "hijas shower (Pandora)",                                    [](const C& c){ return c.nDauShw; },     8, 0, 8},
      {"dauHits",     "hits de las hijas",                                         [](const C& c){ return c.dauHits; },     50, 0, 500},
      {"dauNearStop", "hijas a < 10 cm de la parada",                              [](const C& c){ return c.dauNearStop; }, 6, 0, 6},
      {"dauMaxMichel","CNN Michel máxima de las hijas",                            [](const C& c){ return c.dauMaxMichel; },50, 0, 1},
      {"cosBeam",     "cos(inicio -> parada, eje del haz)",                        [](const C& c){ return c.cosBeam; },     40, -1, 1},
      {"startZ",      "z del inicio (cm)",                                         [](const C& c){ return c.startZ; },      70, 0, 700},
      {"stopZ",       "z de la parada (cm)",                                       [](const C& c){ return c.stopZ; },       70, 0, 700},
    };
    const char* r[3] = {"5", "10", "20"};
    for(int i = 0; i < 3; i++){
      v.push_back({Form("arN%s", r[i]),        Form("hits a < %s cm de la parada", r[i]),               [i](const C& c){ return c.arN[i]; },        50, 0, i == 2 ? 300. : 100.});
      v.push_back({Form("arNMichel%s", r[i]),  Form("hits Michel (CNN > 0,5) a < %s cm de la parada", r[i]), [i](const C& c){ return c.arNMichel[i]; }, 20, 0, i == 2 ? 40. : 20.});
      v.push_back({Form("arMichel%s", r[i]),   Form("CNN Michel media a < %s cm de la parada", r[i]),   [i](const C& c){ return c.arMichel[i]; },   50, 0, 1});
      v.push_back({Form("arCharge%s", r[i]),   Form("carga a < %s cm de la parada", r[i]),              [i](const C& c){ return c.arCharge[i]; },   50, 0, i == 2 ? 100000. : 40000.});
    }
    return v;
  }

  // AUC por suma de rangos (empates a la mitad), con pesos. Centinelas: el valor más bajo.
  struct AUC { double auc; double validS, validB; };
  AUC ComputeAUC(const std::vector<std::pair<double,double>>& s, const std::vector<std::pair<double,double>>& b){
    std::vector<std::tuple<double,double,int>> all;   // valor, peso, 0 señal / 1 fondo
    double wS = 0, wB = 0, vS = 0, vB = 0;
    auto norm = [](double x){ return x < -900 || x == -1 ? -1e30 : x; };
    for(auto& p : s){ all.emplace_back(norm(p.first), p.second, 0); wS += p.second; if(norm(p.first) > -1e29)vS += p.second; }
    for(auto& p : b){ all.emplace_back(norm(p.first), p.second, 1); wB += p.second; if(norm(p.first) > -1e29)vB += p.second; }
    if(wS <= 0 || wB <= 0)return {0.5, 0, 0};
    std::sort(all.begin(), all.end(), [](const auto& a, const auto& c){ return std::get<0>(a) < std::get<0>(c); });
    // P(señal > fondo) + 0,5 P(igual)
    double below = 0, sum = 0;
    for(size_t i = 0; i < all.size();){
      size_t j = i; double tS = 0, tB = 0;
      while(j < all.size() && std::get<0>(all[j]) == std::get<0>(all[i])){ (std::get<2>(all[j]) ? tB : tS) += std::get<1>(all[j]); j++; }
      sum += tS*(below + 0.5*tB);
      below += tB; i = j;
    }
    return {sum/(wS*wB), vS/wS, vB/wB};
  }

  std::string Fmt(const AUC& a){
    const bool up = a.auc >= 0.5;
    return Form("%.3f %s", up ? a.auc : 1 - a.auc, up ? ">" : "<");
  }

  // ROOT (TLatex) no entiende UTF-8 en títulos y leyendas
  std::string Latex(std::string t){
    const std::pair<const char*, const char*> m[] = {{"á","#acute{a}"},{"é","#acute{e}"},{"í","#acute{i}"},{"ó","#acute{o}"},{"ú","#acute{u}"},{"ñ","#tilde{n}"},{"→","#rightarrow"}};
    for(auto& kv : m)for(size_t pos; (pos = t.find(kv.first)) != std::string::npos;)t.replace(pos, std::string(kv.first).size(), kv.second);
    return t;
  }

  bool IsA(const C& c){ return c.tKOrigin == 0 || c.tKOrigin == 1; }   // señal de la medida A
}

int mu_cand_newvars(const char* in, const char* outdir = ".", long nEvents = 0, long nEventsNormal = 0){
  using namespace mucut; using namespace newvars;
  gStyle->SetOptStat(0);
  std::vector<C> cands = Load(in);
  if(cands.empty())return 1;
  // Escala del fondo (solo para los números por evento). Por defecto la de la muestra entera (73.070 / 65.977)
  const double scaleB = (nEvents > 0 && nEventsNormal > 0) ? double(nEvents)/nEventsNormal : 73070./65977.;

  const std::vector<Var> vars = Vars();
  const std::vector<Step> steps[2] = {OutSliceSteps(), InSliceSteps()};
  const char* popName[2] = {"fuera de la beam slice", "en la beam slice"};
  const char* stageName[2] = {"suelta (módulo + [205, 240])", "cadena corta"};

  auto stagePass = [&](const C& c, int pop, int stage){
    if(c.inSlice != (bool)pop)return false;
    if(stage == 0)return c.pRange >= 205 && c.pRange <= 240;
    return Pass(steps[pop], c);
  };

  // Recuentos
  std::cout << std::endl << "Candidatos (señal = mu+ de K+ en reposo, sin escalar; fondo escalado x" << scaleB << "):" << std::endl;
  for(int pop = 1; pop >= 0; pop--)for(int st = 0; st < 2; st++){
    double S = 0, SA = 0, B[kNGroups] = {0};
    std::map<int, int> orig;
    for(const C& c : cands){
      if(c.group < 0 || !stagePass(c, pop, st))continue;
      if(c.group == 0){ S++; if(IsA(c))SA++; orig[c.tKOrigin]++; }
      else B[c.group] += scaleB;
    }
    double Bt = 0; for(int g = 1; g < kNGroups; g++)Bt += B[g];
    std::cout << Form("  %-24s %-30s señal %5.0f (medida A %5.0f)  fondo %8.0f  [", popName[pop], stageName[st], S, SA, Bt);
    for(int g = 1; g < kNGroups; g++)std::cout << Form(" %s %.0f", kGroupName[g], B[g]);
    std::cout << " ]" << std::endl << "      origen de la K de la señal (tKOrigin):";
    for(auto& kv : orig)std::cout << Form(" %d: %d", kv.first, kv.second);
    std::cout << std::endl;
  }
  std::cout << "  tKOrigin: 0 primera interacción, 1 primera con K+ -> K+, 2 secundaria, 3 K+ del haz, 4 cadena cortada," << std::endl
            << "            6 otra raíz, -999 rama ausente" << std::endl;

  // AUCs
  std::cout << std::endl << "Separación (max(AUC, 1 - AUC) y sentido de la señal; 0,5 = nada). S/B: señal frente a todo el fondo;"
            << std::endl << "A/otros: señal de la medida A frente a señal de otros orígenes (suelta). válidos: fracción sin centinela (S, B)." << std::endl;
  for(int pop = 1; pop >= 0; pop--){
    std::cout << std::endl << "== " << popName[pop] << " ==" << std::endl
              << Form("  %-14s %-12s %-12s %-12s %-12s %-12s %-12s   %s", "variable", "S/B suelta", "S/B cadena", "S/pi suelta", "S/p suelta", "S/mu cós.", "A/otros", "válidos (suelta)") << std::endl;
    for(const Var& v : vars){
      std::vector<std::pair<double,double>> s[2], b[2], bpi, bp, bmu, sA, sO;
      for(const C& c : cands){
        if(c.group < 0)continue;
        for(int st = 0; st < 2; st++){
          if(!stagePass(c, pop, st))continue;
          const double x = v.f(c);
          if(c.group == 0){ s[st].push_back({x, 1}); if(st == 0)(IsA(c) ? sA : sO).push_back({x, 1}); }
          else {
            b[st].push_back({x, scaleB});
            if(st == 0){ if(c.group == 3)bpi.push_back({x, 1}); if(c.group == 2)bp.push_back({x, 1}); if(c.group == 1)bmu.push_back({x, 1}); }
          }
        }
      }
      const AUC a0 = ComputeAUC(s[0], b[0]), a1 = ComputeAUC(s[1], b[1]);
      std::cout << Form("  %-14s %-12s %-12s %-12s %-12s %-12s %-12s   %.2f, %.2f", v.name.c_str(), Fmt(a0).c_str(), Fmt(a1).c_str(),
                        Fmt(ComputeAUC(s[0], bpi)).c_str(), Fmt(ComputeAUC(s[0], bp)).c_str(), Fmt(ComputeAUC(s[0], bmu)).c_str(),
                        Fmt(ComputeAUC(sA, sO)).c_str(), a0.validS, a0.validB) << std::endl;
    }
  }

  // Diverters: z de los extremos en el pico de ~570 MeV/c frente a la zona de al lado
  std::cout << std::endl << "z de los extremos (cortes del módulo, todos los candidatos con truth de fondo), por 20 cm:" << std::endl;
  for(int pop = 1; pop >= 0; pop--){
    for(int w = 0; w < 2; w++){
      const double lo = w == 0 ? 540 : 420, hi = w == 0 ? 600 : 520;
      TH1D hz("hz", "", 35, 0, 700);
      long n = 0, nNear = 0;
      for(const C& c : cands){
        if(c.group < 1 || c.inSlice != (bool)pop || c.pRange < lo || c.pRange > hi)continue;
        n++;
        const double zs[2] = {c.startZ, c.stopZ};
        bool near = false;
        for(double z : zs){ if(z > -900)hz.Fill(z); if(std::abs(z - 220) < 15 || std::abs(z - 440) < 15)near = true; }
        if(near)nNear++;
      }
      std::cout << Form("  %-24s p de rango [%3.0f, %3.0f]: %ld candidatos, con un extremo a < 15 cm de z = 220 o 440: %.1f %%",
                        popName[pop], lo, hi, n, n ? 100.*nNear/n : 0.) << std::endl << "    ";
      for(int i = 1; i <= hz.GetNbinsX(); i++)std::cout << Form("%.0f ", hz.GetBinContent(i));
      std::cout << std::endl;
    }
  }

  // Figuras: formas normalizadas de señal y de cada grupo de fondo
  gSystem->mkdir(outdir, true);
  const std::string pdf = std::string(outdir) + "/mu_cand_newvars.pdf";
  TCanvas cv("cv", "", 1200, 900);
  cv.Print((pdf + "[").c_str());
  const int col[kNGroups] = {kRed + 1, kGray + 2, kBlue + 1, kGreen + 2, kOrange + 1, kMagenta + 1};
  for(const Var& v : vars){
    cv.Clear(); cv.Divide(2, 2);
    std::vector<std::unique_ptr<TH1D>> keep;
    for(int pop = 1, pad = 1; pop >= 0; pop--)for(int st = 0; st < 2; st++, pad++){
      cv.cd(pad);
      std::vector<TH1D*> h(kNGroups);
      for(int g = 0; g < kNGroups; g++){
        keep.emplace_back(new TH1D(Form("h_%s_%d_%d_%d", v.name.c_str(), pop, st, g), Latex(Form("%s, %s;%s;fracci#acute{o}n", popName[pop], stageName[st], v.title.c_str())).c_str(), v.nbins, v.xmin, v.xmax));
        h[g] = keep.back().get(); h[g]->SetLineColor(col[g]); h[g]->SetLineWidth(g == 0 ? 3 : 2);
      }
      for(const C& c : cands){
        if(c.group < 0 || !stagePass(c, pop, st))continue;
        const double x = v.f(c);
        if(x < -900 || x == -1)continue;
        h[c.group]->Fill(std::min(x, v.xmax - 1e-6*(v.xmax - v.xmin)));
      }
      double ymax = 0;
      for(int g = 0; g < kNGroups; g++){ if(h[g]->Integral() > 0)h[g]->Scale(1./h[g]->Integral()); ymax = std::max(ymax, h[g]->GetMaximum()); }
      h[0]->SetMaximum(1.25*ymax); h[0]->Draw("hist");
      for(int g = 1; g < kNGroups; g++)if(h[g]->GetEntries() >= 10)h[g]->Draw("hist same");
      h[0]->Draw("hist same");
      if(pad == 1){
        TLegend* l = new TLegend(0.6, 0.6, 0.89, 0.89);
        for(int g = 0; g < kNGroups; g++)l->AddEntry(h[g], Latex(kGroupName[g]).c_str(), "l");
        l->Draw();
      }
    }
    cv.Print(pdf.c_str());
  }
  cv.Print((pdf + "]").c_str());
  std::cout << std::endl << "Figuras en " << pdf << std::endl;
  return 0;
}
