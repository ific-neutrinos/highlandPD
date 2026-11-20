// Histogramas apilados por categoría de los candidatos a mu de K+ en reposo (árbol "cand" de
// mu_candidates_minitrees.C) y de los cortes de mucut::DefaultSteps() (mu_cand_common.h).
//
//   root -l -b -q 'test/phase2/mu_cand_plots.C+("mu_cand_big19.root", "<directorio de salida>")'
//
// Fondo escalado al número total de eventos (solo eventos normales, como en mu_cand_cutflow.C); señal sin
// escalar, apilada arriba y además superpuesta multiplicada por un factor (línea discontinua) para ver su
// forma. Los valores por debajo del rango (centinelas: -1 sin hits) no se dibujan; el overflow va al
// último bin. Páginas (PDF con todas y un PNG por página):
//   1-3. p de rango con los cortes del módulo salvo la ventana: todos, en la beam slice y fuera.
//        Fuera de la beam slice, fuera de la ventana del módulo (180-300) solo hay resumen: el escalón es
//        del módulo, no de la física.
//   4-16. variables candidatas con los cortes del módulo: apiladas y formas normalizadas por categoría.
//   17... un corte por página: a la izquierda con los cortes anteriores (cadena), a la derecha con todos
//        los demás (N-1). Línea vertical donde se corta. Base: cortes del módulo, ventana incluida
//        (salvo en la página de la ventana estrecha, donde se quita).
#include "mu_cand_common.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TSystem.h"
#include <memory>

namespace muplot {
  using namespace mucut;

  const int kColor[kNGroups] = {kRed + 1, kAzure + 1, kGreen + 2, kOrange + 1, kViolet - 4, kGray + 1};

  // ROOT no dibuja bien los caracteres UTF-8: acentos en TLatex
  std::string Tx(std::string t){
    const char* m[][2] = {{"á", "#acute{a}"}, {"é", "#acute{e}"}, {"í", "#acute{i}"}, {"ó", "#acute{o}"}, {"ú", "#acute{u}"}, {"ñ", "#tilde{n}"}};
    for(auto& r : m)for(size_t pos; (pos = t.find(r[0])) != std::string::npos;)t.replace(pos, std::string(r[0]).size(), r[1]);
    return t;
  }

  struct Drawn {   // objetos de una figura (para que vivan hasta guardar)
    std::vector<std::unique_ptr<TObject>> keep;
  };

  // Rellena y dibuja en el pad actual. sel: qué candidatos entran
  void Draw(const std::vector<C>& cands, const std::function<bool(const C&)>& sel, const std::function<double(const C&)>& var,
            int nbins, double xmin, double xmax, const std::string& title, const std::string& xtitle,
            const std::vector<double>& lines, bool logy, double scaleB, Drawn& d){
    static int id = 0; id++;
    TH1D* h[kNGroups];
    for(int g = 0; g < kNGroups; g++){
      h[g] = new TH1D(Form("h%d_%d", id, g), "", nbins, xmin, xmax);
      h[g]->SetFillColor(kColor[g]); h[g]->SetLineColor(kColor[g]);
      d.keep.emplace_back(h[g]);
    }
    for(const C& c : cands){
      if(c.group < 0 || !sel(c))continue;
      double v = var(c);
      if(v < xmin)continue;   // valores centinela (-1 sin hits, -999): fuera
      v = std::min(xmax - 1e-6*(xmax - xmin), v);   // overflow al último bin
      h[c.group]->Fill(v, c.group == 0 ? 1. : scaleB);
    }
    THStack* st = new THStack(Form("st%d", id), Form("%s;%s;candidatos", Tx(title).c_str(), Tx(xtitle).c_str()));
    d.keep.emplace_back(st);
    for(int g = kNGroups - 1; g >= 0; g--)st->Add(h[g]);
    st->Draw("hist");
    const double ymax = st->GetMaximum();
    // señal superpuesta, por un factor redondo
    TH1D* hs = (TH1D*)h[0]->Clone(Form("hs%d", id)); d.keep.emplace_back(hs);
    double k = 1;
    if(h[0]->GetMaximum() > 0){
      const double raw = 0.5*ymax/h[0]->GetMaximum();
      const double nice[] = {1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000};
      for(double n : nice)if(n <= raw)k = n;
    }
    hs->Scale(k); hs->SetFillStyle(0); hs->SetLineColor(kRed + 2); hs->SetLineWidth(3); hs->SetLineStyle(2);
    hs->Draw("hist same");
    if(logy){ gPad->SetLogy(); st->SetMinimum(0.5); st->SetMaximum(5*ymax); }
    else st->SetMaximum(1.25*ymax);
    gPad->Update();
    const double y1 = logy ? std::pow(10, gPad->GetUymin()) : gPad->GetUymin();
    const double y2 = logy ? std::pow(10, gPad->GetUymax()) : gPad->GetUymax();
    for(double x : lines){
      TLine* l = new TLine(x, y1, x, y2); l->SetLineColor(kBlack); l->SetLineWidth(3); l->Draw(); d.keep.emplace_back(l);
    }
    TLegend* leg = new TLegend(0.60, 0.58, 0.89, 0.89); d.keep.emplace_back(leg);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.035);
    double B = 0; for(int g = 1; g < kNGroups; g++)B += h[g]->Integral();
    for(int g = 0; g < kNGroups; g++)leg->AddEntry(h[g], Form("%s (%.0f)", Tx(kGroupName[g]).c_str(), h[g]->Integral()), "f");
    leg->AddEntry(hs, Form("%s x%g", Tx("señal").c_str(), k), "l");
    leg->Draw();
    TLatex* tx = new TLatex(); d.keep.emplace_back(tx);
    tx->SetNDC(); tx->SetTextSize(0.033);
    tx->DrawLatex(0.13, 0.85, Form("S = %.0f, B = %.0f, S/B = %.3f", h[0]->Integral(), B, B > 0 ? h[0]->Integral()/B : 0.));
  }

  // Formas normalizadas a área 1 por categoría (líneas), para ver si una variable separa
  void DrawShapes(const std::vector<C>& cands, const std::function<bool(const C&)>& sel, const std::function<double(const C&)>& var,
                  int nbins, double xmin, double xmax, const std::string& title, const std::string& xtitle, bool logy, Drawn& d){
    static int id = 0; id++;
    TH1D* h[kNGroups];
    for(int g = 0; g < kNGroups; g++){
      h[g] = new TH1D(Form("n%d_%d", id, g), Form("%s;%s;fracci#acute{o}n", Tx(title).c_str(), Tx(xtitle).c_str()), nbins, xmin, xmax);
      h[g]->SetLineColor(kColor[g]); h[g]->SetLineWidth(g == 0 ? 4 : 2); h[g]->SetFillStyle(0);
      d.keep.emplace_back(h[g]);
    }
    for(const C& c : cands){
      if(c.group < 0 || !sel(c))continue;
      double v = var(c);
      if(v < xmin)continue;
      v = std::min(xmax - 1e-6*(xmax - xmin), v);
      h[c.group]->Fill(v);
    }
    double ymax = 0;
    for(int g = 0; g < kNGroups; g++){ if(h[g]->Integral() > 0)h[g]->Scale(1./h[g]->Integral()); ymax = std::max(ymax, h[g]->GetMaximum()); }
    h[0]->SetMaximum(logy ? 3*ymax : 1.15*ymax);
    if(logy){ h[0]->SetMinimum(1e-4); gPad->SetLogy(); }
    h[0]->Draw("hist");
    for(int g = kNGroups - 1; g >= 1; g--)h[g]->Draw("hist same");
    h[0]->Draw("hist same");
    TLegend* leg = new TLegend(0.65, 0.65, 0.89, 0.89); d.keep.emplace_back(leg);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.035);
    for(int g = 0; g < kNGroups; g++)leg->AddEntry(h[g], Tx(kGroupName[g]).c_str(), "l");
    leg->Draw();
  }
}

int mu_cand_plots(const char* in, const char* outdir, long nEvents = 73070, long nEventsNormal = 65977){
  using namespace mucut; using namespace muplot;
  const std::vector<C> cands = Load(in);
  if(cands.empty())return 1;
  const double scaleB = double(nEvents)/nEventsNormal;
  const std::vector<Step> steps = DefaultSteps();
  gSystem->mkdir(outdir, true);
  gStyle->SetOptStat(0); gStyle->SetTitleFontSize(0.04);
  const std::string pdf = std::string(outdir) + "/mu_cand_plots.pdf";
  int page = 0;
  TCanvas cv("cv", "", 1400, 650);
  cv.Print((pdf + "[").c_str());
  auto Save = [&](const std::string& tag){
    page++;
    cv.Print(pdf.c_str());
    cv.Print(Form("%s/%02d_%s.png", outdir, page, tag.c_str()));
  };
  auto pRange = [](const C& c){ return c.pRange; };

  // 1-3: p de rango inicial
  {
    struct S { const char* tag; const char* title; std::function<bool(const C&)> sel; };
    const S sets[] = {{"prange_todos", "Inicial: cortes del módulo salvo la ventana, todos", [](const C&){ return true; }},
                      {"prange_bs", "Inicial: en la beam slice", [](const C& c){ return c.inSlice; }},
                      {"prange_fuera", "Inicial: fuera de la beam slice (fuera de 180-300 solo resumen)", [](const C& c){ return !c.inSlice; }}};
    for(const S& s : sets){
      cv.Clear(); cv.Divide(2, 1);
      Drawn d;
      cv.cd(1); Draw(cands, s.sel, pRange, 50, 100, 600, s.title, "p por rango (hipótesis #mu, MeV/c)", {180, 300}, false, scaleB, d);
      cv.cd(2); Draw(cands, s.sel, pRange, 50, 100, 600, "(escala log)", "p por rango (hipótesis #mu, MeV/c)", {180, 300}, true, scaleB, d);
      Save(s.tag);
    }
  }

  // Variables candidatas (sin cortar): con los cortes del módulo (ventana incluida)
  {
    struct V { const char* tag; const char* name; std::function<double(const C&)> var; int nb; double x0, x1; bool logy; };
    const V vars[] = {
      {"chi2mu",      "#chi^{2}_{#mu}/ndf (parada en el extremo de mayor dE/dx)", [](const C& c){ return c.chi2muStop; }, 60, 0, 60, false},
      {"nhijas",      "n#acute{u}mero de hijas (Pandora)", [](const C& c){ return double(c.nDaughters); }, 8, -0.5, 7.5, true},
      {"caras_min",   "distancia m#acute{i}nima de los dos extremos a las caras (cm)", [](const C& c){ return std::min(c.facesStop, c.facesStart); }, 50, 0, 150, false},
      {"caras_inicio","distancia del extremo de inicio (no parada) a las caras (cm)", [](const C& c){ return c.facesStart; }, 50, 0, 150, false},
      {"caras_parada","distancia del extremo de parada a las caras (cm)", [](const C& c){ return c.facesStop; }, 50, 0, 150, false},
      {"eje",         "distancia m#acute{i}nima al eje del haz (cm)", [](const C& c){ return c.axisMin; }, 50, 0, 150, false},
      {"diry",        "componente y de la direcci#acute{o}n hacia la parada", [](const C& c){ return c.dirYStop; }, 40, -1, 1, false},
      {"dedx_parada", "dE/dx medio en los 5 cm del extremo de parada (MeV/cm)", [](const C& c){ return c.dEdxStop; }, 60, 0, 20, false},
      {"dedx_inicio", "dE/dx medio en los 5 cm del extremo de inicio (MeV/cm)", [](const C& c){ return c.dEdxStart; }, 60, 0, 10, false},
      {"michel_propios", "CNN Michel media de los hits propios a < 5 cm de la parada", [](const C& c){ return c.michelStop; }, 50, 0, 1, true},
      {"act_parada",  "part#acute{i}culas con un extremo a < 10 cm de la parada", [](const C& c){ return double(c.actStop); }, 8, -0.5, 7.5, true},
      {"act_michel",  "CNN Michel media de esas part#acute{i}culas (-1: ninguna)", [](const C& c){ return c.actMichelStop; }, 51, -1.04, 1, true},
      {"cnnmichel",   "CNN Michel de la part#acute{i}cula entera", [](const C& c){ return c.cnnMichel; }, 50, 0, 0.5, true},
    };
    auto sel = [](const C& c){ return InWindow(c); };
    for(const V& v : vars){
      cv.Clear(); cv.Divide(2, 1);
      Drawn d;
      cv.cd(1); Draw(cands, sel, v.var, v.nb, v.x0, v.x1, "Cortes del m#acute{o}dulo", v.name, {}, v.logy, scaleB, d);
      cv.cd(2); DrawShapes(cands, sel, v.var, v.nb, v.x0, v.x1, "Formas (#acute{a}rea 1)", v.name, v.logy, d);
      Save(std::string("var_") + v.tag);
    }
  }

  // Un corte por página
  for(int i = 0; i < (int)steps.size(); i++){
    const Step& s = steps[i];
    const bool isWindow = s.name.rfind("p rango", 0) == 0;
    auto base = [isWindow](const C& c){ return isWindow || InWindow(c); };
    auto before = [&, i](const C& c){ return base(c) && Pass(steps, c, -1, i - 1); };
    auto nm1 = [&, i](const C& c){ return base(c) && Pass(steps, c, i); };
    cv.Clear(); cv.Divide(2, 1);
    Drawn d;
    cv.cd(1); Draw(cands, before, s.var, s.nbins, s.xmin, s.xmax, Form("%d. %s: con los cortes anteriores", i + 1, s.name.c_str()), s.xtitle, s.lines, s.logy, scaleB, d);
    cv.cd(2); Draw(cands, nm1, s.var, s.nbins, s.xmin, s.xmax, Form("%d. %s: N-1 (todos los demás)", i + 1, s.name.c_str()), s.xtitle, s.lines, s.logy, scaleB, d);
    Save(Form("corte%d", i + 1));
  }
  cv.Print((pdf + "]").c_str());
  std::cout << page << " páginas en " << pdf << " (y PNG en " << outdir << ")" << std::endl;
  return 0;
}
