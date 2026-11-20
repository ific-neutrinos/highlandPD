// Escaneo de umbrales de la cadena de cortes de los candidatos a mu de K+ en reposo, por población
// (árbol "cand" de mu_candidates_minitrees.C; ver mu_cand_common.h).
//
//   root -l -b -q 'test/phase2/mu_cand_scan.C+("mu_cand_big19.root")'
//
// Para cada población: cortes fijos (los que el N-1 dice que sirven) y rejilla en chi2 mu, cnnMichel y la
// ventana de p de rango. Imprime la tabla en chi2 mu (resto en el valor de referencia), la de cnnMichel y
// las 15 mejores combinaciones por S/sqrt(S+B). Ojo: optimizar en la misma muestra en la que se mide es
// optimista; con tantos puntos, las diferencias de < ~0,2 en S/sqrt(S+B) no significan nada.
#include "mu_cand_common.h"
#include <tuple>

namespace muscan {
  using namespace mucut;

  struct Res { double S, B; double Z() const { return S + B > 0 ? S/std::sqrt(S + B) : 0; } };

  struct Scan {
    const std::vector<C>* cands; double scaleB;
    std::function<bool(const C&)> pop, fixed;
    Res Eval(double chi2mu, double michel, double w1, double w2) const {
      Res r{0, 0};
      for(const C& c : *cands){
        if(c.group < 0 || !pop(c) || !InWindow(c) || !fixed(c))continue;
        if(c.chi2muStop < 0 || c.chi2muStop >= chi2mu || c.cnnMichel <= michel || c.pRange < w1 || c.pRange > w2)continue;
        if(c.group == 0)r.S++; else r.B += scaleB;
      }
      return r;
    }
  };

  void Run(const Scan& sc, const std::string& label){
    std::cout << std::endl << "== " << label << " ==" << std::endl;
    std::cout << "chi2 mu (cnnMichel > 0.02, ventana [200, 245]):" << std::endl;
    for(double a : {3., 4., 5., 6., 7., 8., 10., 12., 15., 1e9}){
      const Res r = sc.Eval(a, 0.02, 200, 245);
      std::cout << Form("  chi2mu < %-6s S %5.0f  B %7.0f  S/B %6.3f  S/√(S+B) %5.2f", a > 1e8 ? "inf" : Form("%g", a), r.S, r.B, r.B > 0 ? r.S/r.B : 0, r.Z()) << std::endl;
    }
    std::cout << "cnnMichel (chi2 mu < 5, ventana [200, 245]):" << std::endl;
    for(double m : {-1., 0.005, 0.01, 0.02, 0.03, 0.05, 0.08}){
      const Res r = sc.Eval(5, m, 200, 245);
      std::cout << Form("  cnnMichel > %-6g S %5.0f  B %7.0f  S/B %6.3f  S/√(S+B) %5.2f", m, r.S, r.B, r.B > 0 ? r.S/r.B : 0, r.Z()) << std::endl;
    }
    std::cout << "ventana (chi2 mu < 5, cnnMichel > 0.02):" << std::endl;
    for(auto w : std::vector<std::pair<double,double>>{{180, 300}, {190, 260}, {195, 250}, {200, 245}, {205, 245}, {205, 240}, {210, 240}}){
      const Res r = sc.Eval(5, 0.02, w.first, w.second);
      std::cout << Form("  [%3.0f, %3.0f]  S %5.0f  B %7.0f  S/B %6.3f  S/√(S+B) %5.2f", w.first, w.second, r.S, r.B, r.B > 0 ? r.S/r.B : 0, r.Z()) << std::endl;
    }
    std::vector<std::tuple<double, double, double, double, double, Res>> all;
    for(double a : {3., 4., 5., 6., 7., 8., 10., 12.})
      for(double m : {-1., 0.01, 0.02, 0.03, 0.05})
        for(auto w : std::vector<std::pair<double,double>>{{195, 250}, {200, 245}, {205, 245}, {205, 240}}){
          const Res r = sc.Eval(a, m, w.first, w.second);
          all.emplace_back(r.Z(), a, m, w.first, w.second, r);
        }
    std::sort(all.begin(), all.end(), [](const auto& x, const auto& y){ return std::get<0>(x) > std::get<0>(y); });
    std::cout << "Mejores 15 por S/√(S+B):" << std::endl;
    for(size_t i = 0; i < 15 && i < all.size(); i++){
      const auto& t = all[i]; const Res& r = std::get<5>(t);
      std::cout << Form("  chi2mu < %-4g Michel > %-5g [%3.0f, %3.0f]  S %5.0f  B %6.0f  S/B %6.3f  S/√(S+B) %5.2f",
                        std::get<1>(t), std::get<2>(t), std::get<3>(t), std::get<4>(t), r.S, r.B, r.B > 0 ? r.S/r.B : 0, r.Z()) << std::endl;
    }
  }
}

int mu_cand_scan(const char* in, long nEvents = 73070, long nEventsNormal = 65977){
  using namespace muscan;
  const std::vector<C> cands = Load(in);
  if(cands.empty())return 1;
  const double scaleB = double(nEvents)/nEventsNormal;
  // Cortes fijos: los que sobreviven al N-1 de la cadena larga
  Scan sIn{&cands, scaleB, [](const C& c){ return c.inSlice; },
          [](const C& c){ return std::min(c.facesStop, c.facesStart) > 35 && (c.vtxMin < 0 || c.vtxMin > 15) && c.dEdxStop > 3; }};
  Run(sIn, "En la beam slice: extremos a > 35 cm de las caras, vértice > 15 cm, dE/dx parada > 3");
  Scan sOut{&cands, scaleB, [](const C& c){ return !c.inSlice; },
           [](const C& c){ return !c.clearCosmic && c.axisMin < 60 && c.dEdxStop > 3 && c.chi2pStop > 80 && c.dEdxStart > 1.8; }};
  Run(sOut, "Fuera: no IsClearCosmic, eje < 60 cm, dE/dx parada > 3, chi2 p > 80, dE/dx inicio > 1.8");
  return 0;
}
