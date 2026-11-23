// Implementación de pimu_endcharge.C (se compila con ACLiC).
//
// ¿Separa el exceso de carga en la parada un pi+ que se para (pi+ -> mu+ nu: el mu+ de 4,1 MeV recorre ~1,4 mm y deja
// su energía en el último hit o el de al lado) de un mu que se para (sin ese depósito)? Los dos dan Michel.
// Trazas de Pandora emparejadas (completitud > kMinEff) con: mu true que se paran (|PDG| 13, momento final < 10 MeV/c) o
// pi+ true que se paran y se desintegran (PDG 211, Decay, momento final < 10 MeV/c). La parada reco es el extremo más
// cercano al final true (a < kMaxEndDist cm). Con los hits de calorimetría del plano 2 ordenados por distancia a la parada
// (residual range respecto a ese extremo):
//   r1   = carga del último hit / mediana de la carga de los hits a 2-6 cm de la parada
//   r2   = (carga de los dos últimos hits) / (2 x esa mediana)
//   rmax = máximo de la carga de los hits a < 2 cm de la parada / esa mediana
//   dEdxLast = dE/dx del último hit
// Para cada variable: mediana y cuantiles 16/84 de cada clase y AUC (pi frente a mu).
#include "kaon_fate_minitrees_impl.C"
#include <algorithm>
#include <map>

namespace pimu {
  using kfate::InputFiles; using kfate::TreeName;
  bool Valid(const Float_t* x){ return x[0] > -900 && x[1] > -900 && x[2] > -900; }
  // completitud mínima y distancia máxima (cm) entre la parada reco y el final true (sin corregir el SCE, que desplaza
  // las posiciones reco hasta ~10 cm)
  const double kMinEff = 0.5, kMaxEndDist = 15.;
  struct Vals { std::vector<double> r1, r2, rmax, dedx; };
  Vals cls[2];   // 0 mu, 1 pi+
  long nCand[2] = {0, 0};

  double Median(std::vector<double> v){ if(v.empty())return -1; std::sort(v.begin(), v.end()); return v[v.size()/2]; }

  void Process(AnaSpillPD* spill){
    const AnaEventInfoPD* info = static_cast<AnaEventInfoPD*>(spill->EventInfo);
    if(!info || !info->IsMC || spill->Bunches.empty())return;
    for(AnaParticleB* pb : static_cast<AnaBunchPD*>(spill->Bunches[0])->Particles){
      const AnaParticlePD* p = static_cast<AnaParticlePD*>(pb);
      if(p->Type != AnaParticlePD::kTrack || p->Hits[2].size() < 10 || !Valid(p->PositionStart) || !Valid(p->PositionEnd))continue;
      const AnaTrueParticleB* t = static_cast<const AnaTrueParticleB*>(p->TrueObject);
      if(!t || p->TrueEff < kMinEff || t->MomentumEnd > 0.010)continue;
      int c = -1;
      if(std::abs(t->PDG) == 13)c = 0;
      else if(t->PDG == 211 && t->ProcessEnd == AnaTrueParticleB::Decay)c = 1;
      if(c < 0)continue;
      // extremo de parada: el más cercano al final true
      const TVector3 te(t->PositionEnd[0], t->PositionEnd[1], t->PositionEnd[2]);
      const TVector3 A(p->PositionStart[0], p->PositionStart[1], p->PositionStart[2]), B(p->PositionEnd[0], p->PositionEnd[1], p->PositionEnd[2]);
      const double dA = (A - te).Mag(), dB = (B - te).Mag();
      if(std::min(dA, dB) > kMaxEndDist)continue;
      const bool stopB = dB <= dA;   // residual range de Pandora: 0 en PositionEnd
      std::vector<std::pair<double, const AnaHitPD*>> hs;
      for(const AnaHitPD& h : p->Hits[2]){
        if(h.ResidualRange < 0 || h.Integral <= 0)continue;
        const double rr = stopB ? h.ResidualRange : p->Length - h.ResidualRange;
        hs.push_back({rr, &h});
      }
      if(hs.size() < 10)continue;
      std::sort(hs.begin(), hs.end(), [](const auto& x, const auto& y){ return x.first < y.first; });
      std::vector<double> ref;
      for(const auto& x : hs)if(x.first > 2 && x.first < 6)ref.push_back(x.second->Integral);
      if(ref.size() < 3)continue;
      const double qref = Median(ref);
      if(qref <= 0)continue;
      double qmax = 0;
      for(const auto& x : hs)if(x.first < 2)qmax = std::max(qmax, double(x.second->Integral));
      nCand[c]++;
      cls[c].r1.push_back(hs[0].second->Integral/qref);
      cls[c].r2.push_back((hs[0].second->Integral + hs[1].second->Integral)/(2*qref));
      cls[c].rmax.push_back(qmax/qref);
      cls[c].dedx.push_back(hs[0].second->dEdx);
    }
  }

  double AUC(std::vector<double> s, std::vector<double> b){   // P(pi > mu)
    if(s.empty() || b.empty())return 0.5;
    std::sort(b.begin(), b.end());
    double sum = 0;
    for(double x : s){
      const double lo = std::lower_bound(b.begin(), b.end(), x) - b.begin(), hi = std::upper_bound(b.begin(), b.end(), x) - b.begin();
      sum += lo + 0.5*(hi - lo);
    }
    return sum/(double(s.size())*b.size());
  }
  void Print(const char* name, const std::vector<double>& mu, const std::vector<double>& pi){
    auto q = [](std::vector<double> v, double f){ if(v.empty())return -1.; std::sort(v.begin(), v.end()); return v[std::min(v.size() - 1, size_t(f*v.size()))]; };
    std::cout << Form("  %-48s mu %5.2f [%5.2f, %5.2f]   pi+ %5.2f [%5.2f, %5.2f]   AUC (pi > mu) %.3f", name,
                      q(mu, .5), q(mu, .16), q(mu, .84), q(pi, .5), q(pi, .16), q(pi, .84), AUC(pi, mu)) << std::endl;
  }
}

int pimu_endcharge_impl(const char* input, Long64_t maxEntries){
  using namespace pimu;
  const std::vector<std::string> files = InputFiles(input);
  if(files.empty()){ std::cerr << "No hay ficheros en " << input << std::endl; return -1; }
  const std::string treeName = TreeName(TString(files[0]).Contains("*") ? std::string("") : files[0]);
  TChain chain(treeName.empty() ? "highlandana/MiniTree" : treeName.c_str());
  for(const std::string& f : files)chain.Add(f.c_str());
  if(chain.LoadTree(0) < 0){ std::cerr << "No hay entradas" << std::endl; return -1; }
  AnaSpillPD* spill = nullptr;
  chain.SetBranchAddress("Spill", &spill);
  LoopMiniTree(chain, spill, maxEntries, files.size(), [](AnaSpillPD* s){ Process(s); });
  std::cout << std::endl << "Trazas que se paran, emparejadas: mu " << nCand[0] << ", pi+ (pi -> mu nu en reposo) " << nCand[1] << std::endl
            << "  variable                                         mediana [16 %, 84 %]" << std::endl;
  Print("carga del último hit / referencia (2-6 cm)", cls[0].r1, cls[1].r1);
  Print("carga de los dos últimos / (2 x referencia)", cls[0].r2, cls[1].r2);
  Print("carga máxima a < 2 cm / referencia", cls[0].rmax, cls[1].rmax);
  Print("dE/dx del último hit (MeV/cm)", cls[0].dedx, cls[1].dedx);
  return 0;
}
