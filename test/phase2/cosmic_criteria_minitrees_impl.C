// Implementación de cosmic_criteria_minitrees.C (se compila con ACLiC).
//
// Criterios candidatos para guardar completos (con hits) los PFParticles de fuera de la beam
// slice, medidos sobre minitrees ya escritos (MC, eventos con beam slice): todos los PFParticles
// de fuera se guardan al menos en resumen, con inicio, fin y longitud, que es lo que usan.
//   (a) cilindro: el segmento inicio-fin del PFParticle pasa a menos de R de la trayectoria de
//       alguna partícula de la beam slice (aproxima el criterio del módulo con los space points)
//   (b) extremos: su inicio o su fin está a menos de R de esa trayectoria
//   (c) rango: momento por rango con hipótesis de muón en una ventana alrededor de 236 MeV/c
//       (mu+ de K+ -> mu+ nu en reposo)
// La trayectoria de una partícula de la beam slice son sus puntos de trayectoria (sin SCE) o, si
// no tiene (showers), el segmento de su inicio a su fin (o a inicio + dirección x longitud).
//
// Para cada criterio:
//   - eficiencia: fracción de los mu de K+ de la cadena del haz emparejados con un PFParticle de
//     fuera que el criterio guardaría completos; y lo mismo para cualquier partícula de la cadena
//     del haz emparejada fuera (más estadística, como referencia)
//   - coste: jerarquías (primario de fuera y sus descendientes; el módulo decide por jerarquía) y
//     PFParticles que hoy van en resumen y pasarían a completos, por evento
#include "kaon_fate_minitrees_impl.C"
#include "pdAnalysisUtils.hxx"
#include "Parameters.hxx"
#include "TVector3.h"
#include <functional>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace ccrit {
  using namespace kfate;

  struct Seg { TVector3 a, b; };

  double PointSeg(const TVector3& p, const Seg& s){
    const TVector3 d = s.b - s.a;
    const double l2 = d.Mag2();
    double t = l2 > 0 ? (p - s.a).Dot(d)/l2 : 0;
    t = std::max(0., std::min(1., t));
    return (p - (s.a + t*d)).Mag();
  }

  // Distancia mínima entre dos segmentos (Ericson, Real-Time Collision Detection, 5.1.9)
  double SegSeg(const Seg& s1, const Seg& s2){
    const TVector3 d1 = s1.b - s1.a, d2 = s2.b - s2.a, r = s1.a - s2.a;
    const double a = d1.Mag2(), e = d2.Mag2(), f = d2.Dot(r);
    double s = 0, t = 0;
    if(a <= 1e-12 && e <= 1e-12)return r.Mag();
    if(a <= 1e-12){ t = std::max(0., std::min(1., f/e)); }
    else{
      const double c = d1.Dot(r);
      if(e <= 1e-12){ s = std::max(0., std::min(1., -c/a)); }
      else{
        const double b = d1.Dot(d2), den = a*e - b*b;
        s = den > 1e-12 ? std::max(0., std::min(1., (b*f - c*e)/den)) : 0;
        t = (b*s + f)/e;
        if(t < 0){ t = 0; s = std::max(0., std::min(1., -c/a)); }
        else if(t > 1){ t = 1; s = std::max(0., std::min(1., (b - c)/a)); }
      }
    }
    return ((s1.a + s*d1) - (s2.a + t*d2)).Mag();
  }

  bool Valid(const Float_t* x){ return x[0] > -900 && x[1] > -900 && x[2] > -900; }

  // Inicio y fin de un PFParticle (los showers no tienen fin: inicio + dirección x longitud)
  bool Ends(const AnaParticlePD* p, TVector3& a, TVector3& b){
    if(!Valid(p->PositionStart))return false;
    a.SetXYZ(p->PositionStart[0], p->PositionStart[1], p->PositionStart[2]);
    if(Valid(p->PositionEnd))b.SetXYZ(p->PositionEnd[0], p->PositionEnd[1], p->PositionEnd[2]);
    else if(p->Length > 0 && Valid(p->DirectionStart))
      b = a + p->Length*TVector3(p->DirectionStart[0], p->DirectionStart[1], p->DirectionStart[2]);
    else b = a;
    return true;
  }

  // Trayectoria de la beam slice como segmentos
  std::vector<Seg> SliceSegments(const std::vector<AnaParticlePD*>& slice){
    std::vector<Seg> segs;
    for(const AnaParticlePD* p : slice){
      std::vector<TVector3> pts;
      for(const AnaTrajectoryPointPD& tp : p->TrjPoints){
        const TVector3& x = tp.Position_NoSCE;
        if(x.X() < -900 || x.Y() < -900 || x.Z() < -900)continue;
        if(!pts.empty() && (x - pts.back()).Mag() < 0.5)continue;   // puntos de menos de 5 mm: uno
        pts.push_back(x);
      }
      if(pts.size() >= 2){
        for(size_t i = 1; i < pts.size(); i++)segs.push_back({pts[i-1], pts[i]});
        continue;
      }
      TVector3 a, b;
      if(Ends(p, a, b))segs.push_back({a, b});
    }
    return segs;
  }

  // momento por rango con hipótesis de muón, en MeV/c (ComputeRangeMomentum devuelve GeV/c)
  double MuRangeP(double length){
    if(length <= 0)return -1;
    const double p = pdAnaUtils::ComputeRangeMomentum(length, 13);
    return p > 0 ? 1000.*p : -1;
  }

  bool IsComplete(const AnaParticlePD* p){
    return !p->TrjPoints.empty() || !p->SpacePoints.empty() || HasHits(p);
  }

  const double kMuP = 236.;               // MeV/c, mu de K -> mu nu en reposo
  const std::vector<double> kR = {5, 10, 20, 30, 50};
  const std::vector<double> kWin = {0.10, 0.20, 0.30};   // ventana relativa en el momento por rango

  // Criterios: nombre y función (dSeg, dEnd, p)
  struct Crit { std::string name; std::function<bool(double,double,double)> pass; };
  std::vector<Crit> Criteria(){
    std::vector<Crit> c;
    for(double r : kR)c.push_back({Form("(a) cilindro R = %g cm", r),  [r](double ds, double, double){ return ds < r; }});
    for(double r : kR)c.push_back({Form("(b) extremos R = %g cm", r),  [r](double, double de, double){ return de < r; }});
    for(double w : kWin)c.push_back({Form("(c) p rango = 236 +- %.0f %%", 100*w),
                                     [w](double, double, double p){ return std::fabs(p - kMuP) < w*kMuP; }});
    c.push_back({"(b) R = 20 cm o (c) +- 20 %", [](double, double de, double p){ return de < 20 || std::fabs(p - kMuP) < 0.2*kMuP; }});
    c.push_back({"(b) R = 50 cm y (c) +- 30 %", [](double, double de, double p){ return de < 50 && std::fabs(p - kMuP) < 0.3*kMuP; }});
    return c;
  }

  struct Target { double dSeg, dEnd, p; bool complete; };

  struct CosmicCriteria {
    std::vector<Crit> crit = Criteria();
    long nEvents = 0;
    std::vector<Target> muons, beamChain;          // emparejados con un PFParticle de fuera
    std::vector<double> muP;                       // momento por rango de todos los mu de K emparejados
    std::vector<double> muPSlice;                  // ... solo los de la beam slice
    std::vector<long> extraHier, extraParts;       // por criterio
    long summaryHier = 0, summaryParts = 0, completeHier = 0;

    CosmicCriteria(){
      extraHier.assign(crit.size(), 0);
      extraParts.assign(crit.size(), 0);
      // ComputeRangeMomentum carga las tablas de PID (pd::services()), que leen parámetros: fuera de
      // un análisis, sin esto ND::params() no deja leerlos
      ND::params().SetReadParamOverrideFilePointPassed();
    }

    void Process(AnaSpillPD* spill){
      const AnaEventInfoPD* info = static_cast<AnaEventInfoPD*>(spill->EventInfo);
      if(!info || !info->IsMC || spill->Bunches.empty())return;
      const AnaBeamPD* beam = static_cast<AnaBeamPD*>(spill->Beam);
      AnaTrueParticleB* beamTrue = (beam && beam->BeamParticle) ? static_cast<AnaTrueParticleB*>(beam->BeamParticle->TrueObject) : nullptr;

      std::vector<AnaParticlePD*> slice, out;
      for(AnaParticleB* pb : static_cast<AnaBunchPD*>(spill->Bunches[0])->Particles){
        AnaParticlePD* p = static_cast<AnaParticlePD*>(pb);
        (p->BeamOrigin ? slice : out).push_back(p);
      }
      if(slice.empty())return;
      nEvents++;
      const std::vector<Seg> segs = SliceSegments(slice);

      // distancias y momento de cada PFParticle de fuera
      std::map<const AnaParticlePD*, Target> val;
      for(AnaParticlePD* p : out){
        Target t{1e9, 1e9, -1, IsComplete(p)};
        TVector3 a, b;
        if(Ends(p, a, b)){
          const Seg s{a, b};
          for(const Seg& g : segs){
            t.dSeg = std::min(t.dSeg, SegSeg(s, g));
            t.dEnd = std::min(t.dEnd, std::min(PointSeg(a, g), PointSeg(b, g)));
          }
        }
        t.p = MuRangeP(p->Length);
        val[p] = t;
      }

      // coste: jerarquías de fuera (por su primario) que hoy van en resumen
      std::map<int, AnaParticlePD*> byUID;
      for(AnaParticlePD* p : out)byUID[p->UniqueID] = p;
      auto top = [&byUID](AnaParticlePD* p){
        for(int i = 0; i < 1000; i++){
          auto it = byUID.find(p->ParentID);
          if(it == byUID.end())return p;
          p = it->second;
        }
        return p;
      };
      std::map<AnaParticlePD*, std::vector<AnaParticlePD*>> hier;
      for(AnaParticlePD* p : out)hier[top(p)].push_back(p);
      for(auto& h : hier){
        bool complete = false;
        for(AnaParticlePD* p : h.second)if(val[p].complete){ complete = true; break; }
        if(complete){ completeHier++; continue; }
        summaryHier++;
        summaryParts += h.second.size();
        for(size_t c = 0; c < crit.size(); c++){
          bool pass = false;
          for(AnaParticlePD* p : h.second){
            const Target& t = val[p];
            if(crit[c].pass(t.dSeg, t.dEnd, t.p)){ pass = true; break; }
          }
          if(pass){ extraHier[c]++; extraParts[c] += h.second.size(); }
        }
      }

      // eficiencia: partículas true de la cadena del haz emparejadas fuera
      if(!beamTrue)return;
      std::map<int, AnaTrueParticleB*> byID;
      for(AnaTrueParticleB* t : spill->TrueParticles)byID[t->ID] = t;
      // por partícula true, el PFParticle emparejado con más hits (de fuera o de la beam slice)
      std::map<int, AnaParticlePD*> bestOut, bestSlice;
      for(AnaParticlePD* p : out){
        const AnaTrueObjectC* to = p->TrueObject;
        if(!to)continue;
        AnaParticlePD*& b = bestOut[to->ID];
        if(!b || p->NHits > b->NHits)b = p;
      }
      for(AnaParticlePD* p : slice){
        const AnaTrueObjectC* to = p->TrueObject;
        if(!to)continue;
        AnaParticlePD*& b = bestSlice[to->ID];
        if(!b || p->NHits > b->NHits)b = p;
      }
      for(AnaTrueParticleB* t : spill->TrueParticles){
        if(t->ID == beamTrue->ID || !DescendsFrom(t->ID, beamTrue->ID, byID))continue;
        const bool inSlice = bestSlice.count(t->ID);
        auto ito = bestOut.find(t->ID);
        // mu de un K+ de la cadena del haz que se desintegra
        bool muFromK = false;
        if(t->PDG == -13){
          auto ip = byID.find(t->ParentID);
          muFromK = ip != byID.end() && ip->second->PDG == 321 && ip->second->ProcessEnd == AnaTrueParticleB::Decay;
        }
        if(muFromK){
          if(inSlice && bestSlice[t->ID]->Length > 0)muPSlice.push_back(MuRangeP(bestSlice[t->ID]->Length));
          if(!inSlice && ito != bestOut.end())muP.push_back(val[ito->second].p);
          else if(inSlice && bestSlice[t->ID]->Length > 0)muP.push_back(muPSlice.back());
        }
        if(inSlice || ito == bestOut.end())continue;   // solo las que están fuera
        beamChain.push_back(val[ito->second]);
        if(muFromK)muons.push_back(val[ito->second]);
      }
    }

    static std::string Eff(const std::vector<Target>& v, const Crit& c, bool onlyNew){
      long n = 0, pass = 0;
      for(const Target& t : v){
        if(onlyNew && t.complete)continue;
        n++;
        if(c.pass(t.dSeg, t.dEnd, t.p))pass++;
      }
      return Frac(pass, n);
    }

    void Print() const {
      std::cout << std::endl << "Eventos de MC con beam slice: " << nEvents << std::endl;
      std::cout << "Jerarquías de fuera por evento: completas " << Form("%.1f", double(completeHier)/nEvents)
                << ", en resumen " << Form("%.1f", double(summaryHier)/nEvents)
                << " (" << Form("%.1f", double(summaryParts)/nEvents) << " PFParticles)" << std::endl;
      long muNow = 0, bcNow = 0;
      for(const Target& t : muons)if(t.complete)muNow++;
      for(const Target& t : beamChain)if(t.complete)bcNow++;
      std::cout << "mu de K+ emparejados fuera de la beam slice: " << muons.size() << " (hoy completos: " << muNow << ")" << std::endl;
      std::cout << "partículas de la cadena del haz emparejadas fuera: " << beamChain.size() << " (hoy completas: " << bcNow << ")" << std::endl;

      std::vector<double> all = muP, sl = muPSlice;
      std::sort(all.begin(), all.end());
      std::sort(sl.begin(), sl.end());
      auto q = [](const std::vector<double>& v, double f){ return v.empty() ? -1. : v[std::min(v.size()-1, size_t(f*v.size()))]; };
      std::cout << std::endl << "Momento por rango (muón) de los PFParticles emparejados con un mu de K+ (MeV/c):" << std::endl;
      std::cout << Form("  todos (%zu): cuantiles 10/25/50/75/90 %% = %.0f / %.0f / %.0f / %.0f / %.0f", all.size(),
                        q(all,.1), q(all,.25), q(all,.5), q(all,.75), q(all,.9)) << std::endl;
      std::cout << Form("  beam slice (%zu): cuantiles 10/25/50/75/90 %% = %.0f / %.0f / %.0f / %.0f / %.0f", sl.size(),
                        q(sl,.1), q(sl,.25), q(sl,.5), q(sl,.75), q(sl,.9)) << std::endl;
      std::cout << "  valores (todos):";
      for(double p : all)std::cout << " " << Form("%.0f", p);
      std::cout << std::endl;

      std::cout << std::endl << "Por criterio: eficiencia sobre los mu de K+ de fuera y sobre la cadena del haz de fuera" << std::endl
                << "(todas / solo las que hoy van en resumen); coste: jerarquías y PFParticles de fuera que pasarían" << std::endl
                << "de resumen a completos, por evento" << std::endl;
      std::cout << std::left << std::setw(34) << "  criterio" << std::setw(28) << "mu de K+ (todos / nuevos)"
                << std::setw(34) << "cadena del haz (todas / nuevas)" << "jerarquías  PFParticles" << std::right << std::endl;
      for(size_t c = 0; c < crit.size(); c++){
        std::cout << "  " << std::left << std::setw(32) << crit[c].name
                  << std::setw(28) << (Eff(muons, crit[c], false) + " / " + Eff(muons, crit[c], true))
                  << std::setw(34) << (Eff(beamChain, crit[c], false) + " / " + Eff(beamChain, crit[c], true))
                  << std::right << Form("%10.2f %12.2f", double(extraHier[c])/nEvents, double(extraParts[c])/nEvents) << std::endl;
      }
      std::cout << std::endl;
    }
  };
}

int cosmic_criteria_minitrees_impl(const char* input, Long64_t maxEntries){
  using namespace ccrit;
  const std::vector<std::string> files = InputFiles(input);
  if(files.empty()){ std::cerr << "No hay ficheros en " << input << std::endl; return -1; }
  const std::string treeName = TreeName(TString(files[0]).Contains("*") ? std::string("") : files[0]);
  TChain chain(treeName.empty() ? "highlandana/MiniTree" : treeName.c_str());
  for(const std::string& f : files)chain.Add(f.c_str());
  const Long64_t nEntries = chain.GetEntries();
  if(nEntries <= 0){ std::cerr << "No hay entradas (¿árbol " << chain.GetName() << "?)" << std::endl; return -1; }

  AnaSpillPD* spill = nullptr;
  chain.SetBranchAddress("Spill", &spill);
  CosmicCriteria cc;
  const Long64_t nMax = (maxEntries > 0 && maxEntries < nEntries) ? maxEntries : nEntries;
  for(Long64_t i = 0; i < nMax; i++){
    if(chain.GetEntry(i) <= 0 || !spill)continue;
    spill->CopyVectorsIntoArrays();
    spill->RedoLinks();
    cc.Process(spill);
  }
  cc.Print();
  return 0;
}
