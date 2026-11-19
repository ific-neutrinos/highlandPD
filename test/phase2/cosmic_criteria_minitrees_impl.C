// Implementación de cosmic_criteria_minitrees.C (se compila con ACLiC).
//
// Criterios candidatos para guardar completos (con hits) los PFParticles de fuera de la beam
// slice, medidos sobre minitrees ya escritos (MC, eventos con beam slice): todos los PFParticles
// de fuera se guardan al menos en resumen, con inicio, fin, longitud, slice y puntuación de
// Pandora, que es lo que usan los criterios.
//
// Criterios por PFParticle (pasa la jerarquía entera si pasa alguno de sus PFParticles, como en el
// módulo):
//   (a) cilindro: el segmento inicio-fin pasa a menos de R de la trayectoria de alguna partícula
//       de la beam slice (aproxima el criterio del módulo con los space points)
//   (b) extremos: su inicio o su fin está a menos de R de esa trayectoria
//   (c) rango: momento por rango con hipótesis de muón en una ventana alrededor de 236 MeV/c
//       (mu+ de K+ -> mu+ nu en reposo)
// Criterios por slice de Pandora (pasan todas las jerarquías del slice):
//   (d) puntuación: la puntuación de haz del slice (máximo de BeamCosmicScore de sus PFParticles,
//       el TestBeamScore de Pandora) por encima de un umbral. Los umbrales son percentiles de la
//       puntuación de los slices sin nada de la cadena del haz: el percentil 90 deja pasar el 10 %
//       de esos slices
//   (e) slice cercano: algún PFParticle del slice cumple (b)
// La trayectoria de una partícula de la beam slice son sus puntos de trayectoria (sin SCE) o, si
// no tiene (showers), el segmento de su inicio a su fin (o a inicio + dirección x longitud).
//
// Objetivos (para la eficiencia): las partículas true de la cadena del haz (el haz y sus
// descendientes) que no están emparejadas con ninguna partícula de la beam slice pero sí con un
// PFParticle de fuera (el de más hits), y entre ellas los mu de K+ que se desintegran. Eficiencia
// en número y ponderada por los hits del PFParticle; "nuevos" son los que hoy están en una
// jerarquía en resumen.
// Coste: jerarquías y PFParticles que hoy van en resumen y pasarían a completos, por evento.
//
// Eventos con volcado completo (EventInfo.FullDump, en MC los que tienen un K+ secundario): ahí
// todas las jerarquías están completas, así que "hoy en resumen" se emula con el criterio del
// módulo en CosmicsMode 2 (algún punto de la jerarquía a menos de CosmicsRadius de un extremo de
// la beam slice). Los puntos son los space points y, en las trazas (sin space points guardados),
// los puntos de trayectoria sin SCE: es lo que más se parece a los space points del módulo. En los
// eventos normales se compara la emulación con lo que hizo el módulo, para validarla. Los dos
// grupos se imprimen por separado: los eventos con volcado están sesgados (tienen un K+), pero
// en ellos todo está emparejado con hits.
//
// Guarda en memoria un resumen de todos los PFParticles de fuera (~40 B cada uno): para muestras
// muy grandes, usar maxEntries o "which" (1 solo eventos normales, 2 solo con volcado completo).
#include "kaon_fate_minitrees_impl.C"
#include "pdAnalysisUtils.hxx"
#include "Parameters.hxx"
#include "TVector3.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
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

  bool IsComplete(const AnaParticlePD* p){
    return !p->TrjPoints.empty() || !p->SpacePoints.empty() || HasHits(p);
  }

  // momento por rango con hipótesis de muón, en MeV/c (ComputeRangeMomentum devuelve GeV/c)
  double MuRangeP(double length){
    if(length <= 0)return -1;
    const double p = pdAnaUtils::ComputeRangeMomentum(length, 13);
    return p > 0 ? 1000.*p : -1;
  }

  const double kMuP = 236.;   // MeV/c, mu de K -> mu nu en reposo

  // Lo que se guarda de cada PFParticle de fuera, de cada jerarquía y de cada slice
  struct PartV {
    double dSeg = 1e9, dEnd = 1e9, p = -1;
    int nhits = 0;
    bool target = false, muK = false;   // objetivo de la cadena del haz; mu de K+
  };
  struct Hier {
    std::vector<PartV> parts;
    bool complete = false;   // guardada completa en el minitree
    bool emulated = false;   // lo estaría con CosmicsMode 2 (criterio del módulo, emulado)
  };
  struct Slice { double score = -999; bool hasTarget = false; std::vector<Hier> hiers; };
  struct Event { std::vector<Slice> slices; bool fullDump = false; };

  // Extremos de la beam slice, como GetBeamSliceEnds del módulo
  std::vector<TVector3> SliceEnds(const std::vector<AnaParticlePD*>& slice){
    std::vector<TVector3> ends;
    for(const AnaParticlePD* p : slice){
      if(p->PositionStart[0] > -900)ends.emplace_back(p->PositionStart[0], p->PositionStart[1], p->PositionStart[2]);
      if(p->PositionEnd[0]   > -900)ends.emplace_back(p->PositionEnd[0],   p->PositionEnd[1],   p->PositionEnd[2]);
    }
    return ends;
  }

  // Criterio de CosmicsMode 2 del módulo (IsNearBeamSlice) con los puntos guardados de un PFParticle
  bool NearEnds(const AnaParticlePD* p, const std::vector<TVector3>& ends, double r){
    const double r2 = r*r;
    auto near = [&](double x, double y, double z){
      if(x < -900 || y < -900 || z < -900)return false;
      for(const TVector3& e : ends){
        const double dx = x-e.X(), dy = y-e.Y(), dz = z-e.Z();
        if(dx*dx + dy*dy + dz*dz < r2)return true;
      }
      return false;
    };
    for(const AnaSpacePointPD& sp : p->SpacePoints)if(near(sp.Position[0], sp.Position[1], sp.Position[2]))return true;
    for(const AnaTrajectoryPointPD& tp : p->TrjPoints)if(near(tp.Position_NoSCE.X(), tp.Position_NoSCE.Y(), tp.Position_NoSCE.Z()))return true;
    return false;
  }
  const double kCosmicsRadius = 50.;   // cm, CosmicsRadius del fcl

  // Criterio: decide para una jerarquía dentro de su slice
  struct Crit { std::string name; std::function<bool(const Slice&, const Hier&)> pass; };

  bool AnyPart(const Hier& h, const std::function<bool(const PartV&)>& f){
    for(const PartV& p : h.parts)if(f(p))return true;
    return false;
  }
  bool AnyPart(const Slice& s, const std::function<bool(const PartV&)>& f){
    for(const Hier& h : s.hiers)if(AnyPart(h, f))return true;
    return false;
  }

  struct CosmicCriteria {
    int which = 0;                       // 0 todos los eventos, 1 solo normales, 2 solo con volcado completo
    std::vector<Event> events;
    long emulAgree[2][2] = {{0}};        // eventos normales: [completa en el minitree][emulada]
    std::vector<double> muP, muPSlice;   // momento por rango de los mu de K+ emparejados (todos; beam slice)
    std::vector<double> scoreBkg, scoreTarget;
    std::vector<int> nTargetSlices;      // por evento: slices de fuera con algún objetivo
    std::vector<double> topSliceFrac;    // por evento: fracción de los hits de objetivos en su slice principal

    CosmicCriteria(){
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
      if((which == 1 && info->FullDump) || (which == 2 && !info->FullDump))return;
      const std::vector<Seg> segs = SliceSegments(slice);
      const std::vector<TVector3> ends = SliceEnds(slice);

      // objetivos: partículas de la cadena del haz sin emparejar en la beam slice -> PFParticle de fuera con más hits
      std::map<const AnaParticlePD*, std::pair<bool,bool>> targets;   // (objetivo, mu de K)
      if(beamTrue){
        std::map<int, AnaTrueParticleB*> byID;
        for(AnaTrueParticleB* t : spill->TrueParticles)byID[t->ID] = t;
        std::map<int, AnaParticlePD*> bestOut, bestSlice;
        for(AnaParticlePD* p : out)if(p->TrueObject){
          AnaParticlePD*& b = bestOut[p->TrueObject->ID];
          if(!b || p->NHits > b->NHits)b = p;
        }
        for(AnaParticlePD* p : slice)if(p->TrueObject){
          AnaParticlePD*& b = bestSlice[p->TrueObject->ID];
          if(!b || p->NHits > b->NHits)b = p;
        }
        for(AnaTrueParticleB* t : spill->TrueParticles){
          if(!DescendsFrom(t->ID, beamTrue->ID, byID))continue;
          bool muFromK = false;
          if(t->PDG == -13){
            auto ip = byID.find(t->ParentID);
            muFromK = ip != byID.end() && ip->second->PDG == 321 && ip->second->ProcessEnd == AnaTrueParticleB::Decay;
          }
          auto is = bestSlice.find(t->ID);
          auto io = bestOut.find(t->ID);
          if(muFromK){
            if(is != bestSlice.end()){ muPSlice.push_back(MuRangeP(is->second->Length)); muP.push_back(muPSlice.back()); }
            else if(io != bestOut.end())muP.push_back(MuRangeP(io->second->Length));
          }
          if(is != bestSlice.end() || io == bestOut.end())continue;
          targets[io->second] = {true, muFromK};
        }
      }

      // jerarquías de fuera (por su primario) y slices (el del primario; sin slice, uno por jerarquía)
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
      std::map<AnaParticlePD*, std::vector<AnaParticlePD*>> hierParts;
      for(AnaParticlePD* p : out)hierParts[top(p)].push_back(p);

      std::map<long, Slice> slices;
      std::map<long, double> targetHits;
      for(auto& hp : hierParts){
        const long key = hp.first->SliceID >= 0 ? (long)hp.first->SliceID : 1000000L + hp.first->UniqueID;
        Slice& s = slices[key];
        Hier h;
        for(AnaParticlePD* p : hp.second){
          PartV v;
          TVector3 a, b;
          if(Ends(p, a, b)){
            const Seg sg{a, b};
            for(const Seg& g : segs){
              v.dSeg = std::min(v.dSeg, SegSeg(sg, g));
              v.dEnd = std::min(v.dEnd, std::min(PointSeg(a, g), PointSeg(b, g)));
            }
          }
          v.p = MuRangeP(p->Length);
          v.nhits = p->NHits;
          auto it = targets.find(p);
          if(it != targets.end()){ v.target = true; v.muK = it->second.second; s.hasTarget = true; targetHits[key] += p->NHits; }
          if(IsComplete(p))h.complete = true;
          if(!h.emulated && NearEnds(p, ends, kCosmicsRadius))h.emulated = true;
          if(p->BeamCosmicScore > -900)s.score = std::max(s.score, (double)p->BeamCosmicScore);
          h.parts.push_back(v);
        }
        if(!info->FullDump)emulAgree[h.complete][h.emulated]++;
        s.hiers.push_back(h);
      }

      Event ev;
      ev.fullDump = info->FullDump;
      for(auto& ks : slices){
        (ks.second.hasTarget ? scoreTarget : scoreBkg).push_back(ks.second.score);
        ev.slices.push_back(ks.second);
      }
      events.push_back(ev);
      if(!targetHits.empty()){
        double sum = 0, mx = 0;
        for(auto& th : targetHits){ sum += th.second; mx = std::max(mx, th.second); }
        nTargetSlices.push_back(targetHits.size());
        if(sum > 0)topSliceFrac.push_back(mx/sum);
      }
    }

    static double Quantile(std::vector<double> v, double f){
      if(v.empty())return -999;
      std::sort(v.begin(), v.end());
      return v[std::min(v.size()-1, size_t(f*v.size()))];
    }

    std::vector<Crit> Criteria() const {
      std::vector<Crit> c;
      for(double r : {10., 20., 50.})
        c.push_back({Form("(a) cilindro R = %g cm", r), [r](const Slice&, const Hier& h){ return AnyPart(h, [r](const PartV& p){ return p.dSeg < r; }); }});
      for(double r : {10., 20., 50.})
        c.push_back({Form("(b) extremos R = %g cm", r), [r](const Slice&, const Hier& h){ return AnyPart(h, [r](const PartV& p){ return p.dEnd < r; }); }});
      c.push_back({"(c) p rango = 236 +- 20 %", [](const Slice&, const Hier& h){ return AnyPart(h, [](const PartV& p){ return std::fabs(p.p - kMuP) < 0.2*kMuP; }); }});
      for(double q : {0.5, 0.8, 0.9, 0.95, 0.99}){
        const double thr = Quantile(scoreBkg, q);
        c.push_back({Form("(d) puntuación > p%.0f (%.3g)", 100*q, thr), [thr](const Slice& s, const Hier&){ return s.score > thr; }});
      }
      for(double r : {20., 50.})
        c.push_back({Form("(e) slice cercano R = %g cm", r), [r](const Slice& s, const Hier&){ return AnyPart(s, [r](const PartV& p){ return p.dEnd < r; }); }});
      const double thr95 = Quantile(scoreBkg, 0.95);
      c.push_back({"(d) p95 o (e) R = 50 cm", [thr95](const Slice& s, const Hier&){
        return s.score > thr95 || AnyPart(s, [](const PartV& p){ return p.dEnd < 50; }); }});
      return c;
    }

    void Print() const {
      long nFD = 0;
      for(const Event& e : events)if(e.fullDump)nFD++;
      // en resumen no hay puntos guardados: la emulación solo puede comprobarse en las completas
      const long nC = emulAgree[1][0] + emulAgree[1][1];
      if(nC){
        std::cout << std::endl << "Validación de la emulación de CosmicsMode 2 (eventos normales): de las " << nC
                  << " jerarquías completas en el minitree, emuladas completas " << Frac(emulAgree[1][1], nC)
                  << " (las demás: CosmicsMode 3 o puntos no guardados; en resumen no hay puntos con que comprobar)" << std::endl;
      }
      if(long(events.size()) > nFD)PrintGroup("EVENTOS NORMALES (completas: las del minitree)", false);
      if(nFD)PrintGroup("EVENTOS CON VOLCADO COMPLETO (completas: emuladas con CosmicsMode 2; sesgados: tienen un K+)", true);
    }

    void PrintGroup(const std::string& title, bool fullDump) const {
      std::vector<Event> sel;
      for(const Event& e : events)if(e.fullDump == fullDump){
        Event c = e;
        if(fullDump)for(Slice& s : c.slices)for(Hier& h : s.hiers)h.complete = h.emulated;
        sel.push_back(c);
      }
      std::cout << std::endl << "===== " << title << " =====" << std::endl;
      PrintEvents(sel);
    }

    void PrintEvents(const std::vector<Event>& events) const {
      const double nEv = events.size();
      long complete = 0, summary = 0, summaryParts = 0, nSlices = 0;
      long nT = 0, nTNow = 0, nMu = 0, nMuNow = 0;
      double hT = 0;
      for(const Event& e : events)for(const Slice& s : e.slices){
        nSlices++;
        for(const Hier& h : s.hiers){
          if(h.complete)complete++;
          else{ summary++; summaryParts += h.parts.size(); }
          for(const PartV& p : h.parts){
            if(p.target){ nT++; hT += p.nhits; if(h.complete)nTNow++; }
            if(p.muK){ nMu++; if(h.complete)nMuNow++; }
          }
        }
      }
      std::cout << std::endl << "Eventos de MC con beam slice: " << events.size() << std::endl;
      std::cout << Form("Por evento: %.1f slices de fuera; jerarquías completas %.1f, en resumen %.1f (%.1f PFParticles)",
                        nSlices/nEv, complete/nEv, summary/nEv, summaryParts/nEv) << std::endl;
      std::cout << "Objetivos (cadena del haz emparejada fuera): " << nT << " (hoy completos: " << nTNow
                << ", " << Form("%.0f", hT/std::max(1L, nT)) << " hits de media)" << std::endl;
      std::cout << "  de ellos mu de K+: " << nMu << " (hoy completos: " << nMuNow << ")" << std::endl;
      if(!nTargetSlices.empty()){
        double m = 0, f = 0;
        for(int n : nTargetSlices)m += n;
        for(double x : topSliceFrac)f += x;
        std::cout << Form("  eventos con objetivos: %zu; slices de fuera con objetivos por evento %.2f; fracción de sus hits en el slice principal %.0f %%",
                          nTargetSlices.size(), m/nTargetSlices.size(), 100*f/std::max<size_t>(1, topSliceFrac.size())) << std::endl;
      }
      std::cout << std::endl << "Puntuación de haz de los slices de fuera (cuantiles 10/50/90/99 %):" << std::endl;
      std::cout << Form("  sin objetivos (%zu): %.3g / %.3g / %.3g / %.3g", scoreBkg.size(), Quantile(scoreBkg,.1), Quantile(scoreBkg,.5), Quantile(scoreBkg,.9), Quantile(scoreBkg,.99)) << std::endl;
      std::cout << Form("  con objetivos (%zu): %.3g / %.3g / %.3g / %.3g", scoreTarget.size(), Quantile(scoreTarget,.1), Quantile(scoreTarget,.5), Quantile(scoreTarget,.9), Quantile(scoreTarget,.99)) << std::endl;

      std::cout << std::endl << "Momento por rango (muón) de los PFParticles emparejados con un mu de K+ (MeV/c, cuantiles 10/25/50/75/90 %):" << std::endl;
      std::cout << Form("  todos (%zu): %.0f / %.0f / %.0f / %.0f / %.0f", muP.size(), Quantile(muP,.1), Quantile(muP,.25), Quantile(muP,.5), Quantile(muP,.75), Quantile(muP,.9)) << std::endl;
      std::cout << Form("  beam slice (%zu): %.0f / %.0f / %.0f / %.0f / %.0f", muPSlice.size(), Quantile(muPSlice,.1), Quantile(muPSlice,.25), Quantile(muPSlice,.5), Quantile(muPSlice,.75), Quantile(muPSlice,.9)) << std::endl;

      std::cout << std::endl << "Por criterio. Eficiencia sobre los objetivos que hoy van en resumen (en número y ponderada por hits)" << std::endl
                << "y sobre los mu de K+ (todos / hoy en resumen); coste: jerarquías y PFParticles por evento que pasarían" << std::endl
                << "de resumen a completos" << std::endl;
      std::cout << "  " << std::left << std::setw(36) << "criterio" << std::setw(24) << "objetivos nuevos" << std::setw(12) << "por hits"
                << std::setw(14) << "mu de K+" << std::right << "  jerarquías  PFParticles" << std::endl;
      for(const Crit& c : Criteria()){
        long tN = 0, tP = 0, mA = 0, mAP = 0, mN = 0, mNP = 0, xH = 0, xP = 0;
        double hN = 0, hP = 0;
        for(const Event& e : events)for(const Slice& s : e.slices)for(const Hier& h : s.hiers){
          const bool pass = c.pass(s, h);
          if(pass && !h.complete){ xH++; xP += h.parts.size(); }
          for(const PartV& p : h.parts){
            if(p.target && !h.complete){ tN++; hN += p.nhits; if(pass){ tP++; hP += p.nhits; } }
            if(p.muK){ mA++; if(pass || h.complete)mAP++; if(!h.complete){ mN++; if(pass)mNP++; } }
          }
        }
        std::cout << "  " << std::left << std::setw(36) << c.name << std::setw(24) << Frac(tP, tN)
                  << std::setw(12) << Form("%.1f %%", hN > 0 ? 100*hP/hN : 0.)
                  << std::setw(14) << Form("%ld/%ld  %ld/%ld", mAP, mA, mNP, mN) << std::right
                  << Form("%12.2f %12.2f", xH/nEv, xP/nEv) << std::endl;
      }
      std::cout << "  (mu de K+: con el criterio, cuántos quedan completos de todos, y cuántos se recuperan de los que hoy van en resumen)" << std::endl << std::endl;
    }
  };
}

int cosmic_criteria_minitrees_impl(const char* input, Long64_t maxEntries, int which = 0){
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
  cc.which = which;
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
