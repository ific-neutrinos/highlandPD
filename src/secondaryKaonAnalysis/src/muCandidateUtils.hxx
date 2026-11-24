#ifndef muCandidateUtils_h
#define muCandidateUtils_h
// Candidatos a mu+ de K+ -> mu+ nu en reposo: variables de cada candidato (PFParticle con traza, o sin traza con su
// traza forzada) a partir de un minitree del módulo CreateHighlandMiniTree (highlandAnalysis). Una sola definición
// para el análisis (secondaryKaonProductionXS) y las macros de test/phase2 (mu_candidates_minitrees.C): mismos
// números en los dos sitios. Solo cabecera (funciones inline) para que la usen también las macros con ACLiC.
//
// Uso: MakeContext una vez por evento (partículas del bunch, hits no asociados, truth, haz) y Compute para cada
// partícula; Compute devuelve false si la partícula no es candidato (preselección: longitud y extremos válidos, p de
// rango de mu en [kPMin, kPMax] MeV/c). Las variables y su significado, en Candidate.
//
// Eventos con volcado completo (AnaEventInfoPD::FullDump, MC con un K+ en la cadena del haz): el módulo guarda
// completos todos los PFParticles, más de lo que guardaría en un evento normal (y en datos). Para tratarlos como uno
// normal sin perder nada, MakeNormalView da la "vista normal" del evento (lo que el módulo habría guardado), sin
// tocar el spill; la vista completa es la lista original. Ver NormalSaveMode.
#include "pdDataClasses.hxx"
#include "pdAnalysisUtils.hxx"
#include "pdServiceManager.hxx"
#include "PIDTables.hxx"
#include "TProfile.h"
#include "TVector3.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVectorD.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <memory>
#include <set>
#include <vector>

namespace muCandUtils {

  // Como el módulo (CreateHighlandMiniTree.fcl)
  const double kActive[6] = {-360., 360., 0., 607.5, 0., 695.};
  const TVector3 kAxisStart[2] = {TVector3(-28.6, 424.6, 0.96), TVector3(-30.1, 423.0, 0.0)};   // datos, MC
  const TVector3 kAxisDir = TVector3(-0.18, -0.20, 0.96).Unit();
  const double kAxisLength = 700.;
  const double kPMin = 100., kPMax = 600.;      // MeV/c, candidato amplio
  const double kWinMin = 180., kWinMax = 300.;  // MeV/c, ventana del módulo (MuCandidateLength)
  const double kNear = 10.;                     // cm, actividad en un extremo
  const double kBragg = 5.;                     // cm, dE/dx en el extremo

  inline bool Valid(const Float_t* x){ return x[0] > -900 && x[1] > -900 && x[2] > -900; }
  inline TVector3 V(const Float_t* x){ return TVector3(x[0], x[1], x[2]); }

  inline double DistToFaces(const TVector3& p){
    return std::min({p.X()-kActive[0], kActive[1]-p.X(), p.Y()-kActive[2], kActive[3]-p.Y(), p.Z()-kActive[4], kActive[5]-p.Z()});
  }
  inline double DistToSegment(const TVector3& p, const TVector3& a, const TVector3& b){
    const TVector3 ab = b - a;
    const double l2 = ab.Mag2();
    double t = l2 > 0 ? (p - a).Dot(ab)/l2 : 0.;
    t = std::max(0., std::min(1., t));
    return (p - (a + t*ab)).Mag();
  }
  inline double DistToAxis(const TVector3& p, int isMC){
    return DistToSegment(p, kAxisStart[isMC], kAxisStart[isMC] + kAxisLength*kAxisDir);
  }
  inline double MuRangeP(double length){   // MeV/c
    return length > 0 ? 1000.*pdAnaUtils::ComputeRangeMomentum(length, 13) : -1;
  }
  inline bool HasHits(const AnaParticlePD* p){ return !p->Hits[2].empty() || p->NHits > 0; }

  // Los dos extremos en la escala del rango residual (RR) de la calorimetría (plano 2): B, el final de Pandora, en
  // RR = 0; A, el inicio, en RRMax, el RR del hit más lejano. No se usa p->Length: es la longitud reconstruida sin
  // corregir de SCE, y el RR de los hits está corregido (y la calorimetría puede no tener hits en los primeros cm de
  // la traza): RRMax - Length va de -14 a +14 cm (percentiles 5 y 95, trazas de 30-90 cm), y con Length la ventana
  // del extremo A quedaba vacía en el 14 % de las trazas (corregido el 24-nov). -1 sin hits con RR.
  inline double RRMax(const AnaParticlePD* p){
    double m = -1;
    for(const AnaHitPD& h : p->Hits[2])if(h.ResidualRange >= 0 && h.ResidualRange > m)m = h.ResidualRange;
    return m;
  }

  // Hits con dE/dx utilizable (los demás se ignoran en todas las variables de dE/dx y chi2)
  inline bool GooddEdx(const AnaHitPD& h){ return h.dEdx > 0 && h.dEdx <= 100 && h.ResidualRange >= 0; }

  // dE/dx medio de los hits del plano 2 a < kBragg cm de cada extremo (por residual range: RR
  // pequeño = cerca de PositionEnd, en la convención de la calorimetría de Pandora; el extremo A, en RRMax). -1 sin hits
  inline void EndsdEdx(const AnaParticlePD* p, double& dEdxA, double& dEdxB, int& nA, int& nB){
    double sA = 0, sB = 0; nA = nB = 0;
    const double L = RRMax(p);
    for(const AnaHitPD& h : p->Hits[2]){
      if(!GooddEdx(h))continue;
      if(h.ResidualRange < kBragg){ sB += h.dEdx; nB++; }
      if(h.ResidualRange > L - kBragg){ sA += h.dEdx; nA++; }
    }
    dEdxA = nA ? sA/nA : -1; dEdxB = nB ? sB/nB : -1;
  }

  inline double Chi2Trimmed(const AnaParticlePD* p, int pdg, bool stopAtB, double trim);

  // chi2/ndf de una hipótesis con la orientación de Pandora (fwd: se para en PositionEnd) y al revés
  // (rev: RR' = RRMax - RR, se para en PositionStart). -1 si no hay. Con Chi2Trimmed (sin recorte), no con
  // pdAnaUtils::Chi2PID: mismo filtro de hits (GooddEdx) y misma referencia del extremo A que el resto de variables
  inline void Chi2Both(const AnaParticlePD* p, int pdg, double& fwd, double& rev){
    fwd = Chi2Trimmed(p, pdg, true, 0);
    rev = Chi2Trimmed(p, pdg, false, 0);
  }

  // chi2/ndf de la hipótesis pdg con la parada en B (stopAtB, el final de Pandora) o en A, quitando los hits a menos
  // de trim cm de la parada y restando trim al rango residual de los demás: la parada se mueve trim cm hacia dentro.
  // Para trazas con hits de otra partícula pegados al final (el Michel) o que se pasan de la parada. Fórmula de
  // pdAnaUtils::Chi2PID (plano 2, sin el primer ni el último hit en el orden guardado, plantilla de dE/dx frente al
  // rango residual, resolución en dE/dx de la calorimetría), pero solo con hits GooddEdx (Chi2PID solo quita dE/dx >
  // 1000 y -999) y la orientación contraria referida a RRMax. Con trim = 0, el chi2 de Chi2Both. -1 sin hits.
  inline double Chi2Trimmed(const AnaParticlePD* p, int pdg, bool stopAtB, double trim){
    const std::vector<AnaHitPD>& hits = p->Hits[2];
    if(hits.size() < 3)return -1;
    TProfile* profile = const_cast<TProfile*>(pd::services().PID().GetdEdxTemplate(pdg));
    if(!profile)return -1;
    const double L = stopAtB ? 0 : RRMax(p);
    double chi2 = 0; int npt = 0;
    for(size_t i = 1; i + 1 < hits.size(); i++){
      const AnaHitPD& h = hits[i];
      if(!GooddEdx(h))continue;
      double rr = stopAtB ? h.ResidualRange : std::max(0., L - h.ResidualRange);
      if(trim > 0){ if(rr < trim)continue; rr -= trim; }
      const int bin = profile->FindBin(rr);
      if(bin < 1 || bin > profile->GetNbinsX())continue;
      double t = profile->GetBinContent(bin), te = profile->GetBinError(bin);
      if(t < 1e-6)t = (profile->GetBinContent(bin - 1) + profile->GetBinContent(bin + 1))/2.;
      if(te < 1e-6)te = (profile->GetBinError(bin - 1) + profile->GetBinError(bin + 1))/2.;
      const double res = (0.04231 + 0.0001783*h.dEdx*h.dEdx)*h.dEdx;
      chi2 += std::pow(h.dEdx - t, 2)/(te*te + res*res);
      npt++;
    }
    return npt > 0 ? chi2/npt : -1;
  }

  // Mejor chi2/ndf (el menor) recortando trim = 0, kTrimStep, ..., kTrimMax cm de la parada (Chi2Trimmed); en bestTrim,
  // el recorte que lo da (-1 sin chi2)
  const double kTrimMax = 10., kTrimStep = 1.;
  inline double Chi2BestTrim(const AnaParticlePD* p, int pdg, bool stopAtB, double& bestTrim){
    double best = -1; bestTrim = -1;
    for(double trim = 0; trim <= kTrimMax + 1e-6; trim += kTrimStep){
      const double c = Chi2Trimmed(p, pdg, stopAtB, trim);
      if(c >= 0 && (best < 0 || c < best)){ best = c; bestTrim = trim; }
    }
    return best;
  }

  // Puntuación CNN de Michel media de los hits propios (plano 2) a < kBragg cm de cada extremo
  inline void OwnMichel(const AnaParticlePD* p, double& mA, double& mB){
    double sA = 0, sB = 0; int nA = 0, nB = 0;
    const double L = RRMax(p);
    for(const AnaHitPD& h : p->Hits[2]){
      if(h.ResidualRange < 0 || h.CNN[2] < 0)continue;
      if(h.ResidualRange < kBragg){ sB += h.CNN[2]; nB++; }
      if(h.ResidualRange > L - kBragg){ sA += h.CNN[2]; nA++; }
    }
    mA = nA ? sA/nA : -1; mB = nB ? sB/nB : -1;
  }

  // Hits propios (plano 2) con CNN de Michel > 0,5 a < kBragg cm del extremo atB (B = final de Pandora, RR pequeño;
  // A = inicio, RR cerca de RRMax): el Michel que Pandora mete en la traza del mu no aparece en los hits de alrededor
  inline int OwnMichelHits(const AnaParticlePD* p, bool atB){
    const double L = RRMax(p);
    int n = 0;
    for(const AnaHitPD& h : p->Hits[2]){
      if(h.ResidualRange < 0 || h.CNN[2] <= 0.5)continue;
      if(atB ? h.ResidualRange < kBragg : h.ResidualRange > L - kBragg)n++;
    }
    return n;
  }

  // Pico de Bragg interior: máximo de la media de 5 hits consecutivos (plano 2, por residual range)
  // con todos a > kBragg cm de los dos extremos. Trazas K+mu fundidas: el pico de la K a ~54 cm del
  // final del mu
  inline void InteriorPeak(const AnaParticlePD* p, double& peak, double& rrPeak){
    peak = -1; rrPeak = -1;
    std::vector<std::pair<double,double>> h;   // (RR, dE/dx)
    for(const AnaHitPD& x : p->Hits[2])if(GooddEdx(x))h.push_back({x.ResidualRange, x.dEdx});
    std::sort(h.begin(), h.end());
    const int w = 5;
    const double L = RRMax(p);
    for(int i = 0; i + w <= int(h.size()); i++){
      if(h[i].first < kBragg || h[i+w-1].first > L - kBragg)continue;
      double m = 0; for(int k = 0; k < w; k++)m += h[i+k].second; m /= w;
      if(m > peak){ peak = m; rrPeak = h[i+w/2].first; }
    }
  }

  struct Activity { int nParts = 0, nHits = 0; double michel = -1; };
  // Otros PFParticles con algún extremo a < kNear cm de x: cuántos, sus hits del plano 2 y la
  // puntuación CNN de Michel media de esos hits
  inline Activity ActivityAt(const TVector3& x, const AnaParticlePD* self, const std::vector<AnaParticlePD*>& parts){
    Activity a; double sm = 0; int nm = 0;
    for(const AnaParticlePD* q : parts){
      if(q == self)continue;
      bool near = false;
      if(Valid(q->PositionStart) && (V(q->PositionStart) - x).Mag() < kNear)near = true;
      if(Valid(q->PositionEnd)   && (V(q->PositionEnd)   - x).Mag() < kNear)near = true;
      if(!near)continue;
      a.nParts++;
      a.nHits += q->Hits[2].size();
      for(const AnaHitPD& h : q->Hits[2])if(h.CNN[2] >= 0){ sm += h.CNN[2]; nm++; }
    }
    if(nm)a.michel = sm/nm;
    return a;
  }

  // ---- Variables añadidas el 20-nov ----
  const double kWirePitch = 0.479;   // cm, plano de colección
  const double kTickCm = 0.0785;     // cm por tick (0,5 us x ~1,57 mm/us)
  const double kRadii[3] = {5., 10., 20.};

  // Hit de calorimetría del plano 2 más cercano a un extremo (por residual range): atB = PositionEnd
  inline const AnaHitPD* EndHit(const AnaParticlePD* p, bool atB){
    const AnaHitPD* best = nullptr;
    for(const AnaHitPD& h : p->Hits[2]){
      if(h.ResidualRange < 0)continue;
      if(!best || (atB ? h.ResidualRange < best->ResidualRange : h.ResidualRange > best->ResidualRange))best = &h;
    }
    return best;
  }

  // Hits del plano 2 alrededor de un hit (distancia 2D hilo-tiempo en la misma TPC), para kRadii: todos los
  // PFParticles (salvo los hits de calorimetría del propio candidato) y los hits no asociados
  struct Around { int n[3] = {0, 0, 0}, nMichel[3] = {0, 0, 0}; double michel[3] = {-1, -1, -1}, charge[3] = {0, 0, 0}; };
  inline Around HitsAround(const AnaHitPD* h0, const AnaParticlePD* self, const std::vector<AnaParticlePD*>& parts, const std::vector<AnaHitPD>& unassoc){
    Around a;
    if(!h0)return a;
    double sm[3] = {0, 0, 0}; int nm[3] = {0, 0, 0};
    auto add = [&](const AnaHitPD& h){
      if(h.WireID.Plane != 2 || h.WireID.TPC != h0->WireID.TPC)return;
      const double dw = (h.WireID.Wire - h0->WireID.Wire)*kWirePitch, dt = (h.PeakTime - h0->PeakTime)*kTickCm;
      if(std::abs(dw) >= kRadii[2] || std::abs(dt) >= kRadii[2])return;
      const double d = std::hypot(dw, dt);
      for(int r = 0; r < 3; r++){
        if(d >= kRadii[r])continue;
        a.n[r]++; a.charge[r] += h.Integral;
        if(h.CNN[2] >= 0){ sm[r] += h.CNN[2]; nm[r]++; if(h.CNN[2] > 0.5)a.nMichel[r]++; }
      }
    };
    for(const AnaParticlePD* q : parts){
      if(q != self)for(const AnaHitPD& h : q->Hits[2])add(h);
      for(const AnaHitPD& h : q->HitsNoCalo[2])add(h);
    }
    for(const AnaHitPD& h : unassoc)add(h);
    for(int r = 0; r < 3; r++)if(nm[r])a.michel[r] = sm[r]/nm[r];
    return a;
  }

  // Media truncada del dE/dx (plano 2): se ordenan y se quitan el 10 % más bajo y el 40 % más alto.
  // body: sin los kBragg cm de cada extremo. -1 con menos de 5 hits
  inline double TruncMean(const AnaParticlePD* p, bool body){
    std::vector<double> v;
    const double L = RRMax(p);
    for(const AnaHitPD& h : p->Hits[2]){
      if(!GooddEdx(h))continue;
      if(body && (h.ResidualRange < kBragg || h.ResidualRange > L - kBragg))continue;
      v.push_back(h.dEdx);
    }
    if(v.size() < 5)return -1;
    std::sort(v.begin(), v.end());
    const size_t i0 = size_t(0.1*v.size()), i1 = size_t(0.6*v.size());
    double sum = 0; for(size_t i = i0; i < i1; i++)sum += v[i];
    return i1 > i0 ? sum/(i1 - i0) : -1;
  }

  // Lo que acaba en el inicio del mu (la K): otros PFParticles con un extremo a < kNear cm de x; el de mayor
  // dE/dx (kBragg cm) en ese extremo
  struct KAtStart { int n = 0, isTrack = -1; double dEdx = -1, length = -1, pRangeK = -1; };
  inline KAtStart KAt(const TVector3& x, const AnaParticlePD* self, const std::vector<AnaParticlePD*>& parts){
    KAtStart k;
    for(const AnaParticlePD* q : parts){
      if(q == self)continue;
      const bool nearA = Valid(q->PositionStart) && (V(q->PositionStart) - x).Mag() < kNear;
      const bool nearB = Valid(q->PositionEnd) && (V(q->PositionEnd) - x).Mag() < kNear;
      if(!nearA && !nearB)continue;
      k.n++;
      double dA, dB; int nA, nB; EndsdEdx(q, dA, dB, nA, nB);
      const double d = std::max(nearA ? dA : -1., nearB ? dB : -1.);
      if(d > k.dEdx){
        k.dEdx = d; k.length = q->Length; k.isTrack = q->Type == AnaParticlePD::kTrack;
        // ComputeRangeMomentum solo tiene mu y p: escalado con la masa (mismo beta*gamma), R_K(p) = (mK/mmu) R_mu(p mmu/mK)
        const double r = pdAnaUtils::kKaonMass/pdAnaUtils::kMuonMass;
        k.pRangeK = q->Length > 0 ? 1000.*r*pdAnaUtils::ComputeRangeMomentum(q->Length/r, 13) : -1;
      }
    }
    return k;
  }

  // ---- Conexión con el vértice primario (23-nov) ----
  // Si la K sale del vértice del haz V y se para en el inicio S del mu, entre V y S hay una traza (casi recta) que
  // acaba con su pico de Bragg en S. Nube de puntos 3D del evento: hits de calorimetría del plano 2
  // (Position_NoSCE, en el mismo sistema que los extremos de las trazas) y space points de los PFParticles sin
  // traza. Corredor de radio kVcRadius alrededor del segmento V -> S; se ignoran los kVcSkipV cm junto a V (ahí
  // se juntan todas las trazas del vértice).
  const double kVcRadius = 2., kVcSkipV = 3., kVcBragg = 5.;
  struct VcPoint { TVector3 x; Int_t owner; double dEdx; };
  inline std::vector<VcPoint> VcCloud(const std::vector<AnaParticlePD*>& parts){
    std::vector<VcPoint> cloud;
    for(const AnaParticlePD* q : parts){
      for(const AnaHitPD& h : q->Hits[2])
        if(h.Position_NoSCE.X() > -900 && h.Position_NoSCE.Mag2() > 0)cloud.push_back({h.Position_NoSCE, q->UniqueID, (h.dEdx > 0 && h.dEdx < 100) ? h.dEdx : -1.});
      for(const AnaSpacePointPD& sp : q->SpacePoints)
        if(Valid(sp.Position))cloud.push_back({V(sp.Position), q->UniqueID, -1.});
    }
    return cloud;
  }
  struct VtxConn { double dist = -1, cover = -1, maxGap = -1, bragg = -1; int owners = -1, nPts = -1; };
  inline VtxConn VertexConnection(const TVector3& vtx, const TVector3& start, Int_t self, const std::vector<VcPoint>& cloud){
    VtxConn v;
    const TVector3 d = start - vtx;
    const double L = d.Mag();
    v.dist = L;
    if(L < kVcSkipV + 2.)return v;   // demasiado cerca: sin corredor que mirar (dist sí vale)
    const TVector3 u = d.Unit();
    const int nb = int(L - kVcSkipV);   // pasos de 1 cm desde kVcSkipV hasta S
    std::vector<char> hit(nb, 0);
    std::set<Int_t> owners;
    double sb = 0; int nbg = 0; v.nPts = 0;
    for(const VcPoint& p : cloud){
      if(p.owner == self)continue;
      const TVector3 w = p.x - vtx;
      const double t = w.Dot(u);
      if(t < kVcSkipV || t > L)continue;
      if((w - t*u).Mag() > kVcRadius)continue;
      v.nPts++;
      hit[std::min(nb - 1, int(t - kVcSkipV))] = 1;
      owners.insert(p.owner);
      if(t > L - kVcBragg && p.dEdx > 0){ sb += p.dEdx; nbg++; }
    }
    int nOn = 0, gap = 0, maxGap = 0;
    for(char h : hit){ if(h){ nOn++; gap = 0; } else maxGap = std::max(maxGap, ++gap); }
    v.cover = double(nOn)/nb; v.maxGap = maxGap; v.owners = owners.size(); v.bragg = nbg ? sb/nbg : -1;
    return v;
  }

  // Conexión topológica (23-nov): la K de 1 GeV/c recorre > 1 m y se curva, así que el corredor recto V -> S falla
  // aunque Pandora la haya reconstruido. Se busca el PFParticle (no el candidato) con un extremo cerca de V y el
  // otro cerca de un extremo E del mu, minimizando |extremo - V| + |otro extremo - E|. Se hace con los dos
  // extremos del mu (las K cortas ponen su pico de Bragg en el inicio del mu y la orientación por dE/dx falla).
  struct Link { double dV = -1, dE = -1, length = -1, dEdxE = -1; int isTrack = -1; Int_t uid = -1; };
  inline Link BestLink(const TVector3& vtx, const TVector3& e, const AnaParticlePD* self, const std::vector<AnaParticlePD*>& parts){
    Link l; double best = 1e9;
    for(const AnaParticlePD* q : parts){
      if(q == self || !Valid(q->PositionStart) || !Valid(q->PositionEnd))continue;
      const TVector3 qa = V(q->PositionStart), qb = V(q->PositionEnd);
      for(int o = 0; o < 2; o++){   // o = 0: qa en V y qb en E; o = 1 al revés
        const TVector3& nv = o ? qb : qa; const TVector3& ne = o ? qa : qb;
        const double sc = (nv - vtx).Mag() + (ne - e).Mag();
        if(sc >= best)continue;
        best = sc; l.dV = (nv - vtx).Mag(); l.dE = (ne - e).Mag(); l.length = q->Length; l.isTrack = q->Type == AnaParticlePD::kTrack; l.uid = q->UniqueID;
        double dA, dB; int nA, nB; EndsdEdx(q, dA, dB, nA, nB); l.dEdxE = o ? dA : dB;   // dE/dx en el extremo que toca al mu
      }
    }
    return l;
  }

  // ---- Longitud de los candidatos tipo shower desde sus space points (23-nov) ----
  // La traza forzada de los mu que Pandora llama shower no da la longitud (forzada/true 0,25-1,6). Con los space
  // points del PFParticle (uno por hit 2D): eje principal (PCA), rehecho con los puntos a < kShCore cm del eje (el
  // Michel no lo tuerce); longitud del núcleo entre los cuantiles 2 % y 98 % a lo largo del eje; y longitud hasta el
  // pico de Bragg: perfil de carga del plano de colección en pasos de 1 cm, la parada es el extremo con más carga en
  // sus últimos 5 cm, el pico el máximo de la media de 3 cm en los últimos 15 cm, y la longitud va del otro extremo
  // al pico (sin el Michel). lin = mayor valor propio / suma (1 = recta; para separar de cascadas EM).
  const double kShCore = 3.;
  struct ShowerLen { int nSP = -1; double lin = -1, lenExt = -1, lenCore = -1, lenBragg = -1, braggRatio = -1; };
  inline bool PCAxis(const std::vector<TVector3>& pts, TVector3& c, TVector3& u, double& lin){
    if(pts.size() < 3)return false;
    c = TVector3(0, 0, 0); for(const TVector3& x : pts)c += x; c *= 1./pts.size();
    TMatrixDSym m(3);
    for(const TVector3& x : pts){ const TVector3 d = x - c; const double v[3] = {d.X(), d.Y(), d.Z()}; for(int i = 0; i < 3; i++)for(int j = 0; j < 3; j++)m(i, j) += v[i]*v[j]; }
    TMatrixDSymEigen e(m);
    const TVectorD& ev = e.GetEigenValues();   // en orden decreciente
    const TMatrixD& vec = e.GetEigenVectors();
    u = TVector3(vec(0, 0), vec(1, 0), vec(2, 0)).Unit();
    const double sum = ev(0) + ev(1) + ev(2);
    lin = sum > 0 ? ev(0)/sum : -1;
    return true;
  }
  inline double QuantileOf(std::vector<double> v, double f){ std::sort(v.begin(), v.end()); return v[std::min(v.size() - 1, size_t(f*v.size()))]; }
  inline ShowerLen ShowerLength(const AnaParticlePD* p){
    ShowerLen s;
    s.nSP = p->SpacePoints.size();
    if(s.nSP < 5)return s;
    std::map<Int_t, double> q2;   // UniqueID -> carga, plano de colección
    for(const AnaHitPD& h : p->HitsNoCalo[2])q2[h.UniqueID] = h.Integral;
    for(const AnaHitPD& h : p->Hits[2])q2[h.UniqueID] = h.Integral;
    std::vector<TVector3> pts; std::vector<double> q;
    for(const AnaSpacePointPD& sp : p->SpacePoints){
      if(!Valid(sp.Position))continue;
      pts.push_back(V(sp.Position));
      auto it = q2.find(sp.HitID); q.push_back(it == q2.end() ? -1. : it->second);
    }
    TVector3 c, u; double lin;
    if(!PCAxis(pts, c, u, lin))return s;
    s.lin = lin;
    std::vector<TVector3> core;
    for(const TVector3& x : pts){ const TVector3 d = x - c; if((d - d.Dot(u)*u).Mag() < kShCore)core.push_back(x); }
    if(core.size() >= 5){ TVector3 c2, u2; double l2; if(PCAxis(core, c2, u2, l2)){ c = c2; u = u2; } }
    std::vector<double> sAll, sCore; std::vector<std::pair<double, double>> prof;
    for(size_t i = 0; i < pts.size(); i++){
      const TVector3 d = pts[i] - c; const double t = d.Dot(u);
      sAll.push_back(t);
      if((d - t*u).Mag() < kShCore){ sCore.push_back(t); if(q[i] > 0)prof.push_back({t, q[i]}); }
    }
    if(sCore.size() < 5)return s;
    s.lenExt = QuantileOf(sAll, 0.98) - QuantileOf(sAll, 0.02);
    const double lo = QuantileOf(sCore, 0.02), hi = QuantileOf(sCore, 0.98);
    s.lenCore = hi - lo;
    if(prof.size() < 10 || s.lenCore < 6)return s;
    const int nb = int(std::ceil(s.lenCore));
    std::vector<double> qb(nb, 0.);
    for(const auto& x : prof){ const int b = int(x.first - lo); if(b >= 0 && b < nb)qb[b] += x.second; }
    auto meanBins = [&](int a, int b){ double sm = 0; int n = 0; for(int i = std::max(0, a); i < std::min(nb, b); i++){ sm += qb[i]; n++; } return n ? sm/n : 0.; };
    const bool stopHi = meanBins(nb - 5, nb) >= meanBins(0, 5);
    std::vector<double> nz; for(double x : qb)if(x > 0)nz.push_back(x);
    const double med = nz.empty() ? 0. : QuantileOf(nz, 0.5);
    int best = -1; double bq = -1;
    const int i0 = stopHi ? std::max(1, nb - 15) : 1, i1 = stopHi ? nb - 1 : std::min(nb - 1, 15);
    for(int i = i0; i < i1; i++){ const double m3 = (qb[i-1] + qb[i] + qb[i+1])/3.; if(m3 > bq){ bq = m3; best = i; } }
    if(best < 0)return s;
    const double sb = lo + best + 0.5;
    s.lenBragg = stopHi ? sb - lo : hi - sb;
    s.braggRatio = med > 0 ? bq/med : -1;
    return s;
  }

  // Hijas de Pandora del candidato
  struct Dau { int nTrk = 0, nShw = 0, hits = 0, nearStop = 0; double maxMichel = -1; };
  inline Dau Daughters(const AnaParticlePD* p, const TVector3& stop, const std::map<UInt_t, const AnaParticlePD*>& byUID){
    Dau d;
    for(UInt_t id : p->DaughtersIDs){
      auto it = byUID.find(id);
      if(it == byUID.end())continue;
      const AnaParticlePD* q = it->second;
      (q->Type == AnaParticlePD::kTrack ? d.nTrk : d.nShw)++;
      d.hits += q->Hits[2].size() + q->HitsNoCalo[2].size();
      d.maxMichel = std::max(d.maxMichel, double(q->CNNscore[2]));
      double dist = 1e9;
      if(Valid(q->PositionStart))dist = std::min(dist, (V(q->PositionStart) - stop).Mag());
      if(Valid(q->PositionEnd))dist = std::min(dist, (V(q->PositionEnd) - stop).Mag());
      if(dist < kNear)d.nearStop++;
    }
    return d;
  }

  // Signal muon (nota técnica, sec. "The signal muon"): mu+ de K+ -> mu+ nu de una K+ en reposo. En la truth: mu+
  // creado en una desintegración, no cósmico, cuya madre es una K+ que acaba en Decay con momento final cero
  // (desintegración en reposo de Geant4). El canal mu nu se identifica por el momento del mu, que en reposo es
  // 235,5 MeV/c (K -> pi0 mu nu da como mucho ~215 MeV/c): no es un corte de la definición.
  inline bool IsSignalMuon(const AnaTrueParticleB* mu, const std::map<int, AnaTrueParticleB*>& byID){
    if(!mu || mu->PDG != -13 || mu->ProcessStart != AnaTrueParticleB::Decay || static_cast<const AnaTrueParticlePD*>(mu)->Origin == 2)return false;
    auto itK = byID.find(mu->ParentID);
    if(itK == byID.end())return false;
    const AnaTrueParticleB* k = itK->second;
    if(k->PDG != 321 || k->ProcessEnd != AnaTrueParticleB::Decay || k->MomentumEnd > 1e-6)return false;   // en reposo
    return mu->Momentum > 0.230 && mu->Momentum < 0.241;   // canal mu nu (GeV/c)
  }

  // Señal de cada medida, en capas sobre el signal muon: C cualquier signal muon; B la K viene de la cascada del haz
  // que dispara (KOrigin 0-3); A la K se produce en la primera interacción (KOrigin 0-1)
  enum Measurement { kMeasC = 0, kMeasB = 1, kMeasA = 2 };
  inline bool IsSignalOf(int measurement, bool signalMuon, int kOrigin){
    if(!signalMuon)return false;
    if(measurement == kMeasB)return kOrigin >= 0 && kOrigin <= 3;
    if(measurement == kMeasA)return kOrigin == 0 || kOrigin == 1;
    return true;
  }

  // Origen de la K de un mu+ (como kaon_origin_minitrees.C): 0 primera interacción del haz que dispara,
  // 1 primera interacción con K+ -> K+ intermedias, 2 interacción secundaria, 3 la K+ del haz tras K+ -> K+,
  // 4 cadena cortada (falta un antecesor de la K), 6 otra raíz (no el haz que dispara), -1 no es un mu+ de
  // desintegración no cósmico con madre K+ (si la madre no está guardada no es señal: -1)
  inline int KOrigin(const AnaTrueParticleB* mu, const AnaTrueParticleB* beamTrue, const std::map<int, AnaTrueParticleB*>& byID){
    if(!mu || mu->PDG != -13 || mu->ProcessStart != AnaTrueParticleB::Decay || static_cast<const AnaTrueParticlePD*>(mu)->Origin == 2)return -1;
    auto itK = byID.find(mu->ParentID);
    if(itK == byID.end())return -1;
    if(itK->second->PDG != 321)return -1;
    std::vector<const AnaTrueParticleB*> chain{itK->second};
    const AnaTrueParticleB* cur = itK->second;
    for(int depth = 0; depth < 1000 && cur->ParentID != 0; depth++){
      auto it = byID.find(cur->ParentID);
      if(it == byID.end())return 4;
      cur = it->second; chain.push_back(cur);
    }
    if(!beamTrue || cur->ID != beamTrue->ID)return 6;
    size_t i0 = 0;
    while(i0 + 1 < chain.size() && chain[i0 + 1]->PDG == 321)i0++;
    if(i0 + 1 == chain.size())return 3;
    if(i0 + 2 == chain.size())return i0 == 0 ? 0 : 1;
    return 2;
  }


  // ¿Desciende la partícula id de ancestor? (por ParentID; como kfate::DescendsFrom de las macros)
  inline bool IsDescendantOf(int id, int ancestor, const std::map<int, AnaTrueParticleB*>& byID){
    for(int depth = 0; depth < 1000; depth++){
      if(id == ancestor)return true;
      auto it = byID.find(id);
      if(it == byID.end() || it->second->ParentID == 0)return false;
      id = it->second->ParentID;
    }
    return false;
  }

  // Variables de un candidato. "A" = PositionStart y "B" = PositionEnd de Pandora (o de la traza forzada); la parada es
  // el extremo de mayor dE/dx (stopAtB). -1 o -999: no se puede calcular. t*: truth (0 o -1 sin truth).
  struct Candidate {
    const AnaParticlePD* part = nullptr;   // no se guarda en los árboles
    Int_t hasSlice, inSlice, complete, isPrimary, nDaughters;
    Float_t length, pRange, chi2p, chi2mu, ndf, cnnTrk, cnnEM, cnnMichel;
    Float_t facesA, facesB, axisA, axisB, sliceA, sliceB, vtxA, vtxB;
    Float_t dEdxA, dEdxB; Int_t nHitsA, nHitsB, nHits2;
    Int_t actPartsA, actPartsB, actHitsA, actHitsB; Float_t actMichelA, actMichelB;
    Int_t clearCosmic, hasT0, cosmicCriteria; Float_t beamCosmicScore;
    Float_t chi2muF, chi2muR, chi2pF, chi2pR, michelOwnA, michelOwnB, intPeak, intPeakRR, dirYStop; Int_t stopAtB;
    // añadidas el 24-nov: chi2 de mu recortando la parada (Chi2BestTrim) con la parada elegida (trim, trimCm) y el
    // mejor de las dos orientaciones (trimBoth); -1 sin chi2
    Float_t chi2muTrim, chi2muTrimCm, chi2muTrimBoth;
    // añadidas el 24-nov: Michel en los hits propios de la parada (CNN media y hits con CNN > 0,5 a < kBragg cm)
    Float_t michelOwnStop; Int_t nMichelOwnStop;
    // truth
    Int_t tPDG, tParentPDG, tProcStart, tOrigin, tFromBeam, tSignal, tGen, tFromBeamK; Float_t tKMomEnd, tPur, tEff, tLength;
    Int_t tSignalMuon = 0;   // signal muon (IsSignalMuon)
    // añadidas el 20-nov
    Int_t isShower;
    Float_t startX, startY, startZ, stopX, stopY, stopZ, cosBeam, trMean, trMeanBody;
    Int_t arN[3], arNMichel[3]; Float_t arMichel[3], arCharge[3];
    Int_t asN5, asN10; Float_t asCharge5, asCharge10;
    Int_t kN, kIsTrack; Float_t kDEdx, kLength, kPRange;
    Int_t nDauTrk, nDauShw, dauHits, dauNearStop; Float_t dauMaxMichel;
    Int_t tKOrigin; Float_t tStartZ, tEndZ, tMom;
    // añadidas el 23-nov: conexión con el vértice primario (vc*, -1 sin vértice) y error del vértice (truth)
    Float_t vcDist, vcCover, vcMaxGap, vcBragg; Int_t vcOwners, vcNPts; Float_t tVtxErr;
    // conexión topológica (lk*): con el inicio (S, orientación por dE/dx) y con la parada (P); distancia de cada extremo
    // del mu a V; truth: ¿el PFParticle del enlace con S es la K madre del mu? (1 sí, 0 no, -1 sin enlace o sin truth)
    Float_t lkSdV, lkSdE, lkSLen, lkSdEdx; Int_t lkSIsTrack; Float_t lkPdV, lkPdE, lkPLen, lkPdEdx; Int_t lkPIsTrack;
    Float_t dStartV, dStopV; Int_t tLkSIsK, tLkPIsK;
    // MC, añadidas el 23-nov: lo mismo con el vértice perfecto (EventContext::tvtx; -1 sin él) y la posición true de la
    // interacción del haz que dispara (final de la partícula true del haz; -999 sin ella)
    Float_t lktSdV, lktSdE, lktPdV, lktPdE, dStartTV, dStopTV; Float_t tBeamEndX, tBeamEndY, tBeamEndZ;
    // candidatos tipo shower: longitud desde los space points (sh*; -1 en las trazas o sin space points)
    Int_t shNSP; Float_t shLin, shLenExt, shLenCore, shLenBragg, shBraggRatio, shPRangeCore, shPRangeBragg;
  };

  // Lo que se calcula una vez por evento
  struct EventContext {
    std::vector<AnaParticlePD*> parts, slice;   // todas las partículas del bunch y las de la beam slice
    bool hasSlice = false;
    std::map<UInt_t, const AnaParticlePD*> byUID;
    const std::vector<AnaHitPD>* unassoc = nullptr;
    TVector3 vtx; bool hasVtx = false;           // vértice: final del primario de la beam slice
    std::vector<VcPoint> cloud;                  // nube de puntos 3D (solo con vértice)
    std::map<int, AnaTrueParticleB*> byID;
    const AnaTrueParticleB* beamTrue = nullptr;  // partícula true del haz que dispara
    int isMC = 1;
    // MC: "vértice perfecto" en el sistema de la reconstrucción: el último hit (plano 2, con posición 3D) que la truth
    // atribuye a la partícula del haz, a lo largo del eje del haz. Es donde el haz deja de depositar carga, con las
    // mismas distorsiones (SCE) que las trazas: sirve para saber qué daría un buscador de vértice perfecto.
    TVector3 tvtx; bool hasTVtx = false;
  };

  // ---- Vista normal de los eventos con volcado completo ----
  //
  // Cómo guardaría el módulo el PFParticle p en un evento normal, con la configuración de la producción (CosmicsMode
  // 2, MuCandidates, MuCandidatesWithoutBeamSlice; CreateHighlandMiniTree.fcl): la beam slice, completa; con beam
  // slice en el evento, las jerarquías de cósmicos completas si están cerca de la beam slice (bit 1 de CosmicCriteria)
  // o tienen un candidato a mu (bit 4), y en resumen las demás; sin beam slice, solo las jerarquías con un candidato a
  // mu, completas. CosmicCriteria es de la jerarquía entera (las hijas lo heredan) y el módulo lo rellena también en
  // los eventos con volcado. Sin CosmicCriteria (-999 fuera de la beam slice: minitrees anteriores o CosmicsMode < 2)
  // no se sabe: se deja completa.
  enum NormalSave { kNotSaved = 0, kSummary = 1, kComplete = 2 };
  inline int NormalSaveMode(const AnaParticlePD* p, bool hasSlice){
    if(p->BeamOrigin)return kComplete;
    const int crit = p->CosmicCriteria;
    if(crit < 0)return kComplete;
    if(hasSlice)return (crit & (AnaParticlePD::kCosmicNearBeamSlice | AnaParticlePD::kCosmicMuCandidate)) ? kComplete : kSummary;
    return (crit & AnaParticlePD::kCosmicMuCandidate) ? kComplete : kNotSaved;
  }

  // Copia de p como la guarda el módulo en resumen (fFillDetail = false en CreateHighlandMiniTree): sin hits (ni
  // HitsNoCalo), puntos de trayectoria, space points ni traza forzada, y sin pureza ni completitud. Conserva el resto
  // (NHits, extremos, longitud, chi2, CNN, TrueObject...).
  inline std::shared_ptr<AnaParticlePD> SummaryCopy(const AnaParticlePD* p){
    std::shared_ptr<AnaParticlePD> s(const_cast<AnaParticlePD*>(p)->Clone());
    for(int i = 0; i < 3; i++){ s->Hits[i].clear(); s->HitsNoCalo[i].clear(); }
    s->TrjPoints.clear(); s->SpacePoints.clear(); s->ForcedTrack.clear();
    s->TruePur = -999; s->TrueEff = -999;   // kFloatUnassigned (no está en ninguna cabecera)
    return s;
  }

  // Vista normal de un evento con volcado completo: las partículas que el módulo guardaría completas, tal cual; las
  // que guardaría en resumen, sustituidas por su SummaryCopy (que se guardan en copies, para que vivan lo que la
  // vista); las que no guardaría, fuera. El orden es el de parts. En un evento normal, la vista es parts.
  inline void MakeNormalView(const std::vector<AnaParticlePD*>& parts, bool fullDump, std::vector<AnaParticlePD*>& view,
                             std::vector<std::shared_ptr<AnaParticlePD>>& copies){
    view.clear();
    if(!fullDump){ view = parts; return; }
    bool hasSlice = false;
    for(const AnaParticlePD* p : parts)if(p->BeamOrigin){ hasSlice = true; break; }
    for(AnaParticlePD* p : parts){
      const int mode = NormalSaveMode(p, hasSlice);
      if(mode == kComplete)view.push_back(p);
      else if(mode == kSummary){ copies.push_back(SummaryCopy(p)); view.push_back(copies.back().get()); }
    }
  }

  inline void MakeContext(EventContext& ctx, const std::vector<AnaParticlePD*>& parts, const std::vector<AnaHitPD>& unassoc,
                          const std::vector<AnaTrueParticleB*>& trueParts, const AnaTrueParticleB* beamTrue, int isMC){
    ctx = EventContext();
    ctx.parts = parts;
    for(AnaParticlePD* p : parts)if(p->BeamOrigin)ctx.slice.push_back(p);
    ctx.hasSlice = !ctx.slice.empty();
    for(const AnaParticlePD* p : parts)ctx.byUID[p->UniqueID] = p;
    ctx.unassoc = &unassoc;
    if(ctx.hasSlice && Valid(ctx.slice[0]->PositionEnd)){ ctx.vtx = V(ctx.slice[0]->PositionEnd); ctx.hasVtx = true; }
    if(ctx.hasVtx)ctx.cloud = VcCloud(parts);
    for(AnaTrueParticleB* t : trueParts)ctx.byID[t->ID] = t;
    ctx.beamTrue = beamTrue;
    ctx.isMC = isMC;
    if(beamTrue){
      double best = -1e9;
      for(const AnaParticlePD* p : parts)
        for(const AnaHitPD& h : p->Hits[2]){
          if(h.TrueID != beamTrue->ID || h.Position_NoSCE.X() < -900 || h.Position_NoSCE.Mag2() == 0)continue;
          const double t = h.Position_NoSCE.Dot(kAxisDir);
          if(t > best){ best = t; ctx.tvtx = h.Position_NoSCE; ctx.hasTVtx = true; }
        }
    }
  }

  // Distancia a la beam slice (segmentos de sus partículas); -1 sin beam slice
  inline double DistToSlice(const TVector3& x, const EventContext& ctx){
    double d = 1e9;
    for(const AnaParticlePD* q : ctx.slice){
      if(Valid(q->PositionStart) && Valid(q->PositionEnd))d = std::min(d, DistToSegment(x, V(q->PositionStart), V(q->PositionEnd)));
      else if(Valid(q->PositionStart))d = std::min(d, (x - V(q->PositionStart)).Mag());
    }
    return d < 1e8 ? d : -1.;
  }

  // Señal antigua (tSignal): mu+ cuya madre es una K+ que acaba en Decay y desciende del haz que dispara; kMomEnd en MeV/c
  inline bool IsOldSignal(const AnaTrueParticleB* t, const EventContext& ctx, double& kMomEnd){
    kMomEnd = -1;
    if(!t || t->PDG != -13 || !ctx.beamTrue)return false;
    auto ip = ctx.byID.find(t->ParentID);
    if(ip == ctx.byID.end() || ip->second->PDG != 321 || ip->second->ProcessEnd != AnaTrueParticleB::Decay)return false;
    if(!IsDescendantOf(ip->second->ID, ctx.beamTrue->ID, ctx.byID))return false;
    kMomEnd = 1000.*ip->second->MomentumEnd;
    return true;
  }

  // Variables de la partícula p como candidato; false si no pasa la preselección
  inline bool Compute(const AnaParticlePD* p, const EventContext& ctx, Candidate& c){
    const std::vector<AnaParticlePD*>& parts = ctx.parts;
    const std::vector<AnaHitPD>& unassoc = *ctx.unassoc;
    const auto& byID = ctx.byID; const auto& byUID = ctx.byUID;
    const AnaTrueParticleB* beamTrue = ctx.beamTrue;
    const bool hasVtx = ctx.hasVtx; const TVector3& vtx = ctx.vtx;
    const int isMC = ctx.isMC;
    // trazas, o PFParticles sin traza con traza forzada (candidatos tipo shower)
    const bool isTrk = p->Type == AnaParticlePD::kTrack;
    const AnaForcedTrackPD* ft = (!isTrk && !p->ForcedTrack.empty()) ? &p->ForcedTrack[0] : nullptr;
    if(!isTrk && !ft)return false;
    const double L = isTrk ? p->Length : ft->Length;
    const Float_t* pa = isTrk ? p->PositionStart : ft->PositionStart;
    const Float_t* pb = isTrk ? p->PositionEnd : ft->PositionEnd;
    if(L <= 0 || !Valid(pa) || !Valid(pb))return false;
    const double pr = MuRangeP(L);
    if(pr < kPMin || pr > kPMax)return false;
    const TVector3 a = V(pa), b = V(pb);
    c = Candidate();
    c.part = p;
    c.isShower = !isTrk;
    c.hasSlice = ctx.hasSlice; c.inSlice = p->BeamOrigin; c.complete = HasHits(p);
    c.isPrimary = p->isPandora; c.nDaughters = p->DaughtersIDs.size();
    c.length = L; c.pRange = pr;
    c.chi2p = isTrk ? p->Chi2Proton : ft->Chi2Proton; c.chi2mu = isTrk ? p->Chi2Muon : ft->Chi2Muon; c.ndf = isTrk ? p->Chi2ndf : ft->Chi2ndf;
    c.cnnTrk = p->CNNscore[0]; c.cnnEM = p->CNNscore[1]; c.cnnMichel = p->CNNscore[2];
    c.facesA = DistToFaces(a); c.facesB = DistToFaces(b);
    c.axisA = DistToAxis(a, isMC); c.axisB = DistToAxis(b, isMC);
    c.sliceA = p->BeamOrigin ? 0 : DistToSlice(a, ctx); c.sliceB = p->BeamOrigin ? 0 : DistToSlice(b, ctx);
    c.vtxA = hasVtx ? (a - vtx).Mag() : -1; c.vtxB = hasVtx ? (b - vtx).Mag() : -1;
    double dA, dB; EndsdEdx(p, dA, dB, c.nHitsA, c.nHitsB); c.dEdxA = dA; c.dEdxB = dB;
    c.nHits2 = p->Hits[2].size();
    const Activity actA = ActivityAt(a, p, parts), actB = ActivityAt(b, p, parts);
    c.actPartsA = actA.nParts; c.actPartsB = actB.nParts; c.actHitsA = actA.nHits; c.actHitsB = actB.nHits;
    c.actMichelA = actA.michel; c.actMichelB = actB.michel;
    { double f, r; Chi2Both(p, 13, f, r); c.chi2muF = f; c.chi2muR = r; Chi2Both(p, 2212, f, r); c.chi2pF = f; c.chi2pR = r; }
    { double mA, mB; OwnMichel(p, mA, mB); c.michelOwnA = mA; c.michelOwnB = mB; }
    { double pk, rr; InteriorPeak(p, pk, rr); c.intPeak = pk; c.intPeakRR = rr; }
    // extremo de parada: el de mayor dE/dx (pico de Bragg); sin dE/dx en los dos, la orientación
    // con mejor chi2 de mu; sin nada, la de Pandora
    if(c.dEdxA > 0 && c.dEdxB > 0)c.stopAtB = c.dEdxB >= c.dEdxA;
    else c.stopAtB = !(c.chi2muR >= 0 && (c.chi2muF < 0 || c.chi2muR < c.chi2muF));
    c.dirYStop = c.stopAtB ? (b - a).Unit().Y() : (a - b).Unit().Y();
    { double cm, cmO; c.chi2muTrim = Chi2BestTrim(p, 13, c.stopAtB, cm); c.chi2muTrimCm = cm;
      const double other = Chi2BestTrim(p, 13, !c.stopAtB, cmO);
      c.chi2muTrimBoth = (other >= 0 && (c.chi2muTrim < 0 || other < c.chi2muTrim)) ? other : c.chi2muTrim; }
    c.michelOwnStop = c.stopAtB ? c.michelOwnB : c.michelOwnA;
    c.nMichelOwnStop = OwnMichelHits(p, c.stopAtB);
    c.clearCosmic = p->IsClearCosmic; c.hasT0 = p->T0 > -900; c.cosmicCriteria = p->CosmicCriteria; c.beamCosmicScore = p->BeamCosmicScore;
    const AnaTrueParticlePD* t = static_cast<const AnaTrueParticlePD*>(p->TrueObject);
    c.tPDG = t ? t->PDG : 0;
    { auto ip = t ? byID.find(t->ParentID) : byID.end(); c.tParentPDG = ip != byID.end() ? ip->second->PDG : 0; }   // ParentPDG no lo rellena el módulo
    c.tFromBeamK = t && beamTrue && beamTrue->PDG == 321 && t->PDG == -13 && t->ParentID == beamTrue->ID; c.tProcStart = t ? int(t->ProcessStart) : -1;
    c.tOrigin = t ? t->Origin : -1; c.tGen = t ? t->Generation : -999;
    c.tFromBeam = (t && beamTrue) ? IsDescendantOf(t->ID, beamTrue->ID, byID) : 0;
    double km; c.tSignal = IsOldSignal(t, ctx, km);
    // momento final de la K madre para cualquier mu cuya madre sea una K+ (no solo la señal antigua)
    c.tKMomEnd = -1;
    if(t){ auto ik = ctx.byID.find(t->ParentID); if(ik != ctx.byID.end() && ik->second->PDG == 321)c.tKMomEnd = 1000.*ik->second->MomentumEnd; }
    c.tSignalMuon = IsSignalMuon(t, ctx.byID);
    c.tPur = p->TruePur; c.tEff = p->TrueEff; c.tLength = t ? (V(t->PositionEnd) - V(t->Position)).Mag() : -1;   // distancia inicio-fin
    // añadidas el 20-nov
    const TVector3 stop = c.stopAtB ? b : a, start = c.stopAtB ? a : b;
    c.startX = start.X(); c.startY = start.Y(); c.startZ = start.Z(); c.stopX = stop.X(); c.stopY = stop.Y(); c.stopZ = stop.Z();
    c.cosBeam = (stop - start).Unit().Dot(kAxisDir);
    c.trMean = TruncMean(p, false); c.trMeanBody = TruncMean(p, true);
    { const Around ar = HitsAround(EndHit(p, c.stopAtB), p, parts, unassoc);
      for(int r = 0; r < 3; r++){ c.arN[r] = ar.n[r]; c.arNMichel[r] = ar.nMichel[r]; c.arMichel[r] = ar.michel[r]; c.arCharge[r] = ar.charge[r]; }
      if(!EndHit(p, c.stopAtB))for(int r = 0; r < 3; r++){ c.arN[r] = c.arNMichel[r] = -1; c.arCharge[r] = -1; }
      const AnaHitPD* hs = EndHit(p, !c.stopAtB);
      const Around as = HitsAround(hs, p, parts, unassoc);
      c.asN5 = hs ? as.n[0] : -1; c.asN10 = hs ? as.n[1] : -1; c.asCharge5 = hs ? as.charge[0] : -1; c.asCharge10 = hs ? as.charge[1] : -1; }
    { const KAtStart k = KAt(start, p, parts); c.kN = k.n; c.kIsTrack = k.isTrack; c.kDEdx = k.dEdx; c.kLength = k.length; c.kPRange = k.pRangeK; }
    { const Dau d = Daughters(p, stop, byUID); c.nDauTrk = d.nTrk; c.nDauShw = d.nShw; c.dauHits = d.hits; c.dauNearStop = d.nearStop; c.dauMaxMichel = d.maxMichel; }
    c.tKOrigin = KOrigin(t, beamTrue, byID);
    c.tStartZ = t ? t->Position[2] : -999; c.tEndZ = t ? t->PositionEnd[2] : -999; c.tMom = t ? 1000.*t->Momentum : -1;
    // añadidas el 23-nov
    { VtxConn vc; if(hasVtx)vc = VertexConnection(vtx, start, p->UniqueID, ctx.cloud);
      c.vcDist = vc.dist; c.vcCover = vc.cover; c.vcMaxGap = vc.maxGap; c.vcBragg = vc.bragg; c.vcOwners = vc.owners; c.vcNPts = vc.nPts;
      c.tVtxErr = (hasVtx && beamTrue) ? (vtx - V(beamTrue->PositionEnd)).Mag() : -1; }
    { Link ls, lp; if(hasVtx){ ls = BestLink(vtx, start, p, parts); lp = BestLink(vtx, stop, p, parts); }
      c.lkSdV = ls.dV; c.lkSdE = ls.dE; c.lkSLen = ls.length; c.lkSdEdx = ls.dEdxE; c.lkSIsTrack = ls.isTrack;
      c.lkPdV = lp.dV; c.lkPdE = lp.dE; c.lkPLen = lp.length; c.lkPdEdx = lp.dEdxE; c.lkPIsTrack = lp.isTrack;
      c.dStartV = hasVtx ? (start - vtx).Mag() : -1; c.dStopV = hasVtx ? (stop - vtx).Mag() : -1;
      // truth del enlace: ¿su partícula true es la madre del mu?
      auto isK = [&](const Link& l) -> int {
        if(l.uid < 0 || !t)return -1;
        auto iq = byUID.find(l.uid); if(iq == byUID.end() || !iq->second->TrueObject)return -1;
        return static_cast<AnaTrueParticleB*>(iq->second->TrueObject)->ID == t->ParentID ? 1 : 0;
      };
      c.tLkSIsK = isK(ls); c.tLkPIsK = isK(lp); }
    { Link ls, lp; if(ctx.hasTVtx){ ls = BestLink(ctx.tvtx, start, p, parts); lp = BestLink(ctx.tvtx, stop, p, parts); }
      c.lktSdV = ls.dV; c.lktSdE = ls.dE; c.lktPdV = lp.dV; c.lktPdE = lp.dE;
      c.dStartTV = ctx.hasTVtx ? (start - ctx.tvtx).Mag() : -1; c.dStopTV = ctx.hasTVtx ? (stop - ctx.tvtx).Mag() : -1;
      c.tBeamEndX = beamTrue ? beamTrue->PositionEnd[0] : -999; c.tBeamEndY = beamTrue ? beamTrue->PositionEnd[1] : -999;
      c.tBeamEndZ = beamTrue ? beamTrue->PositionEnd[2] : -999; }
    { ShowerLen sl; if(!isTrk)sl = ShowerLength(p);
      c.shNSP = sl.nSP; c.shLin = sl.lin; c.shLenExt = sl.lenExt; c.shLenCore = sl.lenCore; c.shLenBragg = sl.lenBragg; c.shBraggRatio = sl.braggRatio;
      c.shPRangeCore = sl.lenCore > 0 ? MuRangeP(sl.lenCore) : -1; c.shPRangeBragg = sl.lenBragg > 0 ? MuRangeP(sl.lenBragg) : -1; }
    return true;
  }
}
#endif
