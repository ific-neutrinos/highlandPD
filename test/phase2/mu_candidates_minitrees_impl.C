// Implementación de mu_candidates_minitrees.C (se compila con ACLiC).
//
// Candidato amplio: PFParticle con traza (Type kTrack), longitud > 0, extremos válidos y momento por
// rango de mu en [100, 600] MeV/c. Todos los PFParticles del evento (beam slice, fuera completos y
// fuera en resumen: estos sin hits, así que sin calorimetría ni actividad).
//
// Extremos: "A" = PositionStart, "B" = PositionEnd tal como los da Pandora (la orientación puede
// estar al revés: para un mu de K parada, el extremo con pico de Bragg es donde se para y el otro
// es donde estaba la K). Para cada extremo: distancia a las caras del volumen activo, al eje del haz,
// a la beam slice, al final del primario; dE/dx medio de los hits del plano 2 a < 5 cm; actividad
// (otros PFParticles con un extremo a < 10 cm, sus hits y su puntuación CNN de Michel media).
//
// Truth (TrueObject: la partícula true que más hits aporta): PDG, PDG de la madre, proceso de
// creación, Origin (2 cósmico, 4 haz), si desciende de la partícula del haz y si es señal: mu+ cuya
// madre es un K+ que acaba en Decay y desciende del haz. kMomEnd: momento final de esa K (MeV/c;
// ~0 = en reposo).
//
// Añadido el 20-nov (las ramas anteriores no cambian para las trazas):
//   - Candidatos tipo shower (isShower): PFParticles sin traza con traza forzada (pandora2Track); longitud,
//     extremos y chi2 de la traza forzada; sin hits de calorimetría (dE/dx, chi2 por orientación y
//     actividad en hits a -1).
//   - Posiciones de los extremos de inicio y de parada (start*, stop*; extremo de parada como stopAtB), coseno
//     con el eje del haz (cosBeam, dirección inicio -> parada) y media truncada del dE/dx (trMean; trMeanBody
//     sin los 5 cm de cada extremo).
//   - Hits del plano 2 alrededor del extremo de parada a < 5, 10 y 20 cm (distancia 2D hilo-tiempo en la
//     misma TPC, todos los PFParticles salvo los hits de calorimetría del candidato, más los no asociados):
//     número, CNN de Michel media, número con CNN de Michel > 0,5 y carga (ar*). Alrededor del inicio: número
//     y carga a < 5 y 10 cm (as*).
//   - Lo que acaba en el inicio (la K): otros PFParticles con un extremo a < 10 cm; el de mayor dE/dx en ese
//     extremo, su longitud y su p de rango con hipótesis de K (k*).
//   - Hijas de Pandora: trazas, showers, hits, Michel máxima y cuántas empiezan o acaban a < 10 cm de la
//     parada (dau*).
//   - Truth: tKOrigin (origen de la K del mu+, ver KOrigin y kaon_origin_minitrees.C), z true de inicio y fin
//     y momento true.
//   - Árbol sig: además de los mu+ de antes (isSignalOld), todos los mu+ de desintegración no cósmicos con
//     p true en [230, 241] MeV/c (K en reposo), con tKOrigin; tipo del PFParticle emparejado y p de rango de
//     la traza forzada si es un shower.
// Añadido el 23-nov: conexión con el vértice primario (vc*: longitud V -> inicio, fracción del segmento con puntos 3D
// a < 2 cm, tramo vacío más largo, dE/dx en los últimos 5 cm antes del inicio, cuántos PFParticles; ver
// VertexConnection) y tVtxErr (distancia del vértice reco al final true del haz). Conexión topológica (lk*: el
// PFParticle que une V con el inicio o con la parada del mu; ver BestLink) y distancias de los extremos a V.
// Longitud de los candidatos tipo shower desde sus space points (sh*: eje principal, núcleo, hasta el pico de
// Bragg, linealidad; ver ShowerLength).
#include "kaon_fate_minitrees_impl.C"
#include "pdAnalysisUtils.hxx"
#include "Parameters.hxx"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include "TVectorD.h"
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <vector>

namespace mucand {
  using kfate::InputFiles; using kfate::TreeName; using kfate::DescendsFrom;

  // Como el módulo (CreateHighlandMiniTree.fcl)
  const double kActive[6] = {-360., 360., 0., 607.5, 0., 695.};
  const TVector3 kAxisStart[2] = {TVector3(-28.6, 424.6, 0.96), TVector3(-30.1, 423.0, 0.0)};   // datos, MC
  const TVector3 kAxisDir = TVector3(-0.18, -0.20, 0.96).Unit();
  const double kAxisLength = 700.;
  const double kPMin = 100., kPMax = 600.;      // MeV/c, candidato amplio
  const double kWinMin = 180., kWinMax = 300.;  // MeV/c, ventana del módulo (MuCandidateLength)
  const double kNear = 10.;                     // cm, actividad en un extremo
  const double kBragg = 5.;                     // cm, dE/dx en el extremo

  bool Valid(const Float_t* x){ return x[0] > -900 && x[1] > -900 && x[2] > -900; }
  TVector3 V(const Float_t* x){ return TVector3(x[0], x[1], x[2]); }

  double DistToFaces(const TVector3& p){
    return std::min({p.X()-kActive[0], kActive[1]-p.X(), p.Y()-kActive[2], kActive[3]-p.Y(), p.Z()-kActive[4], kActive[5]-p.Z()});
  }
  double DistToSegment(const TVector3& p, const TVector3& a, const TVector3& b){
    const TVector3 ab = b - a;
    const double l2 = ab.Mag2();
    double t = l2 > 0 ? (p - a).Dot(ab)/l2 : 0.;
    t = std::max(0., std::min(1., t));
    return (p - (a + t*ab)).Mag();
  }
  double DistToAxis(const TVector3& p, int isMC){
    return DistToSegment(p, kAxisStart[isMC], kAxisStart[isMC] + kAxisLength*kAxisDir);
  }
  double MuRangeP(double length){   // MeV/c
    return length > 0 ? 1000.*pdAnaUtils::ComputeRangeMomentum(length, 13) : -1;
  }
  bool HasHits(const AnaParticlePD* p){ return !p->Hits[2].empty() || p->NHits > 0; }

  // dE/dx medio de los hits del plano 2 a < kBragg cm de cada extremo (por residual range: RR
  // pequeño = cerca de PositionEnd, en la convención de la calorimetría de Pandora). -1 sin hits
  void EndsdEdx(const AnaParticlePD* p, double& dEdxA, double& dEdxB, int& nA, int& nB){
    double sA = 0, sB = 0; nA = nB = 0;
    const double L = p->Length;
    for(const AnaHitPD& h : p->Hits[2]){
      if(h.dEdx <= 0 || h.dEdx > 100 || h.ResidualRange < 0)continue;
      if(h.ResidualRange < kBragg){ sB += h.dEdx; nB++; }
      if(h.ResidualRange > L - kBragg){ sA += h.dEdx; nA++; }
    }
    dEdxA = nA ? sA/nA : -1; dEdxB = nB ? sB/nB : -1;
  }

  // chi2 de una hipótesis con la orientación de Pandora (fwd: se para en PositionEnd) y al revés
  // (rev: residual range RR' = L - RR, se para en PositionStart). first = chi2/ndf (-1 si no hay)
  void Chi2Both(const AnaParticlePD* p, int pdg, double& fwd, double& rev){
    fwd = rev = -1;
    if(p->Hits[2].size() < 3)return;
    const std::pair<double,int> f = pdAnaUtils::Chi2PID(*p, pdg, 1e9);
    if(f.second > 0)fwd = f.first/f.second;
    AnaParticlePD* q = static_cast<AnaParticlePD*>(const_cast<AnaParticlePD*>(p)->Clone());
    for(AnaHitPD& h : q->Hits[2])if(h.ResidualRange >= 0)h.ResidualRange = std::max(0.f, Float_t(p->Length - h.ResidualRange));
    const std::pair<double,int> r = pdAnaUtils::Chi2PID(*q, pdg, 1e9);
    if(r.second > 0)rev = r.first/r.second;
    delete q;
  }

  // Puntuación CNN de Michel media de los hits propios (plano 2) a < kBragg cm de cada extremo
  void OwnMichel(const AnaParticlePD* p, double& mA, double& mB){
    double sA = 0, sB = 0; int nA = 0, nB = 0;
    for(const AnaHitPD& h : p->Hits[2]){
      if(h.ResidualRange < 0 || h.CNN[2] < 0)continue;
      if(h.ResidualRange < kBragg){ sB += h.CNN[2]; nB++; }
      if(h.ResidualRange > p->Length - kBragg){ sA += h.CNN[2]; nA++; }
    }
    mA = nA ? sA/nA : -1; mB = nB ? sB/nB : -1;
  }

  // Pico de Bragg interior: máximo de la media de 5 hits consecutivos (plano 2, por residual range)
  // con todos a > kBragg cm de los dos extremos. Trazas K+mu fundidas: el pico de la K a ~54 cm del
  // final del mu
  void InteriorPeak(const AnaParticlePD* p, double& peak, double& rrPeak){
    peak = -1; rrPeak = -1;
    std::vector<std::pair<double,double>> h;   // (RR, dE/dx)
    for(const AnaHitPD& x : p->Hits[2])if(x.ResidualRange >= 0 && x.dEdx > 0 && x.dEdx < 100)h.push_back({x.ResidualRange, x.dEdx});
    std::sort(h.begin(), h.end());
    const int w = 5;
    for(int i = 0; i + w <= int(h.size()); i++){
      if(h[i].first < kBragg || h[i+w-1].first > p->Length - kBragg)continue;
      double m = 0; for(int k = 0; k < w; k++)m += h[i+k].second; m /= w;
      if(m > peak){ peak = m; rrPeak = h[i+w/2].first; }
    }
  }

  struct Activity { int nParts = 0, nHits = 0; double michel = -1; };
  // Otros PFParticles con algún extremo a < kNear cm de x: cuántos, sus hits del plano 2 y la
  // puntuación CNN de Michel media de esos hits
  Activity ActivityAt(const TVector3& x, const AnaParticlePD* self, const std::vector<AnaParticlePD*>& parts){
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
  const AnaHitPD* EndHit(const AnaParticlePD* p, bool atB){
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
  Around HitsAround(const AnaHitPD* h0, const AnaParticlePD* self, const std::vector<AnaParticlePD*>& parts, const std::vector<AnaHitPD>& unassoc){
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
  double TruncMean(const AnaParticlePD* p, bool body){
    std::vector<double> v;
    for(const AnaHitPD& h : p->Hits[2]){
      if(h.dEdx <= 0 || h.dEdx > 100 || h.ResidualRange < 0)continue;
      if(body && (h.ResidualRange < kBragg || h.ResidualRange > p->Length - kBragg))continue;
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
  KAtStart KAt(const TVector3& x, const AnaParticlePD* self, const std::vector<AnaParticlePD*>& parts){
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
  std::vector<VcPoint> VcCloud(const std::vector<AnaParticlePD*>& parts){
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
  VtxConn VertexConnection(const TVector3& vtx, const TVector3& start, Int_t self, const std::vector<VcPoint>& cloud){
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
  Link BestLink(const TVector3& vtx, const TVector3& e, const AnaParticlePD* self, const std::vector<AnaParticlePD*>& parts){
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
  bool PCAxis(const std::vector<TVector3>& pts, TVector3& c, TVector3& u, double& lin){
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
  double QuantileOf(std::vector<double> v, double f){ std::sort(v.begin(), v.end()); return v[std::min(v.size() - 1, size_t(f*v.size()))]; }
  ShowerLen ShowerLength(const AnaParticlePD* p){
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
  Dau Daughters(const AnaParticlePD* p, const TVector3& stop, const std::map<UInt_t, const AnaParticlePD*>& byUID){
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

  // Origen de la K de un mu+ (como kaon_origin_minitrees.C): 0 primera interacción del haz que dispara,
  // 1 primera interacción con K+ -> K+ intermedias, 2 interacción secundaria, 3 la K+ del haz tras K+ -> K+,
  // 4 cadena cortada (falta un antecesor de la K), 5 madre del mu no guardada, 6 otra raíz (no el haz que
  // dispara), -1 no es un mu+ de desintegración no cósmico con madre K+ (o madre no guardada)
  int KOrigin(const AnaTrueParticleB* mu, const AnaTrueParticleB* beamTrue, const std::map<int, AnaTrueParticleB*>& byID){
    if(!mu || mu->PDG != -13 || mu->ProcessStart != AnaTrueParticleB::Decay || static_cast<const AnaTrueParticlePD*>(mu)->Origin == 2)return -1;
    auto itK = byID.find(mu->ParentID);
    if(itK == byID.end())return 5;
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

  // Variables del árbol de candidatos (una entrada por candidato)
  struct Cand {
    Int_t entry, fullDump, hasSlice, inSlice, complete, isPrimary, nDaughters;
    Float_t length, pRange, chi2p, chi2mu, ndf, cnnTrk, cnnEM, cnnMichel;
    Float_t facesA, facesB, axisA, axisB, sliceA, sliceB, vtxA, vtxB;
    Float_t dEdxA, dEdxB; Int_t nHitsA, nHitsB, nHits2;
    Int_t actPartsA, actPartsB, actHitsA, actHitsB; Float_t actMichelA, actMichelB;
    Int_t clearCosmic, hasT0, cosmicCriteria; Float_t beamCosmicScore;
    Float_t chi2muF, chi2muR, chi2pF, chi2pR, michelOwnA, michelOwnB, intPeak, intPeakRR, dirYStop; Int_t stopAtB;
    // truth
    Int_t tPDG, tParentPDG, tProcStart, tOrigin, tFromBeam, tSignal, tGen, tFromBeamK; Float_t tKMomEnd, tPur, tEff, tLength;
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
    // candidatos tipo shower: longitud desde los space points (sh*; -1 en las trazas o sin space points)
    Int_t shNSP; Float_t shLin, shLenExt, shLenCore, shLenBragg, shBraggRatio, shPRangeCore, shPRangeBragg;
  };

  // mu+ true de K+ de la cadena del haz: qué le pasa en la reco
  struct Sig {
    Int_t entry, fullDump, hasSlice, matched, isTrack, inSlice, complete, inWindow, contained; Float_t kMomEnd, pRange, trueLength, axisMin, nearestToEnd;
    // añadidas el 20-nov
    Int_t isSignalOld, tKOrigin, matchedType; Float_t muMom, tStartZ, tEndZ, pRangeShower;
  };

  struct MuCandidates {
    TFile* fout = nullptr;
    TTree* tc = nullptr; TTree* ts = nullptr;
    Cand c; Sig s;
    Long64_t entry = 0;
    long nEvents[2] = {0, 0}, nEventsSlice[2] = {0, 0};   // [normal, volcado]

    void Book(const char* out){
      fout = TFile::Open(out, "RECREATE");
      tc = new TTree("cand", "candidatos amplios a mu de K+ en reposo");
      tc->Branch("entry", &c.entry); tc->Branch("fullDump", &c.fullDump); tc->Branch("hasSlice", &c.hasSlice);
      tc->Branch("inSlice", &c.inSlice); tc->Branch("complete", &c.complete); tc->Branch("isPrimary", &c.isPrimary); tc->Branch("nDaughters", &c.nDaughters);
      tc->Branch("length", &c.length); tc->Branch("pRange", &c.pRange); tc->Branch("chi2p", &c.chi2p); tc->Branch("chi2mu", &c.chi2mu); tc->Branch("ndf", &c.ndf);
      tc->Branch("cnnTrk", &c.cnnTrk); tc->Branch("cnnEM", &c.cnnEM); tc->Branch("cnnMichel", &c.cnnMichel);
      tc->Branch("facesA", &c.facesA); tc->Branch("facesB", &c.facesB); tc->Branch("axisA", &c.axisA); tc->Branch("axisB", &c.axisB);
      tc->Branch("sliceA", &c.sliceA); tc->Branch("sliceB", &c.sliceB); tc->Branch("vtxA", &c.vtxA); tc->Branch("vtxB", &c.vtxB);
      tc->Branch("dEdxA", &c.dEdxA); tc->Branch("dEdxB", &c.dEdxB); tc->Branch("nHitsA", &c.nHitsA); tc->Branch("nHitsB", &c.nHitsB); tc->Branch("nHits2", &c.nHits2);
      tc->Branch("actPartsA", &c.actPartsA); tc->Branch("actPartsB", &c.actPartsB); tc->Branch("actHitsA", &c.actHitsA); tc->Branch("actHitsB", &c.actHitsB);
      tc->Branch("actMichelA", &c.actMichelA); tc->Branch("actMichelB", &c.actMichelB);
      tc->Branch("chi2muF", &c.chi2muF); tc->Branch("chi2muR", &c.chi2muR); tc->Branch("chi2pF", &c.chi2pF); tc->Branch("chi2pR", &c.chi2pR);
      tc->Branch("michelOwnA", &c.michelOwnA); tc->Branch("michelOwnB", &c.michelOwnB); tc->Branch("intPeak", &c.intPeak); tc->Branch("intPeakRR", &c.intPeakRR);
      tc->Branch("dirYStop", &c.dirYStop); tc->Branch("stopAtB", &c.stopAtB);
      tc->Branch("clearCosmic", &c.clearCosmic); tc->Branch("hasT0", &c.hasT0); tc->Branch("cosmicCriteria", &c.cosmicCriteria); tc->Branch("beamCosmicScore", &c.beamCosmicScore);
      tc->Branch("tPDG", &c.tPDG); tc->Branch("tParentPDG", &c.tParentPDG); tc->Branch("tProcStart", &c.tProcStart); tc->Branch("tOrigin", &c.tOrigin);
      tc->Branch("tFromBeam", &c.tFromBeam); tc->Branch("tSignal", &c.tSignal); tc->Branch("tGen", &c.tGen); tc->Branch("tFromBeamK", &c.tFromBeamK); tc->Branch("tKMomEnd", &c.tKMomEnd);
      tc->Branch("tPur", &c.tPur); tc->Branch("tEff", &c.tEff); tc->Branch("tLength", &c.tLength);
      tc->Branch("isShower", &c.isShower);
      tc->Branch("startX", &c.startX); tc->Branch("startY", &c.startY); tc->Branch("startZ", &c.startZ);
      tc->Branch("stopX", &c.stopX); tc->Branch("stopY", &c.stopY); tc->Branch("stopZ", &c.stopZ);
      tc->Branch("cosBeam", &c.cosBeam); tc->Branch("trMean", &c.trMean); tc->Branch("trMeanBody", &c.trMeanBody);
      tc->Branch("arN", c.arN, "arN[3]/I"); tc->Branch("arNMichel", c.arNMichel, "arNMichel[3]/I");
      tc->Branch("arMichel", c.arMichel, "arMichel[3]/F"); tc->Branch("arCharge", c.arCharge, "arCharge[3]/F");
      tc->Branch("asN5", &c.asN5); tc->Branch("asN10", &c.asN10); tc->Branch("asCharge5", &c.asCharge5); tc->Branch("asCharge10", &c.asCharge10);
      tc->Branch("kN", &c.kN); tc->Branch("kIsTrack", &c.kIsTrack); tc->Branch("kDEdx", &c.kDEdx); tc->Branch("kLength", &c.kLength); tc->Branch("kPRange", &c.kPRange);
      tc->Branch("nDauTrk", &c.nDauTrk); tc->Branch("nDauShw", &c.nDauShw); tc->Branch("dauHits", &c.dauHits); tc->Branch("dauNearStop", &c.dauNearStop);
      tc->Branch("dauMaxMichel", &c.dauMaxMichel);
      tc->Branch("tKOrigin", &c.tKOrigin); tc->Branch("tStartZ", &c.tStartZ); tc->Branch("tEndZ", &c.tEndZ); tc->Branch("tMom", &c.tMom);
      tc->Branch("vcDist", &c.vcDist); tc->Branch("vcCover", &c.vcCover); tc->Branch("vcMaxGap", &c.vcMaxGap); tc->Branch("vcBragg", &c.vcBragg);
      tc->Branch("vcOwners", &c.vcOwners); tc->Branch("vcNPts", &c.vcNPts); tc->Branch("tVtxErr", &c.tVtxErr);
      tc->Branch("lkSdV", &c.lkSdV); tc->Branch("lkSdE", &c.lkSdE); tc->Branch("lkSLen", &c.lkSLen); tc->Branch("lkSdEdx", &c.lkSdEdx); tc->Branch("lkSIsTrack", &c.lkSIsTrack);
      tc->Branch("lkPdV", &c.lkPdV); tc->Branch("lkPdE", &c.lkPdE); tc->Branch("lkPLen", &c.lkPLen); tc->Branch("lkPdEdx", &c.lkPdEdx); tc->Branch("lkPIsTrack", &c.lkPIsTrack);
      tc->Branch("dStartV", &c.dStartV); tc->Branch("dStopV", &c.dStopV); tc->Branch("tLkSIsK", &c.tLkSIsK); tc->Branch("tLkPIsK", &c.tLkPIsK);
      tc->Branch("shNSP", &c.shNSP); tc->Branch("shLin", &c.shLin); tc->Branch("shLenExt", &c.shLenExt); tc->Branch("shLenCore", &c.shLenCore);
      tc->Branch("shLenBragg", &c.shLenBragg); tc->Branch("shBraggRatio", &c.shBraggRatio); tc->Branch("shPRangeCore", &c.shPRangeCore); tc->Branch("shPRangeBragg", &c.shPRangeBragg);
      ts = new TTree("sig", "mu+ true de K+ de la cadena del haz");
      ts->Branch("entry", &s.entry); ts->Branch("fullDump", &s.fullDump); ts->Branch("hasSlice", &s.hasSlice); ts->Branch("matched", &s.matched);
      ts->Branch("isTrack", &s.isTrack); ts->Branch("inSlice", &s.inSlice); ts->Branch("complete", &s.complete); ts->Branch("inWindow", &s.inWindow);
      ts->Branch("contained", &s.contained); ts->Branch("kMomEnd", &s.kMomEnd); ts->Branch("pRange", &s.pRange); ts->Branch("trueLength", &s.trueLength); ts->Branch("axisMin", &s.axisMin); ts->Branch("nearestToEnd", &s.nearestToEnd);
      ts->Branch("isSignalOld", &s.isSignalOld); ts->Branch("tKOrigin", &s.tKOrigin); ts->Branch("matchedType", &s.matchedType);
      ts->Branch("muMom", &s.muMom); ts->Branch("tStartZ", &s.tStartZ); ts->Branch("tEndZ", &s.tEndZ); ts->Branch("pRangeShower", &s.pRangeShower);
    }

    void Process(AnaSpillPD* spill){
      entry++;
      const AnaEventInfoPD* info = static_cast<AnaEventInfoPD*>(spill->EventInfo);
      if(!info || !info->IsMC || spill->Bunches.empty())return;
      const int isMC = 1, fd = info->FullDump ? 1 : 0;
      const AnaBeamPD* beam = static_cast<AnaBeamPD*>(spill->Beam);
      const AnaTrueParticleB* beamTrue = (beam && beam->BeamParticle) ? static_cast<AnaTrueParticleB*>(beam->BeamParticle->TrueObject) : nullptr;

      std::vector<AnaParticlePD*> parts, slice;
      for(AnaParticleB* pb : static_cast<AnaBunchPD*>(spill->Bunches[0])->Particles){
        AnaParticlePD* p = static_cast<AnaParticlePD*>(pb);
        parts.push_back(p);
        if(p->BeamOrigin)slice.push_back(p);
      }
      const bool hasSlice = !slice.empty();
      std::map<UInt_t, const AnaParticlePD*> byUID;
      for(const AnaParticlePD* p : parts)byUID[p->UniqueID] = p;
      const std::vector<AnaHitPD>& unassoc = static_cast<AnaBunchPD*>(spill->Bunches[0])->UnassociatedHits;
      nEvents[fd]++; if(hasSlice)nEventsSlice[fd]++;
      TVector3 vtx; bool hasVtx = false;
      if(hasSlice && Valid(slice[0]->PositionEnd)){ vtx = V(slice[0]->PositionEnd); hasVtx = true; }
      const std::vector<VcPoint> cloud = hasVtx ? VcCloud(parts) : std::vector<VcPoint>();
      auto distSlice = [&slice](const TVector3& x){
        double d = 1e9;
        for(const AnaParticlePD* q : slice){
          if(Valid(q->PositionStart) && Valid(q->PositionEnd))d = std::min(d, DistToSegment(x, V(q->PositionStart), V(q->PositionEnd)));
          else if(Valid(q->PositionStart))d = std::min(d, (x - V(q->PositionStart)).Mag());
        }
        return d < 1e8 ? d : -1.;
      };

      std::map<int, AnaTrueParticleB*> byID;
      for(AnaTrueParticleB* t : spill->TrueParticles)byID[t->ID] = t;
      auto isSignal = [&](const AnaTrueParticleB* t, double& kMomEnd){
        kMomEnd = -1;
        if(!t || t->PDG != -13 || !beamTrue)return false;
        auto ip = byID.find(t->ParentID);
        if(ip == byID.end() || ip->second->PDG != 321 || ip->second->ProcessEnd != AnaTrueParticleB::Decay)return false;
        if(!DescendsFrom(ip->second->ID, beamTrue->ID, byID))return false;
        kMomEnd = 1000.*ip->second->MomentumEnd;
        return true;
      };

      // candidatos
      for(const AnaParticlePD* p : parts){
        // trazas, o PFParticles sin traza con traza forzada (candidatos tipo shower)
        const bool isTrk = p->Type == AnaParticlePD::kTrack;
        const AnaForcedTrackPD* ft = (!isTrk && !p->ForcedTrack.empty()) ? &p->ForcedTrack[0] : nullptr;
        if(!isTrk && !ft)continue;
        const double L = isTrk ? p->Length : ft->Length;
        const Float_t* pa = isTrk ? p->PositionStart : ft->PositionStart;
        const Float_t* pb = isTrk ? p->PositionEnd : ft->PositionEnd;
        if(L <= 0 || !Valid(pa) || !Valid(pb))continue;
        const double pr = MuRangeP(L);
        if(pr < kPMin || pr > kPMax)continue;
        const TVector3 a = V(pa), b = V(pb);
        c = Cand();
        c.isShower = !isTrk;
        c.entry = entry-1; c.fullDump = fd; c.hasSlice = hasSlice; c.inSlice = p->BeamOrigin; c.complete = HasHits(p);
        c.isPrimary = p->isPandora; c.nDaughters = p->DaughtersIDs.size();
        c.length = L; c.pRange = pr;
        c.chi2p = isTrk ? p->Chi2Proton : ft->Chi2Proton; c.chi2mu = isTrk ? p->Chi2Muon : ft->Chi2Muon; c.ndf = isTrk ? p->Chi2ndf : ft->Chi2ndf;
        c.cnnTrk = p->CNNscore[0]; c.cnnEM = p->CNNscore[1]; c.cnnMichel = p->CNNscore[2];
        c.facesA = DistToFaces(a); c.facesB = DistToFaces(b);
        c.axisA = DistToAxis(a, isMC); c.axisB = DistToAxis(b, isMC);
        c.sliceA = p->BeamOrigin ? 0 : distSlice(a); c.sliceB = p->BeamOrigin ? 0 : distSlice(b);
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
        c.clearCosmic = p->IsClearCosmic; c.hasT0 = p->T0 > -900; c.cosmicCriteria = p->CosmicCriteria; c.beamCosmicScore = p->BeamCosmicScore;
        const AnaTrueParticlePD* t = static_cast<const AnaTrueParticlePD*>(p->TrueObject);
        c.tPDG = t ? t->PDG : 0;
        { auto ip = t ? byID.find(t->ParentID) : byID.end(); c.tParentPDG = ip != byID.end() ? ip->second->PDG : 0; }   // ParentPDG no lo rellena el módulo
        c.tFromBeamK = t && beamTrue && beamTrue->PDG == 321 && t->PDG == -13 && t->ParentID == beamTrue->ID; c.tProcStart = t ? int(t->ProcessStart) : -1;
        c.tOrigin = t ? t->Origin : -1; c.tGen = t ? t->Generation : -999;
        c.tFromBeam = (t && beamTrue) ? DescendsFrom(t->ID, beamTrue->ID, byID) : 0;
        double km; c.tSignal = isSignal(t, km); c.tKMomEnd = km;
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
        { VtxConn vc; if(hasVtx)vc = VertexConnection(vtx, start, p->UniqueID, cloud);
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
        { ShowerLen sl; if(!isTrk)sl = ShowerLength(p);
          c.shNSP = sl.nSP; c.shLin = sl.lin; c.shLenExt = sl.lenExt; c.shLenCore = sl.lenCore; c.shLenBragg = sl.lenBragg; c.shBraggRatio = sl.braggRatio;
          c.shPRangeCore = sl.lenCore > 0 ? MuRangeP(sl.lenCore) : -1; c.shPRangeBragg = sl.lenBragg > 0 ? MuRangeP(sl.lenBragg) : -1; }
        tc->Fill();
      }

      // señal true: para cada mu+ de K+ de la cadena del haz, el PFParticle emparejado con más hits
      for(AnaTrueParticleB* t : spill->TrueParticles){
        double km;
        const bool old = isSignal(t, km);
        // añadido el 20-nov: todos los mu+ de desintegración no cósmicos con p true de K en reposo
        const bool atRest = t->PDG == -13 && t->ProcessStart == AnaTrueParticleB::Decay && static_cast<const AnaTrueParticlePD*>(t)->Origin != 2 &&
                            t->Momentum > 0.230 && t->Momentum < 0.241;   // GeV/c
        if(!old && !atRest)continue;
        if(!old){ auto ip = byID.find(t->ParentID); if(ip != byID.end() && ip->second->PDG == 321)km = 1000.*ip->second->MomentumEnd; }
        s = Sig();
        s.isSignalOld = old; s.tKOrigin = KOrigin(t, beamTrue, byID); s.muMom = 1000.*t->Momentum;
        s.tStartZ = t->Position[2]; s.tEndZ = t->PositionEnd[2];
        s.entry = entry-1; s.fullDump = fd; s.hasSlice = hasSlice; s.kMomEnd = km; s.trueLength = (V(t->PositionEnd) - V(t->Position)).Mag();   // distancia inicio-fin
        s.axisMin = std::min(DistToAxis(V(t->Position), isMC), DistToAxis(V(t->PositionEnd), isMC));
        s.contained = DistToFaces(V(t->Position)) > 0 && DistToFaces(V(t->PositionEnd)) > 0;
        const AnaParticlePD* best = nullptr;
        for(const AnaParticlePD* p : parts)if(p->TrueObject && p->TrueObject->ID == t->ID && (!best || p->NHits > best->NHits))best = p;
        // ¿hay algo reconstruido donde se para el mu (el Michel)? Extremo más cercano de otro PFParticle
        s.nearestToEnd = -1;
        { const TVector3 e = V(t->PositionEnd); double d = 1e9;
          for(const AnaParticlePD* p : parts){
            if(p == best)continue;
            if(Valid(p->PositionStart))d = std::min(d, (V(p->PositionStart) - e).Mag());
            if(Valid(p->PositionEnd))d = std::min(d, (V(p->PositionEnd) - e).Mag());
          }
          if(d < 1e8)s.nearestToEnd = d; }
        s.matched = best != nullptr;
        s.matchedType = best ? int(best->Type) : -1;
        s.pRangeShower = (best && best->Type != AnaParticlePD::kTrack && !best->ForcedTrack.empty()) ? MuRangeP(best->ForcedTrack[0].Length) : -1;
        if(best){
          s.isTrack = best->Type == AnaParticlePD::kTrack; s.inSlice = best->BeamOrigin; s.complete = HasHits(best);
          s.pRange = MuRangeP(best->Length);
          s.inWindow = s.isTrack && s.pRange >= kWinMin && s.pRange <= kWinMax;
        }
        else s.pRange = -1;
        ts->Fill();
      }
    }

    // Categoría true de un candidato, para el resumen
    static std::string Category(const Cand& x){
      if(x.tPDG == 0)return "sin truth";
      if(x.tSignal)return x.tKMomEnd < 50 ? "SEÑAL mu+ de K+ (K en reposo)" : "SEÑAL mu+ de K+ (K en vuelo)";
      if(x.tFromBeamK)return "mu+ de la K+ del haz";
      const char* org = x.tOrigin == 2 ? "cósmico" : (x.tFromBeam ? "haz" : "otro");
      switch(std::abs(x.tPDG)){
        case 13: return Form("mu%s (%s%s)", x.tPDG < 0 ? "+" : "-", org, x.tOrigin == 2 ? "" : Form(", madre %d", x.tParentPDG));
        case 211: return Form("pi%s (%s)", x.tPDG > 0 ? "+" : "-", org);
        case 2212: return Form("p (%s)", org);
        case 321: return Form("K%s (%s)", x.tPDG > 0 ? "+" : "-", org);
        case 11: case 22: return Form("EM: e o gamma (%s)", org);
        default: return Form("PDG %d (%s)", x.tPDG, org);
      }
    }
    // Grupo grueso para comparar variables: señal en reposo, mu cósmico, p, pi, EM, resto
    static int Group(const Cand& x){
      if(x.tSignal)return x.tKMomEnd < 50 ? 0 : 6;
      if(x.tFromBeamK)return 7;
      if(std::abs(x.tPDG) == 13 && x.tOrigin == 2)return 1;
      if(std::abs(x.tPDG) == 2212)return 2;
      if(std::abs(x.tPDG) == 211)return 3;
      if(std::abs(x.tPDG) == 11 || x.tPDG == 22)return 4;
      return 5;
    }
    static constexpr int kNGroups = 8;
    static const char* GroupName(int g){
      static const char* n[kNGroups] = {"señal (K en reposo)", "mu cósmico", "protón", "pi+-", "EM (e, gamma)", "resto", "señal (K en vuelo)", "mu+ de K+ del haz"};
      return n[g];
    }
    // Pasa los cortes del módulo (ventana, se para, cerca del eje)
    static bool ModuleCuts(const Cand& x){
      return x.pRange >= kWinMin && x.pRange <= kWinMax && std::max(x.facesA, x.facesB) >= 10 && std::min(x.axisA, x.axisB) <= 150;
    }

    void SetAddresses(TTree* t, Cand& x) const {
      for(auto* br : *t->GetListOfBranches()){
        const std::string n = br->GetName();
        void* a = nullptr;
#define MUCAND_BR(v) if(n == #v)a = &x.v;
        MUCAND_BR(entry) MUCAND_BR(fullDump) MUCAND_BR(hasSlice) MUCAND_BR(inSlice) MUCAND_BR(complete) MUCAND_BR(isPrimary) MUCAND_BR(nDaughters)
        MUCAND_BR(length) MUCAND_BR(pRange) MUCAND_BR(chi2p) MUCAND_BR(chi2mu) MUCAND_BR(ndf) MUCAND_BR(cnnTrk) MUCAND_BR(cnnEM) MUCAND_BR(cnnMichel)
        MUCAND_BR(facesA) MUCAND_BR(facesB) MUCAND_BR(axisA) MUCAND_BR(axisB) MUCAND_BR(sliceA) MUCAND_BR(sliceB) MUCAND_BR(vtxA) MUCAND_BR(vtxB)
        MUCAND_BR(dEdxA) MUCAND_BR(dEdxB) MUCAND_BR(nHitsA) MUCAND_BR(nHitsB) MUCAND_BR(nHits2) MUCAND_BR(actPartsA) MUCAND_BR(actPartsB)
        MUCAND_BR(actHitsA) MUCAND_BR(actHitsB) MUCAND_BR(actMichelA) MUCAND_BR(actMichelB) MUCAND_BR(clearCosmic) MUCAND_BR(hasT0)
        MUCAND_BR(cosmicCriteria) MUCAND_BR(beamCosmicScore) MUCAND_BR(tPDG) MUCAND_BR(tParentPDG) MUCAND_BR(tProcStart) MUCAND_BR(tOrigin)
        MUCAND_BR(chi2muF) MUCAND_BR(chi2muR) MUCAND_BR(chi2pF) MUCAND_BR(chi2pR) MUCAND_BR(michelOwnA) MUCAND_BR(michelOwnB)
        MUCAND_BR(intPeak) MUCAND_BR(intPeakRR) MUCAND_BR(dirYStop) MUCAND_BR(stopAtB)
        MUCAND_BR(tFromBeam) MUCAND_BR(tSignal) MUCAND_BR(tGen) MUCAND_BR(tFromBeamK) MUCAND_BR(tKMomEnd) MUCAND_BR(tPur) MUCAND_BR(tEff) MUCAND_BR(tLength)
        MUCAND_BR(isShower)
#undef MUCAND_BR
        if(a)t->SetBranchAddress(n.c_str(), a);
      }
    }

    static double Q(std::vector<double> v, double f){
      if(v.empty())return -999;
      std::sort(v.begin(), v.end());
      return v[std::min(v.size()-1, size_t(f*v.size()))];
    }

    void Print() const {
      std::cout << std::endl << "Eventos: normales " << nEvents[0] << " (con beam slice " << nEventsSlice[0] << "), con volcado "
                << nEvents[1] << " (con beam slice " << nEventsSlice[1] << ")" << std::endl
                << "(fondo: eventos normales; señal: todos, casi siempre con volcado, que se hace en los eventos con un K+ secundario)" << std::endl;
      Cand x; SetAddresses(tc, x);
      std::map<std::string, long> cat[2];
      long total = 0;
      // variables por grupo (candidatos con los cortes del módulo; fondo solo de eventos normales)
      enum { vMaxdEdx, vMindEdx, vChi2mu, vChi2p, vMichelOwn, vMichelStop, vActStart, vIntPeak, vDirY, vAxis, vVtx, vSlice, vClear, vT0, vInSlice, vReversed, kNV };
      static const char* vname[kNV] = {"dE/dx extremo de parada (MeV/cm)", "dE/dx otro extremo", "chi2 mu / ndf (parada por dE/dx)", "chi2 p / ndf (misma orient.)",
                                       "CNN Michel hits propios (parada)", "CNN Michel activ. extremo parada", "partículas en el otro extremo",
                                       "pico de dE/dx interior", "dir. y hacia la parada", "dist. mín. al eje (cm)",
                                       "dist. mín. al vértice del haz (cm)", "dist. mín. a la beam slice (cm)", "IsClearCosmic", "tiene T0", "en la beam slice",
                                       "Pandora al revés"};
      std::vector<double> var[kNGroups][kNV];
      long nGroup[kNGroups] = {0};
      for(Long64_t i = 0; i < tc->GetEntries(); i++){
        tc->GetEntry(i);
        if(x.isShower || !ModuleCuts(x))continue;   // el resumen, solo con trazas (como antes del 20-nov)
        const int g = Group(x);
        const bool isSig = g == 0 || g == 6;
        if(!x.fullDump){ cat[x.inSlice][Category(x)]++; total++; }
        if(!isSig && x.fullDump)continue;   // el fondo, solo de eventos normales
        nGroup[g]++;
        // extremo de parada: el de mayor dE/dx (ver Process)
        const bool stopB = x.stopAtB;
        const double dStop = stopB ? x.dEdxB : x.dEdxA, dOther = stopB ? x.dEdxA : x.dEdxB;
        if(dStop > 0)var[g][vMaxdEdx].push_back(dStop);
        if(dOther > 0)var[g][vMindEdx].push_back(dOther);
        const double c2mu = stopB ? x.chi2muF : x.chi2muR, c2p = stopB ? x.chi2pF : x.chi2pR;
        if(c2mu >= 0)var[g][vChi2mu].push_back(c2mu);
        if(c2p >= 0)var[g][vChi2p].push_back(c2p);
        const double mo = stopB ? x.michelOwnB : x.michelOwnA;
        if(mo >= 0)var[g][vMichelOwn].push_back(mo);
        if(x.intPeak >= 0)var[g][vIntPeak].push_back(x.intPeak);
        var[g][vDirY].push_back(x.dirYStop);
        if(x.chi2muF >= 0)var[g][vReversed].push_back(!x.stopAtB);
        const double mStop = stopB ? x.actMichelB : x.actMichelA;
        var[g][vMichelStop].push_back(mStop);   // -1: sin actividad
        var[g][vActStart].push_back(stopB ? x.actPartsA : x.actPartsB);
        var[g][vAxis].push_back(std::min(x.axisA, x.axisB));
        if(x.vtxA >= 0)var[g][vVtx].push_back(std::min(x.vtxA, x.vtxB));
        if(!x.inSlice && x.sliceA >= 0)var[g][vSlice].push_back(std::min(x.sliceA, x.sliceB));
        var[g][vClear].push_back(x.clearCosmic); var[g][vT0].push_back(x.hasT0); var[g][vInSlice].push_back(x.inSlice);
      }
      std::cout << std::endl << "Candidatos con los cortes del módulo (p de rango en [180, 300] MeV/c, un extremo a > 10 cm de las caras y uno a < 150 cm del eje), eventos normales: "
                << total << " (" << Form("%.2f", nEvents[0] ? double(total)/nEvents[0] : 0.) << " por evento)" << std::endl;
      for(int in = 1; in >= 0; in--){
        long n = 0; for(auto& kv : cat[in])n += kv.second;
        std::cout << (in ? "  en la beam slice: " : "  fuera de la beam slice: ") << n << std::endl;
        std::vector<std::pair<long, std::string>> v;
        for(auto& kv : cat[in])v.push_back({kv.second, kv.first});
        std::sort(v.rbegin(), v.rend());
        for(size_t k = 0; k < v.size() && k < 15; k++)
          std::cout << Form("    %-44s %8ld (%5.1f %%)", v[k].second.c_str(), v[k].first, n ? 100.*v[k].first/n : 0.) << std::endl;
      }
      std::cout << std::endl << "Variables por grupo (mediana [16 %, 84 %]; las 0/1, fracción), candidatos con los cortes del módulo:" << std::endl;
      std::cout << Form("  %-36s", "");
      for(int g = 0; g < kNGroups; g++)if(nGroup[g])std::cout << Form("%22s", GroupName(g));
      std::cout << std::endl << Form("  %-36s", "candidatos");
      for(int g = 0; g < kNGroups; g++)if(nGroup[g])std::cout << Form("%22ld", nGroup[g]);
      std::cout << std::endl;
      for(int v = 0; v < kNV; v++){
        std::cout << Form("  %-36s", vname[v]);
        const bool frac = v == vClear || v == vT0 || v == vInSlice || v == vReversed;
        for(int g = 0; g < kNGroups; g++){
          if(!nGroup[g])continue;
          const std::vector<double>& d = var[g][v];
          if(d.empty()){ std::cout << Form("%22s", "-"); continue; }
          if(frac){ double m = 0; for(double y : d)m += y; std::cout << Form("%22.2f", m/d.size()); }
          else std::cout << Form("%10.3g [%4.3g,%4.3g]", Q(d,.5), Q(d,.16), Q(d,.84));
        }
        std::cout << std::endl;
      }

      // señal true
      Sig y; TTree* u = ts;
      u->SetBranchAddress("fullDump", &y.fullDump); u->SetBranchAddress("hasSlice", &y.hasSlice); u->SetBranchAddress("matched", &y.matched);
      u->SetBranchAddress("isTrack", &y.isTrack); u->SetBranchAddress("inSlice", &y.inSlice); u->SetBranchAddress("complete", &y.complete);
      u->SetBranchAddress("inWindow", &y.inWindow); u->SetBranchAddress("contained", &y.contained); u->SetBranchAddress("kMomEnd", &y.kMomEnd);
      u->SetBranchAddress("nearestToEnd", &y.nearestToEnd); u->SetBranchAddress("isSignalOld", &y.isSignalOld);
      long n = 0, rest = 0, contained = 0, matched = 0, track = 0, window = 0, windowRest = 0;
      std::vector<double> nearest;
      for(Long64_t i = 0; i < u->GetEntries(); i++){
        u->GetEntry(i);
        if(!y.isSignalOld)continue;   // el resumen, con la señal de antes del 20-nov
        n++; if(y.kMomEnd < 50)rest++;
        if(y.contained)contained++;
        if(y.matched)matched++;
        if(y.matched && y.isTrack)track++;
        if(y.inWindow){ window++; if(y.kMomEnd < 50)windowRest++; }
        if(y.contained && y.nearestToEnd >= 0)nearest.push_back(y.nearestToEnd);
      }
      std::cout << std::endl << "mu+ true de K+ de la cadena del haz (todos los eventos): " << n << std::endl;
      std::cout << Form("  K en reposo (p final < 50 MeV/c)       %ld (%.1f %%)", rest, n ? 100.*rest/n : 0.) << std::endl;
      std::cout << Form("  mu contenido en el volumen activo      %ld (%.1f %%)", contained, n ? 100.*contained/n : 0.) << std::endl;
      std::cout << Form("  con PFParticle emparejado              %ld (%.1f %%)", matched, n ? 100.*matched/n : 0.) << std::endl;
      std::cout << Form("    traza                                %ld (%.1f %%)", track, n ? 100.*track/n : 0.) << std::endl;
      std::cout << Form("    traza con p de rango en [180, 300]   %ld (%.1f %%), de K en reposo %ld", window, n ? 100.*window/n : 0., windowRest) << std::endl;
      std::cout << Form("  mu contenidos: distancia del final true al extremo de otro PFParticle (cm), cuantiles 16/50/84: %.1f / %.1f / %.1f; < 10 cm: %.0f %%",
                        Q(nearest,.16), Q(nearest,.5), Q(nearest,.84),
                        nearest.empty() ? 0. : 100.*std::count_if(nearest.begin(), nearest.end(), [](double d){ return d < 10; })/nearest.size()) << std::endl;
      std::cout << std::endl << "Árboles cand y sig en " << fout->GetName() << std::endl;
    }

    void Write(){ fout->cd(); tc->Write(); ts->Write(); fout->Close(); }
  };
}

int mu_candidates_minitrees_impl(const char* input, const char* output, Long64_t maxEntries){
  using namespace mucand;
  ND::params().SetReadParamOverrideFilePointPassed();   // ComputeRangeMomentum lee parámetros
  const std::vector<std::string> files = InputFiles(input);
  if(files.empty()){ std::cerr << "No hay ficheros en " << input << std::endl; return -1; }
  const std::string treeName = TreeName(TString(files[0]).Contains("*") ? std::string("") : files[0]);
  TChain chain(treeName.empty() ? "highlandana/MiniTree" : treeName.c_str());
  for(const std::string& f : files)chain.Add(f.c_str());
  std::cout << "Ficheros: " << files.size() << std::endl;
  if(chain.LoadTree(0) < 0){ std::cerr << "No hay entradas (¿árbol " << chain.GetName() << "?)" << std::endl; return -1; }
  AnaSpillPD* spill = nullptr;
  chain.SetBranchAddress("Spill", &spill);
  MuCandidates mc;
  mc.Book(output);
  LoopMiniTree(chain, spill, maxEntries, files.size(), [&mc](AnaSpillPD* s){ mc.Process(s); });
  mc.Print();
  mc.Write();
  return 0;
}
