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
#include "kaon_fate_minitrees_impl.C"
#include "pdAnalysisUtils.hxx"
#include "Parameters.hxx"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
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
  };

  // mu+ true de K+ de la cadena del haz: qué le pasa en la reco
  struct Sig {
    Int_t entry, fullDump, hasSlice, matched, isTrack, inSlice, complete, inWindow, contained; Float_t kMomEnd, pRange, trueLength, axisMin, nearestToEnd;
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
      ts = new TTree("sig", "mu+ true de K+ de la cadena del haz");
      ts->Branch("entry", &s.entry); ts->Branch("fullDump", &s.fullDump); ts->Branch("hasSlice", &s.hasSlice); ts->Branch("matched", &s.matched);
      ts->Branch("isTrack", &s.isTrack); ts->Branch("inSlice", &s.inSlice); ts->Branch("complete", &s.complete); ts->Branch("inWindow", &s.inWindow);
      ts->Branch("contained", &s.contained); ts->Branch("kMomEnd", &s.kMomEnd); ts->Branch("pRange", &s.pRange); ts->Branch("trueLength", &s.trueLength); ts->Branch("axisMin", &s.axisMin); ts->Branch("nearestToEnd", &s.nearestToEnd);
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
      nEvents[fd]++; if(hasSlice)nEventsSlice[fd]++;
      TVector3 vtx; bool hasVtx = false;
      if(hasSlice && Valid(slice[0]->PositionEnd)){ vtx = V(slice[0]->PositionEnd); hasVtx = true; }
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
        if(p->Type != AnaParticlePD::kTrack || p->Length <= 0 || !Valid(p->PositionStart) || !Valid(p->PositionEnd))continue;
        const double pr = MuRangeP(p->Length);
        if(pr < kPMin || pr > kPMax)continue;
        const TVector3 a = V(p->PositionStart), b = V(p->PositionEnd);
        c = Cand();
        c.entry = entry-1; c.fullDump = fd; c.hasSlice = hasSlice; c.inSlice = p->BeamOrigin; c.complete = HasHits(p);
        c.isPrimary = p->isPandora; c.nDaughters = p->DaughtersIDs.size();
        c.length = p->Length; c.pRange = pr; c.chi2p = p->Chi2Proton; c.chi2mu = p->Chi2Muon; c.ndf = p->Chi2ndf;
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
        tc->Fill();
      }

      // señal true: para cada mu+ de K+ de la cadena del haz, el PFParticle emparejado con más hits
      for(AnaTrueParticleB* t : spill->TrueParticles){
        double km;
        if(!isSignal(t, km))continue;
        s = Sig();
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
        if(!ModuleCuts(x))continue;
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
      u->SetBranchAddress("nearestToEnd", &y.nearestToEnd);
      long n = 0, rest = 0, contained = 0, matched = 0, track = 0, window = 0, windowRest = 0;
      std::vector<double> nearest;
      for(Long64_t i = 0; i < u->GetEntries(); i++){
        u->GetEntry(i);
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
