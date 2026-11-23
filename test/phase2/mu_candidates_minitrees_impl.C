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
#include "muCandidateUtils.hxx"
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
  // Variables de los candidatos: muCandidateUtils.hxx (secondaryKaonAnalysis), compartidas con el análisis
  // secondaryKaonProductionXS (23-nov).
  using namespace muCandUtils;

  // Variables del árbol de candidatos (una entrada por candidato): las de muCandUtils::Candidate y la entrada
  struct Cand : public Candidate {
    Int_t entry, fullDump;
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

      std::vector<AnaParticlePD*> parts;
      for(AnaParticleB* pb : static_cast<AnaBunchPD*>(spill->Bunches[0])->Particles)parts.push_back(static_cast<AnaParticlePD*>(pb));
      EventContext ctx;
      MakeContext(ctx, parts, static_cast<AnaBunchPD*>(spill->Bunches[0])->UnassociatedHits, spill->TrueParticles, beamTrue, isMC);
      const bool hasSlice = ctx.hasSlice;
      const auto& byID = ctx.byID;
      nEvents[fd]++; if(hasSlice)nEventsSlice[fd]++;
      auto isSignal = [&ctx](const AnaTrueParticleB* t, double& kMomEnd){ return IsOldSignal(t, ctx, kMomEnd); };

      // candidatos
      for(const AnaParticlePD* p : parts){
        c = Cand();
        if(!Compute(p, ctx, c))continue;
        c.entry = entry-1; c.fullDump = fd;
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
