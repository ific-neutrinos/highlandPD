// Implementación de kmu_pairing_minitrees.C (se compila con ACLiC).
//
// secondaryKaon exige que el mu sea la única hija del candidato K en Pandora (candidato: no
// primario, con una sola hija). Aquí se mide cuánto se pierde por eso y si se recupera buscando el
// mu cerca del final de la K aunque no cuelgue de ella. Solo MC y eventos con beam slice.
//
// 1. Relación K-mu. Para cada K+ secundario de la cadena del haz que se desintegra con un mu+ y
//    cuya K está emparejada con una partícula de la beam slice (TrueObject, como kaon_fate): qué es
//    el mu respecto a esa K en la reconstrucción.
//      - hija única de la K (lo que pide la selección)
//      - hija de la K, pero la K tiene más hijas
//      - en la beam slice sin ser hija de la K (madre de la K: la dirección está invertida; otra)
//      - fuera de la beam slice (completo o en resumen)
//      - dentro de la traza de la K (sus hits en la K: K y mu unidos en una traza)
//      - dentro de otra partícula, no reconstruido o sin hits
//    Para las que no son hija: distancia del final de la K al extremo más cercano del mu (con el
//    inicio, como la selección, o con el fin si el mu está al revés), y cuántas pasan d < 7,7 cm y
//    la ventana del momento por rango del mu.
//    Si la K emparejada es el primario de la beam slice, la selección la salta (isPandora): se
//    cuenta aparte.
//
// 2. Selección aproximada. Solo los cortes geométricos y de momento del mu de secondaryKaon:
//    mu traza (Type 2), momento por rango (muón) en [0,221, 0,245] GeV/c y distancia del final de la
//    K al mu < 7,7 cm. Sin los cortes de chi2, ángulo y FV, ni las correcciones (posiciones tal
//    como están en el minitree): sirve para comparar reglas de emparejamiento, no para eficiencias
//    absolutas.
//      "hija"      : K no primaria con una sola hija; mu = esa hija, distancia a su inicio
//      "vecino"    : K no primaria de la beam slice con cualquier número de hijas; mu = cualquier
//                    otra partícula guardada (de la beam slice o de fuera), con su inicio o su fin
//                    cerca del final de la K, sin ser la madre de la K
//      "vecino, solo slice": igual, con el mu solo de la beam slice
//    Señal: K emparejada con un K+ secundario que se desintegra y mu emparejado con su mu+. Se
//    cuentan las K de señal encontradas (cada K+ true una vez) y las parejas de fondo por evento.
//    Además, "hija, todas": la regla actual sin restringir la K a la beam slice. Con los minitrees
//    del módulo los PFParticles de fuera también están en el evento, y GetKaonsAction de
//    secondaryKaon no mira BeamOrigin: esta fila dice cuántas K de fuera entrarían.
#include "kaon_fate_minitrees_impl.C"
#include "pdAnalysisUtils.hxx"
#include "Parameters.hxx"
#include "TVector3.h"
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <vector>

namespace kmu {
  using namespace kfate;

  const double kDistCut = 7.7;                  // cm, KaonMuonDistanceCut
  const double kMomLow = 0.221, kMomHigh = 0.245;   // GeV/c, MuonRangeMomCut

  bool Valid(const Float_t* x){ return x[0] > -900 && x[1] > -900 && x[2] > -900; }
  TVector3 V(const Float_t* x){ return TVector3(x[0], x[1], x[2]); }

  double RangeMom(const AnaParticlePD* p){
    return p->Length > 0 ? pdAnaUtils::ComputeRangeMomentum(p->Length, 13) : -1;
  }
  bool MuMomOK(const AnaParticlePD* p){
    const double m = RangeMom(p);
    return m > kMomLow && m < kMomHigh;
  }

  // Distancia del final de la K al mu: con su inicio (como la selección) y con el extremo más cercano
  double DistStart(const AnaParticlePD* k, const AnaParticlePD* mu){
    if(!Valid(k->PositionEnd) || !Valid(mu->PositionStart))return 1e9;
    return (V(k->PositionEnd) - V(mu->PositionStart)).Mag();
  }
  double DistClosest(const AnaParticlePD* k, const AnaParticlePD* mu, bool* reversed = nullptr){
    const double ds = DistStart(k, mu);
    const double de = (Valid(k->PositionEnd) && Valid(mu->PositionEnd)) ? (V(k->PositionEnd) - V(mu->PositionEnd)).Mag() : 1e9;
    if(reversed)*reversed = de < ds;
    return std::min(ds, de);
  }

  enum Rel {
    kOnlyDaughter = 0, kDaughterMany, kParentOfK, kSliceOther, kRelOutComplete, kRelOutSummary,
    kInK, kInOther, kNotReco, kNRel
  };
  const char* RelName(int r){
    static const char* names[kNRel] = {
      "hija única de la K (la selección la ve)", "hija de la K, que tiene más hijas",
      "en la beam slice: madre de la K (invertida)", "en la beam slice: otra relación",
      "fuera de la beam slice: completo", "fuera de la beam slice: en resumen",
      "dentro de la traza de la K (hits)", "dentro de otra partícula (hits)",
      "no reconstruido o sin hits" };
    return names[r];
  }

  struct RelCounts {
    long n[kNRel] = {0};
    long nNear[kNRel] = {0}, nNearMom[kNRel] = {0}, nReversed[kNRel] = {0}, nWithDist[kNRel] = {0};
    std::vector<double> dist[kNRel];
  };

  // Reglas de la selección aproximada
  enum Rule { kHija = 0, kVecino, kVecinoSlice, kHijaTodas, kNRules };
  const char* RuleName(int r){
    static const char* names[kNRules] = { "hija", "vecino", "vecino, solo slice", "hija, todas (K también de fuera)" };
    return names[r];
  }

  struct KMuPairing {
    long nEvents = 0, nEventsSlice = 0, nFullDump = 0;
    long nPairs = 0, nKPrimary = 0;
    RelCounts rel;
    // selección aproximada
    long nSignalK = 0;                         // K+ true con mu+ y K emparejada en la beam slice (no primaria)
    long sigFound[kNRules] = {0};              // de esas, encontradas con la regla
    long bkgPairs[kNRules] = {0}, bkgEvents[kNRules] = {0}, outsideK[kNRules] = {0};

    KMuPairing(){ ND::params().SetReadParamOverrideFilePointPassed(); }

    void Process(AnaSpillPD* spill){
      const AnaEventInfoPD* info = static_cast<AnaEventInfoPD*>(spill->EventInfo);
      if(!info || !info->IsMC || spill->Bunches.empty())return;
      nEvents++;
      const AnaBeamPD* beam = static_cast<AnaBeamPD*>(spill->Beam);
      AnaTrueParticleB* beamTrue = (beam && beam->BeamParticle) ? static_cast<AnaTrueParticleB*>(beam->BeamParticle->TrueObject) : nullptr;

      std::vector<AnaParticlePD*> parts;
      bool hasSlice = false;
      for(AnaParticleB* pb : static_cast<AnaBunchPD*>(spill->Bunches[0])->Particles){
        parts.push_back(static_cast<AnaParticlePD*>(pb));
        if(parts.back()->BeamOrigin)hasSlice = true;
      }
      if(!hasSlice)return;
      nEventsSlice++;
      if(info->FullDump)nFullDump++;

      std::map<int, AnaParticlePD*> byUID;
      for(AnaParticlePD* p : parts)byUID[p->UniqueID] = p;
      std::map<int, AnaTrueParticleB*> byID;
      for(AnaTrueParticleB* t : spill->TrueParticles)byID[t->ID] = t;

      // pares de señal: true K ID -> true mu ID
      std::map<int, int> muOfK;
      if(beamTrue){
        for(AnaTrueParticleB* tb : spill->TrueParticles){
          if(tb->PDG != 321 || tb->ID == beamTrue->ID || tb->ProcessEnd != AnaTrueParticleB::Decay)continue;
          if(!DescendsFrom(tb->ID, beamTrue->ID, byID))continue;
          for(int d : static_cast<AnaTrueParticle*>(tb)->Daughters){
            auto it = byID.find(d);
            if(it != byID.end() && it->second->PDG == -13){ muOfK[tb->ID] = d; break; }
          }
        }
      }

      // 1. relación K-mu
      std::set<int> signalK;   // true K con K reco elegible (en la beam slice, no primaria)
      for(auto& km : muOfK){
        const AnaTrueParticlePD* k  = static_cast<AnaTrueParticlePD*>(byID[km.first]);
        const AnaTrueParticlePD* mu = static_cast<AnaTrueParticlePD*>(byID[km.second]);
        const Result rk = Classify(k, parts, true);
        if(!rk.part || !rk.part->BeamOrigin || rk.fate == kSliceInOther)continue;
        nPairs++;
        if(rk.part->isPandora){ nKPrimary++; continue; }
        signalK.insert(k->ID);
        const AnaParticlePD* kr = rk.part;
        const Result rm = Classify(mu, parts, true);
        int r = kNotReco;
        const AnaParticlePD* mr = rm.part;
        if(rm.fate <= kSliceDeeper || rm.fate == kfate::kOutComplete || rm.fate == kfate::kOutSummary){
          if(mr->ParentID == kr->UniqueID)r = kr->DaughtersIDs.size() == 1 ? kOnlyDaughter : kDaughterMany;
          else if(mr->BeamOrigin)r = (kr->ParentID == mr->UniqueID) ? kParentOfK : kSliceOther;
          else r = (rm.fate == kfate::kOutComplete) ? kRelOutComplete : kRelOutSummary;
        }
        else if(rm.fate == kSliceInOther || rm.fate == kOutInOther)r = (mr == kr) ? kInK : kInOther;
        rel.n[r]++;
        if(r != kInK && r != kInOther && r != kNotReco && mr){
          bool rev = false;
          const double d = DistClosest(kr, mr, &rev);
          if(d < 1e8){
            rel.nWithDist[r]++;
            rel.dist[r].push_back(d);
            if(rev)rel.nReversed[r]++;
            if(d < kDistCut){ rel.nNear[r]++; if(MuMomOK(mr) && mr->Type == AnaParticlePD::kTrack)rel.nNearMom[r]++; }
          }
        }
      }
      nSignalK += signalK.size();

      // 2. selección aproximada
      for(int rule = 0; rule < kNRules; rule++){
        std::set<int> found;
        long bkg = 0;
        for(AnaParticlePD* k : parts){
          if(k->isPandora || k->ParentID == -1)continue;   // como GetKaonsAction (primarios fuera)
          if(rule != kHijaTodas && !k->BeamOrigin)continue;
          if(!Valid(k->PositionEnd))continue;
          std::vector<AnaParticlePD*> mus;
          if(rule == kHija || rule == kHijaTodas){
            if(k->DaughtersIDs.size() != 1)continue;
            auto it = byUID.find(k->DaughtersIDs[0]);
            if(it == byUID.end())continue;
            AnaParticlePD* mu = it->second;
            if(mu->Type == AnaParticlePD::kTrack && MuMomOK(mu) && DistStart(k, mu) < kDistCut)mus.push_back(mu);
          }
          else{
            for(AnaParticlePD* mu : parts){
              if(mu == k || mu->UniqueID == k->ParentID)continue;
              if(rule == kVecinoSlice && !mu->BeamOrigin)continue;
              if(mu->Type == AnaParticlePD::kTrack && MuMomOK(mu) && DistClosest(k, mu) < kDistCut)mus.push_back(mu);
            }
          }
          for(AnaParticlePD* mu : mus){
            const int kt = k->TrueObject ? k->TrueObject->ID : -999999;
            const int mt = mu->TrueObject ? mu->TrueObject->ID : -999999;
            auto it = muOfK.find(kt);
            if(it != muOfK.end() && it->second == mt && signalK.count(kt))found.insert(kt);
            else{ bkg++; if(!k->BeamOrigin)outsideK[rule]++; }
          }
        }
        sigFound[rule] += found.size();
        bkgPairs[rule] += bkg;
        if(bkg)bkgEvents[rule]++;
      }
    }

    static double Median(std::vector<double> v){
      if(v.empty())return -1;
      std::sort(v.begin(), v.end());
      return v[v.size()/2];
    }

    void Print() const {
      std::cout << std::endl << "Eventos de MC: " << nEvents << ", con beam slice: " << nEventsSlice
                << " (con volcado completo: " << nFullDump << ")" << std::endl << std::endl;
      std::cout << "1. K+ secundarios con mu+ y la K emparejada en la beam slice: " << nPairs << std::endl;
      std::cout << "   la K es el primario de la beam slice (la selección la salta): " << Frac(nKPrimary, nPairs) << std::endl;
      const long n = nPairs - nKPrimary;
      std::cout << "   resto: " << n << ". Dónde está el mu respecto a la K:" << std::endl;
      for(int r = 0; r < kNRel; r++){
        std::string label = std::string("     ") + RelName(r);
        int len = 0;
        for(unsigned char c : label)if((c & 0xC0) != 0x80)len++;
        std::cout << label << std::string(len < 52 ? 52-len : 1, ' ') << Frac(rel.n[r], n);
        if(rel.nWithDist[r])
          std::cout << Form("   d mediana %.1f cm; d < %.1f: %ld (con p y traza: %ld); al revés: %ld",
                            Median(rel.dist[r]), kDistCut, rel.nNear[r], rel.nNearMom[r], rel.nReversed[r]);
        std::cout << std::endl;
      }
      std::cout << std::endl << "2. Selección aproximada (mu traza, p rango del mu en [" << kMomLow << ", " << kMomHigh
                << "] GeV/c, d < " << kDistCut << " cm; sin chi2, ángulo, FV ni correcciones)" << std::endl;
      std::cout << "   K de señal posibles (K no primaria en la beam slice, con mu+): " << nSignalK << std::endl;
      std::cout << Form("   %-34s %22s %16s %16s %14s", "regla", "señal encontrada", "fondo/evento", "eventos c/fondo", "fondo K fuera") << std::endl;
      for(int r = 0; r < kNRules; r++)
        std::cout << Form("   %-34s %22s %16.3f %16s %14ld", RuleName(r), Frac(sigFound[r], nSignalK).c_str(),
                          nEventsSlice ? double(bkgPairs[r])/nEventsSlice : 0., Frac(bkgEvents[r], nEventsSlice).c_str(), outsideK[r]) << std::endl;
      std::cout << "   (los eventos con volcado completo tienen un K+ y no son representativos del fondo: usar which = 1 para el fondo)" << std::endl << std::endl;
    }
  };
}

int kmu_pairing_minitrees_impl(const char* input, Long64_t maxEntries, int which = 0){
  using namespace kmu;
  const std::vector<std::string> files = InputFiles(input);
  if(files.empty()){ std::cerr << "No hay ficheros en " << input << std::endl; return -1; }
  const std::string treeName = TreeName(TString(files[0]).Contains("*") ? std::string("") : files[0]);
  TChain chain(treeName.empty() ? "highlandana/MiniTree" : treeName.c_str());
  for(const std::string& f : files)chain.Add(f.c_str());
  if(chain.LoadTree(0) < 0){ std::cerr << "No hay entradas (¿árbol " << chain.GetName() << "?)" << std::endl; return -1; }

  AnaSpillPD* spill = nullptr;
  chain.SetBranchAddress("Spill", &spill);
  KMuPairing km;
  LoopMiniTree(chain, spill, maxEntries, files.size(), [&km, which](AnaSpillPD* s){
    const AnaEventInfoPD* info = static_cast<AnaEventInfoPD*>(s->EventInfo);
    if(info && ((which == 1 && info->FullDump) || (which == 2 && !info->FullDump)))return;
    km.Process(s);
  });
  std::cout << std::endl << "Ficheros: " << files.size() << (which == 1 ? " (solo eventos normales)" : which == 2 ? " (solo volcado completo)" : "") << std::endl;
  km.Print();
  return 0;
}
