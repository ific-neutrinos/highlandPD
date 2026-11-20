// Implementación de kaon_fate_minitrees.C (se compila con ACLiC).
//
// Para cada K+ true de la cadena del haz (la partícula del haz si es un K+ y los K+ que salen de
// ella o de sus descendientes) y para el mu+ de su desintegración: dónde ha acabado en la
// reconstrucción guardada en el minitree.
//
// Cómo se decide, en este orden:
//   1. "evento sin beam slice": el módulo no guarda cósmicos en esos eventos (CosmicsWithoutBeamSlice
//      false), así que no se puede saber más.
//   2. Emparejamiento de la reco (TrueObject = la partícula true que más hits aporta): si alguna
//      partícula de la beam slice está emparejada con ella, "beam slice" (primario, hija del
//      primario o generación >= 2); si no, si lo está alguna de fuera, "fuera de la beam slice"
//      (cósmico completo o en resumen, y si Pandora lo marca como IsClearCosmic).
//   3. Hits: si ninguna está emparejada pero sus hits del plano 2 (TrueID) están en alguna
//      partícula guardada con hits, "dentro de otra partícula" (de la beam slice o de fuera): por
//      ejemplo K+ y mu+ unidos en una sola traza emparejada con el otro. Los cósmicos en resumen
//      no guardan hits: ahí este caso no se ve y sale como "no reconstruida".
//   4. "sin hits en el TPC" si no deja hits (NHitsTrue = 0 en los tres planos); si no,
//      "no reconstruida": deja hits pero no está en ninguna partícula guardada.
#include "pdDataClasses.hxx"
#include "minitree_loop.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace kfate {

  // Ficheros de entrada: un .txt con uno por línea, o un fichero o patrón de ROOT
  std::vector<std::string> InputFiles(const std::string& input){
    std::vector<std::string> files;
    if(TString(input).EndsWith(".txt")){
      std::ifstream in(input);
      std::string line;
      while(std::getline(in, line)){
        line = TString(line).Strip(TString::kBoth).Data();
        if(!line.empty() && line[0] != '#')files.push_back(line);
      }
    }
    else files.push_back(input);
    return files;
  }

  std::string TreeName(const std::string& file){
    if(file.empty())return "";   // patrón: el árbol por defecto (sin TFile::Open(""), que se queja)
    std::unique_ptr<TFile> f(TFile::Open(file.c_str()));
    if(!f || f->IsZombie())return "";
    if(f->Get("highlandana/MiniTree"))return "highlandana/MiniTree";
    if(f->Get("MiniTree"))return "MiniTree";
    return "";
  }

  std::string Frac(long n, long total){
    return Form("%ld (%.1f %%)", n, total > 0 ? 100.*n/total : 0.);
  }

  bool DescendsFrom(int id, int ancestor, const std::map<int, AnaTrueParticleB*>& byID){
    for(int depth = 0; depth < 1000; depth++){
      if(id == ancestor)return true;
      auto it = byID.find(id);
      if(it == byID.end() || it->second->ParentID == 0)return false;
      id = it->second->ParentID;
    }
    return false;
  }

  // Categorías, en el orden en que se imprimen
  enum Fate {
    kSlicePrimary = 0, kSliceDaughter, kSliceDeeper, kSliceInOther,
    kOutComplete, kOutSummary, kOutClearCosmic, kOutInOther,
    kNoSlice, kNotReco, kNoHits, kNFates
  };
  const char* FateName(int f){
    static const char* names[kNFates] = {
      "beam slice: primario", "beam slice: hija del primario", "beam slice: generación >= 2",
      "beam slice: dentro de otra partícula (hits)",
      "fuera: cósmico completo", "fuera: cósmico en resumen", "fuera (completo o en resumen): IsClearCosmic",
      "fuera: dentro de otra partícula (hits)",
      "evento sin beam slice", "no reconstruida (deja hits)", "sin hits en el TPC" };
    return names[f];
  }
  // Grupos para la tabla K+ x mu+
  enum Coarse { kCSlice = 0, kCOut, kCNoSlice, kCNotReco, kCNoHits, kNCoarse };
  const char* CoarseName(int c){
    static const char* names[kNCoarse] = { "beam slice", "fuera", "sin slice", "no reco", "sin hits" };
    return names[c];
  }
  int ToCoarse(int f){
    if(f <= kSliceInOther)return kCSlice;
    if(f <= kOutInOther)return kCOut;
    if(f == kNoSlice)return kCNoSlice;
    if(f == kNotReco)return kCNotReco;
    return kCNoHits;
  }

  struct Result {
    int fate = kNotReco;
    bool clearCosmic = false;
    AnaParticlePD* part = nullptr;   // partícula reconstruida (emparejada o con más hits suyos)
    double hitFrac = -1;             // hits suyos del plano 2 en esa partícula / NHitsTrue[2]
  };

  int NHitsOf(const AnaParticlePD* p, int trueID){
    int n = 0;
    for(const AnaHitPD& h : p->Hits[2])     if(h.TrueID == trueID)n++;
    for(const AnaHitPD& h : p->HitsNoCalo[2])if(h.TrueID == trueID)n++;
    return n;
  }
  bool HasHits(const AnaParticlePD* p){
    for(int pl = 0; pl < 3; pl++)if(!p->Hits[pl].empty() || !p->HitsNoCalo[pl].empty())return true;
    return false;
  }

  Result Classify(const AnaTrueParticlePD* t, const std::vector<AnaParticlePD*>& parts, bool hasSlice){
    Result r;
    const int nTrue = t->NHitsTrue[0] + t->NHitsTrue[1] + t->NHitsTrue[2];
    if(!hasSlice){ r.fate = kNoSlice; return r; }

    // 2. emparejamiento: mejor la de la beam slice; entre varias, la de más hits
    AnaParticlePD* bestSlice = nullptr;
    AnaParticlePD* bestOut = nullptr;
    for(AnaParticlePD* p : parts){
      const AnaTrueObjectC* to = p->TrueObject;
      if(!to || to->ID != t->ID)continue;
      AnaParticlePD*& best = p->BeamOrigin ? bestSlice : bestOut;
      if(!best || p->NHits > best->NHits)best = p;
    }
    auto frac = [&](AnaParticlePD* p){ return t->NHitsTrue[2] > 0 && HasHits(p) ? double(NHitsOf(p, t->ID))/t->NHitsTrue[2] : -1.; };
    if(bestSlice){
      r.part = bestSlice;
      r.fate = bestSlice->isPandora ? kSlicePrimary : (bestSlice->Generation == 1 ? kSliceDaughter : kSliceDeeper);
      r.hitFrac = frac(bestSlice);
      return r;
    }
    if(bestOut){
      r.part = bestOut;
      r.fate = HasHits(bestOut) ? kOutComplete : kOutSummary;
      r.clearCosmic = bestOut->IsClearCosmic;
      r.hitFrac = frac(bestOut);
      return r;
    }

    // 3. sus hits dentro de otra partícula
    AnaParticlePD* bestHits = nullptr;
    int nBest = 0;
    for(AnaParticlePD* p : parts){
      const int n = NHitsOf(p, t->ID);
      if(n > nBest || (n == nBest && n > 0 && p->BeamOrigin && !bestHits->BeamOrigin)){ nBest = n; bestHits = p; }
    }
    if(bestHits){
      r.part = bestHits;
      r.fate = bestHits->BeamOrigin ? kSliceInOther : kOutInOther;
      r.clearCosmic = bestHits->IsClearCosmic;
      r.hitFrac = t->NHitsTrue[2] > 0 ? double(nBest)/t->NHitsTrue[2] : -1;
      return r;
    }

    // 4.
    r.fate = nTrue > 0 ? kNotReco : kNoHits;
    return r;
  }

  struct Counts {
    long all = 0;
    long fate[kNFates] = {0};
    double hitFracSum[kNFates] = {0};
    long hitFracN[kNFates] = {0};
    void Add(const Result& r){
      all++;
      fate[r.fate]++;
      if(r.clearCosmic && (r.fate == kOutComplete || r.fate == kOutSummary))fate[kOutClearCosmic]++;
      if(r.hitFrac >= 0){ hitFracSum[r.fate] += r.hitFrac; hitFracN[r.fate]++; }
    }
    void Print(const std::string& title) const {
      std::cout << title << ": " << all << std::endl;
      for(int f = 0; f < kNFates; f++){
        std::string label = std::string("  ") + FateName(f);
        int len = 0;
        for(unsigned char c : label)if((c & 0xC0) != 0x80)len++;
        std::cout << label << std::string(len < 50 ? 50-len : 1, ' ') << Frac(fate[f], all);
        if(hitFracN[f])std::cout << Form("   hits del plano 2 en la partícula: %.0f %%", 100.*hitFracSum[f]/hitFracN[f]);
        std::cout << std::endl;
      }
    }
  };

  // Todo lo que se cuenta, para poder llamarlo spill a spill (también desde el test)
  struct KaonFate {
    int kaonPDG = 321;
    long nEntries = 0, nMC = 0, nBeamTruth = 0, nWithSlice = 0, nNoHitCounts = 0;
    long nEventsWithKaon = 0;
    std::map<std::string, long> endProcess;   // proceso final de los K+ secundarios
    Counts beamKaon, secKaon, secKaonDecay, muon;
    // eventos sin beam slice con todos los PFParticles guardados (EventInfo.FullDump): dónde están
    Counts secKaonNoSlice, muonNoSlice;
    long joint[kNCoarse][kNCoarse] = {{0}};   // K+ (fila) x mu+ (columna), K+ secundarios con mu+
    long muDaughterOfKaonReco = 0, bothInSlice = 0;
    long nKaonsNoMuon = 0;

    void Process(AnaSpillPD* spill){
      nEntries++;
      const AnaEventInfoPD* info = static_cast<AnaEventInfoPD*>(spill->EventInfo);
      if(!info || !info->IsMC)return;
      nMC++;
      const AnaBeamPD* beam = static_cast<AnaBeamPD*>(spill->Beam);
      AnaTrueParticleB* beamTrue = (beam && beam->BeamParticle) ? static_cast<AnaTrueParticleB*>(beam->BeamParticle->TrueObject) : nullptr;
      if(!beamTrue)return;
      nBeamTruth++;

      std::map<int, AnaTrueParticleB*> byID;
      bool anyHitCounts = false;
      for(AnaTrueParticleB* t : spill->TrueParticles){
        byID[t->ID] = t;
        const AnaTrueParticlePD* tpd = static_cast<AnaTrueParticlePD*>(t);
        if(tpd->NHitsTrue[0] + tpd->NHitsTrue[1] + tpd->NHitsTrue[2] > 0)anyHitCounts = true;
      }
      if(!anyHitCounts)nNoHitCounts++;   // sin SaveTrueHitCounts todo sale "sin hits en el TPC"

      std::vector<AnaParticlePD*> parts;
      bool hasSlice = false;
      if(!spill->Bunches.empty()){
        for(AnaParticleB* p : static_cast<AnaBunchPD*>(spill->Bunches[0])->Particles){
          parts.push_back(static_cast<AnaParticlePD*>(p));
          if(parts.back()->BeamOrigin)hasSlice = true;
        }
      }
      if(hasSlice)nWithSlice++;

      const int muPDG = kaonPDG > 0 ? -13 : 13;
      bool anyKaon = false;
      for(AnaTrueParticleB* tb : spill->TrueParticles){
        if(tb->PDG != kaonPDG || !DescendsFrom(tb->ID, beamTrue->ID, byID))continue;
        anyKaon = true;
        const AnaTrueParticlePD* k = static_cast<AnaTrueParticlePD*>(tb);
        const Result rk = Classify(k, parts, hasSlice);
        if(k->ID == beamTrue->ID){ beamKaon.Add(rk); continue; }
        secKaon.Add(rk);
        const bool dumpNoSlice = !hasSlice && info->FullDump;
        if(dumpNoSlice)secKaonNoSlice.Add(Classify(k, parts, true));
        endProcess[k->ConvertProcess(k->ProcessEnd)]++;
        if(k->ProcessEnd != AnaTrueParticleB::Decay)continue;
        secKaonDecay.Add(rk);

        // mu+ de la desintegración (si hay más de uno, el primero: no debería)
        const AnaTrueParticlePD* mu = nullptr;
        for(int d : k->Daughters){
          auto it = byID.find(d);
          if(it != byID.end() && it->second->PDG == muPDG){ mu = static_cast<AnaTrueParticlePD*>(it->second); break; }
        }
        if(!mu){ nKaonsNoMuon++; continue; }
        const Result rm = Classify(mu, parts, hasSlice);
        muon.Add(rm);
        if(dumpNoSlice)muonNoSlice.Add(Classify(mu, parts, true));
        joint[ToCoarse(rk.fate)][ToCoarse(rm.fate)]++;
        if(rk.part && rm.part && rk.part->BeamOrigin && rm.part->BeamOrigin){
          bothInSlice++;
          if(rm.part->ParentID == rk.part->UniqueID)muDaughterOfKaonReco++;
        }
      }
      if(anyKaon)nEventsWithKaon++;
    }

    void Print() const {
      std::cout << std::endl << "Entradas: " << nEntries << " (MC: " << nMC << ", con partícula true del haz: " << nBeamTruth
                << ", con beam slice: " << nWithSlice << ")" << std::endl;
      if(nNoHitCounts)std::cout << "AVISO: " << nNoHitCounts << " eventos sin NHitsTrue (¿SaveTrueHitCounts false?): ahí todo sale como sin hits" << std::endl;
      std::cout << "Eventos con algún K (PDG " << kaonPDG << ") en la cadena del haz: " << nEventsWithKaon << std::endl << std::endl;
      beamKaon.Print("K de la partícula del haz");
      std::cout << std::endl;
      secKaon.Print("K secundarios (de la cadena del haz)");
      std::cout << "  Proceso final:" << std::endl;
      for(auto& kv : endProcess)std::cout << "    " << kv.first << ": " << Frac(kv.second, secKaon.all) << std::endl;
      std::cout << std::endl;
      secKaonDecay.Print("K secundarios que se desintegran");
      std::cout << "  sin mu en la truth (otros modos: pi+ pi0, ...): " << Frac(nKaonsNoMuon, secKaonDecay.all) << std::endl << std::endl;
      muon.Print("mu de esas desintegraciones");
      std::cout << std::endl << "K secundario (filas) frente a su mu (columnas):" << std::endl << "  " << std::string(14, ' ');
      for(int c = 0; c < kNCoarse; c++)std::cout << Form("%12s", CoarseName(c));
      std::cout << std::endl;
      for(int r = 0; r < kNCoarse; r++){
        std::cout << Form("  %-14s", CoarseName(r));
        for(int c = 0; c < kNCoarse; c++)std::cout << Form("%12ld", joint[r][c]);
        std::cout << std::endl;
      }
      if(secKaonNoSlice.all || muonNoSlice.all){
        std::cout << std::endl << "Eventos sin beam slice con volcado completo (FullDump): dónde están" << std::endl;
        secKaonNoSlice.Print("K secundarios");
        muonNoSlice.Print("mu de sus desintegraciones");
      }
      std::cout << "  K y mu en la beam slice: " << bothInSlice << "; de ellos, el mu es hija de la K en Pandora: "
                << Frac(muDaughterOfKaonReco, bothInSlice) << std::endl << std::endl;
    }
  };
}

int kaon_fate_minitrees_impl(const char* input, int kaonPDG, Long64_t maxEntries){
  using namespace kfate;
  const std::vector<std::string> files = InputFiles(input);
  if(files.empty()){ std::cerr << "No hay ficheros en " << input << std::endl; return -1; }
  const std::string treeName = TreeName(TString(files[0]).Contains("*") ? std::string("") : files[0]);
  TChain chain(treeName.empty() ? "highlandana/MiniTree" : treeName.c_str());
  for(const std::string& f : files)chain.Add(f.c_str());
  if(chain.LoadTree(0) < 0){ std::cerr << "No hay entradas (¿árbol " << chain.GetName() << "?)" << std::endl; return -1; }

  AnaSpillPD* spill = nullptr;
  chain.SetBranchAddress("Spill", &spill);
  KaonFate kf;
  kf.kaonPDG = kaonPDG;
  LoopMiniTree(chain, spill, maxEntries, files.size(), [&kf](AnaSpillPD* s){ kf.Process(s); });
  std::cout << std::endl << "Ficheros: " << files.size() << std::endl;
  kf.Print();
  return 0;
}
