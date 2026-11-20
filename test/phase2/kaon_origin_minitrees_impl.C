// Implementación de kaon_origin_minitrees.C (se compila con ACLiC).
//
// De dónde vienen los mu+ de K+ en reposo (mu+ creados en Decay con p true en [230, 241] MeV/c, no
// cósmicos): se sube por la cadena de madres de la truth hasta la raíz y se clasifica:
//   - la partícula del haz que dispara (beam->BeamParticle): K+ de la primera interacción (madre de la K =
//     haz) o de interacciones secundarias (la K desciende del haz por otras partículas);
//   - otra partícula primaria del haz (Generation 0, Origin 4, no la que dispara: pileup del haz);
//   - cadena cortada: alguna madre no está guardada en el minitree (se da el PDG de la última guardada y,
//     si la hay, ParentPDG de esa, que el módulo rellena aunque la madre no se guarde).
// Para cada clase: PDG de la madre de la K, proceso de creación de la K y PDG de la raíz.
#include "kaon_fate_minitrees_impl.C"
#include <map>

namespace korigin {
  using kfate::InputFiles; using kfate::TreeName;

  const char* Proc(int p){
    static const char* n[] = {"Unknown", "primary", "Decay", "kaonplusInelastic", "kaonminusInelastic", "neutronInelastic", "hadElastic",
                              "nCapture", "protonInelastic", "piplusInelastic", "piminusInelastic", "hBertiniCaptureAtRest", "CoulombScat",
                              "kaon0LInelastic", "electronNuclear", "muMinusCaptureAtRest", "dInelastic", "kaon0SInelastic", "positronNuclear",
                              "lambdaInelastic", "tInelastic"};
    return p >= 0 && p <= 20 ? n[p] : "?";
  }

  // Clases gruesas de posición, tiempo (respecto al haz que dispara) y momento, para el detalle
  const char* ZBin(double z){ return z < -50 ? "z < -50" : (z < 0 ? "-50 < z < 0" : (z < 50 ? "0-50" : (z < 220 ? "50-220" : "> 220"))); }
  const char* TBin(double t, const AnaTrueParticleB* beam){
    if(!beam)return "sin haz";
    const double dt = t - beam->Position[3];   // ns
    return std::abs(dt) < 10 ? "|dt| < 10 ns" : (std::abs(dt) < 1000 ? "|dt| < 1 us" : "|dt| > 1 us");
  }
  const char* PBin(double p){ return p < 0.3 ? "< 0,3 GeV/c" : (p < 1 ? "0,3-1 GeV/c" : "> 1 GeV/c"); }

  struct Origin {
    long nEvents = 0, nNoBeamTruth = 0, nMu = 0;
    std::map<std::string, long> cls;                       // clase
    std::map<std::string, std::map<std::string, long>> detail;   // clase -> detalle
    std::map<std::string, std::map<std::string, long>> kProc;    // clase -> proceso de creación de la K

    void Add(const std::string& c, const std::string& d, const std::string& kp){ cls[c]++; detail[c][d]++; kProc[c][kp]++; }

    void Process(AnaSpillPD* spill){
      nEvents++;
      const AnaEventInfoPD* info = static_cast<AnaEventInfoPD*>(spill->EventInfo);
      if(!info || !info->IsMC)return;
      const AnaBeamPD* beam = static_cast<AnaBeamPD*>(spill->Beam);
      const AnaTrueParticleB* beamTrue = (beam && beam->BeamParticle) ? static_cast<AnaTrueParticleB*>(beam->BeamParticle->TrueObject) : nullptr;
      if(!beamTrue)nNoBeamTruth++;
      std::map<int, AnaTrueParticleB*> byID;
      std::map<int, const AnaTrueParticleB*> listedBy;   // ID de hija -> partícula guardada que la lista
      for(AnaTrueParticleB* t : spill->TrueParticles){
        byID[t->ID] = t;
        for(int d : static_cast<AnaTrueParticlePD*>(t)->Daughters)listedBy[d] = t;
      }
      // Dónde se corta la cadena: quién lista como hija a la partícula que falta
      auto Gap = [&](int missingID) -> std::string {
        auto it = listedBy.find(missingID);
        if(it == listedBy.end())return "nadie la lista como hija";
        const AnaTrueParticleB* p = it->second;
        const bool inChain = static_cast<const AnaTrueParticlePD*>(p)->Generation >= 0;
        return Form("hija no guardada de %d (%s, generación %d)", p->PDG, inChain ? "en la cadena del haz" : "suelta",
                    static_cast<const AnaTrueParticlePD*>(p)->Generation);
      };

      for(AnaTrueParticleB* t : spill->TrueParticles){
        const AnaTrueParticlePD* mu = static_cast<AnaTrueParticlePD*>(t);
        if(mu->PDG != -13 || mu->ProcessStart != AnaTrueParticleB::Decay || mu->Origin == 2)continue;
        if(mu->Momentum < 0.230 || mu->Momentum > 0.241)continue;   // GeV/c en la truth del módulo
        nMu++;
        // madre del mu
        auto itK = byID.find(mu->ParentID);
        if(itK == byID.end()){
          Add("madre del mu no guardada", Gap(mu->ParentID), "-");
          continue;
        }
        const AnaTrueParticlePD* k = static_cast<AnaTrueParticlePD*>(itK->second);
        if(k->PDG != 321){ Add("madre del mu no es K+", Form("PDG madre %d", k->PDG), Proc(k->ProcessStart)); continue; }
        // subir desde la K
        std::vector<const AnaTrueParticlePD*> chain{k};
        const AnaTrueParticlePD* cur = k;
        bool cut = false;
        for(int depth = 0; depth < 1000 && cur->ParentID != 0; depth++){
          auto it = byID.find(cur->ParentID);
          if(it == byID.end()){ cut = true; break; }
          cur = static_cast<AnaTrueParticlePD*>(it->second);
          chain.push_back(cur);
        }
        const std::string kp = Proc(k->ProcessStart);
        const int kParentPDG = chain.size() > 1 ? chain[1]->PDG : k->ParentPDG;
        if(cut){
          Add("cadena cortada", Form("última guardada %d; su madre: %s", cur->PDG, Gap(cur->ParentID).c_str()), kp);
          continue;
        }
        if(beamTrue && cur->ID == beamTrue->ID){
          // primera K de la cadena: se suben las K+ -> K+ (kaonplusInelastic con una K+ que sobrevive)
          size_t i0 = 0;
          while(i0 + 1 < chain.size() && chain[i0 + 1]->PDG == 321)i0++;
          if(i0 + 1 == chain.size()){ Add("K+ del haz (tras K+ -> K+): no es producción", Form("%zu K+ -> K+", i0), kp); continue; }
          if(chain.size() == 2)Add("haz que dispara: primera interacción", Form("haz %d", cur->PDG), kp);
          else if(i0 + 2 == chain.size())Add("haz que dispara: primera interacción, con K+ -> K+ intermedias", Form("%zu K+ -> K+", i0), kp);
          else Add("haz que dispara: interacción secundaria", Form("madre de la primera K %d (generación %d)", chain[i0 + 1]->PDG, chain[i0]->Generation), kp);
          continue;
        }
        if(cur->Generation == 0 && cur->Origin == 4){
          Add(chain.size() == 2 ? "otra partícula del haz (pileup): primera interacción" : "otra partícula del haz (pileup): secundaria",
              Form("raíz %d, madre de la K %d", cur->PDG, kParentPDG), kp);
          continue;
        }
        Add("otra raíz", Form("raíz %d (Generation %d, Origin %d), madre de la K %d", cur->PDG, cur->Generation, cur->Origin, kParentPDG), kp);
      }
    }

    static void Top(const std::map<std::string, long>& m, long tot, size_t n){
      std::vector<std::pair<long, std::string>> v; for(auto& kv : m)v.push_back({kv.second, kv.first});
      std::sort(v.rbegin(), v.rend());
      for(size_t i = 0; i < std::min(n, v.size()); i++)std::cout << Form("      %-90s %6ld (%5.1f %%)", v[i].second.c_str(), v[i].first, 100.*v[i].first/tot) << std::endl;
    }

    void Print() const {
      std::cout << std::endl << "Entradas " << nEvents << " (sin truth del haz " << nNoBeamTruth << ")" << std::endl
                << "mu+ de desintegración, no cósmicos, p true en [230, 241] MeV/c: " << nMu << std::endl;
      std::vector<std::pair<long, std::string>> v; for(auto& kv : cls)v.push_back({kv.second, kv.first});
      std::sort(v.rbegin(), v.rend());
      for(auto& x : v){
        std::cout << Form("  %-55s %6ld (%5.1f %%)", x.second.c_str(), x.first, 100.*x.first/nMu) << std::endl;
        std::cout << "    detalle:" << std::endl; Top(detail.at(x.second), x.first, 8);
        std::cout << "    proceso de creación de la K:" << std::endl; Top(kProc.at(x.second), x.first, 5);
      }
    }
  };
}

int kaon_origin_minitrees_impl(const char* input, Long64_t maxEntries){
  using namespace korigin;
  const std::vector<std::string> files = InputFiles(input);
  if(files.empty()){ std::cerr << "No hay ficheros en " << input << std::endl; return -1; }
  const std::string treeName = TreeName(TString(files[0]).Contains("*") ? std::string("") : files[0]);
  TChain chain(treeName.empty() ? "highlandana/MiniTree" : treeName.c_str());
  for(const std::string& f : files)chain.Add(f.c_str());
  if(chain.LoadTree(0) < 0){ std::cerr << "No hay entradas (¿árbol " << chain.GetName() << "?)" << std::endl; return -1; }
  AnaSpillPD* spill = nullptr;
  chain.SetBranchAddress("Spill", &spill);
  Origin o;
  LoopMiniTree(chain, spill, maxEntries, files.size(), [&o](AnaSpillPD* s){ o.Process(s); });
  std::cout << std::endl << "Ficheros: " << files.size() << std::endl;
  o.Print();
  return 0;
}
