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
// Al final, la composición para la medida A por especie del haz (Composition, 23-nov).
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

  // Composición para la medida A (23-nov): por especie del haz true, cuántos mu+ de K+ en reposo son de la
  // primera interacción (A), de K secundarias o de la K+ del haz, y de qué madre y con qué momento salen las K
  // secundarias (para asignar una incertidumbre de modelo por componente).
  const char* BeamName(int pdg){
    switch(pdg){ case 211: return "pi+"; case 2212: return "p"; case 321: return "K+"; case -13: return "mu+"; case -11: return "e+"; default: return "otro"; }
  }
  const char* MomBin(double p){   // GeV/c
    return p < 0 ? "sin momento" : (p < 0.5 ? "< 0,5" : (p < 1 ? "0,5-1" : (p < 2 ? "1-2" : (p < 4 ? "2-4" : "> 4"))));
  }
  const char* MotherName(int pdg){
    switch(pdg){ case 211: return "pi+"; case -211: return "pi-"; case 2212: return "p"; case 2112: return "n"; case 130: return "K0L"; case 310: return "K0S";
                 case 321: return "K+"; case -321: return "K-"; case 3122: return "Lambda"; case 3222: return "Sigma+"; case 3112: return "Sigma-"; default: return "otra"; }
  }
  struct Composition {
    std::map<std::string, std::map<std::string, long>> n;             // haz -> clase -> mu
    std::map<std::string, std::map<std::string, long>> secMother;     // haz -> madre -> mu (secundarias)
    std::map<std::string, std::map<std::string, long>> secMotherMom;  // madre -> bin de momento -> mu (todas las especies del haz)
    std::map<std::string, std::vector<double>> beamMomA;              // haz -> momento del haz al final (A)
    long secAtMotherEnd = 0, secMidMother = 0, secZeroEndMom = 0;
    // Opción a de la definición de A (decidida el 23-nov): subir también por la partícula del haz que sobrevive a
    // una interacción casi elástica (misma especie que el haz). Para cada secundaria cuya cadena hasta el haz es
    // solo de partículas de la especie del haz, la fracción mínima de momento conservada en cada paso
    // (p inicial de la hija / p final de la madre).
    std::map<std::string, std::vector<double>> leadingRatio;   // haz -> razones (solo cadenas de la especie del haz)

    void AddA(int beamPDG, const AnaTrueParticleB* beam){
      n[BeamName(beamPDG)]["A (primera interacción)"]++;
      beamMomA[BeamName(beamPDG)].push_back(beam->MomentumEnd);
    }
    void AddBeamK(int beamPDG){ n[BeamName(beamPDG)]["K+ del haz"]++; }
    void AddLeading(int beamPDG, double ratio){ leadingRatio[BeamName(beamPDG)].push_back(ratio); }
    void AddSecondary(int beamPDG, const AnaTrueParticleB* k, const AnaTrueParticleB* mother){
      n[BeamName(beamPDG)]["K secundaria"]++;
      secMother[BeamName(beamPDG)][MotherName(mother->PDG)]++;
      // ¿la K nace en el final de la madre (interacción que acaba con ella)? Entonces su momento en la
      // interacción es MomentumEnd; si no, no lo sabemos con la truth guardada
      double d2 = 0; for(int i = 0; i < 3; i++)d2 += std::pow(k->Position[i] - mother->PositionEnd[i], 2);
      double pm = -1;
      if(d2 < 1.){ secAtMotherEnd++; pm = mother->MomentumEnd; if(pm <= 0)secZeroEndMom++; }
      else secMidMother++;
      secMotherMom[MotherName(mother->PDG)][MomBin(pm)]++;
    }
    static double Median(std::vector<double> v){ if(v.empty())return -1; std::sort(v.begin(), v.end()); return v[v.size()/2]; }
    void Print() const {
      std::cout << std::endl << "Composición para la medida A (mu+ de K+ en reposo con la cadena hasta el haz que dispara), por especie del haz true:" << std::endl;
      for(auto& kv : n){
        long tot = 0; for(auto& c : kv.second)tot += c.second;
        std::cout << Form("  haz %-5s %6ld:", kv.first.c_str(), tot);
        for(auto& c : kv.second)std::cout << Form("  %s %ld (%.1f %%)", c.first.c_str(), c.second, 100.*c.second/tot);
        auto itA = beamMomA.find(kv.first);
        if(itA != beamMomA.end())std::cout << Form("  [A: p del haz al final, mediana %.2f GeV/c]", Median(itA->second));
        std::cout << std::endl;
        auto itS = secMother.find(kv.first);
        if(itS != secMother.end()){
          std::cout << "      madre de la primera K (secundarias):";
          std::vector<std::pair<long, std::string>> v; for(auto& m : itS->second)v.push_back({m.second, m.first});
          std::sort(v.rbegin(), v.rend());
          for(auto& m : v)std::cout << Form(" %s %ld", m.second.c_str(), m.first);
          std::cout << std::endl;
        }
      }
      std::cout << "  Opción a: secundarias cuya madre desciende del haz solo por partículas de su especie (haz que sobrevive);" << std::endl
                << "  pasarían a A si conservan en cada paso al menos una fracción r del momento:" << std::endl;
      for(auto& kv : leadingRatio){
        long nSec = 0; auto itn = n.find(kv.first);
        if(itn != n.end()){ auto it = itn->second.find("K secundaria"); if(it != itn->second.end())nSec = it->second; }
        std::cout << Form("    haz %-5s secundarias %ld, con cadena de la especie del haz %zu; pasan a A con r >=", kv.first.c_str(), nSec, kv.second.size());
        for(double r : {0.5, 0.7, 0.8, 0.9, 0.95}){
          long m = 0; for(double x : kv.second)if(x >= r)m++;
          std::cout << Form("  %.2f: %ld", r, m);
        }
        std::cout << std::endl;
      }
      std::cout << "  Secundarias: la K nace en el final de la madre " << secAtMotherEnd << " (de ellas con MomentumEnd <= 0: " << secZeroEndMom
                << "), a mitad de la madre " << secMidMother << " (sin momento)" << std::endl
                << "  Momento de la madre en la interacción (GeV/c), todas las especies del haz:" << std::endl;
      const char* bins[] = {"< 0,5", "0,5-1", "1-2", "2-4", "> 4", "sin momento"};
      for(auto& kv : secMotherMom){
        long tot = 0; for(auto& c : kv.second)tot += c.second;
        std::cout << Form("    %-7s %5ld:", kv.first.c_str(), tot);
        for(const char* b : bins){ auto it = kv.second.find(b); std::cout << Form("  %s %ld", b, it == kv.second.end() ? 0L : it->second); }
        std::cout << std::endl;
      }
    }
  };

  struct Origin {
    Composition comp;
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
          if(i0 + 1 == chain.size()){ Add("K+ del haz (tras K+ -> K+): no es producción", Form("%zu K+ -> K+", i0), kp); comp.AddBeamK(cur->PDG); continue; }
          if(chain.size() == 2)Add("haz que dispara: primera interacción", Form("haz %d", cur->PDG), kp);
          else if(i0 + 2 == chain.size())Add("haz que dispara: primera interacción, con K+ -> K+ intermedias", Form("%zu K+ -> K+", i0), kp);
          else Add("haz que dispara: interacción secundaria", Form("madre de la primera K %d (generación %d)", chain[i0 + 1]->PDG, chain[i0]->Generation), kp);
          if(i0 + 2 == chain.size())comp.AddA(cur->PDG, cur);
          else {
            comp.AddSecondary(cur->PDG, chain[i0], chain[i0 + 1]);
            // opción a: ¿la cadena de la madre de la primera K hasta el haz es solo de la especie del haz?
            bool same = true; double rmin = 1e9;
            for(size_t j = i0 + 1; j + 1 < chain.size(); j++){
              if(chain[j]->PDG != cur->PDG){ same = false; break; }
              const double pEnd = chain[j + 1]->MomentumEnd;
              rmin = std::min(rmin, pEnd > 0 ? chain[j]->Momentum/pEnd : -1.);
            }
            if(same)comp.AddLeading(cur->PDG, rmin);
          }
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
      comp.Print();
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
