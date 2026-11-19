// Implementación de count_module_minitrees.C (se compila con ACLiC).
#include "pdDataClasses.hxx"
#include "minitree_loop.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include <algorithm>
#include <cmath>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace {

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

  // El árbol es highlandana/MiniTree (salida de lar) o MiniTree
  std::string TreeName(const std::string& file){
    if(file.empty())return "";   // patrón: el árbol por defecto (sin TFile::Open(""), que se queja)
    std::unique_ptr<TFile> f(TFile::Open(file.c_str()));
    if(!f || f->IsZombie())return "";
    if(f->Get("highlandana/MiniTree"))return "highlandana/MiniTree";
    if(f->Get("MiniTree"))return "MiniTree";
    return "";
  }

  // Contador con su total, para imprimir "n (x %)"
  std::string Frac(long n, long total){
    return Form("%ld (%.1f %%)", n, total > 0 ? 100.*n/total : 0.);
  }

  // setw cuenta bytes: rellenar a mano contando caracteres UTF-8
  void Line(const std::string& label, const std::string& value){
    int len = 0;
    for(unsigned char c : label)if((c & 0xC0) != 0x80)len++;
    std::cout << "  " << label << std::string(len < 58 ? 58-len : 1, ' ') << value << std::endl;
  }

  // ¿Es "id" la partícula del haz o una descendiente suya? (las partículas true guardadas son la
  // cadena del haz más las que emparejó la reco, con sus ParentID)
  bool DescendsFrom(int id, int ancestor, const std::map<int, AnaTrueParticleB*>& byID){
    for(int depth = 0; depth < 1000; depth++){
      if(id == ancestor)return true;
      auto it = byID.find(id);
      if(it == byID.end() || it->second->ParentID == 0)return false;
      id = it->second->ParentID;
    }
    return false;
  }

  // Las "nmax" claves más frecuentes de un recuento, de más a menos
  template<class K> void PrintTop(const std::map<K, long>& counts, long total, size_t nmax, const std::string& indent){
    std::vector<std::pair<long, K>> v;
    for(auto& kv : counts)v.emplace_back(kv.second, kv.first);
    std::sort(v.begin(), v.end(), [](const std::pair<long, K>& a, const std::pair<long, K>& b){ return a.first > b.first; });
    long rest = 0;
    for(size_t i = 0; i < v.size(); i++){
      if(i < nmax){ std::ostringstream os; os << indent << v[i].second; Line(os.str(), Frac(v[i].first, total)); }
      else rest += v[i].first;
    }
    if(rest)Line(indent + "resto", Frac(rest, total));
  }

  // ¿Guardado completo? (hits, puntos de trayectoria o space points)
  bool IsComplete(const AnaParticlePD* p){
    if(!p->TrjPoints.empty() || !p->SpacePoints.empty())return true;
    for(int pl = 0; pl < 3; pl++)if(!p->Hits[pl].empty() || !p->HitsNoCalo[pl].empty())return true;
    return false;
  }

  // ¿Es un mu+ de K+ -> mu+ nu?
  bool IsMuFromK(const AnaTrueParticleB* t, const std::map<int, AnaTrueParticleB*>& byID){
    if(!t || t->PDG != -13)return false;
    auto ip = byID.find(t->ParentID);
    return ip != byID.end() && ip->second->PDG == 321 && ip->second->ProcessEnd == AnaTrueParticleB::Decay;
  }

  // Jerarquías de los PFParticles de fuera de la beam slice, por su primario
  std::map<AnaParticlePD*, std::vector<AnaParticlePD*>> Hierarchies(const std::vector<AnaParticlePD*>& out){
    std::map<int, AnaParticlePD*> byUID;
    for(AnaParticlePD* p : out)byUID[p->UniqueID] = p;
    std::map<AnaParticlePD*, std::vector<AnaParticlePD*>> h;
    for(AnaParticlePD* p : out){
      AnaParticlePD* t = p;
      for(int i = 0; i < 1000; i++){
        auto it = byUID.find(t->ParentID);
        if(it == byUID.end())break;
        t = it->second;
      }
      h[t].push_back(p);
    }
    return h;
  }

  struct Mean{
    double sum = 0; long n = 0;
    void Add(double x){ if(x < -900)return; sum += x; n++; }   // sin -999 (sin rellenar)
    std::string Str(const char* fmt = "%.2f") const { return n ? Form(fmt, sum/n) : std::string("-"); }
  };

  // Cósmicos según AnaParticlePD::CosmicCriteria (módulo desde el 19-nov-2026)
  //  - eventos normales con beam slice: jerarquías completas por evento, y cuántas solo por el
  //    candidato a mu (bit 4 sin bit 1)
  //  - eventos normales sin beam slice (solo se guardan las jerarquías con candidato a mu): cuántos
  //    tienen alguna, cuántas y de qué tamaño; incoherencias (sin bit 4 o en resumen); en MC, pureza
  //  - MC con volcado completo sin beam slice: lo que guardaría un evento normal (bit 4) con los mu
  //    de K+ de la cadena del haz, y si su K está en la misma jerarquía
  // Los eventos sin beam slice y sin cósmicos guardados solo se reconocen como del módulo nuevo si
  // tienen NBeamRegionSpacePoints rellenado (EmptyEventFinder, del mismo día)
  struct CritCounter {
    long events = 0;                          // eventos normales del módulo nuevo
    Mean compH, compP, muOnlyH, muOnlyP;      // con beam slice
    long nsNormal = 0, nsWith = 0, nsBadCrit = 0, nsIncomplete = 0;   // sin beam slice
    Mean nsH, nsP;
    long nsHier = 0, nsHierMuK = 0;
    long fdNS = 0, fdMu = 0, fdMuNoReco = 0, fdMuBit4 = 0, fdMuBit4K = 0, fdBit4H = 0, fdBit4HMuK = 0;
    Mean fdBit4P;

    void Process(AnaSpillPD* spill, AnaTrueParticleB* beamTrue){
      const AnaEventInfoPD* info = static_cast<AnaEventInfoPD*>(spill->EventInfo);
      if(!info || spill->Bunches.empty())return;
      bool hasSlice = false, filled = false;
      std::vector<AnaParticlePD*> out;
      for(AnaParticleB* p : static_cast<AnaBunchPD*>(spill->Bunches[0])->Particles){
        AnaParticlePD* part = static_cast<AnaParticlePD*>(p);
        if(part->BeamOrigin){ hasSlice = true; continue; }
        out.push_back(part);
        if(part->CosmicCriteria > -900)filled = true;
      }
      if(!filled && !(!hasSlice && out.empty() && info->NBeamRegionSpacePoints >= 0))return;
      const bool isMC = info->IsMC;
      const auto hiers = Hierarchies(out);
      std::map<int, AnaTrueParticleB*> byID;
      if(isMC)for(AnaTrueParticleB* t : spill->TrueParticles)byID[t->ID] = t;
      auto hasMuK = [&byID](const std::vector<AnaParticlePD*>& ps){
        for(const AnaParticlePD* p : ps)if(IsMuFromK(static_cast<AnaTrueParticleB*>(p->TrueObject), byID))return true;
        return false;
      };
      auto bit = [](int c, int b){ return c > -900 && (c & b); };

      if(!info->FullDump){
        events++;
        int nH = 0, nP = 0, nMuH = 0, nMuP = 0;
        for(const auto& h : hiers){
          const int c = h.first->CosmicCriteria;
          bool complete = false;
          for(const AnaParticlePD* p : h.second)complete = complete || IsComplete(p);
          if(hasSlice){
            if(!complete)continue;
            nH++; nP += h.second.size();
            if(bit(c, AnaParticlePD::kCosmicMuCandidate) && !bit(c, AnaParticlePD::kCosmicNearBeamSlice)){ nMuH++; nMuP += h.second.size(); }
          }
          else{
            nH++; nP += h.second.size();
            if(!bit(c, AnaParticlePD::kCosmicMuCandidate))nsBadCrit++;
            if(!complete)nsIncomplete++;
            if(isMC){ nsHier++; if(hasMuK(h.second))nsHierMuK++; }
          }
        }
        if(hasSlice){ compH.Add(nH); compP.Add(nP); muOnlyH.Add(nMuH); muOnlyP.Add(nMuP); }
        else{ nsNormal++; if(nH)nsWith++; nsH.Add(nH); nsP.Add(nP); }
        return;
      }
      if(hasSlice || !isMC || !beamTrue)return;

      // volcado completo sin beam slice: todo está guardado; el bit 4 dice qué guardaría un evento normal
      fdNS++;
      std::map<int, const AnaParticlePD*> best;   // ID true -> PFParticle con más hits
      std::map<const AnaParticlePD*, const AnaParticlePD*> topOf;
      for(const auto& h : hiers){
        if(bit(h.first->CosmicCriteria, AnaParticlePD::kCosmicMuCandidate)){
          fdBit4H++; fdBit4P.Add(h.second.size());
          if(hasMuK(h.second))fdBit4HMuK++;
        }
        for(const AnaParticlePD* p : h.second){
          topOf[p] = h.first;
          if(!p->TrueObject)continue;
          const AnaParticlePD*& b = best[p->TrueObject->ID];
          if(!b || p->NHits > b->NHits)b = p;
        }
      }
      for(AnaTrueParticleB* t : spill->TrueParticles){
        if(!IsMuFromK(t, byID) || !DescendsFrom(t->ID, beamTrue->ID, byID))continue;
        fdMu++;
        auto im = best.find(t->ID);
        if(im == best.end()){ fdMuNoReco++; continue; }
        if(!bit(topOf[im->second]->CosmicCriteria, AnaParticlePD::kCosmicMuCandidate))continue;
        fdMuBit4++;
        auto ik = best.find(t->ParentID);   // la K, ¿en la misma jerarquía?
        if(ik != best.end() && topOf[ik->second] == topOf[im->second])fdMuBit4K++;
      }
    }

    void Print() const {
      if(!events && !fdNS)return;
      std::cout << std::endl << "Cósmicos por criterio del módulo (CosmicCriteria; eventos sin volcado completo)" << std::endl;
      Line("con beam slice: jerarquías completas por evento", compH.Str() + " (" + compP.Str() + " PFParticles)");
      Line("  de ellas solo por el candidato a mu (bit 4 sin bit 1)", muOnlyH.Str() + " (" + muOnlyP.Str() + " PFParticles)");
      Line("sin beam slice: eventos", Form("%ld", nsNormal));
      Line("  con alguna jerarquía de candidato a mu guardada", Frac(nsWith, nsNormal));
      Line("  jerarquías por evento", nsH.Str() + " (" + nsP.Str() + " PFParticles)");
      if(nsBadCrit || nsIncomplete){
        std::cout << "  INCOHERENCIAS (sin beam slice solo deberían ir jerarquías completas con el bit 4):" << std::endl;
        Line("    sin el bit 4", Form("%ld", nsBadCrit));
        Line("    en resumen", Form("%ld", nsIncomplete));
      }
      if(nsHier)Line("  MC: jerarquías con un mu de K+ (pureza)", Frac(nsHierMuK, nsHier));
      if(fdNS){
        std::cout << "  MC con volcado completo y sin beam slice (" << fdNS << " eventos): lo que guardaría un evento normal" << std::endl;
        Line("    mu de K+ de la cadena del haz", Form("%ld", fdMu));
        Line("    sin PFParticle", Frac(fdMuNoReco, fdMu));
        Line("    en una jerarquía con el bit 4 (se guardaría)", Frac(fdMuBit4, fdMu));
        Line("      con su K en la misma jerarquía", Frac(fdMuBit4K, fdMuBit4));
        Line("    jerarquías con el bit 4 por evento", std::string(Form("%.2f", 1.*fdBit4H/fdNS)) + " (" + fdBit4P.Str() + " PFParticles cada una)");
        Line("      con un mu de K+ (pureza, sesgada: hay un K+)", Frac(fdBit4HMuK, fdBit4H));
      }
    }
  };
}

int count_module_minitrees_impl(const char* input, Long64_t maxEntries){

  const std::vector<std::string> files = InputFiles(input);
  if(files.empty()){ std::cerr << "No hay ficheros en " << input << std::endl; return -1; }
  const std::string treeName = TreeName(TString(files[0]).Contains("*") ? std::string("") : files[0]);
  TChain chain(treeName.empty() ? "highlandana/MiniTree" : treeName.c_str());
  for(const std::string& f : files)chain.Add(f.c_str());
  if(chain.LoadTree(0) < 0){ std::cerr << "No hay entradas (¿árbol " << chain.GetName() << "?)" << std::endl; return -1; }

  AnaSpillPD* spill = nullptr;
  chain.SetBranchAddress("Spill", &spill);

  // contadores
  long n = 0, nMC = 0, badDQ = 0, noBunch = 0, withBunch = 0, noSlice = 0, withSlice = 0, missingReco = 0;
  // empty event finder (EventInfo.NBeamRegionSpacePoints, -999 en minitrees anteriores)
  long efFilled = 0, efNoReco3d = 0, efEmpty = 0;
  long efBins[4] = {0};                      // space points en la caja: 0, 1-9, 10-49, >= 50
  long efVsSlice[2][2] = {{0}};              // eventos con bunch: [EmptyEvent][con beam slice]
  long efVsBeamEnd[2][2] = {{0}};            // MC con partícula true del haz: [EmptyEvent][acaba antes del TPC]
  long beamWithMom = 0, beamWithTrack = 0;
  long badDQWithBunch = 0, flagNoSliceButParticles = 0, flagSliceButEmpty = 0;   // incoherencias
  std::map<int, long> beamPDGs;              // datos: hipótesis de la instrumentación del haz (trigger bueno)
  long nPID = 0;
  Mean beamMom, nSlice, nCosmics, nCosmicsWithHits;
  // MC
  long noBeamTruth = 0, beamEndsBeforeTPC = 0, noSliceEndsBeforeTPC = 0, noSliceWithTruth = 0;
  long matched = 0, mismatched = 0, primNoTruth = 0, primDescendant = 0, primCosmic = 0, primOther = 0;
  long beamInSliceNotPrimary = 0;
  struct PerPDG{ long all = 0, noSlice = 0, matched = 0, mismatched = 0, descendant = 0; };
  std::map<int, PerPDG> byPDG;               // por PDG de la partícula true del haz
  std::map<std::string, long> mismatchPDGs;  // "PDG del haz -> PDG emparejado con el primario"
  Mean purMatched, effMatched, purMismatched, effMismatched;
  CritCounter crit;                          // CosmicCriteria (módulo desde el 19-nov-2026)

  // sin GetEntries (abre todos los ficheros antes de empezar) y con el avance por std::cerr, como LoopMiniTree
  const auto t0 = std::chrono::steady_clock::now();
  for(Long64_t i = 0; maxEntries <= 0 || i < maxEntries; i++){
    if(chain.LoadTree(i) < 0)break;
    if(chain.GetEntry(i) <= 0 || !spill)continue;
    spill->CopyVectorsIntoArrays();
    spill->RedoLinks();
    n++;
    if(n % 1000 == 0)
      std::cerr << Form("[avance] %ld entradas, fichero %d de %zu, %.0f s", n, chain.GetTreeNumber()+1, files.size(),
                        std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count()) << std::endl;

    const AnaEventInfoPD* info = static_cast<AnaEventInfoPD*>(spill->EventInfo);
    const bool isMC = info && info->IsMC;
    if(isMC)nMC++;
    if(info && info->MissingReco)missingReco++;
    const bool efOK = info && info->NBeamRegionSpacePoints >= 0;
    if(info && info->NBeamRegionSpacePoints == -1)efNoReco3d++;
    if(efOK){
      efFilled++;
      if(info->EmptyEvent)efEmpty++;
      const int nsp = info->NBeamRegionSpacePoints;
      efBins[nsp == 0 ? 0 : nsp < 10 ? 1 : nsp < 50 ? 2 : 3]++;
    }
    if(spill->DataQuality && !spill->DataQuality->GoodDaq)badDQ++;

    // instrumentación del haz (todas las entradas tienen Beam)
    const AnaBeamPD* beam = static_cast<AnaBeamPD*>(spill->Beam);
    if(beam){
      // nMomenta y nTracks son size_t: -999 (sin evento del haz) sale como un número enorme
      if(beam->nMomenta > 0 && beam->nMomenta < 1000){
        beamWithMom++;
        if(beam->BeamMomentum > 0)beamMom.Add(beam->BeamMomentum);
      }
      if(beam->nTracks  > 0 && beam->nTracks  < 1000)beamWithTrack++;
      // solo con trigger bueno: con trigger malo los Cherenkov están a 0 y a 6-7 GeV salen como protón
      if(!isMC && (!spill->DataQuality || spill->DataQuality->GoodDaq)){
        nPID++;
        if(beam->PDGs.empty())beamPDGs[0]++;
        for(int pdg : beam->PDGs)beamPDGs[pdg]++;
      }
    }

    // partícula true del haz (la de la instrumentación del haz en MC)
    AnaTrueParticleB* beamTrue = nullptr;
    if(isMC && beam && beam->BeamParticle)beamTrue = static_cast<AnaTrueParticleB*>(beam->BeamParticle->TrueObject);
    if(isMC){
      if(!beamTrue)noBeamTruth++;
      else{
        byPDG[beamTrue->PDG].all++;
        if(beamTrue->PositionEnd[2] < 0)beamEndsBeforeTPC++;   // la cara frontal del TPC está en z = 0
        if(efOK)efVsBeamEnd[info->EmptyEvent][beamTrue->PositionEnd[2] < 0]++;
      }
    }

    if(spill->Bunches.empty()){ noBunch++; continue; }
    withBunch++;
    if(spill->DataQuality && !spill->DataQuality->GoodDaq)badDQWithBunch++;
    const AnaBunchPD* bunch = static_cast<AnaBunchPD*>(spill->Bunches[0]);

    // beam slice (BeamOrigin) y cósmicos guardados
    std::vector<AnaParticlePD*> slice;
    int ncos = 0, ncosHits = 0;
    for(AnaParticleB* p : bunch->Particles){
      AnaParticlePD* part = static_cast<AnaParticlePD*>(p);
      if(part->BeamOrigin)slice.push_back(part);
      else{
        ncos++;
        if(!part->Hits[2].empty() || !part->HitsNoCalo[2].empty())ncosHits++;
      }
    }
    const bool hasSlice = !slice.empty();
    crit.Process(spill, beamTrue);
    if(info && info->HasPandora && !hasSlice)flagSliceButEmpty++;
    if(info && !info->HasPandora && hasSlice)flagNoSliceButParticles++;
    if(efOK)efVsSlice[info->EmptyEvent][hasSlice]++;
    if(!hasSlice){
      noSlice++;
      if(isMC && beamTrue){
        noSliceWithTruth++;
        byPDG[beamTrue->PDG].noSlice++;
        if(beamTrue->PositionEnd[2] < 0)noSliceEndsBeforeTPC++;
      }
      continue;
    }
    withSlice++;
    nSlice.Add(slice.size());
    nCosmics.Add(ncos);
    nCosmicsWithHits.Add(ncosHits);

    // MC: ¿el primario de la beam slice es la partícula true del haz?
    if(!isMC || !beamTrue || slice.empty())continue;
    std::map<int, AnaTrueParticleB*> byID;
    for(AnaTrueParticleB* t : spill->TrueParticles)byID[t->ID] = t;

    AnaParticlePD* primary = slice[0];   // el módulo guarda primero el primario y luego sus hijas
    const AnaTrueParticleB* primTrue = static_cast<AnaTrueParticleB*>(primary->TrueObject);
    if(primTrue && primTrue->ID == beamTrue->ID){
      matched++;
      byPDG[beamTrue->PDG].matched++;
      purMatched.Add(primary->TruePur);
      effMatched.Add(primary->TrueEff);
      continue;
    }
    mismatched++;
    byPDG[beamTrue->PDG].mismatched++;
    purMismatched.Add(primary->TruePur);
    effMismatched.Add(primary->TrueEff);
    if(!primTrue)primNoTruth++;
    else if(DescendsFrom(primTrue->ID, beamTrue->ID, byID)){ primDescendant++; byPDG[beamTrue->PDG].descendant++; }
    else if(static_cast<const AnaTrueParticlePD*>(primTrue)->Origin == 2)primCosmic++;
    else primOther++;
    mismatchPDGs[Form("%d -> %s", beamTrue->PDG, primTrue ? Form("%d", primTrue->PDG) : "sin truth")]++;
    // ¿está la partícula del haz en otra partícula de la beam slice (traza rota)?
    for(size_t k = 1; k < slice.size(); k++){
      const AnaTrueParticleB* t = static_cast<AnaTrueParticleB*>(slice[k]->TrueObject);
      if(t && t->ID == beamTrue->ID){ beamInSliceNotPrimary++; break; }
    }
  }

  // resumen
  std::cout << std::endl << "Ficheros: " << files.size() << ", entradas leídas: " << n
            << " (" << nMC << " de MC)" << std::endl << std::endl;

  std::cout << "Eventos" << std::endl;
  Line("DQ malo (GoodDaq = false, sin bunch)", Frac(badDQ, n));
  Line("sin bunch", Frac(noBunch, n));
  Line("con bunch, sin beam slice", Frac(noSlice, n));
  Line("con bunch y beam slice", Frac(withSlice, n));
  Line("sin PFParticles o productos de Pandora (MissingReco)", Frac(missingReco, n));
  if(badDQWithBunch || flagSliceButEmpty || flagNoSliceButParticles){
    std::cout << "  INCOHERENCIAS (no debería haber ninguna con minitrees del módulo):" << std::endl;
    Line("    DQ malo con bunch", Frac(badDQWithBunch, n));
    Line("    HasPandora = true sin partículas de la beam slice", Frac(flagSliceButEmpty, n));
    Line("    HasPandora = false con partículas de la beam slice", Frac(flagNoSliceButParticles, n));
  }
  std::cout << std::endl << "Empty event finder (space points de reco3d en la entrada del haz; EmptyEvent: < 10)" << std::endl;
  if(!efFilled)Line("sin rellenar (minitrees anteriores al 19-nov-2026 o sin reco3d)", Frac(n - efFilled, n));
  else{
    Line("rellenado", Frac(efFilled, n) + (efNoReco3d ? "   sin reco3d: " + Frac(efNoReco3d, n) : std::string("")));
    Line("vacíos (EmptyEvent)", Frac(efEmpty, efFilled));
    Line("space points en la caja: 0 / 1-9 / 10-49 / >= 50", Form("%ld / %ld / %ld / %ld", efBins[0], efBins[1], efBins[2], efBins[3]));
    const long withB = efVsSlice[0][0] + efVsSlice[0][1] + efVsSlice[1][0] + efVsSlice[1][1];
    if(withB){
      std::cout << "  Eventos con bunch (" << withB << "): EmptyEvent frente a beam slice" << std::endl;
      Line("    vacío, sin beam slice", Frac(efVsSlice[1][0], withB));
      Line("    vacío, con beam slice", Frac(efVsSlice[1][1], withB));
      Line("    no vacío, sin beam slice (¿fallo de Pandora?)", Frac(efVsSlice[0][0], withB));
      Line("    no vacío, con beam slice", Frac(efVsSlice[0][1], withB));
    }
    const long withT = efVsBeamEnd[0][0] + efVsBeamEnd[0][1] + efVsBeamEnd[1][0] + efVsBeamEnd[1][1];
    if(withT){
      std::cout << "  MC (" << withT << "): EmptyEvent frente a la partícula true del haz" << std::endl;
      Line("    vacío, acaba antes del TPC", Frac(efVsBeamEnd[1][1], withT));
      Line("    vacío, entra en el TPC", Frac(efVsBeamEnd[1][0], withT));
      Line("    no vacío, acaba antes del TPC", Frac(efVsBeamEnd[0][1], withT));
      Line("    no vacío, entra en el TPC", Frac(efVsBeamEnd[0][0], withT));
    }
  }

  std::cout << std::endl << "Instrumentación del haz" << std::endl;
  Line("con momento reconstruido", Frac(beamWithMom, n) + "   <p> = " + beamMom.Str("%.3f") + " GeV");
  Line("con traza del haz", Frac(beamWithTrack, n));
  if(!beamPDGs.empty()){
    std::cout << "  Hipótesis de PDG de la instrumentación, eventos con trigger bueno (0 = ninguna; a 6-7 GeV" << std::endl
              << "  e/mu/pi salen juntas):" << std::endl;
    PrintTop(beamPDGs, nPID, 10, "    ");
  }
  std::cout << std::endl << "Beam slice y cósmicos (eventos con beam slice)" << std::endl;
  Line("partículas de la beam slice por evento", nSlice.Str());
  Line("cósmicos guardados por evento", nCosmics.Str());
  Line("  de ellos con hits en el plano 2 (completos)", nCosmicsWithHits.Str());

  crit.Print();

  if(nMC > 0){
    std::cout << std::endl << "Truth (MC)" << std::endl;
    Line("sin partícula true del haz", Frac(noBeamTruth, nMC));
    Line("la partícula del haz acaba antes del TPC (z < 0)", Frac(beamEndsBeforeTPC, nMC));
    Line("sin beam slice, con truth del haz", Frac(noSliceWithTruth, nMC));
    Line("  de ellos con el haz acabando antes del TPC", Frac(noSliceEndsBeforeTPC, noSliceWithTruth));
    std::cout << "  Por PDG de la partícula true del haz (% sobre cada fila; primario sobre los eventos con beam slice):" << std::endl;
    std::cout << "        PDG  eventos   sin beam slice   primario = haz   primario = descendiente" << std::endl;
    std::vector<std::pair<long, int>> order;
    for(auto& kv : byPDG)order.emplace_back(kv.second.all, kv.first);
    std::sort(order.rbegin(), order.rend());
    for(auto& o : order){
      const PerPDG& c = byPDG[o.second];
      const long withSlice = c.matched + c.mismatched;
      std::cout << Form("    %7d  %7ld   %14s   %14s   %22s", o.second, c.all,
                        Frac(c.noSlice, c.all).c_str(), Frac(c.matched, withSlice).c_str(),
                        Frac(c.descendant, withSlice).c_str()) << std::endl;
    }

    const long checked = matched + mismatched;
    std::cout << std::endl << "Primario de la beam slice frente a la partícula true del haz (MC con beam slice)" << std::endl;
    Line("es la partícula del haz", Frac(matched, checked) + "   pureza " + purMatched.Str() + ", completitud " + effMatched.Str());
    Line("no lo es", Frac(mismatched, checked) + "   pureza " + purMismatched.Str() + ", completitud " + effMismatched.Str());
    Line("  sin partícula true emparejada", Frac(primNoTruth, checked));
    Line("  una descendiente del haz", Frac(primDescendant, checked));
    Line("  un cósmico", Frac(primCosmic, checked));
    Line("  otra", Frac(primOther, checked));
    Line("  el haz es otra partícula de la beam slice", Frac(beamInSliceNotPrimary, checked));
    if(!mismatchPDGs.empty()){
      std::cout << "  PDG del haz -> PDG emparejado con el primario:" << std::endl;
      PrintTop(mismatchPDGs, mismatched, 10, "    ");
    }
  }
  std::cout << std::endl;
  return 0;
}
