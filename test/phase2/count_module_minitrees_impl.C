// Implementación de count_module_minitrees.C (se compila con ACLiC).
#include "pdDataClasses.hxx"
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

  struct Mean{
    double sum = 0; long n = 0;
    void Add(double x){ if(x < -900)return; sum += x; n++; }   // sin -999 (sin rellenar)
    std::string Str(const char* fmt = "%.2f") const { return n ? Form(fmt, sum/n) : std::string("-"); }
  };
}

int count_module_minitrees_impl(const char* input, Long64_t maxEntries){

  const std::vector<std::string> files = InputFiles(input);
  if(files.empty()){ std::cerr << "No hay ficheros en " << input << std::endl; return -1; }
  const std::string treeName = TreeName(TString(files[0]).Contains("*") ? std::string("") : files[0]);
  TChain chain(treeName.empty() ? "highlandana/MiniTree" : treeName.c_str());
  for(const std::string& f : files)chain.Add(f.c_str());
  const Long64_t nEntries = chain.GetEntries();
  if(nEntries <= 0){ std::cerr << "No hay entradas (¿árbol " << chain.GetName() << "?)" << std::endl; return -1; }

  AnaSpillPD* spill = nullptr;
  chain.SetBranchAddress("Spill", &spill);

  // contadores
  long n = 0, nMC = 0, badDQ = 0, noBunch = 0, withBunch = 0, noSlice = 0, withSlice = 0, emptyEvent = 0;
  long beamWithMom = 0, beamWithTrack = 0;
  long badDQWithBunch = 0, flagNoSliceButParticles = 0, flagSliceButEmpty = 0;   // incoherencias
  std::map<int, long> beamPDGs;              // datos: hipótesis de la instrumentación del haz
  Mean beamMom, nSlice, nCosmics, nCosmicsWithHits;
  // MC
  long noBeamTruth = 0, beamEndsBeforeTPC = 0, noSliceEndsBeforeTPC = 0, noSliceWithTruth = 0;
  long matched = 0, mismatched = 0, primNoTruth = 0, primDescendant = 0, primCosmic = 0, primOther = 0;
  long beamInSliceNotPrimary = 0;
  std::map<int, long> truePDGs;
  std::map<std::string, long> mismatchPDGs;  // "PDG del haz -> PDG emparejado con el primario"
  Mean purMatched, effMatched, purMismatched, effMismatched;

  const Long64_t nMax = (maxEntries > 0 && maxEntries < nEntries) ? maxEntries : nEntries;
  for(Long64_t i = 0; i < nMax; i++){
    if(chain.GetEntry(i) <= 0 || !spill)continue;
    spill->CopyVectorsIntoArrays();
    spill->RedoLinks();
    n++;

    const AnaEventInfoPD* info = static_cast<AnaEventInfoPD*>(spill->EventInfo);
    const bool isMC = info && info->IsMC;
    if(isMC)nMC++;
    if(info && info->EmptyEvent)emptyEvent++;
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
      if(!isMC){
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
        truePDGs[beamTrue->PDG]++;
        if(beamTrue->PositionEnd[2] < 0)beamEndsBeforeTPC++;   // la cara frontal del TPC está en z = 0
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
    if(info && info->HasPandora && !hasSlice)flagSliceButEmpty++;
    if(info && !info->HasPandora && hasSlice)flagNoSliceButParticles++;
    if(!hasSlice){
      noSlice++;
      if(isMC && beamTrue){
        noSliceWithTruth++;
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
      purMatched.Add(primary->TruePur);
      effMatched.Add(primary->TrueEff);
      continue;
    }
    mismatched++;
    purMismatched.Add(primary->TruePur);
    effMismatched.Add(primary->TrueEff);
    if(!primTrue)primNoTruth++;
    else if(DescendsFrom(primTrue->ID, beamTrue->ID, byID))primDescendant++;
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
  std::cout << std::endl << "Ficheros: " << files.size() << ", entradas leídas: " << n << " de " << nEntries
            << " (" << nMC << " de MC)" << std::endl << std::endl;

  std::cout << "Eventos" << std::endl;
  Line("DQ malo (GoodDaq = false, sin bunch)", Frac(badDQ, n));
  Line("sin bunch", Frac(noBunch, n));
  Line("con bunch, sin beam slice", Frac(noSlice, n));
  Line("con bunch y beam slice", Frac(withSlice, n));
  Line("sin PFParticles (EmptyEvent)", Frac(emptyEvent, n));
  if(badDQWithBunch || flagSliceButEmpty || flagNoSliceButParticles){
    std::cout << "  INCOHERENCIAS (no debería haber ninguna con minitrees del módulo):" << std::endl;
    Line("    DQ malo con bunch", Frac(badDQWithBunch, n));
    Line("    HasPandora = true sin partículas de la beam slice", Frac(flagSliceButEmpty, n));
    Line("    HasPandora = false con partículas de la beam slice", Frac(flagNoSliceButParticles, n));
  }
  std::cout << std::endl << "Instrumentación del haz" << std::endl;
  Line("con momento reconstruido", Frac(beamWithMom, n) + "   <p> = " + beamMom.Str("%.3f") + " GeV");
  Line("con traza del haz", Frac(beamWithTrack, n));
  if(!beamPDGs.empty()){
    std::cout << "  Hipótesis de PDG de la instrumentación (0 = ninguna; un evento puede tener varias):" << std::endl;
    PrintTop(beamPDGs, n, 10, "    ");
  }
  std::cout << std::endl << "Beam slice y cósmicos (eventos con beam slice)" << std::endl;
  Line("partículas de la beam slice por evento", nSlice.Str());
  Line("cósmicos guardados por evento", nCosmics.Str());
  Line("  de ellos con hits en el plano 2 (completos)", nCosmicsWithHits.Str());

  if(nMC > 0){
    std::cout << std::endl << "Truth (MC)" << std::endl;
    Line("sin partícula true del haz", Frac(noBeamTruth, nMC));
    Line("la partícula del haz acaba antes del TPC (z < 0)", Frac(beamEndsBeforeTPC, nMC));
    Line("sin beam slice, con truth del haz", Frac(noSliceWithTruth, nMC));
    Line("  de ellos con el haz acabando antes del TPC", Frac(noSliceEndsBeforeTPC, noSliceWithTruth));
    std::cout << "  PDG de la partícula true del haz:" << std::endl;
    PrintTop(truePDGs, nMC - noBeamTruth, 10, "    ");

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
