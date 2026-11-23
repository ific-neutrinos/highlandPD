// Bucle común de las macros de test/phase2 sobre minitrees del módulo.
//
// No llama a TChain::GetEntries(): con miles de ficheros por xrootd abre todos antes de empezar
// (minutos u horas). Lee hasta que LoadTree dice que no hay más entradas. Cada "every" entradas
// escribe el avance por std::cerr (sin buffer: se ve en el log de nohup).
#ifndef MINITREE_LOOP_H
#define MINITREE_LOOP_H

#include "pdDataClasses.hxx"
#include "TChain.h"
#include "TNamed.h"
#include "TSystem.h"
#include "TDatime.h"
#include "TFile.h"
#include <chrono>
#include <iostream>

// Versión del código con el que se corre: commit de highlandPD ($HIGHLANDPDPATH), con "+cambios" si hay ficheros
// seguidos por git modificados sin commitear. Para saber con qué versión se hizo cada log y cada fichero de salida.
inline std::string CodeVersion(){
  const char* path = gSystem->Getenv("HIGHLANDPDPATH");
  if(!path)return "desconocida (sin HIGHLANDPDPATH)";
  TString hash = gSystem->GetFromPipe(Form("git -C %s rev-parse --short HEAD 2>/dev/null", path));
  if(hash.IsNull())return Form("desconocida (%s no es un repositorio git)", path);
  TString dirty = gSystem->GetFromPipe(Form("git -C %s status --porcelain --untracked-files=no 2>/dev/null | head -1", path));
  return std::string(hash.Data()) + (dirty.IsNull() ? "" : "+cambios");
}

// Guarda en el fichero de salida un TNamed "version" con el commit, la fecha y lo que se corrió
inline void WriteCodeVersion(TFile* f, const std::string& what){
  if(!f)return;
  TDirectory* old = gDirectory;
  f->cd();
  TNamed v("version", Form("highlandPD %s; %s; %s", CodeVersion().c_str(), TDatime().AsSQLString(), what.c_str()));
  v.Write();
  if(old)old->cd();
}

// process(spill) para cada entrada leída, ya con CopyVectorsIntoArrays y RedoLinks. Devuelve el
// número de entradas leídas.
template<class F>
Long64_t LoopMiniTree(TChain& chain, AnaSpillPD*& spill, Long64_t maxEntries, size_t nFiles, F process, Long64_t every = 1000){
  std::cerr << "[versión] highlandPD " << CodeVersion() << std::endl;
  const auto t0 = std::chrono::steady_clock::now();
  auto elapsed = [&t0](){ return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); };
  Long64_t n = 0;
  for(Long64_t i = 0; maxEntries <= 0 || i < maxEntries; i++){
    if(chain.LoadTree(i) < 0)break;   // no hay más entradas (o no se puede abrir el fichero siguiente)
    if(chain.GetEntry(i) <= 0 || !spill)continue;
    spill->CopyVectorsIntoArrays();
    spill->RedoLinks();
    process(spill);
    n++;
    if(every > 0 && n % every == 0)
      std::cerr << Form("[avance] %lld entradas, fichero %d de %zu, %.0f s", n, chain.GetTreeNumber()+1, nFiles, elapsed()) << std::endl;
  }
  std::cerr << Form("[avance] fin: %lld entradas en %.0f s", n, elapsed()) << std::endl;
  return n;
}

#endif
