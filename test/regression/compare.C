// Compara dos ficheros de salida de highland árbol a árbol, rama a rama y
// entrada a entrada. Pensado para el test de regresión: dos ejecuciones con
// el mismo input y la misma configuración deberían dar ficheros idénticos.
//
// Uso:
//   root -l -b -q 'test/regression/compare.C+("ref/mc.root","new/mc.root")'
//
// Por defecto se salta el árbol "config", que guarda rutas (HIGHLANDPATH,
// fichero de entrada...) que cambian entre instalaciones aunque los
// resultados sean los mismos. Para incluirlo: compare.C+("a","b",false).
//
// Para cada hoja compara los bytes del valor (igualdad exacta) y, si es
// numérica, la diferencia absoluta máxima. Devuelve el número de hojas con
// diferencias (0 = ficheros equivalentes).

#include "TFile.h"
#include "TKey.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TObjArray.h"
#include "TMath.h"
#include <cstring>
#include <iostream>
#include <set>
#include <string>
#include <vector>

namespace {

struct LeafPair {
  TLeaf* a;
  TLeaf* b;
  Long64_t ndiff = 0;     // entradas con alguna diferencia
  double   maxdiff = 0;   // diferencia absoluta máxima (si es numérica)
  Long64_t firstdiff = -1;
};

int CompareTree(TTree* ta, TTree* tb){
  const std::string name = ta->GetName();
  int nbad = 0;

  if(ta->GetEntries() != tb->GetEntries()){
    std::cout << "  [" << name << "] número de entradas distinto: "
              << ta->GetEntries() << " vs " << tb->GetEntries() << std::endl;
    nbad++;
  }

  // Emparejar hojas por nombre completo
  std::vector<LeafPair> pairs;
  std::set<std::string> namesA;
  TObjArray* leavesA = ta->GetListOfLeaves();
  for(int i = 0; i < leavesA->GetEntriesFast(); i++){
    TLeaf* la = static_cast<TLeaf*>(leavesA->At(i));
    namesA.insert(la->GetFullName().Data());
    TLeaf* lb = tb->GetLeaf(la->GetFullName());
    if(!lb){
      std::cout << "  [" << name << "] hoja solo en el primero: " << la->GetFullName() << std::endl;
      nbad++;
      continue;
    }
    if(std::string(la->GetTypeName()) != lb->GetTypeName()){
      std::cout << "  [" << name << "] tipo distinto en " << la->GetFullName() << ": "
                << la->GetTypeName() << " vs " << lb->GetTypeName() << std::endl;
      nbad++;
      continue;
    }
    pairs.push_back({la, lb});
  }
  TObjArray* leavesB = tb->GetListOfLeaves();
  for(int i = 0; i < leavesB->GetEntriesFast(); i++){
    TLeaf* lb = static_cast<TLeaf*>(leavesB->At(i));
    if(!namesA.count(lb->GetFullName().Data())){
      std::cout << "  [" << name << "] hoja solo en el segundo: " << lb->GetFullName() << std::endl;
      nbad++;
    }
  }

  const Long64_t n = TMath::Min(ta->GetEntries(), tb->GetEntries());
  for(Long64_t e = 0; e < n; e++){
    ta->GetEntry(e);
    tb->GetEntry(e);
    for(auto& p : pairs){
      const int lenA = p.a->GetLen();
      const int lenB = p.b->GetLen();
      bool differs = (lenA != lenB);
      if(!differs){
        const size_t nbytes = (size_t)lenA * p.a->GetLenType();
        void* va = p.a->GetValuePointer();
        void* vb = p.b->GetValuePointer();
        if(va && vb && nbytes > 0 && std::memcmp(va, vb, nbytes) != 0){
          differs = true;
          // Diferencia numérica máxima (no aplica a cadenas de caracteres)
          if(std::string(p.a->GetTypeName()) != "Char_t"){
            for(int k = 0; k < lenA; k++){
              double d = TMath::Abs(p.a->GetValue(k) - p.b->GetValue(k));
              if(d > p.maxdiff) p.maxdiff = d;
            }
          }
        }
      }
      if(differs){
        if(p.firstdiff < 0) p.firstdiff = e;
        p.ndiff++;
      }
    }
  }

  for(auto& p : pairs){
    if(p.ndiff == 0) continue;
    std::cout << "  [" << name << "] " << p.a->GetFullName() << ": " << p.ndiff << "/" << n
              << " entradas distintas, primera " << p.firstdiff
              << ", |dif| máx " << p.maxdiff << std::endl;
    nbad++;
  }
  if(nbad == 0)
    std::cout << "  [" << name << "] idéntico (" << n << " entradas, "
              << pairs.size() << " hojas)" << std::endl;
  return nbad;
}

} // namespace

int compare(const char* fileA, const char* fileB, bool skipConfig = true){
  TFile* fa = TFile::Open(fileA);
  TFile* fb = TFile::Open(fileB);
  if(!fa || fa->IsZombie() || !fb || fb->IsZombie()){
    std::cout << "No se pueden abrir los ficheros" << std::endl;
    return -1;
  }
  std::cout << "Comparando " << fileA << " con " << fileB << std::endl;

  int nbad = 0;
  std::set<std::string> seen;
  for(TObject* obj : *fa->GetListOfKeys()){
    TKey* key = static_cast<TKey*>(obj);
    if(!TClass::GetClass(key->GetClassName())->InheritsFrom(TTree::Class())) continue;
    if(!seen.insert(key->GetName()).second) continue;   // solo el ciclo más reciente
    if(skipConfig && std::string(key->GetName()) == "config"){
      std::cout << "  [config] no se compara (rutas de la instalación)" << std::endl;
      continue;
    }
    TTree* ta = static_cast<TTree*>(fa->Get(key->GetName()));
    TTree* tb = static_cast<TTree*>(fb->Get(key->GetName()));
    if(!tb){
      std::cout << "  [" << key->GetName() << "] árbol solo en el primero" << std::endl;
      nbad++;
      continue;
    }
    nbad += CompareTree(ta, tb);
  }
  for(TObject* obj : *fb->GetListOfKeys()){
    TKey* key = static_cast<TKey*>(obj);
    if(!TClass::GetClass(key->GetClassName())->InheritsFrom(TTree::Class())) continue;
    if(!seen.count(key->GetName()) && !(skipConfig && std::string(key->GetName()) == "config")){
      std::cout << "  [" << key->GetName() << "] árbol solo en el segundo" << std::endl;
      nbad++;
    }
  }

  std::cout << (nbad == 0 ? "RESULTADO: idénticos" : "RESULTADO: hay diferencias")
            << " (" << nbad << ")" << std::endl;
  return nbad;
}
