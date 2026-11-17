// Desglose del tamaño de un minitree escrito por CreateHighlandMiniTree (highlandAnalysis):
// cuánto ocupa cada componente (partículas de la beam slice y cósmicos, con sus hits por plano,
// HitsNoCalo, puntos de trayectoria, space points y trazas forzadas; hits no asociados;
// partículas true y sus puntos de trayectoria; el resto del spill).
//
// Cada componente se serializa sin comprimir; el tamaño en disco se estima con el factor de
// compresión global del árbol (no es exacto por componente, pero sirve para ver qué pesa).
//
// Uso, con el entorno de highland y highlandPD cargado ($HIGHLANDPDPATH, $HIGHLANDPATH):
//   root -l -b -q 'test/phase2/size_breakdown_minitree.C("minitree.root")'

int size_breakdown_minitree(const char* file, Long64_t maxEntries = -1){
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPDPATH")));
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPATH")));
  if(gSystem->Load("libhighlandPD") < 0){ std::cerr << "No se puede cargar libhighlandPD" << std::endl; return -1; }
  TString dir = gSystem->DirName(gInterpreter->GetCurrentMacroName());
  gROOT->ProcessLine(Form(".L %s/size_breakdown_minitree_impl.C+", dir.Data()));
  return gROOT->ProcessLine(Form("size_breakdown_minitree_impl(\"%s\",%lld)", file, maxEntries));
}
