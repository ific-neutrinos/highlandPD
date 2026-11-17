// Compara dos minitrees escritos por CreateHighlandMiniTree (highlandAnalysis), por
// ejemplo antes y después de un cambio que no debería alterar la salida.
//
// Cada spill de los dos ficheros se lee y se vuelve a serializar: si los bytes coinciden,
// el contenido es idéntico (todas las clases y miembros, sin listarlos). En las entradas que
// difieren dice qué cambia en las partículas, para localizarlo.
//
// Uso, con el entorno de highland y highlandPD cargado ($HIGHLANDPDPATH, $HIGHLANDPATH):
//   root -l -b -q 'test/phase2/compare_module_minitrees.C("antes.root","despues.root")'
// Devuelve 0 si son idénticos.

int compare_module_minitrees(const char* fileA, const char* fileB, int maxPrint = 10){
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPDPATH")));
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPATH")));
  if(gSystem->Load("libhighlandPD") < 0){ std::cerr << "No se puede cargar libhighlandPD" << std::endl; return -1; }
  TString dir = gSystem->DirName(gInterpreter->GetCurrentMacroName());
  gROOT->ProcessLine(Form(".L %s/compare_module_minitrees_impl.C+", dir.Data()));
  return gROOT->ProcessLine(Form("compare_module_minitrees_impl(\"%s\",\"%s\",%d)", fileA, fileB, maxPrint));
}
