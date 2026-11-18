// Recuentos sobre minitrees escritos por CreateHighlandMiniTree (highlandAnalysis), para validar
// una muestra de la grid: eventos marcados (DQ malo, sin beam slice, sin truth del haz), la
// instrumentación del haz, el tamaño de la beam slice y de los cósmicos guardados y, en MC, si el
// primario de la beam slice es la partícula true del haz (y, si no, qué es).
//
// La entrada es un .txt con un fichero por línea (rutas o URLs xrootd; '#' para comentarios) o un
// fichero o patrón de ROOT (p. ej. "dir/minitree_*.root").
//
// Uso, con el entorno de highland y highlandPD cargado ($HIGHLANDPDPATH, $HIGHLANDPATH):
//   root -l -b -q 'test/phase2/count_module_minitrees.C("lista.txt")'
// Con maxEntries > 0 se para tras ese número de entradas.

int count_module_minitrees(const char* input, Long64_t maxEntries = -1){
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPDPATH")));
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPATH")));
  if(gSystem->Load("libhighlandPD") < 0){ std::cerr << "No se puede cargar libhighlandPD" << std::endl; return -1; }
  TString dir = gSystem->DirName(gInterpreter->GetCurrentMacroName());
  gROOT->ProcessLine(Form(".L %s/count_module_minitrees_impl.C+", dir.Data()));
  return gROOT->ProcessLine(Form("count_module_minitrees_impl(\"%s\",%lld)", input, maxEntries));
}
