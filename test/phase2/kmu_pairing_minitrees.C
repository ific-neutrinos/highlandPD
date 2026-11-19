// Pérdida de secondaryKaon por exigir que el mu sea la única hija de la K en Pandora, y si se
// recupera buscando el mu cerca del final de la K, sobre minitrees de CreateHighlandMiniTree (MC).
// Detalles en kmu_pairing_minitrees_impl.C.
//
// La entrada es un .txt con un fichero por línea o un fichero o patrón de ROOT.
// Uso, con el entorno de highland y highlandPD cargado ($HIGHLANDPDPATH, $HIGHLANDPATH):
//   root -l -b -q 'test/phase2/kmu_pairing_minitrees.C("dir/minitree_*.root")'
// which: 0 todos los eventos, 1 solo normales (para el fondo), 2 solo con volcado completo.

int kmu_pairing_minitrees(const char* input, Long64_t maxEntries = -1, int which = 0){
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPDPATH")));
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPATH")));
  if(gSystem->Load("libhighlandPD") < 0){ std::cerr << "No se puede cargar libhighlandPD" << std::endl; return -1; }
  TString dir = gSystem->DirName(gInterpreter->GetCurrentMacroName());
  gROOT->ProcessLine(Form(".L %s/kmu_pairing_minitrees_impl.C+", dir.Data()));
  return gROOT->ProcessLine(Form("kmu_pairing_minitrees_impl(\"%s\",%lld,%d)", input, maxEntries, which));
}
