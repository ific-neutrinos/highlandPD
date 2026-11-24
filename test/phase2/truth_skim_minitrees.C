// Truth skim: árbol plano "skim" (una entrada por evento) con toda la truth no cósmica del minitree, la truth del
// haz que dispara (con su trayectoria), la instrumentación del haz y la reconstrucción de la beam slice (con los
// hits del plano de colección del primario). Se hace una vez en FNAL y los estudios de truth (origen de la K,
// aceptancia, pileup, ventana del mu, pion del haz que sobrevive) se hacen en local sobre él. Ver
// truth_skim_minitrees_impl.C.
//
// Uso, con el entorno de highland y highlandPD cargado:
//   root -l -b -q 'test/phase2/truth_skim_minitrees.C("lista.txt", "truth_skim.root")'
// La entrada es un .txt con un fichero por línea (rutas o URLs xrootd) o un fichero o patrón de ROOT. Con
// maxEntries > 0 se para tras ese número de entradas; con keepCosmics = 1 guarda también las partículas true de
// los cósmicos (Origin 2).

int truth_skim_minitrees(const char* input, const char* output = "truth_skim.root", Long64_t maxEntries = -1, int keepCosmics = 0){
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPDPATH")));
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPATH")));
  if(gSystem->Load("libhighlandPD") < 0){ std::cerr << "No se puede cargar libhighlandPD" << std::endl; return -1; }
  TString dir = gSystem->DirName(gInterpreter->GetCurrentMacroName());
  gROOT->ProcessLine(Form(".L %s/truth_skim_minitrees_impl.C+", dir.Data()));
  return gROOT->ProcessLine(Form("truth_skim_minitrees_impl(\"%s\",\"%s\",%lld,%d)", input, output, maxEntries, keepCosmics));
}
