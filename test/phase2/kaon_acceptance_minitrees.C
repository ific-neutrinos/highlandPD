// Aceptancia física (solo truth) de la señal mu+ de K+ en reposo para la medida A: de las K+ producidas en la primera
// interacción inelástica del haz que dispara, en qué (T_K, theta_K) acaban dando un mu+ de K+ -> mu+ nu en reposo y
// contenido; destino del resto; K+ por interacción; K de la opción a (nivel, fracción de momento). Ver
// kaon_acceptance_minitrees_impl.C.
//
// Uso, con el entorno de highland y highlandPD cargado:
//   root -l -b -q 'test/phase2/kaon_acceptance_minitrees.C("lista.txt", "kaon_acceptance.root")'
// La entrada es un .txt con un fichero por línea (rutas o URLs xrootd) o un fichero o patrón de ROOT. Con
// maxEntries > 0 se para tras ese número de entradas.

int kaon_acceptance_minitrees(const char* input, const char* output = "kaon_acceptance.root", Long64_t maxEntries = -1){
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPDPATH")));
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPATH")));
  if(gSystem->Load("libhighlandPD") < 0){ std::cerr << "No se puede cargar libhighlandPD" << std::endl; return -1; }
  TString dir = gSystem->DirName(gInterpreter->GetCurrentMacroName());
  gROOT->ProcessLine(Form(".L %s/kaon_acceptance_minitrees_impl.C+", dir.Data()));
  return gROOT->ProcessLine(Form("kaon_acceptance_minitrees_impl(\"%s\",\"%s\",%lld)", input, output, maxEntries));
}
