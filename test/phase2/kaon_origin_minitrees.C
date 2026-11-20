// De dónde vienen los mu+ de K+ en reposo en la truth de los minitrees del módulo (MC): K+ de la primera
// interacción de la partícula del haz, de interacciones secundarias, de otras partículas del haz (pileup) o
// de cadenas que el minitree no guarda enteras. Sirve para definir la señal del análisis de producción
// (sigma de la primera interacción) y saber qué K+ son fondo. Ver kaon_origin_minitrees_impl.C.
//
// La entrada es un .txt con un fichero por línea (rutas o URLs xrootd; '#' para comentarios) o un fichero o
// patrón de ROOT. Con el entorno de highland y highlandPD cargado:
//   root -l -b -q 'test/phase2/kaon_origin_minitrees.C("lista.txt")'

int kaon_origin_minitrees(const char* input, Long64_t maxEntries = -1){
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPDPATH")));
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPATH")));
  if(gSystem->Load("libhighlandPD") < 0){ std::cerr << "No se puede cargar libhighlandPD" << std::endl; return -1; }
  TString dir = gSystem->DirName(gInterpreter->GetCurrentMacroName());
  gROOT->ProcessLine(Form(".L %s/kaon_origin_minitrees_impl.C+", dir.Data()));
  return gROOT->ProcessLine(Form("kaon_origin_minitrees_impl(\"%s\",%lld)", input, maxEntries));
}
