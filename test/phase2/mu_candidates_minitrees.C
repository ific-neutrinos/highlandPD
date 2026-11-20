// Candidatos a mu de K+ -> mu+ nu en reposo (p = 236 MeV/c, ~54 cm en LAr), sin reconstruir la K,
// sobre minitrees del módulo (MC): qué son de verdad las trazas que se parecen a ese mu y qué
// variables los separan. Base para la selección del análisis de producción
// (secondaryKaonProductionXS).
//
// Escribe un árbol plano "cand" (una entrada por traza candidata amplia: momento por rango de mu en
// [100, 600] MeV/c, para tener bandas laterales) con geometría, calorimetría, actividad en cada
// extremo, etiquetas de Pandora y truth, y otro "sig" con cada mu+ true de K+ de la cadena del
// haz y lo que le pasa en la reco. Imprime un resumen.
//
// La entrada es un .txt con un fichero por línea (rutas o URLs xrootd; '#' para comentarios) o un
// fichero o patrón de ROOT.
//
// Uso, con el entorno de highland y highlandPD cargado ($HIGHLANDPDPATH, $HIGHLANDPATH):
//   root -l -b -q 'test/phase2/mu_candidates_minitrees.C("lista.txt","mu_candidates.root")'
// Con maxEntries > 0 se para tras ese número de entradas.

int mu_candidates_minitrees(const char* input, const char* output = "mu_candidates.root", Long64_t maxEntries = -1){
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPDPATH")));
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPATH")));
  if(gSystem->Load("libhighlandPD") < 0){ std::cerr << "No se puede cargar libhighlandPD" << std::endl; return -1; }
  TString dir = gSystem->DirName(gInterpreter->GetCurrentMacroName());
  gROOT->ProcessLine(Form(".L %s/mu_candidates_minitrees_impl.C+", dir.Data()));
  return gROOT->ProcessLine(Form("mu_candidates_minitrees_impl(\"%s\",\"%s\",%lld)", input, output, maxEntries));
}
