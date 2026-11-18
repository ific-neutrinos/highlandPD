// Criterios candidatos para guardar completos los PFParticles de fuera de la beam slice (los que
// hoy van en resumen), medidos sobre minitrees del módulo (MC): (a) cilindro alrededor de la
// trayectoria de la beam slice, (b) inicio o fin cerca de esa trayectoria, (c) momento por rango
// (muón) cerca de 236 MeV/c, y combinaciones. Para cada uno, la eficiencia sobre los mu de K+ de la
// cadena del haz que están fuera de la beam slice y el coste (PFParticles que pasarían a
// completos por evento). Detalles en cosmic_criteria_minitrees_impl.C.
//
// La entrada es un .txt con un fichero por línea o un fichero o patrón de ROOT.
// Uso, con el entorno de highland y highlandPD cargado ($HIGHLANDPDPATH, $HIGHLANDPATH):
//   root -l -b -q 'test/phase2/cosmic_criteria_minitrees.C("dir/minitree_*.root")'

int cosmic_criteria_minitrees(const char* input, Long64_t maxEntries = -1){
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPDPATH")));
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPATH")));
  if(gSystem->Load("libhighlandPD") < 0){ std::cerr << "No se puede cargar libhighlandPD" << std::endl; return -1; }
  TString dir = gSystem->DirName(gInterpreter->GetCurrentMacroName());
  gROOT->ProcessLine(Form(".L %s/cosmic_criteria_minitrees_impl.C+", dir.Data()));
  return gROOT->ProcessLine(Form("cosmic_criteria_minitrees_impl(\"%s\",%lld)", input, maxEntries));
}
