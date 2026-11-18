// Qué le pasa en la reconstrucción a cada K+ true de la cadena del haz y al mu+ de su
// desintegración, sobre minitrees de CreateHighlandMiniTree (MC): en la beam slice (primario,
// hija, generación >= 2), en un PFParticle fuera de ella (cósmico completo o en resumen,
// IsClearCosmic), dentro de otra partícula (por los hits), en un evento sin beam slice, sin
// reconstruir o sin hits en el TPC. Detalles de la clasificación en kaon_fate_minitrees_impl.C.
//
// Para decidir qué información de los cósmicos hace falta guardar (plan de la muestra grande de MC).
//
// La entrada es un .txt con un fichero por línea (rutas o URLs xrootd; '#' para comentarios) o un
// fichero o patrón de ROOT (p. ej. "dir/minitree_*.root").
//
// Uso, con el entorno de highland y highlandPD cargado ($HIGHLANDPDPATH, $HIGHLANDPATH):
//   root -l -b -q 'test/phase2/kaon_fate_minitrees.C("lista.txt")'
//   root -l -b -q 'test/phase2/kaon_fate_minitrees.C("lista.txt", -321)'    // K-
// Con maxEntries > 0 se para tras ese número de entradas.

int kaon_fate_minitrees(const char* input, int kaonPDG = 321, Long64_t maxEntries = -1){
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPDPATH")));
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPATH")));
  if(gSystem->Load("libhighlandPD") < 0){ std::cerr << "No se puede cargar libhighlandPD" << std::endl; return -1; }
  TString dir = gSystem->DirName(gInterpreter->GetCurrentMacroName());
  gROOT->ProcessLine(Form(".L %s/kaon_fate_minitrees_impl.C+", dir.Data()));
  return gROOT->ProcessLine(Form("kaon_fate_minitrees_impl(\"%s\",%d,%lld)", input, kaonPDG, maxEntries));
}
