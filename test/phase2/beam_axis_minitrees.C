// Eje del haz en el TPC (BeamVolumeStart{Data,MC} y BeamVolumeDirection{Data,MC} del fcl de
// CreateHighlandMiniTree, highlandAnalysis), medido sobre minitrees del módulo, por momento nominal:
//   - traza de la instrumentación del haz (Beam->BeamParticle, extremo final);
//   - primario de la beam slice: en MC cuando es la partícula true del haz; en MC y en datos, con
//     un corte parecido al de calidad (empieza a z < 100 cm y va en la dirección de la
//     instrumentación, cos > 0,9);
//   - en MC, la trayectoria true del haz al cruzar z = 0 (cara frontal del TPC).
// Para cada una: punto de cruce con z = 0 (x, y) y dirección (cuantiles 16/50/84 %), ángulo con
// la dirección del fcl y distancia del punto de inicio de los primarios al eje del fcl.
//
// La entrada es un .txt con un fichero por línea (rutas o URLs xrootd; '#' para comentarios) o un
// fichero o patrón de ROOT.
//
// Uso, con el entorno de highland y highlandPD cargado ($HIGHLANDPDPATH, $HIGHLANDPATH):
//   root -l -b -q 'test/phase2/beam_axis_minitrees.C("lista.txt")'
// Con maxEntries > 0 se para tras ese número de entradas.

int beam_axis_minitrees(const char* input, Long64_t maxEntries = -1){
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPDPATH")));
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPATH")));
  if(gSystem->Load("libhighlandPD") < 0){ std::cerr << "No se puede cargar libhighlandPD" << std::endl; return -1; }
  TString dir = gSystem->DirName(gInterpreter->GetCurrentMacroName());
  gROOT->ProcessLine(Form(".L %s/beam_axis_minitrees_impl.C+", dir.Data()));
  return gROOT->ProcessLine(Form("beam_axis_minitrees_impl(\"%s\",%lld)", input, maxEntries));
}
