// Test de beam_axis_minitrees_impl.C con eventos construidos a mano (geometría conocida).
// Uso, con el entorno de highland y highlandPD cargado:
//   root -l -b -q test/phase2/test_beam_axis.C
// Imprime "TEST OK" o lo que no cuadra (test_beam_axis_impl.C incluye la macro).

int test_beam_axis(){
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPDPATH")));
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPATH")));
  if(gSystem->Load("libhighlandPD") < 0) return -1;
  TString dir = gSystem->DirName(gInterpreter->GetCurrentMacroName());
  gROOT->ProcessLine(Form(".L %s/test_beam_axis_impl.C+", dir.Data()));
  return gROOT->ProcessLine("test_beam_axis_impl()");
}
