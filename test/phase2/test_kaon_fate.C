// Test de kaon_fate_minitrees_impl.C con eventos construidos a mano (un caso por categoría).
// Uso, con el entorno de highland y highlandPD cargado:
//   root -l -b -q test/phase2/test_kaon_fate.C
// Imprime "TEST OK" o los recuentos que no cuadran (test_kaon_fate_impl.C incluye la macro).

int test_kaon_fate(){
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPDPATH")));
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPATH")));
  if(gSystem->Load("libhighlandPD") < 0) return -1;
  TString dir = gSystem->DirName(gInterpreter->GetCurrentMacroName());
  gROOT->ProcessLine(Form(".L %s/test_kaon_fate_impl.C+", dir.Data()));
  return gROOT->ProcessLine("test_kaon_fate_impl()");
}
