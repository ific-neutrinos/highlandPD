// Carga highlandPD, escribe el minitree sintético y lo comprueba.
// Uso, con el entorno de highland y highlandPD cargado:
//   root -l -b -q 'test/phase2/run_test_minitree.C("<minitree.root>","<salida.root>",200)'
int run_test_minitree(const char* in, const char* out, int n = 200){
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPDPATH")));
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPATH")));
  if(gSystem->Load("libhighlandPD") < 0) return -1;
  TString dir = gSystem->DirName(gInterpreter->GetCurrentMacroName());
  gROOT->ProcessLine(Form(".L %s/make_test_minitree.C+", dir.Data()));
  if(gROOT->ProcessLine(Form("make_test_minitree(\"%s\",\"%s\",%d)", in, out, n)) != 0) return 1;
  return gROOT->ProcessLine(Form("check_test_minitree(\"%s\")", out));
}
