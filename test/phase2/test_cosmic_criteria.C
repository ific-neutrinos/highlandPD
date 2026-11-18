// Test de cosmic_criteria_minitrees_impl.C: geometría y un evento hecho a mano.
// Uso, con el entorno de highland y highlandPD cargado:
//   root -l -b -q test/phase2/test_cosmic_criteria.C
int test_cosmic_criteria(){
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPDPATH")));
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPATH")));
  if(gSystem->Load("libhighlandPD") < 0) return -1;
  TString dir = gSystem->DirName(gInterpreter->GetCurrentMacroName());
  gROOT->ProcessLine(Form(".L %s/test_cosmic_criteria_impl.C+", dir.Data()));
  return gROOT->ProcessLine("test_cosmic_criteria_impl()");
}
