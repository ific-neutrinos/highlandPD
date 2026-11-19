// Test de CritCounter (count_module_minitrees_impl.C: cósmicos según CosmicCriteria) con eventos hechos a mano.
// Uso, con el entorno de highland y highlandPD cargado:
//   root -l -b -q test/phase2/test_count_module_crit.C
// Imprime "TEST OK" o los recuentos que no cuadran.

int test_count_module_crit(){
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPDPATH")));
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPATH")));
  if(gSystem->Load("libhighlandPD") < 0) return -1;
  TString dir = gSystem->DirName(gInterpreter->GetCurrentMacroName());
  gROOT->ProcessLine(Form(".L %s/test_count_module_crit_impl.C+", dir.Data()));
  return gROOT->ProcessLine("test_count_module_crit_impl()");
}
