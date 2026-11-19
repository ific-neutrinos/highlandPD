// Test de kmu_pairing_minitrees_impl.C con eventos construidos a mano.
// Uso, con el entorno de highland y highlandPD cargado:
//   root -l -b -q test/phase2/test_kmu_pairing.C
// Imprime "TEST OK" o los recuentos que no cuadran.

int test_kmu_pairing(){
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPDPATH")));
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPATH")));
  if(gSystem->Load("libhighlandPD") < 0) return -1;
  TString dir = gSystem->DirName(gInterpreter->GetCurrentMacroName());
  gROOT->ProcessLine(Form(".L %s/test_kmu_pairing_impl.C+", dir.Data()));
  return gROOT->ProcessLine("test_kmu_pairing_impl()");
}
