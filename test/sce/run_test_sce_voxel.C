// Carga highlandPD y compila/ejecuta test_sce_voxel.C.
// Uso, con el entorno de highland y highlandPD cargado:
//   root -l -b -q test/sce/run_test_sce_voxel.C
int run_test_sce_voxel(int npoints = 5000){
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPDPATH")));
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPATH")));
  if(gSystem->Load("libhighlandPD") < 0) return -1;
  TString dir = gSystem->DirName(gInterpreter->GetCurrentMacroName());
  gROOT->ProcessLine(Form(".L %s/test_sce_voxel.C+", dir.Data()));
  return gROOT->ProcessLine(Form("test_sce_voxel(%d)", npoints));
}
