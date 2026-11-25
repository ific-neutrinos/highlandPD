// Contenido de los minitrees de una producción nueva: hits de la traza forzada y tiempo de las partículas true.
//   root -l -b -q 'test/phase2/check_production_minitree.C("mc_xxx.txt")'   (lista con URLs xrootd, receta 4)
//   root -l -b -q 'test/phase2/check_production_minitree.C("dir/minitree_*.root")'   (patrón local)
int check_production_minitree(const char* pattern, double lmin = 7.52, double lmax = 219.10){
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPDPATH")));
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPATH")));
  if(gSystem->Load("libhighlandPD") < 0){ std::cerr << "No se puede cargar libhighlandPD" << std::endl; return -1; }
  TString dir = gSystem->DirName(gInterpreter->GetCurrentMacroName());
  gROOT->ProcessLine(Form(".L %s/check_production_minitree_impl.C+", dir.Data()));
  return gROOT->ProcessLine(Form("check_production_minitree_impl(\"%s\",%g,%g)", pattern, lmin, lmax));
}
