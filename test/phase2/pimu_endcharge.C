// ¿Separa el exceso de carga en la parada (pi+ -> mu+ nu en reposo) un pi+ que se para de un mu? Ver pimu_endcharge_impl.C.
//
//   root -l -b -q 'test/phase2/pimu_endcharge.C("lista.txt")'


int pimu_endcharge(const char* input, Long64_t maxEntries = -1){
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPDPATH")));
  gSystem->AddIncludePath(Form("-I%s/include", gSystem->Getenv("HIGHLANDPATH")));
  if(gSystem->Load("libhighlandPD") < 0){ std::cerr << "No se puede cargar libhighlandPD" << std::endl; return -1; }
  TString dir = gSystem->DirName(gInterpreter->GetCurrentMacroName());
  gROOT->ProcessLine(Form(".L %s/pimu_endcharge_impl.C+", dir.Data()));
  return gROOT->ProcessLine(Form("pimu_endcharge_impl(\"%s\",%lld)", input, maxEntries));
}
