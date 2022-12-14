void selfsysttest(){

  gStyle->SetOptStat(0);
  
  CoherentFit* cf_mc = new CoherentFit("../../../../test.root");
  cf_mc->CreateCoherentSamples(84);
  
  cf_mc->GenerateTrueMCHistograms(1,31,2,84,2);
  cf_mc->SetBackgroundModel(CoherentSample::BackgroundModelEnum::kQuadraticWidths);
  cf_mc->SequentialCoherentFit();
  cf_mc->InitializeHistogramsForSystematicErrors();
  cf_mc->ComputeSelfSystematicError(false,false,false,false,false,false,false,true);
}

