// microtree: micro-tree de sistemáticos de entrada
void selfsysttest(const std::string& microtree){

  gStyle->SetOptStat(0);

  CoherentFit* cf_mc = new CoherentFit(microtree,true);
  cf_mc->CreateCoherentSamples(999);
  cf_mc->GenerateTrueMCHistograms(1,61,2,999,1,20);
  CoherentSample* sb_sample = cf_mc->GetSignalPlusBackgroundSample();
  CoherentSample* ts_sample = cf_mc->GetTrueSignalSample();
  CoherentSample* tb_sample = cf_mc->GetTrueBackgroundSample();

  ts_sample->SequentialCoherentFit();
  tb_sample->SequentialCoherentFit();
  sb_sample->GetInitialParValuesForCoherentFit();
  
  cf_mc->ComputeSelfSystematicError();
}

