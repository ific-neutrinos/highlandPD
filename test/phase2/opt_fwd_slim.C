// Árboles preseleccionados por población para opt_fwd.C (26-nov-2026): la preselección suelta de opt_bdt.C, gapmin
// como rama y, fuera, el fondo submuestreado 1 de 3 (entry%3 == 0; opt_fwd lo pesa x3).
//
//   root -l -b -q 'test/phase2/opt_fwd_slim.C("opt_flat2_v100to105d01.root", "opt_fwd")'
void opt_fwd_slim(const char* in, const char* dir){
  gSystem->mkdir(dir, true);
  ROOT::RDataFrame df("flat", in);
  auto d = df.Define("gapmin", "std::min(std::min(std::abs(zstart-230.f),std::abs(zstart-462.f)),std::min(std::abs(zstop-230.f),std::abs(zstop-462.f)))");
  d.Filter("pop==0 && chi2mu<20 && faces>10").Snapshot("flat", Form("%s/slim_0.root", dir));
  d.Filter("pop==1 && axis<150 && chi2mu<20 && (cat<=1 || entry%3==0)").Snapshot("flat", Form("%s/slim_1.root", dir));
  for(int p = 0; p < 2; p++){ TFile f(Form("%s/slim_%d.root", dir, p)); printf("slim_%d: %lld candidatos\n", p, ((TTree*)f.Get("flat"))->GetEntries()); }
}
