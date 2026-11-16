{

  // Directory with the micro-trees, from the environment variable MICROTREEDIR
  const char* envdir = gSystem->Getenv("MICROTREEDIR");
  if(!envdir){ std::cerr << "load.C: set MICROTREEDIR to the micro-tree directory" << std::endl; return; }
  std::string dir = std::string(envdir)+"/";

  std::string fdata = dir+"pionana_data_new.root";
  std::string fmc   = dir+"6GeV_prod4_2.root";
  std::string fmcs  = dir+"pionana_mc_syst.root";

  //TFile *_file0 = TFile::Open(fdata.c_str());
  TFile *_file1 = TFile::Open(fmc.c_str());
  //TFile *_file2 = TFile::Open(fmcs.c_str());

  //TTree* d  = (TTree*)_file0->Get("ana");
  TTree* mc = (TTree*)_file1->Get("ana");
  //TTree* mct= (TTree*)_file1->Get("truth");
  //TTree* mcs= (TTree*)_file2->Get("ana");

  DrawingTools draw(fmc.c_str());
}
