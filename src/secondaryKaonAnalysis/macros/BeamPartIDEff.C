void BeamPartIDEff(){

  //list file
  std::string list = "efflist.txt";
  FILE* lfile = fopen(list.c_str(),"r");
  if(!lfile)std::exit(1);

  //efficiency variables
  int PLikeBeamEvents               = 0;
  int PiLikeBeamEvents              = 0;
  int KLikeBeamEvents               = 0;
  int PLikeBeamEventsReconstructed  = 0;
  int PiLikeBeamEventsReconstructed = 0;
  int KLikeBeamEventsReconstructed  = 0;

  //tree variables
  int tPLikeBeamEvents               = 0;
  int tPiLikeBeamEvents              = 0;
  int tKLikeBeamEvents               = 0;
  int tPLikeBeamEventsReconstructed  = 0;
  int tPiLikeBeamEventsReconstructed = 0;
  int tKLikeBeamEventsReconstructed  = 0;

  //loop over list
  char filename[500];
  TFile* tfile;
  TTree* tree;
  while(fscanf(lfile,"%s",filename) == 1){
    //open root file and get tree
    std::string sfilename(filename);
    std::cout << "running over file " << sfilename << std::endl;
    tfile = TFile::Open(sfilename.c_str());
    tree = (TTree*)tfile->Get("highlandana/EffTree");
    
    //set tree variables
    tree->SetBranchAddress("PLikeBeamEvents" ,              &tPLikeBeamEvents);
    tree->SetBranchAddress("PiLikeBeamEvents",              &tPiLikeBeamEvents);
    tree->SetBranchAddress("KLikeBeamEvents" ,              &tKLikeBeamEvents);
    tree->SetBranchAddress("PLikeBeamEventsReconstructed" , &tPLikeBeamEventsReconstructed);
    tree->SetBranchAddress("PiLikeBeamEventsReconstructed", &tPiLikeBeamEventsReconstructed);
    tree->SetBranchAddress("KLikeBeamEventsReconstructed" , &tKLikeBeamEventsReconstructed);

    //get values and update efficiencies
    tree->GetEntry(0);
    PLikeBeamEvents               += tPLikeBeamEvents;             
    PiLikeBeamEvents              += tPiLikeBeamEvents;            
    KLikeBeamEvents               += tKLikeBeamEvents;             
    PLikeBeamEventsReconstructed  += tPLikeBeamEventsReconstructed;
    PiLikeBeamEventsReconstructed += tPiLikeBeamEventsReconstructed;
    KLikeBeamEventsReconstructed  += tKLikeBeamEventsReconstructed;

    //close root file
    tfile->Close();
  }

  std::cout << "Proton Eff " << (double)PLikeBeamEventsReconstructed/PLikeBeamEvents << std::endl; 
  std::cout << "Pion Eff   " << (double)PiLikeBeamEventsReconstructed/PiLikeBeamEvents << std::endl; 
  std::cout << "Kaon Eff   " << (double)KLikeBeamEventsReconstructed/KLikeBeamEvents << std::endl; 

}
