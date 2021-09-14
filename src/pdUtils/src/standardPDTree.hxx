#ifndef standardPDTree_h
#define standardPDTree_h

#include "OutputManager.hxx"
#include "baseAnalysis.hxx"
#include "pdDataClasses.hxx"

namespace standardPDTree{

  // Methods to add to the output tree the standard sets of variables
  void AddStandardVariables_CountersTrue(OutputManager& output);
  
  void AddStandardVariables_BeamReco(OutputManager& output);
  void AddStandardVariables_BeamTrue(OutputManager& output);

  void AddStandardVariables_CandidateReco(OutputManager& output);
  void AddStandardVariables_CandidateTrue(OutputManager& output);
  void AddStandardVariables_CandidateHitsReco(OutputManager& output);

  void AddStandardVariables_AllParticlesReco(OutputManager& output, UInt_t nmax);
  void AddStandardVariables_AllParticlesTrue(OutputManager& output, UInt_t nmax);

  void AddStandardVariables_CandidateDaughtersReco(OutputManager& output, UInt_t nmax);
  void AddStandardVariables_CandidateDaughtersTrue(OutputManager& output, UInt_t nmax);
  void AddStandardVariables_CandidateDaughtersHitsReco(OutputManager& output, UInt_t nmax, UInt_t nmaxhitsperplane = NMAXHITSPERPLANE);

  void AddStandardVariables_CosmicsReco(OutputManager& output, UInt_t nmax = NMAXCOSMICS);
  void AddStandardVariables_CosmicsTrue(OutputManager& output, UInt_t nmax = NMAXCOSMICS);

  // Methods to fill the standard sets of variables in the output tree
  void FillStandardVariables_CountersTrue(OutputManager& output, PDCounters& counters);
  
  void FillStandardVariables_BeamTrue(OutputManager& output, AnaBeamB* beamB);
  void FillStandardVariables_BeamReco(OutputManager& output, AnaBeamB* beamB);

  void FillStandardVariables_AllParticlesReco(OutputManager& output, AnaParticlePD* part);
  void FillStandardVariables_AllParticlesTrue(OutputManager& output, AnaParticlePD* part);

  void FillStandardVariables_CandidateTrue(OutputManager& output, AnaParticlePD* part);
  void FillStandardVariables_CandidateReco(OutputManager& output, AnaParticlePD* part);
  void FillStandardVariables_CandidateHitsReco(OutputManager& output, AnaParticlePD* part);

  void FillStandardVariables_CandidateDaughterReco(OutputManager& output, AnaParticlePD* part);
  void FillStandardVariables_CandidateDaughterTrue(OutputManager& output, AnaParticlePD* part);  
  void FillStandardVariables_CandidateDaughterHitsReco(OutputManager& output, AnaParticlePD* part, UInt_t nmaxsavedhits);

  void FillStandardVariables_CosmicsReco(OutputManager& output, AnaParticlePD* part);
  void FillStandardVariables_CosmicsTrue(OutputManager& output, AnaParticlePD* part);

  // Enum with unique indexes for output tree variables  
  enum enumStandardMicroTrees_standardPDTree{

    // selected track true info
    seltrk_truemom = baseAnalysis::enumStandardMicroTreesLast_baseAnalysis+1,
    seltrk_trueendmom,
    seltrk_truepdg,
    seltrk_trueproc,
    seltrk_trueendproc,
    seltrk_truendau,
    seltrk_truepos,
    seltrk_trueendpos,
    seltrk_truedir,
    seltrk_trueeff,
    seltrk_truepur,

    // selected track reco info
    seltrk_ndau,
    seltrk_pos,
    seltrk_endpos,
    seltrk_dir,
    seltrk_enddir,
    seltrk_costheta,
    seltrk_length,
    seltrk_mom_muon,
    seltrk_mom_prot,
    seltrk_nhits,
    seltrk_dqdx,
    seltrk_dedx,
    seltrk_dedx_raw,
    seltrk_CNNscore,
    seltrk_chi2_prot,
    seltrk_chi2_muon,
    seltrk_chi2_ndf,
    //seltrk_calo,
    //seltrk_nhitsperplane,
    //seltrk_mom_muon_alt,
    //seltrk_mom_prot_alt,
    //seltrk_length_alt,
    //seltrk_length_raw,
    //seltrk_dedx_exp,
    //seltrk_pid,

    seltrk_hit_dedx,
    seltrk_hit_dedx_raw,
    seltrk_hit_dedx_cal,
    seltrk_hit_dqdx,
    seltrk_hit_dqdx_raw,
    seltrk_hit_dqdx_noSCE,
    seltrk_hit_ch,
    seltrk_hit_t0,
    seltrk_hit_x,
    seltrk_hit_y,
    seltrk_hit_z,
    seltrk_hit_x_raw,
    seltrk_hit_y_raw,
    seltrk_hit_z_raw,
    seltrk_hit_resrange,
    seltrk_hit_cnn,

    // selected track daughters reco info
    seltrk_dau_ndau,
    seltrk_dau_pos,
    seltrk_dau_dir,
    seltrk_dau_endpos,
    seltrk_dau_enddir,
    seltrk_dau_length,
    seltrk_dau_mom_muon,
    seltrk_dau_mom_prot,
    seltrk_dau_CNNscore,
    seltrk_dau_chi2_prot,
    seltrk_dau_chi2_muon,
    seltrk_dau_chi2_ndf,

    seltrk_dau_nhits,
    seltrk_dau_hit_x,
    seltrk_dau_hit_y,
    seltrk_dau_hit_z,
    seltrk_dau_hit_dedx,
    seltrk_dau_hit_dedx_cal,
    seltrk_dau_hit_dqdx_raw,    
    seltrk_dau_hit_resrange,
    seltrk_dau_hit_cnn,
    seltrk_dau_hit_ch,
    //seltrk_dau_nhits2,

    // selected track daughters true info
    seltrk_dau_truendau,
    seltrk_dau_truepdg,
    seltrk_dau_truepos,
    seltrk_dau_trueendpos,
    seltrk_dau_trueproc,
    seltrk_dau_trueendproc,
    seltrk_dau_truemom,
    seltrk_dau_trueendmom,

    // all particles in the event, reco info
    ntracks,
    trk_ndau,
    trk_pos,
    trk_dir,
    trk_endpos,
    trk_enddir,
    trk_length,
    trk_mom_muon,
    trk_mom_prot,
    //trk_dedx,

    // all particles in the event, true info    
    trk_truendau,
    trk_truepdg,
    trk_truepos,
    trk_trueendpos,
    trk_trueproc,
    trk_trueendproc,
    trk_truemom,
    trk_trueendmom,

    // beam true info
    beam_truepos,
    beam_trueendpos,
    beam_truemom,
    beam_truemom_tpc,
    beam_truedir,
    beam_truepdg,
    beam_trueendproc,

    // beam reco info
    beam_endpos,          
    beam_enddir,          
    beam_mom,                                 
    beam_mom_raw,                                 
    beam_mom_tpc,                                 
    beam_mom_tpc_raw,                                 
    beam_trigger,      
    beam_tof,
    beam_time,                                 
    beam_ckov_status,  
    beam_ckov_time,    
    beam_ckov_pressure,
    beam_npdgs,
    beam_pdgs,

    //cosmic info
    ncosmics,
    cosmics_generation,
    cosmics_ndau,
    cosmics_pos,
    cosmics_dir,
    cosmics_endpos,
    cosmics_enddir,
    cosmics_length,
    cosmics_mom_muon,
    cosmics_mom_prot,
    cosmics_type,
    cosmics_CNNscore,
    cosmics_chi2_prot,
    cosmics_chi2_muon,
    cosmics_chi2_ndf,
    cosmics_averagedEdx,
    cosmics_nhits,

    //cosmics true info
    cosmics_truendau,
    cosmics_truegeneration,
    cosmics_truepdg,
    cosmics_trueorigin,
    cosmics_truepos,
    cosmics_trueendpos,
    cosmics_trueproc,
    cosmics_trueendproc,
    cosmics_truemom,
    cosmics_trueendmom,
    
    // global counters
    truebeamdau_npi0,
    truebeamdau_npiplus,
    truebeamdau_npiminus,
    truebeamdau_nproton,
    truebeamdau_nneutron,
    truebeamdau_nnucleus,
        
    enumStandardMicroTreesLast_standardPDTree
  };
}



#endif
