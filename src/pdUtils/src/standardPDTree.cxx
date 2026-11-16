#include "standardPDTree.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
void standardPDTree::AddStandardVariables_EventInfo(OutputManager& output){
//********************************************************************
  
  AddVarI(output, beam_quality,        "beam quality cuts");
  AddVarF(output, beam_nominal_mom,    "beam nominal momentum");
}

//********************************************************************
void standardPDTree::AddStandardVariables_BeamInstrumentationTrue(OutputManager& output){
//********************************************************************

  AddVarI  (output, beam_truepdg    ,  "true beam track true pdg");
  AddVar4VF(output, beam_truepos    ,  "true beam track true start position");
  AddVar3VF(output, beam_truedir    ,  "true beam track true direction");
  AddVarF  (output, beam_truemom    ,  "true beam track true momentum");
  AddVar4VF(output, beam_trueendpos ,  "true beam track true end position");
  AddVarI  (output, beam_trueendproc,  "true beam track true end process");
}
  
//********************************************************************
void standardPDTree::AddStandardVariables_BeamInstrumentationReco(OutputManager& output){
//********************************************************************
  
  AddVar3VF(output, beam_endpos , "beam instrumentation track end position");
  AddVar3VF(output, beam_enddir , "beam instrumentation track end direction");
  AddVarF(  output, beam_mom    , "beam instrumentation track momentum");
  AddVarI(  output, beam_pdg    , "beam instrumentation pdgs");
  AddVarI(  output, beam_ntracks, "beam instrumentation ntracks");
}

//********************************************************************
void standardPDTree::AddStandardVariables_BeamParticleReco(OutputManager& output){
//********************************************************************

  AddVar4VF( output, seltrk_pos            , "beam particle position");
  AddVar3VF( output, seltrk_dir            , "beam particle direction");
  AddVar4VF( output, seltrk_endpos         , "beam particle end position");
  AddVar3VF( output, seltrk_enddir         , "beam particle end direction");
  AddVarI(   output, seltrk_nhits          , "beam particle number of hits");
  AddVarF(   output, seltrk_length         , "beam particle length");
  AddVarF(   output, seltrk_mom_muon       , "beam particle momentum muon");
  AddVarF(   output, seltrk_mom_prot       , "beam particle momentum proton");
  AddVarF(   output, seltrk_csdarange_muon , "beam particle momentum muon");
  AddVarF(   output, seltrk_csdarange_prot , "beam particle momentum proton");
  AddVar3VF( output, seltrk_CNNscore       , "beam particle reconstructed CNN score");
  AddVarF(   output, seltrk_chi2_prot      , "beam particle chi2 prot");
  AddVarF(   output, seltrk_chi2_kaon      , "beam particle chi2 kaon");
  AddVarF(   output, seltrk_chi2_muon      , "beam particle chi2 muon");
  AddVarF(   output, seltrk_chi2_ndf       , "beam particle chi2 ndf");
  AddVarF(   output, seltrk_calE           , "beam particle energy deposited along path");
  AddVarF(   output, seltrk_truncated_dedx , "beam particle energy deposited along path");
}

//********************************************************************
void standardPDTree::AddStandardVariables_BeamParticleHitsReco(OutputManager& output){
//********************************************************************

  AddVarFixMF(output, seltrk_hit_x       , "beam particle hit x"              , 3, NMAXHITSPERPLANE_SELTRK);
  AddVarFixMF(output, seltrk_hit_y       , "beam particle hit y"              , 3, NMAXHITSPERPLANE_SELTRK);
  AddVarFixMF(output, seltrk_hit_z       , "beam particle hit z"              , 3, NMAXHITSPERPLANE_SELTRK);
  AddVarFixMF(output, seltrk_hit_dedx    , "beam particle hit calibrated dEdx", 3, NMAXHITSPERPLANE_SELTRK);
  AddVarFixMF(output, seltrk_hit_dqdx    , "beam particle hit calibrated dQdx", 3, NMAXHITSPERPLANE_SELTRK);
  AddVarFixMF(output, seltrk_hit_resrange, "beam particle hit Residual Range" , 3, NMAXHITSPERPLANE_SELTRK);
}

//********************************************************************
void standardPDTree::AddStandardVariables_BeamParticleTrue(OutputManager& output){
//********************************************************************

  AddVar4VF(output, seltrk_truepos     , "beam particle true position");         
  AddVar3VF(output, seltrk_truedir     , "beam particle true direction");         
  AddVarF(  output, seltrk_truemom     , "beam particle true momentum");
  AddVarI(  output, seltrk_trueproc    , "beam particle true initial process");
  AddVarI(  output, seltrk_truendau    , "beam particle true number of daughters");
  AddVar4VF(output, seltrk_trueendpos  , "beam particle true end position");      
  AddVarF(  output, seltrk_trueendmom  , "beam particle true end momentum");
  AddVarI(  output, seltrk_trueendproc , "beam particle true final process");
  AddVarI(  output, seltrk_truepdg     , "beam particle true pdg code");
  AddVarF(  output, seltrk_trueeff     , "beam particle true-reco matching eff");
  AddVarF(  output, seltrk_truepur     , "beam particle true-reco matching pur");
  AddVarI(  output, seltrk_trueId      , "beam particle true ID object");
  AddVarI(  output, seltrk_true_matched, "beam particle matching");
}

//********************************************************************
void standardPDTree::AddStandardVariables_AllParticlesTrue(OutputManager& output, UInt_t nmax){
//********************************************************************

  AddVarMaxSizeVI (output, trk_truendau      , "all tracks true ndaughters"  , ntracks, nmax);
  AddVarMaxSizeVI (output, trk_truegeneration, "all tracks true generation"  , ntracks, nmax);
  AddVarMaxSizeVI (output, trk_truepdg       , "all tracks true pdg"         , ntracks, nmax);
  AddVarMaxSizeVI (output, trk_trueorigin    , "all tracks true origin"      , ntracks, nmax);
  AddVarMaxSize4MF(output, trk_truepos       , "all tracks true position"    , ntracks, nmax);
  AddVarMaxSize4MF(output, trk_trueendpos    , "all tracks true end position", ntracks, nmax);
  AddVarMaxSizeVI (output, trk_trueproc      , "all tracks true process"     , ntracks, nmax);
  AddVarMaxSizeVI (output, trk_trueendproc   , "all tracks true end process" , ntracks, nmax);
  AddVarMaxSizeVF (output, trk_truemom       , "all tracks true momentum"    , ntracks, nmax);
  AddVarMaxSizeVF (output, trk_trueendmom    , "all tracks true end momentum", ntracks, nmax);
}

//********************************************************************
void standardPDTree::AddStandardVariables_AllParticlesReco(OutputManager& output, UInt_t nmax){
//********************************************************************

  AddVarMaxSizeVI (output, trk_generation, "all tracks generation"       , ntracks, nmax);  
  AddVarMaxSizeVI (output, trk_ndau      , "all tracks #daughters"       , ntracks, nmax);
  AddVarMaxSize4MF(output, trk_pos       , "all tracks position"         , ntracks, nmax); 
  AddVarMaxSize3MF(output, trk_dir       , "all tracks direction"        , ntracks, nmax);
  AddVarMaxSize4MF(output, trk_endpos    , "all tracks position"         , ntracks, nmax); 
  AddVarMaxSize3MF(output, trk_enddir    , "all tracks direction"        , ntracks, nmax);
  AddVarMaxSizeVF (output, trk_length    , "all tracks length"           , ntracks, nmax);
  AddVarMaxSizeVF (output, trk_mom_muon  , "all tracks momentum (muon)"  , ntracks, nmax);
  AddVarMaxSizeVF (output, trk_mom_prot  , "all tracks momentum (proton)", ntracks, nmax);
  AddVarMaxSizeVI (output, trk_type      , "all tracks object type"      , ntracks, nmax);
  AddVarMaxSize3MF(output, trk_CNNscore  , "all tracks CNN score"        , ntracks, nmax);
  AddVarMaxSizeVF (output, trk_chi2_prot , "all tracks chi2 proton"      , ntracks, nmax);
  AddVarMaxSizeVF (output, trk_chi2_kaon , "all tracks chi2 kaon"        , ntracks, nmax);
  AddVarMaxSizeVF (output, trk_chi2_muon , "all tracks chi2 muon"        , ntracks, nmax);
  AddVarMaxSizeVF (output, trk_chi2_ndf  , "all tracks chi2 ndf"         , ntracks, nmax);
  AddVarMaxSizeVI (output, trk_nhits     , "all tracks #hits"            , ntracks, nmax);

}

//********************************************************************
void standardPDTree::AddStandardVariables_BeamParticleDaughtersTrue(OutputManager& output, UInt_t nmax){
//********************************************************************

  AddVarMaxSize4MF(output, seltrk_dau_truepos,     "beam particle daugthers true position",    seltrk_ndau,nmax);
  AddVarMaxSize4MF(output, seltrk_dau_trueendpos,  "beam particle daugthers true end position",seltrk_ndau,nmax);
  AddVarMaxSizeVI(output,  seltrk_dau_trueproc,    "beam particle daugthers true process",     seltrk_ndau,nmax);
  AddVarMaxSizeVI(output,  seltrk_dau_trueendproc, "beam particle daughters true end process", seltrk_ndau,nmax);
  AddVarMaxSizeVF(output,  seltrk_dau_truemom,     "beam particle daughters true momentum",    seltrk_ndau,nmax);
  AddVarMaxSizeVF(output,  seltrk_dau_trueendmom,  "beam particle daughters true end momentum",seltrk_ndau,nmax);
  AddVarMaxSizeVI(output,  seltrk_dau_truepdg,     "beam particle daughters true pdg",         seltrk_ndau,nmax);
  AddVarMaxSizeVI(output,  seltrk_dau_truendau,    "beam particle daughters' true ndaughters", seltrk_ndau,nmax);
}

//********************************************************************
void standardPDTree::AddStandardVariables_BeamParticleDaughtersReco(OutputManager& output, UInt_t nmax){
//********************************************************************
  
  AddVarMaxSize4MF(output, seltrk_dau_pos,            "beam particle daughters position",                 seltrk_ndau,nmax); 
  AddVarMaxSize3MF(output, seltrk_dau_dir,            "beam particle daughters direction",                seltrk_ndau,nmax);
  AddVarMaxSize4MF(output, seltrk_dau_endpos,         "beam particle daughters position",                 seltrk_ndau,nmax); 
  AddVarMaxSize3MF(output, seltrk_dau_enddir,         "beam particle daughters direction",                seltrk_ndau,nmax);
  AddVarMaxSizeVF( output, seltrk_dau_length,         "beam particle daughters length",                   seltrk_ndau,nmax);
  AddVarMaxSizeVF( output, seltrk_dau_mom_muon,       "beam particle daughters momentum (muon)",          seltrk_ndau,nmax);
  AddVarMaxSizeVF( output, seltrk_dau_mom_prot,       "beam particle daughters momentum (proton)",        seltrk_ndau,nmax);
  AddVarMaxSizeVI( output, seltrk_dau_nhits,          "beam particle daughters #hits",                    seltrk_ndau,nmax);
  AddVarMaxSizeVI( output, seltrk_dau_ndau,           "beam particle daughters' daughters"     ,          seltrk_ndau,nmax);
  AddVarMaxSize3MF(output, seltrk_dau_CNNscore,       "beam particle candidate daughters CNN score",      seltrk_ndau,nmax);
  AddVarMaxSizeVF( output, seltrk_dau_chi2_prot,      "beam particle candidate daughters chi2 proton",    seltrk_ndau,nmax);
  AddVarMaxSizeVF( output, seltrk_dau_chi2_kaon,      "beam particle candidate daughters chi2 kaon",      seltrk_ndau,nmax);
  AddVarMaxSizeVF( output, seltrk_dau_chi2_muon,      "beam particle candidate daughters chi2 proton",    seltrk_ndau,nmax);
  AddVarMaxSizeVF( output, seltrk_dau_chi2_ndf,       "beam particle candidate daughters chi2 ndf",       seltrk_ndau,nmax);
  AddVarMaxSizeVF( output, seltrk_dau_calE,           "beam particle candidate daughters deposited E",    seltrk_ndau,nmax);
  AddVarMaxSizeVF( output, seltrk_dau_truncated_dedx, "beam particle candidate daughters truncated dedx", seltrk_ndau,nmax);
}

//********************************************************************
void standardPDTree::AddStandardVariables_BeamParticleDaughtersHitsReco(OutputManager& output, UInt_t nmax, UInt_t nmaxhitsperplane){
//********************************************************************
  
  AddVarMF(output, seltrk_dau_hit_x,        "daughters x per hit",          seltrk_ndau,-nmax,nmaxhitsperplane);
  AddVarMF(output, seltrk_dau_hit_y,        "daughters y per hit",          seltrk_ndau,-nmax,nmaxhitsperplane);
  AddVarMF(output, seltrk_dau_hit_z,        "daughters z per hit",          seltrk_ndau,-nmax,nmaxhitsperplane);

  AddVarMF(output, seltrk_dau_hit_dedx,     "daughters dEdx per hit",       seltrk_ndau,-nmax,nmaxhitsperplane);
  AddVarMF(output, seltrk_dau_hit_resrange, "daughters hit residual range", seltrk_ndau,-nmax,nmaxhitsperplane);
}

//********************************************************************
void standardPDTree::FillStandardVariables_EventInfo(OutputManager& output, AnaEventPD* event){
//********************************************************************
  
  if (!event) return;
  bool quality = pdAnaUtils::GetBeamQuality(event);
  output.FillVar(beam_quality, quality);   

  AnaEventInfoPD* info = static_cast<AnaEventInfoPD*>(event->EventInfo);
  if (!info) return;
  output.FillVar(beam_nominal_mom, info->NominalBeamMom);   
}

//********************************************************************
void standardPDTree::FillStandardVariables_BeamInstrumentationTrue(OutputManager& output, AnaBeamB* beamB){
//********************************************************************

  if (!beamB) return;  
  AnaBeam* beam = static_cast<AnaBeam*>(beamB);
  if (!beam->BeamParticle) return;
  AnaTrueParticlePD* beamTruePart= static_cast<AnaTrueParticlePD*>(beam->BeamParticle->TrueObject);   
  if (!beamTruePart) return;

  output.FillVar(beam_truemom,                    beamTruePart->Momentum);
  output.FillVar(beam_truepdg,                    beamTruePart->PDG);
  output.FillVar(beam_trueendproc,                beamTruePart->ProcessEnd);
  output.FillVectorVarFromArray(beam_truedir,     beamTruePart->Direction,   3);
  output.FillVectorVarFromArray(beam_truepos,     beamTruePart->Position, 4);
  output.FillVectorVarFromArray(beam_trueendpos,  beamTruePart->PositionEnd, 4);
}

//********************************************************************
void standardPDTree::FillStandardVariables_BeamInstrumentationReco(OutputManager& output, AnaBeamB* beamB){
//********************************************************************

  if (!beamB) return;
  AnaBeamPD* beam         = static_cast<AnaBeamPD*>(beamB);

  AnaParticleMomB* beamPart = beam->BeamParticle;    
  if (!beamPart) return;
  
  output.FillVectorVarFromArray(beam_endpos, beamPart->PositionEnd,3);
  output.FillVectorVarFromArray(beam_enddir, beamPart->DirectionEnd,3);
  output.FillVar(               beam_mom   , beamPart->Momentum);   
  
  if((int)beam->PDGs.size()>0)
    output.FillVar(beam_pdg,beam->PDGs[0]);
  
  output.FillVar(beam_ntracks,(int)beam->nTracks);
}

//********************************************************************
void standardPDTree::FillStandardVariables_AllParticlesReco(OutputManager& output, AnaParticlePD* part){
//********************************************************************

  if (!part) return;
  output.FillVectorVar         (trk_generation,       part->Generation       );
  output.FillVectorVar         (trk_ndau,      (Int_t)part->Daughters.size() );
  output.FillMatrixVarFromArray(trk_pos,              part->PositionStart,  4);
  output.FillMatrixVarFromArray(trk_dir,              part->DirectionStart, 3); 
  output.FillMatrixVarFromArray(trk_endpos,           part->PositionEnd,    4);
  output.FillMatrixVarFromArray(trk_enddir,           part->DirectionEnd,   3); 
  output.FillVectorVar         (trk_length,           part->Length           );
  output.FillVectorVar         (trk_mom_prot,         pdAnaUtils::ComputeRangeMomentum(part->Length,2212));
  output.FillVectorVar         (trk_mom_muon,         pdAnaUtils::ComputeRangeMomentum(part->Length,13)  );
  output.FillVectorVar         (trk_type,             part->Type             );
  output.FillMatrixVarFromArray(trk_CNNscore,         part->CNNscore,       3); 
  output.FillVectorVar         (trk_chi2_prot,        part->Chi2Proton       );
  output.FillVectorVar         (trk_chi2_muon,        part->Chi2Muon         );
  output.FillVectorVar         (trk_chi2_ndf,         part->Chi2ndf          );
  output.FillVectorVar         (trk_nhits,            part->NHits            );
}
  
//********************************************************************
void standardPDTree::FillStandardVariables_AllParticlesTrue(OutputManager& output, AnaParticlePD* part){
//********************************************************************

  if (!part) return;
  if (!part->TrueObject) return;
 
  AnaTrueParticlePD* truePart = static_cast<AnaTrueParticlePD*>(part->TrueObject);
  if(!truePart) return;

  output.FillVectorVar         (trk_truendau,    (Int_t)truePart->Daughters.size());
  output.FillVectorVar         (trk_truegeneration,     truePart->Generation      );
  output.FillVectorVar         (trk_truepdg,            truePart->PDG             );
  output.FillVectorVar         (trk_trueorigin,         truePart->Origin          );
  output.FillMatrixVarFromArray(trk_truepos,            truePart->Position,      4);
  output.FillMatrixVarFromArray(trk_trueendpos,         truePart->PositionEnd,   4);
  output.FillVectorVar         (trk_trueproc,    (Int_t)truePart->ProcessStart    );
  output.FillVectorVar         (trk_trueendproc, (Int_t)truePart->ProcessEnd      );
  output.FillVectorVar         (trk_truemom,            truePart->Momentum        );
  output.FillVectorVar         (trk_trueendmom,         truePart->MomentumEnd     );
}

//********************************************************************
void standardPDTree::FillStandardVariables_BeamParticleTrue(OutputManager& output, AnaParticlePD* part){
//********************************************************************

  if (!part) return;
  AnaTrueParticle* truePart = static_cast<AnaTrueParticle*>(part->TrueObject);
  if (!truePart) return;

  output.FillVar               (seltrk_truepdg,             truePart->PDG);
  output.FillVectorVarFromArray(seltrk_truepos,             truePart->Position,    4);
  output.FillVectorVarFromArray(seltrk_trueendpos,          truePart->PositionEnd, 4);
  output.FillVectorVarFromArray(seltrk_truedir,             truePart->Direction,   3);
  output.FillVar               (seltrk_truemom,             truePart->Momentum);
  output.FillVar               (seltrk_trueendmom,          truePart->MomentumEnd);    
  output.FillVar               (seltrk_trueproc,     (Int_t)truePart->ProcessStart);
  output.FillVar               (seltrk_trueendproc,  (Int_t)truePart->ProcessEnd);    
  output.FillVar               (seltrk_truepur,             part->TruePur);
  output.FillVar               (seltrk_trueeff,             part->TrueEff);
  output.FillVar               (seltrk_truendau,     (Int_t)truePart->Daughters.size());
  output.FillVar               (seltrk_trueId,       (Int_t)truePart->ID);
  output.FillVar               (seltrk_true_matched, (Int_t)static_cast<AnaTrueParticlePD*>(truePart)->Matched);
}

//********************************************************************
void standardPDTree::FillStandardVariables_BeamParticleReco(OutputManager& output, AnaParticlePD* part, AnaParticlePD* beamPart){
//********************************************************************

  if (!part) return;

  output.FillVar(               seltrk_mom_prot      , (Float_t)part->RangeMomentum[0]);
  output.FillVar(               seltrk_mom_muon      , (Float_t)part->RangeMomentum[1]);
  output.FillVar(               seltrk_nhits         ,          part->NHits);
  output.FillVar(               seltrk_length        ,          part->Length);
  output.FillVectorVarFromArray(seltrk_pos           ,          part->PositionStart , 4);
  output.FillVectorVarFromArray(seltrk_endpos        ,          part->PositionEnd   , 4);
  output.FillVectorVarFromArray(seltrk_dir           ,          part->DirectionStart, 3);
  output.FillVectorVarFromArray(seltrk_enddir        ,          part->DirectionEnd  , 3);
  output.FillVectorVarFromArray(seltrk_CNNscore      ,          part->CNNscore      , 3);
  output.FillVar(               seltrk_chi2_prot     ,          part->Chi2Proton);
  output.FillVar(               seltrk_chi2_muon     ,          part->Chi2Muon);
  output.FillVar(               seltrk_chi2_ndf      ,          part->Chi2ndf);
  output.FillVar(               seltrk_chi2_kaon     , (Float_t)pdAnaUtils::Chi2PID(*part,321).first);
  output.FillVar(               seltrk_calE          , (Float_t)pdAnaUtils::ComputeDepositedEnergy(part));
  output.FillVar(               seltrk_truncated_dedx, (Float_t)pdAnaUtils::ComputeTruncatedMean(0.16,0.16,part->Hits[2]));

  if(!beamPart) return;
  output.FillVar(               seltrk_csdarange_prot,          pdAnaUtils::ComputeCSDARange(1000*(beamPart->Momentum), 2212));
  output.FillVar(               seltrk_csdarange_muon,          pdAnaUtils::ComputeCSDARange(1000*(beamPart->Momentum), 13));
}

//********************************************************************
void standardPDTree::FillStandardVariables_BeamParticleHitsReco(OutputManager& output, AnaParticlePD* part){
//********************************************************************

  if (!part) return;
  for (int i = 0; i < 3; i++){
    if(part->Hits[i].empty()){
      for (int j = 0; j < (int)NMAXHITSPERPLANE_SELTRK; j++){
        output.FillMatrixVar(seltrk_hit_x       ,(Float_t)-999., i, j);
        output.FillMatrixVar(seltrk_hit_y       ,(Float_t)-999., i, j);
        output.FillMatrixVar(seltrk_hit_z       ,(Float_t)-999., i, j);
        output.FillMatrixVar(seltrk_hit_dedx    ,(Float_t)-999., i, j);
        output.FillMatrixVar(seltrk_hit_dqdx    ,(Float_t)-999., i, j);
        output.FillMatrixVar(seltrk_hit_resrange,(Float_t)-999., i, j);
      }
    }
    else{
      int jmin = std::max<int>(0,(int)part->Hits[i].size()-NMAXHITSPERPLANE_SELTRK);
      int jmax = (int)part->Hits[i].size();
      for (int j = jmin; j < jmax; j++){
        output.FillMatrixVar(seltrk_hit_x    ,(Float_t)part->Hits[i][j].Position.X(), i, j-jmin);
        output.FillMatrixVar(seltrk_hit_y    ,(Float_t)part->Hits[i][j].Position.Y(), i, j-jmin);
        output.FillMatrixVar(seltrk_hit_z    ,(Float_t)part->Hits[i][j].Position.Z(), i, j-jmin);
        output.FillMatrixVar(seltrk_hit_dedx,       part->Hits[i][j].dEdx,              i, j-jmin);
        output.FillMatrixVar(seltrk_hit_dqdx,       part->Hits[i][j].dQdx,              i, j-jmin);
        output.FillMatrixVar(seltrk_hit_resrange,   part->Hits[i][j].ResidualRange,     i, j-jmin);
      }
    }
  }
}
  
//********************************************************************
void standardPDTree::FillStandardVariables_BeamParticleDaughtersReco(OutputManager& output, AnaParticlePD* dau){
//********************************************************************

  if (!dau) return;
  output.FillVectorVar         (seltrk_dau_mom_prot,         dau->RangeMomentum[0]);
  output.FillVectorVar         (seltrk_dau_mom_muon,         dau->RangeMomentum[1]);
  output.FillMatrixVarFromArray(seltrk_dau_pos,              dau->PositionStart, 4);
  output.FillMatrixVarFromArray(seltrk_dau_dir,              dau->DirectionStart,3); 
  output.FillMatrixVarFromArray(seltrk_dau_endpos,           dau->PositionEnd,   4);
  output.FillMatrixVarFromArray(seltrk_dau_enddir,           dau->DirectionEnd,  3); 
  output.FillVectorVar         (seltrk_dau_length,           dau->Length);
  output.FillVectorVar         (seltrk_dau_nhits,            dau->NHits           );
  output.FillVectorVar         (seltrk_dau_ndau,      (Int_t)dau->Daughters.size());
  output.FillMatrixVarFromArray(seltrk_dau_CNNscore,         dau->CNNscore,      3); 
  output.FillVectorVar         (seltrk_dau_chi2_prot,        dau->Chi2Proton      );
  output.FillVectorVar         (seltrk_dau_chi2_kaon,        (Float_t)pdAnaUtils::Chi2PID(*dau,321).first);
  output.FillVectorVar         (seltrk_dau_chi2_muon,        dau->Chi2Muon        );
  output.FillVectorVar         (seltrk_dau_chi2_ndf,         dau->Chi2ndf         );
  output.FillVectorVar         (seltrk_dau_truncated_dedx,   (Float_t)pdAnaUtils::ComputeTruncatedMean(0.16,0.16,dau->Hits[2]));
}

//********************************************************************
void standardPDTree::FillStandardVariables_BeamParticleDaughtersTrue(OutputManager& output, AnaParticlePD* dau){
//********************************************************************

  if (!dau) return;  
  AnaTrueParticle* dauTruePart = static_cast<AnaTrueParticle*>(dau->TrueObject);
  if(!dauTruePart) return;

  output.FillVectorVar(seltrk_dau_truepdg,             dauTruePart->PDG);
  output.FillVectorVar(seltrk_dau_trueproc,     (Int_t)dauTruePart->ProcessStart);
  output.FillVectorVar(seltrk_dau_trueendproc,  (Int_t)dauTruePart->ProcessEnd);
  output.FillMatrixVarFromArray(seltrk_dau_truepos,    dauTruePart->Position,4);
  output.FillMatrixVarFromArray(seltrk_dau_trueendpos, dauTruePart->PositionEnd,4);
  output.FillVectorVar(seltrk_dau_truemom,             dauTruePart->Momentum);
  output.FillVectorVar(seltrk_dau_trueendmom,          dauTruePart->MomentumEnd);
  output.FillVectorVar(seltrk_dau_truendau,     (Int_t)dauTruePart->Daughters.size());
}

//********************************************************************
void standardPDTree::FillStandardVariables_BeamParticleDaughtersHitsReco(OutputManager& output, AnaParticlePD* dau, UInt_t nmaxhitsperplane){
//********************************************************************

  if (!dau) return;
  if(dau->Hits[2].empty()){
    for (int j = 0; j < (int)NMAXHITSPERPLANE; j++){
      output.FillMatrixVar(seltrk_dau_hit_x,        (Float_t)-999., -1, j);
      output.FillMatrixVar(seltrk_dau_hit_y,        (Float_t)-999., -1, j);
      output.FillMatrixVar(seltrk_dau_hit_z,        (Float_t)-999., -1, j);
      output.FillMatrixVar(seltrk_dau_hit_dedx,     (Float_t)-999., -1, j);
      output.FillMatrixVar(seltrk_dau_hit_resrange, (Float_t)-999., -1, j);
    }
  }
  else{
    int jmin = std::max<int>(0,(int)dau->Hits[2].size()-nmaxhitsperplane);
    int jmax = (int)dau->Hits[2].size(); 

    for (int j = jmin; j < jmax; j++){
      output.FillMatrixVar(seltrk_dau_hit_x,         (Float_t)dau->Hits[2][j].Position.X(),   -1, j-jmin);
      output.FillMatrixVar(seltrk_dau_hit_y,         (Float_t)dau->Hits[2][j].Position.Y(),   -1, j-jmin);
      output.FillMatrixVar(seltrk_dau_hit_z,         (Float_t)dau->Hits[2][j].Position.Z(),   -1, j-jmin);
      output.FillMatrixVar(seltrk_dau_hit_dedx,      (Float_t)dau->Hits[2][j].dEdx,           -1, j-jmin);
      output.FillMatrixVar(seltrk_dau_hit_resrange,  (Float_t)dau->Hits[2][j].ResidualRange,  -1, j-jmin);
    }
  }
}
