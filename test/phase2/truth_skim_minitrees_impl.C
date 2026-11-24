// Implementación de truth_skim_minitrees.C (se compila con ACLiC).
//
// Árbol "skim", una entrada por evento del minitree (MC y datos; en datos solo van la instrumentación y la reco).
// Todo en las unidades del minitree: GeV, GeV/c, cm, ns. Las ramas t_* son vectores paralelos (una posición por
// partícula true), las s_* vectores paralelos por PFParticle de la beam slice y las bh_* por hit del plano de
// colección del primario de la beam slice.
//
// Evento:
//   entry (índice en la cadena de entrada), file (índice del fichero en la lista; con run_skim_parallel.sh, los
//   dos son relativos al trozo: el evento lo identifican run, subrun y event), run, subrun, event, isMC,
//   fullDump, emptyEvent, hasPandora, missingReco, nominalP.
// Instrumentación del haz: bi_p (momento del espectrómetro), bi_tof, bi_ck0/bi_ck1 (estado de los Cherenkov),
//   bi_trigger, bi_nPDG y bi_pdg (los PDG compatibles; vector).
// Haz que dispara (truth, beam->BeamParticle->TrueObject): beamIdx (su posición en t_*, -1 si no está), y la
//   trayectoria true guardada (bt_x, bt_y, bt_z, bt_p, bt_proc; vectores; proceso -999 si la producción no lo
//   guardó). bvtx_*: "vértice perfecto" de EventContext (último hit del plano 2 del haz true a lo largo del eje,
//   en el sistema de la reco; hasTVtx = 0 sin él).
// Partículas true (todas las no cósmicas; con keepCosmics también las de Origin 2):
//   t_id, t_parent (ParentID), t_pdg, t_parentPDG (el módulo lo rellena aunque la madre no se guarde), t_origin,
//   t_gen (Generation; -999 las que solo están por estar emparejadas con un PFParticle), t_procStart, t_procEnd
//   (AnaTrueParticleB::ProcessEnum), t_p, t_pEnd, t_x/t_y/t_z/t_t (inicio; t en ns), t_xe/t_ye/t_ze/t_te (final),
//   t_dx/t_dy/t_dz (dirección inicial), t_dxe/t_dye/t_dze (final), t_nDau (hijas en la lista de Geant4, guardadas
//   o no), t_nHits2 (hits del plano de colección donde es la que más deposita), t_eDep (energía en los hits del
//   plano de colección), t_matched (Matched), t_sigMu (IsSignalMuon), t_kOrigin (KOrigin; -1 si no es mu+ de K+).
//   Reco emparejada (el PFParticle con más hits cuyo TrueObject es ella): t_recoN (cuántos PFParticles la tienen
//   como TrueObject), t_recoUID (-1 sin ninguno), t_recoType (0 desconocido, 1 shower, 2 traza), t_recoInSlice,
//   t_recoNHits, t_recoLength.
// Beam slice (PFParticles con BeamOrigin, en el orden del minitree; el primero es el primario): nSlice, s_uid,
//   s_parent (UniqueID de la madre; -1 el primario), s_gen, s_type, s_isBeam (isBeamPart), s_nHits, s_length,
//   s_x/s_y/s_z, s_xe/s_ye/s_ze (inicio y final, sin SCE), s_chi2p, s_chi2mu, s_ndf, s_trkScore, s_nDau, s_trueID,
//   s_truePDG, s_trueEff, s_truePur.
// Hits del plano de colección del primario de la beam slice (s_* posición 0), en el orden guardado: bh_dEdx,
//   bh_rr (rango residual), bh_x/bh_y/bh_z (posición 3D, SCE corregida), bh_trueID.
#include "kaon_fate_minitrees_impl.C"
#include "muCandidateUtils.hxx"
#include "TFile.h"
#include "TTree.h"
#include "Compression.h"
#include <iostream>
#include <map>
#include <vector>

namespace tskim {
  using kfate::InputFiles; using kfate::TreeName;

  struct Row {
    // evento
    Long64_t entry = 0; Int_t file = 0, run = 0, subrun = 0, event = 0, isMC = 0, fullDump = 0, emptyEvent = 0, hasPandora = 0, missingReco = 0;
    Float_t nominalP = 0;
    // instrumentación del haz
    Float_t bi_p = 0, bi_tof = 0; Int_t bi_ck0 = 0, bi_ck1 = 0, bi_trigger = 0, bi_nPDG = 0; std::vector<Int_t> bi_pdg;
    // haz true
    Int_t beamIdx = -1; std::vector<Float_t> bt_x, bt_y, bt_z, bt_p; std::vector<Int_t> bt_proc;
    Int_t hasTVtx = 0; Float_t bvtx_x = -999, bvtx_y = -999, bvtx_z = -999;
    // partículas true
    Int_t nTrue = 0;
    std::vector<Int_t> t_id, t_parent, t_pdg, t_parentPDG, t_origin, t_gen, t_procStart, t_procEnd, t_nDau, t_nHits2, t_matched, t_sigMu, t_kOrigin;
    std::vector<Float_t> t_p, t_pEnd, t_x, t_y, t_z, t_t, t_xe, t_ye, t_ze, t_te, t_dx, t_dy, t_dz, t_dxe, t_dye, t_dze, t_eDep;
    std::vector<Int_t> t_recoN, t_recoUID, t_recoType, t_recoInSlice, t_recoNHits; std::vector<Float_t> t_recoLength;
    // beam slice
    Int_t nSlice = 0;
    std::vector<Int_t> s_uid, s_parent, s_gen, s_type, s_isBeam, s_nHits, s_nDau, s_trueID, s_truePDG;
    std::vector<Float_t> s_length, s_x, s_y, s_z, s_xe, s_ye, s_ze, s_chi2p, s_chi2mu, s_ndf, s_trkScore, s_trueEff, s_truePur;
    // hits del primario
    std::vector<Float_t> bh_dEdx, bh_rr, bh_x, bh_y, bh_z; std::vector<Int_t> bh_trueID;

    void Clear(){ *this = Row(); }
  };

  struct Skim {
    TFile* fout = nullptr; TTree* tr = nullptr; Row r;
    bool keepCosmics = false;
    long nEvents = 0, nMC = 0, nNoBeam = 0, nTrueTot = 0, nSigMu = 0, nSliceEv = 0;
    std::map<int, long> nKOrigin;

    void Open(const char* out, bool cosmics){
      keepCosmics = cosmics;
      fout = TFile::Open(out, "RECREATE");
      fout->SetCompressionSettings(ROOT::CompressionSettings(ROOT::kZSTD, 5));
      tr = new TTree("skim", "Truth skim de los minitrees (una entrada por evento)");
#define SK_BR(x) tr->Branch(#x, &r.x);
      SK_BR(entry) SK_BR(file) SK_BR(run) SK_BR(subrun) SK_BR(event) SK_BR(isMC) SK_BR(fullDump) SK_BR(emptyEvent) SK_BR(hasPandora) SK_BR(missingReco) SK_BR(nominalP)
      SK_BR(bi_p) SK_BR(bi_tof) SK_BR(bi_ck0) SK_BR(bi_ck1) SK_BR(bi_trigger) SK_BR(bi_nPDG) SK_BR(bi_pdg)
      SK_BR(beamIdx) SK_BR(bt_x) SK_BR(bt_y) SK_BR(bt_z) SK_BR(bt_p) SK_BR(bt_proc) SK_BR(hasTVtx) SK_BR(bvtx_x) SK_BR(bvtx_y) SK_BR(bvtx_z)
      SK_BR(nTrue) SK_BR(t_id) SK_BR(t_parent) SK_BR(t_pdg) SK_BR(t_parentPDG) SK_BR(t_origin) SK_BR(t_gen) SK_BR(t_procStart) SK_BR(t_procEnd)
      SK_BR(t_p) SK_BR(t_pEnd) SK_BR(t_x) SK_BR(t_y) SK_BR(t_z) SK_BR(t_t) SK_BR(t_xe) SK_BR(t_ye) SK_BR(t_ze) SK_BR(t_te)
      SK_BR(t_dx) SK_BR(t_dy) SK_BR(t_dz) SK_BR(t_dxe) SK_BR(t_dye) SK_BR(t_dze) SK_BR(t_nDau) SK_BR(t_nHits2) SK_BR(t_eDep) SK_BR(t_matched)
      SK_BR(t_sigMu) SK_BR(t_kOrigin) SK_BR(t_recoN) SK_BR(t_recoUID) SK_BR(t_recoType) SK_BR(t_recoInSlice) SK_BR(t_recoNHits) SK_BR(t_recoLength)
      SK_BR(nSlice) SK_BR(s_uid) SK_BR(s_parent) SK_BR(s_gen) SK_BR(s_type) SK_BR(s_isBeam) SK_BR(s_nHits) SK_BR(s_nDau) SK_BR(s_length)
      SK_BR(s_x) SK_BR(s_y) SK_BR(s_z) SK_BR(s_xe) SK_BR(s_ye) SK_BR(s_ze) SK_BR(s_chi2p) SK_BR(s_chi2mu) SK_BR(s_ndf) SK_BR(s_trkScore)
      SK_BR(s_trueID) SK_BR(s_truePDG) SK_BR(s_trueEff) SK_BR(s_truePur)
      SK_BR(bh_dEdx) SK_BR(bh_rr) SK_BR(bh_x) SK_BR(bh_y) SK_BR(bh_z) SK_BR(bh_trueID)
#undef SK_BR
    }

    void Process(AnaSpillPD* spill, Long64_t entry, int file){
      nEvents++;
      r.Clear();
      r.entry = entry; r.file = file;
      const AnaEventInfoPD* info = static_cast<AnaEventInfoPD*>(spill->EventInfo);
      if(info){
        r.run = info->Run; r.subrun = info->SubRun; r.event = info->Event; r.isMC = info->IsMC;
        r.fullDump = info->FullDump; r.emptyEvent = info->EmptyEvent; r.hasPandora = info->HasPandora; r.missingReco = info->MissingReco;
        r.nominalP = info->NominalBeamMom;
      }
      if(r.isMC)nMC++;

      // instrumentación del haz
      const AnaBeamPD* beam = static_cast<AnaBeamPD*>(spill->Beam);
      if(beam){
        r.bi_p = beam->BeamMomentum; r.bi_tof = beam->TOF; r.bi_ck0 = beam->CerenkovStatus[0]; r.bi_ck1 = beam->CerenkovStatus[1];
        r.bi_trigger = beam->BeamTrigger; r.bi_pdg = beam->PDGs; r.bi_nPDG = r.bi_pdg.size();
      }
      const AnaTrueParticleB* beamTrue = (beam && beam->BeamParticle) ? static_cast<AnaTrueParticleB*>(beam->BeamParticle->TrueObject) : nullptr;

      // reco
      std::vector<AnaParticlePD*> parts;
      const std::vector<AnaHitPD>* unassoc = nullptr;
      static const std::vector<AnaHitPD> noHits;
      if(!spill->Bunches.empty()){
        AnaBunchPD* bunch = static_cast<AnaBunchPD*>(spill->Bunches[0]);
        for(AnaParticleB* pb : bunch->Particles)parts.push_back(static_cast<AnaParticlePD*>(pb));
        unassoc = &bunch->UnassociatedHits;
      }
      muCandUtils::EventContext ctx;
      muCandUtils::MakeContext(ctx, parts, unassoc ? *unassoc : noHits, spill->TrueParticles, beamTrue, r.isMC);
      if(beamTrue){ auto it = ctx.byID.find(beamTrue->ID); if(it != ctx.byID.end())beamTrue = it->second; }   // la de la lista
      else if(r.isMC)nNoBeam++;
      if(ctx.hasTVtx){ r.hasTVtx = 1; r.bvtx_x = ctx.tvtx.X(); r.bvtx_y = ctx.tvtx.Y(); r.bvtx_z = ctx.tvtx.Z(); }

      // reco emparejada con cada partícula true
      struct Match { int n = 0; const AnaParticlePD* best = nullptr; };
      std::map<int, Match> match;
      for(const AnaParticlePD* p : parts){
        if(!p->TrueObject)continue;
        Match& m = match[p->TrueObject->ID];
        m.n++;
        if(!m.best || p->NHits > m.best->NHits)m.best = p;
      }

      // partículas true
      for(AnaTrueParticleB* tb : spill->TrueParticles){
        const AnaTrueParticlePD* t = static_cast<AnaTrueParticlePD*>(tb);
        if(t->Origin == 2 && !keepCosmics)continue;
        if(beamTrue && t->ID == beamTrue->ID)r.beamIdx = r.t_id.size();
        r.t_id.push_back(t->ID); r.t_parent.push_back(t->ParentID); r.t_pdg.push_back(t->PDG); r.t_parentPDG.push_back(t->ParentPDG);
        r.t_origin.push_back(t->Origin); r.t_gen.push_back(t->Generation); r.t_procStart.push_back(t->ProcessStart); r.t_procEnd.push_back(t->ProcessEnd);
        r.t_p.push_back(t->Momentum); r.t_pEnd.push_back(t->MomentumEnd);
        r.t_x.push_back(t->Position[0]); r.t_y.push_back(t->Position[1]); r.t_z.push_back(t->Position[2]); r.t_t.push_back(t->Position[3]);
        r.t_xe.push_back(t->PositionEnd[0]); r.t_ye.push_back(t->PositionEnd[1]); r.t_ze.push_back(t->PositionEnd[2]); r.t_te.push_back(t->PositionEnd[3]);
        r.t_dx.push_back(t->Direction[0]); r.t_dy.push_back(t->Direction[1]); r.t_dz.push_back(t->Direction[2]);
        r.t_dxe.push_back(t->DirectionEnd[0]); r.t_dye.push_back(t->DirectionEnd[1]); r.t_dze.push_back(t->DirectionEnd[2]);
        r.t_nDau.push_back(t->Daughters.size()); r.t_nHits2.push_back(t->NHitsTrue[2]); r.t_eDep.push_back(t->EDepTrue); r.t_matched.push_back(t->Matched);
        const bool sig = muCandUtils::IsSignalMuon(t, ctx.byID);
        const int ko = muCandUtils::KOrigin(t, beamTrue, ctx.byID);
        r.t_sigMu.push_back(sig); r.t_kOrigin.push_back(ko);
        if(sig){ nSigMu++; nKOrigin[ko]++; }
        auto im = match.find(t->ID);
        const AnaParticlePD* b = im == match.end() ? nullptr : im->second.best;
        r.t_recoN.push_back(im == match.end() ? 0 : im->second.n);
        r.t_recoUID.push_back(b ? Int_t(b->UniqueID) : -1); r.t_recoType.push_back(b ? Int_t(b->Type) : -1);
        r.t_recoInSlice.push_back(b ? b->BeamOrigin : 0); r.t_recoNHits.push_back(b ? b->NHits : 0); r.t_recoLength.push_back(b ? b->Length : -1);
      }
      r.nTrue = r.t_id.size(); nTrueTot += r.nTrue;

      // trayectoria true del haz
      if(beamTrue){
        for(const AnaTrueTrajPointPD& pt : static_cast<const AnaTrueParticlePD*>(beamTrue)->TrjPoints){
          r.bt_x.push_back(pt.Position[0]); r.bt_y.push_back(pt.Position[1]); r.bt_z.push_back(pt.Position[2]);
          r.bt_p.push_back(std::sqrt(pt.Momentum[0]*pt.Momentum[0] + pt.Momentum[1]*pt.Momentum[1] + pt.Momentum[2]*pt.Momentum[2]));
          r.bt_proc.push_back(pt.Process);
        }
      }

      // beam slice
      for(const AnaParticlePD* p : ctx.slice){
        r.s_uid.push_back(p->UniqueID); r.s_parent.push_back(p->ParentID); r.s_gen.push_back(p->Generation); r.s_type.push_back(p->Type);
        r.s_isBeam.push_back(p->isBeamPart); r.s_nHits.push_back(p->NHits); r.s_nDau.push_back(p->DaughtersIDs.size()); r.s_length.push_back(p->Length);
        r.s_x.push_back(p->PositionStart[0]); r.s_y.push_back(p->PositionStart[1]); r.s_z.push_back(p->PositionStart[2]);
        r.s_xe.push_back(p->PositionEnd[0]); r.s_ye.push_back(p->PositionEnd[1]); r.s_ze.push_back(p->PositionEnd[2]);
        r.s_chi2p.push_back(p->Chi2Proton); r.s_chi2mu.push_back(p->Chi2Muon); r.s_ndf.push_back(p->Chi2ndf); r.s_trkScore.push_back(p->PandoraTrackScore);
        r.s_trueID.push_back(p->TrueObject ? p->TrueObject->ID : 0);
        r.s_truePDG.push_back(p->TrueObject ? static_cast<const AnaTrueParticleB*>(p->TrueObject)->PDG : 0);
        r.s_trueEff.push_back(p->TrueEff); r.s_truePur.push_back(p->TruePur);
      }
      r.nSlice = r.s_uid.size();
      if(r.nSlice)nSliceEv++;

      // hits del plano de colección del primario
      if(!ctx.slice.empty()){
        for(const AnaHitPD& h : ctx.slice[0]->Hits[2]){
          r.bh_dEdx.push_back(h.dEdx); r.bh_rr.push_back(h.ResidualRange);
          r.bh_x.push_back(h.Position.X()); r.bh_y.push_back(h.Position.Y()); r.bh_z.push_back(h.Position.Z());
          r.bh_trueID.push_back(h.TrueID);
        }
      }
      tr->Fill();
    }

    void Print() const {
      std::cout << std::endl << "Entradas " << nEvents << " (MC " << nMC << ", MC sin truth del haz " << nNoBeam << ", con beam slice " << nSliceEv << ")" << std::endl
                << Form("Partículas true guardadas: %.1f por evento%s", nEvents ? double(nTrueTot)/nEvents : 0., keepCosmics ? " (con cósmicos)" : " (sin cósmicos)") << std::endl
                << "Signal muons (IsSignalMuon): " << nSigMu << "; por KOrigin:";
      for(auto& kv : nKOrigin)std::cout << "  " << kv.first << ": " << kv.second;
      std::cout << std::endl;
    }

    void Close(){ if(fout){ fout->cd(); tr->Write(); WriteCodeVersion(fout, "truth_skim_minitrees"); fout->Close(); } }
  };
}

int truth_skim_minitrees_impl(const char* input, const char* output, Long64_t maxEntries, int keepCosmics){
  using namespace tskim;
  const std::vector<std::string> files = InputFiles(input);
  if(files.empty()){ std::cerr << "No hay ficheros en " << input << std::endl; return -1; }
  const std::string treeName = TreeName(TString(files[0]).Contains("*") ? std::string("") : files[0]);
  TChain chain(treeName.empty() ? "highlandana/MiniTree" : treeName.c_str());
  for(const std::string& f : files)chain.Add(f.c_str());
  if(chain.LoadTree(0) < 0){ std::cerr << "No hay entradas (¿árbol " << chain.GetName() << "?)" << std::endl; return -1; }
  AnaSpillPD* spill = nullptr;
  chain.SetBranchAddress("Spill", &spill);
  Skim s;
  s.Open(output, keepCosmics);
  LoopMiniTree(chain, spill, maxEntries, files.size(), [&s, &chain](AnaSpillPD* sp){ s.Process(sp, chain.GetReadEntry(), chain.GetTreeNumber()); });
  std::cout << std::endl << "Ficheros: " << files.size() << std::endl;
  s.Print();
  s.Close();
  std::cout << "Árbol skim en " << output << std::endl;
  return 0;
}
