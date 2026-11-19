// Minitree sintético para probar la lectura de los minitrees de la Fase 2
// (phase2.md) sin pasar por LArSoft.
//
// Copia las primeras N entradas de un minitree real y:
//   - rellena los miembros nuevos (SliceID, T0, SpacePoints, UniqueID/TrueID de
//     los hits, TrjPoints de la truth, UnassociatedHits) con valores conocidos;
//   - cada 10 entradas añade un evento marcado, como los guarda ahora
//     CreateHighlandMiniTree: con Beam, y sin beam slice (bunch vacío, HasPandora =
//     false) o con DQ malo (GoodDaq = false, sin bunch).
// Luego relee el fichero y comprueba que los valores sobreviven al Clone().
//
// Uso, con el entorno de highland y highlandPD cargado:
//   root -l -b -q 'test/phase2/run_test_minitree.C("<minitree.root>","<salida.root>",200)'
#include "pdDataClasses.hxx"
#include "Header.hxx"
#include "TFile.h"
#include "TTree.h"
#include <iostream>

namespace {
  // Lee una entrada como MiniTreeConverter::GetSpill: AnaSpill nuevo (la rama puede ser
  // AnaSpill o AnaSpillPD: ROOT acepta un puntero a la clase base), vectores a arrays y enlaces
  AnaSpill* ReadSpill(TTree* t, Long64_t i){
    AnaSpill* spill = new AnaSpill();
    t->SetBranchAddress("Spill", &spill);
    t->GetEntry(i);
    spill->CopyVectorsIntoArrays();
    spill->RedoLinks();
    return spill;
  }

  // valores conocidos para comprobar al releer
  const Int_t   kSlice = 7;
  const Float_t kT0    = 123.5;
  const Int_t   kNSP   = 3;
  const Int_t   kNUnas = 5;

  // AnaSpillPD con el contenido de un clon de AnaSpill (el módulo escribe AnaSpillPD); vacía el
  // origen. Como en highland, el clon no es dueño de la truth (la comparte con el spill leído)
  AnaSpillPD* ToSpillPD(AnaSpill* from){
    AnaSpillPD* spill = new AnaSpillPD();
    spill->isClone = true;
    spill->EventInfo     = from->EventInfo;     from->EventInfo = NULL;
    spill->Beam          = from->Beam;          from->Beam = NULL;
    spill->DataQuality   = from->DataQuality;   from->DataQuality = NULL;
    spill->Bunches       = from->Bunches;       from->Bunches.clear();
    spill->TrueParticles = from->TrueParticles; from->TrueParticles.clear();
    spill->TrueVertices  = from->TrueVertices;  from->TrueVertices.clear();
    spill->RedoLinks();
    return spill;
  }

  void FillNew(AnaSpill* spill){
    for(AnaTrueParticleB* tb : spill->TrueParticles){
      AnaTrueParticlePD* tp = static_cast<AnaTrueParticlePD*>(tb);
      tp->TrjPoints.clear();
      for(int i = 0; i < 2; i++){
        AnaTrueTrajPointPD p;
        p.Position[0] = i; p.Momentum[2] = i; p.Process = tp->ID;
        tp->TrjPoints.push_back(p);
      }
      for(int pl = 0; pl < 3; pl++)tp->NHitsTrue[pl] = pl+1;
      tp->EDepTrue = tp->ID % 100;
    }
    for(AnaBunchC* bc : spill->Bunches){
      AnaBunchPD* bunch = static_cast<AnaBunchPD*>(bc);
      bunch->UnassociatedHits.assign(kNUnas, AnaHitPD());
      for(int i = 0; i < kNUnas; i++)bunch->UnassociatedHits[i].UniqueID = 1000+i;
      for(AnaParticleB* pb : bunch->Particles){
        AnaParticlePD* part = static_cast<AnaParticlePD*>(pb);
        part->SliceID = kSlice;
        part->T0 = kT0;
        part->IsClearCosmic = true;
        part->PandoraTrackScore = 0.5;
        part->SpacePoints.assign(kNSP, AnaSpacePointPD());
        for(int i = 0; i < kNSP; i++)part->SpacePoints[i].HitID = part->UniqueID*10+i;
        for(int pl = 0; pl < 3; pl++){
          part->HitsNoCalo[pl].assign(pl+1, AnaHitPD());
          part->HitsNoCalo[pl][pl].UniqueID = part->UniqueID;
        }
        part->BeamCosmicScore = 0.25;
        AnaForcedTrackPD ft;
        ft.ID = part->UniqueID;
        ft.Chi2Muon = 2*part->UniqueID;
        ft.Hits[2].assign(2, AnaHitPD());
        ft.Hits[2][1].dEdx = 2.5;
        part->ForcedTrack.assign(1, ft);
        for(int pl = 0; pl < 3; pl++)
          for(size_t h = 0; h < part->Hits[pl].size(); h++){
            part->Hits[pl][h].UniqueID = h;
            part->Hits[pl][h].TrueID = part->UniqueID;
          }
      }
    }
  }
}

//********************************************************************
int make_test_minitree(const char* in, const char* out, int n = 200){
//********************************************************************

  TFile fin(in);
  TTree* tin = (TTree*)fin.Get("MiniTree");
  if(!tin){std::cerr << "No MiniTree en " << in << std::endl; return 1;}
  TFile fout(out, "RECREATE");
  TTree* tout = new TTree("MiniTree", "event model tree");
  AnaSpillPD* spill = NULL;
  tout->Branch("Spill", "AnaSpillPD", &spill, 64000, 1);

  int nmarked = 0;
  for(int i = 0; i < n && i < tin->GetEntries(); i++){
    AnaSpill* spillIn = ReadSpill(tin, i);
    AnaSpill* tmp = spillIn->Clone();
    spill = ToSpillPD(tmp);
    delete tmp;
    FillNew(spill);
    tout->Fill();
    delete spill;

    if(i % 10 != 0){delete spillIn; continue;}
    // evento marcado, construido como en CreateHighlandMiniTree: siempre con Beam (que cuenta
    // como un spill); sin beam slice (bunch vacío, HasPandora = false) o con DQ malo (sin bunch)
    spill = new AnaSpillPD();
    spill->EventInfo   = spillIn->EventInfo->Clone();
    AnaDataQualityPD* dq = new AnaDataQualityPD();
    for(int apa = 0; apa < 6; apa++)dq->NActiveFembs[apa] = 20-apa;
    dq->UnstableHV = true;
    dq->TimingTimeStamp = 1537315216.123456;
    spill->DataQuality = dq;
    AnaBeamPD* beam = new AnaBeamPD();
    beam->BeamParticle = new AnaParticlePD();
    beam->SpillsSincePreviousSavedSpill = 1;
    beam->POTSincePreviousSavedSpill    = 0;
    spill->Beam = beam;
    AnaEventInfoPD* info = static_cast<AnaEventInfoPD*>(spill->EventInfo);
    info->Event = -info->Event;
    info->HasPandora = false;
    info->FullDump = true;   // miembros de AnaEventInfoPD: se comprueba la ida y vuelta
    info->EmptyEvent = true;
    info->NBeamRegionSpacePoints = 7;
    info->MissingReco = true;
    if((i/10) % 2 == 0){
      AnaBunchPD* bunch = new AnaBunchPD();
      bunch->Bunch = 1;
      bunch->Weight = 1;
      spill->Bunches.push_back(bunch);
    }
    else
      dq->GoodDaq = false;
    tout->Fill();
    delete spill;
    delete spillIn;
    nmarked++;
  }
  spill = NULL;

  // la cabecera es obligatoria para MiniTreeConverter
  TTree* hin = (TTree*)fin.Get("header");
  fout.cd();
  if(hin)hin->CloneTree(-1, "")->Write();
  tout->Write();
  std::cout << "Escritas " << tout->GetEntries() << " entradas (" << nmarked << " marcadas) en " << out << std::endl;
  fout.Close();
  fin.Close();
  return 0;
}

//********************************************************************
int check_test_minitree(const char* file){
//********************************************************************

  TFile f(file);
  TTree* t = (TTree*)f.Get("MiniTree");
  int ndq = 0;
  int nerr = 0, nparts = 0, nsp = 0, ntraj = 0, nunas = 0;
  for(Long64_t i = 0; i < t->GetEntries(); i++){
    AnaSpill* raw = ReadSpill(t, i);
    // lo que ve el análisis pasa por Clone(). El clon comparte la truth con el spill leído,
    // así que este se borra después (como en CreateHighlandMiniTree)
    AnaSpill* spill = raw->Clone();
    for(AnaTrueParticleB* tb : spill->TrueParticles){
      AnaTrueParticlePD* tp = static_cast<AnaTrueParticlePD*>(tb);
      ntraj += tp->TrjPoints.size();
      if(static_cast<AnaEventInfoPD*>(spill->EventInfo)->Event > 0 &&
         (tp->TrjPoints.size() != 2 || tp->TrjPoints[1].Momentum[2] != 1 || tp->TrjPoints[1].Position[0] != 1 ||
          tp->TrjPoints[0].Process != tp->ID || tp->NHitsTrue[2] != 3 || tp->EDepTrue != tp->ID % 100)){
        if(!nerr)std::cout << "Primer error: true " << tp->ID << " de la entrada " << i << " con " << tp->TrjPoints.size() << " puntos" << std::endl;
        nerr++;
      }
    }
    // FullDump solo en los eventos marcados
    const AnaEventInfoPD* einfo = static_cast<AnaEventInfoPD*>(spill->EventInfo);
    if(einfo->FullDump != (einfo->Event < 0) || einfo->EmptyEvent != (einfo->Event < 0) || einfo->MissingReco != (einfo->Event < 0) ||
       (einfo->Event < 0 && einfo->NBeamRegionSpacePoints != 7)){
      if(!nerr)std::cout << "Primer error: FullDump, EmptyEvent, NBeamRegionSpacePoints o MissingReco de la entrada " << i << std::endl;
      nerr++;
    }
    // los eventos marcados llevan AnaDataQualityPD
    if(static_cast<AnaEventInfoPD*>(spill->EventInfo)->Event < 0){
      AnaDataQualityPD* dq = dynamic_cast<AnaDataQualityPD*>(spill->DataQuality);
      if(!dq || dq->NActiveFembs[5] != 15 || !dq->UnstableHV || dq->TimingTimeStamp != 1537315216.123456){
        if(!nerr)std::cout << "Primer error: calidad de datos de la entrada " << i << std::endl;
        nerr++;
      }
      ndq++;
    }
    for(AnaBunchC* bc : spill->Bunches){
      AnaBunchPD* bunch = static_cast<AnaBunchPD*>(bc);
      nunas += bunch->UnassociatedHits.size();
      for(AnaParticleB* pb : bunch->Particles){
        AnaParticlePD* p = static_cast<AnaParticlePD*>(pb);
        nparts++;
        nsp += p->SpacePoints.size();
        if(p->SliceID != kSlice || p->T0 != kT0 || !p->IsClearCosmic || (int)p->SpacePoints.size() != kNSP ||
           p->SpacePoints[kNSP-1].HitID != p->UniqueID*10+kNSP-1 ||
           p->HitsNoCalo[2].size() != 3 || p->HitsNoCalo[2][2].UniqueID != p->UniqueID ||
           p->BeamCosmicScore != 0.25f || p->ForcedTrack.size() != 1 || p->ForcedTrack[0].Chi2Muon != 2*p->UniqueID ||
           p->ForcedTrack[0].Hits[2].size() != 2 || p->ForcedTrack[0].Hits[2][1].dEdx != 2.5f){
          if(!nerr)std::cout << "Primer error: partícula " << p->UniqueID << " de la entrada " << i << std::endl;
          nerr++;
        }
        for(int pl = 0; pl < 3; pl++)
          for(size_t h = 0; h < p->Hits[pl].size(); h++)
            if(p->Hits[pl][h].UniqueID != (Int_t)h || p->Hits[pl][h].TrueID != p->UniqueID){
              if(!nerr)std::cout << "Primer error: hit " << h << " del plano " << pl << " de la partícula " << p->UniqueID << std::endl;
              nerr++;
            }
      }
    }
    delete spill;
    delete raw;
  }
  std::cout << t->GetEntries() << " entradas, " << nparts << " partículas, " << nsp << " space points, "
            << ntraj << " puntos de trayectoria true, " << nunas << " hits no asociados, " << ndq << " AnaDataQualityPD" << std::endl;
  std::cout << (nerr ? "ERROR: " : "OK: ") << nerr << " valores distintos de los escritos" << std::endl;
  return nerr ? 1 : 0;
}
