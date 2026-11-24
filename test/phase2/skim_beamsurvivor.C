// El pion del haz que sobrevive a su primera interacción (opción a de la medida A; nota técnica, cap. de truth), en el
// truth skim. Eventos con un pi+ del haz que dispara que acaba en una interacción inelástica (no desintegración) con el
// vértice en el volumen fiducial de kaon_acceptance (x en [-340, 340], y en [20, 587.5], z en [0, 600] cm).
//
//   root -l -b -q 'test/phase2/skim_beamsurvivor.C+("truth_skim_v100s.root", "skim_beamsurvivor.root")'
//
// Para cada interacción: el pi+ hijo de mayor momento (el "superviviente") y r = p del superviviente / p del haz al
// final; clases de r (sin pi+, < 0,5, 0,5-0,9, 0,9-0,97, > 0,97). Truth, por clase: qué más sale de la interacción
// (hijas con hits propios en el plano de colección, nHits2 >= kVisHits: hadrones cargados, protones, núcleos; neutrones
// y gammas aparte), ángulo entre el haz y el superviviente, energía transferida (E del haz - E del superviviente) y K+
// que produce el superviviente (hijas directas K+: la producción que la opción a pasaría a A). Reco: con qué PFParticle
// está emparejado el superviviente (el primario de la beam slice, es decir, la misma traza que el haz; otra partícula
// de la beam slice; fuera de la beam slice; sin reconstruir) y dónde acaba el primario respecto a la interacción true.
// Salida: resumen por pantalla y árbol "int" (una interacción por entrada).
#include "skim_common.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include <algorithm>
#include <map>
#include <string>

namespace sksurv {
  const double kFV[6] = {-340., 340., 20., 587.5, 0., 600.};
  const int kVisHits = 5;           // hits del plano de colección para considerar visible una hija
  const double kMpi = 0.13957;      // GeV
  bool InFV(double x, double y, double z){ return x > kFV[0] && x < kFV[1] && y > kFV[2] && y < kFV[3] && z > kFV[4] && z < kFV[5]; }
  int RClass(double r){ return r < 0 ? 0 : r < 0.5 ? 1 : r < 0.9 ? 2 : r < 0.97 ? 3 : 4; }
  const char* kRName[5] = {"sin pi+", "r < 0,5", "0,5-0,9", "0,9-0,97", "> 0,97"};
  // topología de lo que acompaña al superviviente (hijas visibles salvo él)
  const char* kTopoName[4] = {"solo el pi", "pi + protones/núcleos", "pi + otros cargados", "nada visible"};
}

int skim_beamsurvivor(const char* input, const char* output = "skim_beamsurvivor.root"){
  using namespace skim; using namespace sksurv;
  SkimReader r(input);
  // también la beam slice (s_*) para la reco
  std::vector<int>* s_uid = nullptr; std::vector<float> *s_xe = nullptr, *s_ye = nullptr, *s_ze = nullptr, *s_length = nullptr;
  r.chain.SetBranchStatus("s_uid", 1); r.chain.SetBranchAddress("s_uid", &s_uid);
  r.chain.SetBranchStatus("s_xe", 1); r.chain.SetBranchAddress("s_xe", &s_xe);
  r.chain.SetBranchStatus("s_ye", 1); r.chain.SetBranchAddress("s_ye", &s_ye);
  r.chain.SetBranchStatus("s_ze", 1); r.chain.SetBranchAddress("s_ze", &s_ze);
  r.chain.SetBranchStatus("s_length", 1); r.chain.SetBranchAddress("s_length", &s_length);
  std::vector<int>* t_recoUID = nullptr;
  r.chain.SetBranchStatus("t_recoUID", 1); r.chain.SetBranchAddress("t_recoUID", &t_recoUID);

  TFile fout(output, "RECREATE");
  Int_t rcls, topo, nVisCh, nVisP, nVisOther, nN, nGamma, nKp, recoCls, hasSlice;
  Float_t rr, pBeam, pLead, angle, dE, dEndBeam, dEndLead;
  TTree tint("int", "interacciones del pi+ del haz en el FV");
  tint.Branch("r", &rr); tint.Branch("rcls", &rcls); tint.Branch("pBeam", &pBeam); tint.Branch("pLead", &pLead); tint.Branch("angle", &angle); tint.Branch("dE", &dE);
  tint.Branch("topo", &topo); tint.Branch("nVisCh", &nVisCh); tint.Branch("nVisP", &nVisP); tint.Branch("nVisOther", &nVisOther); tint.Branch("nN", &nN);
  tint.Branch("nGamma", &nGamma); tint.Branch("nKp", &nKp); tint.Branch("recoCls", &recoCls); tint.Branch("hasSlice", &hasSlice);
  tint.Branch("dEndBeam", &dEndBeam); tint.Branch("dEndLead", &dEndLead);

  long nInt = 0, nCls[5] = {0}, nTopo[5][4] = {{0}}, nReco[5][4] = {{0}}, nK[5] = {0};
  std::map<int, std::vector<double>> vAngle, vdE, vdEndBeam, vdEndLead;
  const char* kRecoName[4] = {"misma traza que el haz", "otra partícula de la beam slice", "fuera de la beam slice", "sin reconstruir"};
  for(Long64_t e = 0; r.Next(e); e++){
    const Event& v = r.ev;
    const int b = v.beamIdx;
    if(!v.isMC || b < 0 || (*v.t_pdg)[b] != 211)continue;
    if((*v.t_procEnd)[b] == kDecay || !InFV((*v.t_xe)[b], (*v.t_ye)[b], (*v.t_ze)[b]))continue;
    nInt++;
    const int bid = (*v.t_id)[b];
    pBeam = (*v.t_pEnd)[b];
    // hijas de la interacción que acaba con el haz
    int lead = -1;
    std::vector<int> dau;
    for(size_t i = 0; i < v.t_id->size(); i++){
      if((*v.t_parent)[i] != bid)continue;
      dau.push_back(i);
      if((*v.t_pdg)[i] == 211 && (lead < 0 || (*v.t_p)[i] > (*v.t_p)[lead]))lead = i;
    }
    rr = (lead >= 0 && pBeam > 0) ? (*v.t_p)[lead]/pBeam : -1;
    rcls = RClass(rr); nCls[rcls]++;
    pLead = lead >= 0 ? (*v.t_p)[lead] : -1;
    nVisCh = nVisP = nVisOther = nN = nGamma = 0;
    for(int i : dau){
      if(i == lead)continue;
      const int pdg = (*v.t_pdg)[i], a = std::abs(pdg);
      if(a == 2112){ nN++; continue; }
      if(pdg == 22){ nGamma++; continue; }
      if((*v.t_nHits2)[i] < kVisHits)continue;
      if(a == 2212 || pdg > 1000000000)nVisP++;
      else if(a == 211 || a == 321 || a == 13 || a == 11)nVisCh++;
      else nVisOther++;
    }
    if(lead < 0)topo = 3;
    else topo = (nVisCh + nVisOther) > 0 ? 2 : nVisP > 0 ? 1 : 0;
    if(lead >= 0 && nVisCh + nVisOther + nVisP == 0)topo = 0;
    nTopo[rcls][topo]++;
    // cinemática
    angle = -1; dE = -1;
    if(lead >= 0){
      TVector3 db((*v.t_dxe)[b], (*v.t_dye)[b], (*v.t_dze)[b]), dl((*v.t_dx)[lead], (*v.t_dy)[lead], (*v.t_dz)[lead]);
      if(db.Mag() > 0.5 && dl.Mag() > 0.5)angle = db.Angle(dl)*180./M_PI;
      dE = 1000.*(std::sqrt(pBeam*pBeam + kMpi*kMpi) - std::sqrt(pLead*pLead + kMpi*kMpi));
      vAngle[rcls].push_back(angle); vdE[rcls].push_back(dE);
    }
    // K+ hijas directas del superviviente
    nKp = 0;
    if(lead >= 0){ const int lid = (*v.t_id)[lead]; for(size_t i = 0; i < v.t_id->size(); i++)if((*v.t_parent)[i] == lid && (*v.t_pdg)[i] == 321)nKp++; }
    nK[rcls] += nKp;
    // reco del superviviente
    hasSlice = !s_uid->empty();
    recoCls = 3;
    if(lead >= 0 && (*v.t_recoN)[lead] > 0){
      const int uid = (*t_recoUID)[lead];
      if(hasSlice && uid == (*s_uid)[0])recoCls = 0;
      else if((*v.t_recoInSlice)[lead])recoCls = 1;
      else recoCls = 2;
    }
    if(lead >= 0)nReco[rcls][recoCls]++;
    // dónde acaba el primario de la beam slice: distancia al final true del haz y al final true del superviviente
    dEndBeam = dEndLead = -1;
    if(hasSlice){
      dEndBeam = Dist((*s_xe)[0], (*s_ye)[0], (*s_ze)[0], (*v.t_xe)[b], (*v.t_ye)[b], (*v.t_ze)[b]);
      if(lead >= 0)dEndLead = Dist((*s_xe)[0], (*s_ye)[0], (*s_ze)[0], (*v.t_xe)[lead], (*v.t_ye)[lead], (*v.t_ze)[lead]);
      vdEndBeam[rcls].push_back(dEndBeam); if(lead >= 0)vdEndLead[rcls].push_back(dEndLead);
    }
    tint.Fill();
  }

  auto Med = [](std::vector<double> x){ if(x.empty())return -1.; std::sort(x.begin(), x.end()); return x[x.size()/2]; };
  printf("\nInteracciones inelásticas del pi+ del haz en el FV: %ld\n", nInt);
  printf("\n  clase de r      interacc.   K+ del superv. |  %s / %s / %s / %s\n", kTopoName[0], kTopoName[1], kTopoName[2], kTopoName[3]);
  for(int c = 0; c < 5; c++){
    printf("  %-12s %10ld %8ld        |", kRName[c], nCls[c], nK[c]);
    for(int t = 0; t < 4; t++)printf("  %6ld (%4.1f %%)", nTopo[c][t], nCls[c] ? 100.*nTopo[c][t]/nCls[c] : 0.);
    printf("\n");
  }
  printf("\n  clase de r     ángulo haz-superv. (med., grados)   E transferida (med., MeV)\n");
  for(int c = 1; c < 5; c++)printf("  %-12s %20.1f %30.0f\n", kRName[c], Med(vAngle[c]), Med(vdE[c]));
  printf("\n  Reco del superviviente:\n  clase de r  ");
  for(int k = 0; k < 4; k++)printf("  %-30s", kRecoName[k]); printf("\n");
  for(int c = 1; c < 5; c++){
    long tot = 0; for(int k = 0; k < 4; k++)tot += nReco[c][k];
    printf("  %-12s", kRName[c]);
    for(int k = 0; k < 4; k++)printf("  %6ld (%4.1f %%)                ", nReco[c][k], tot ? 100.*nReco[c][k]/tot : 0.);
    printf("\n");
  }
  printf("\n  Final del primario de la beam slice (sin SCE) frente al final true (mediana, cm):\n  clase de r    al final del haz   al final del superviviente\n");
  for(int c = 0; c < 5; c++)printf("  %-12s %16.1f %24.1f\n", kRName[c], Med(vdEndBeam[c]), Med(vdEndLead[c]));
  fout.cd(); tint.Write(); fout.Close();
  printf("\nÁrbol int en %s\n", output);
  return 0;
}
