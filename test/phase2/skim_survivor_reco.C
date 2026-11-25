// El pion del haz que sobrevive (opción a de la medida A), en reco: ¿se ve la interacción en la traza del primario de la
// beam slice cuando Pandora la sigue con el superviviente? Continuación de skim_beamsurvivor.C, en el truth skim.
//
//   root -l -b -q 'test/phase2/skim_survivor_reco.C+("truth_skim_v100s.root", "skim_survivor_reco.root")'
//
// Interacciones inelásticas del pi+ del haz que dispara con el vértice en el FV de kaon_acceptance, clases de r como en
// skim_beamsurvivor.C (r = p del pi+ hijo de mayor momento / p del haz al final). Se usan los hits del plano de colección
// del primario de la beam slice (bh_*, posición con SCE corregida, en orden de rango residual): "continúa" si tiene al
// menos kMinHits hits del haz y kMinHits del superviviente (bh_trueID). La transición (la interacción true vista en la
// traza, sin depender del SCE) es el último hit del haz a lo largo de la traza: s_T, con s = RRmax - RR.
// Tres firmas, en la transición frente a la traza del haz lejos de ella (control: el tramo del haz a más de 2L de s_T):
//  1. codo: ángulo entre las direcciones (componente principal) de los hits en [s - L, s] y en [s, s + L], con al
//     menos kMinHits a cada lado; en cada hit de la traza. En la transición, y el máximo del tramo de control;
//  2. dE/dx: máximo del dE/dx de los hits a < kVtxWin cm de s_T, frente al mismo a 30 cm antes (control);
//  3. hijas de Pandora del primario que empiezan lejos de su final (> kMidDist cm, sin SCE): empiezan a mitad de traza,
//     y si su partícula true es hija de la interacción.
// Salida: resumen por pantalla y árbol "sv" (una interacción que continúa por entrada).
#include "skim_common.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include <algorithm>
#include <map>

namespace sksr {
  const double kFV[6] = {-340., 340., 20., 587.5, 0., 600.};
  const int kMinHits = 5;
  const double kL = 10.;        // cm, ventana de cada lado para las direcciones
  const double kVtxWin = 2.;    // cm, ventana de dE/dx alrededor de la transición
  const double kMidDist = 10.;  // cm, una hija empieza "a mitad de traza" si su inicio está a más de esto del final
  const double kThr[6] = {1., 2., 3., 5., 10., 20.};   // grados
  // búsqueda combinada a lo largo de la traza: codo > kCK y dE/dx máximo > kCD en el mismo punto
  const double kCK[3] = {1., 2., 3.}, kCD[4] = {8., 10., 15., 20.};
  const double kFound = 3.;     // cm, "encontrada" si algún punto que pasa está a menos de esto de la transición
  bool InFV(double x, double y, double z){ return x > kFV[0] && x < kFV[1] && y > kFV[2] && y < kFV[3] && z > kFV[4] && z < kFV[5]; }
  int RClass(double r){ return r < 0 ? 0 : r < 0.5 ? 1 : r < 0.9 ? 2 : r < 0.97 ? 3 : 4; }
  const char* kRName[5] = {"sin pi+", "r < 0,5", "0,5-0,9", "0,9-0,97", "> 0,97"};
  struct Hit { double s, dEdx; TVector3 x; int tid; };
  // dirección principal de los puntos (iteración de potencias sobre la covarianza); false con menos de kMinHits
  bool Dir(const std::vector<Hit>& h, double s0, double s1, TVector3& d){
    TVector3 c(0, 0, 0); int n = 0;
    for(const Hit& x : h)if(x.s >= s0 && x.s <= s1){ c += x.x; n++; }
    if(n < kMinHits)return false;
    c *= 1./n;
    double m[3][3] = {{0}};
    for(const Hit& x : h)if(x.s >= s0 && x.s <= s1){ const TVector3 v = x.x - c; for(int i = 0; i < 3; i++)for(int j = 0; j < 3; j++)m[i][j] += v[i]*v[j]; }
    TVector3 v(1, 1, 1);
    for(int it = 0; it < 30; it++){
      TVector3 w(m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2], m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2], m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2]);
      if(w.Mag() <= 0)return false;
      v = w.Unit();
    }
    d = v; return true;
  }
  double Kink(const std::vector<Hit>& h, double s){
    TVector3 a, b;
    if(!Dir(h, s - kL, s, a) || !Dir(h, s, s + kL, b))return -1;
    return std::acos(std::min(1., std::abs(a.Dot(b))))*180./M_PI;
  }
  double MaxdEdx(const std::vector<Hit>& h, double s){ double m = -1; for(const Hit& x : h)if(std::abs(x.s - s) < kVtxWin && x.dEdx > 0 && x.dEdx < 100)m = std::max(m, x.dEdx); return m; }
  double Med(std::vector<double> x){ if(x.empty())return -1; std::sort(x.begin(), x.end()); return x[x.size()/2]; }
}

int skim_survivor_reco(const char* input, const char* output = "skim_survivor_reco.root"){
  using namespace skim; using namespace sksr;
  SkimReader r(input);
  std::vector<int> *s_uid = nullptr, *s_parent = nullptr, *s_trueID = nullptr, *bh_trueID = nullptr;
  std::vector<float> *s_x = nullptr, *s_y = nullptr, *s_z = nullptr, *s_xe = nullptr, *s_ye = nullptr, *s_ze = nullptr, *s_length = nullptr;
  std::vector<float> *bh_dEdx = nullptr, *bh_rr = nullptr, *bh_x = nullptr, *bh_y = nullptr, *bh_z = nullptr;
#define SRB(x) r.chain.SetBranchStatus(#x, 1); r.chain.SetBranchAddress(#x, &x);
  SRB(s_uid) SRB(s_parent) SRB(s_trueID) SRB(s_x) SRB(s_y) SRB(s_z) SRB(s_xe) SRB(s_ye) SRB(s_ze) SRB(s_length)
  SRB(bh_dEdx) SRB(bh_rr) SRB(bh_x) SRB(bh_y) SRB(bh_z) SRB(bh_trueID)
#undef SRB

  TFile fout(output, "RECREATE");
  Float_t rr, angle, kinkT, kinkCtrl, kinkMax, sKinkMax, sT, len, dEdxT, dEdxCtrl; Int_t rcls, nBeamHits, nLeadHits, nMid, nMidTrue, nMidTrueAll;
  TTree tsv("sv", "interacciones del pi+ del haz que el primario sigue con el superviviente");
  tsv.Branch("r", &rr); tsv.Branch("rcls", &rcls); tsv.Branch("angle", &angle); tsv.Branch("kinkT", &kinkT); tsv.Branch("kinkCtrl", &kinkCtrl);
  tsv.Branch("kinkMax", &kinkMax); tsv.Branch("sKinkMax", &sKinkMax); tsv.Branch("sT", &sT); tsv.Branch("len", &len);
  tsv.Branch("dEdxT", &dEdxT); tsv.Branch("dEdxCtrl", &dEdxCtrl); tsv.Branch("nBeamHits", &nBeamHits); tsv.Branch("nLeadHits", &nLeadHits);
  tsv.Branch("nMid", &nMid); tsv.Branch("nMidTrue", &nMidTrue); tsv.Branch("nMidTrueAll", &nMidTrueAll);

  long nInt[5] = {0}, nCont[5] = {0}, nKinkOK[5] = {0}, nCtrlOK[5] = {0};
  long passT[5][6] = {{0}}, passCtrl[5][6] = {{0}}, maxNearT[5][6] = {{0}};
  long midAny[5] = {0}, midTrue[5] = {0};
  long cat[5][4] = {{0}}, leadSlice[5] = {0}, leadOut[5] = {0}; std::map<int, std::vector<double>> vAfter;
  long cFound[5][3][4] = {{{0}}}, cFake[5][3][4] = {{{0}}}, nCtrlTrk[5] = {0}; double ctrlLen[5] = {0};
  std::map<int, std::vector<double>> vKinkT, vKinkCtrl, vdT, vdC, vAng;
  // codo reco frente al ángulo true (superviviente que continúa, todas las clases con r >= 0,5)
  const double angBins[6] = {0, 2, 5, 10, 20, 180}; std::vector<double> vKinkByAng[5], vAngByAng[5];
  for(Long64_t e = 0; r.Next(e); e++){
    const Event& v = r.ev;
    const int b = v.beamIdx;
    if(!v.isMC || b < 0 || (*v.t_pdg)[b] != 211)continue;
    if((*v.t_procEnd)[b] == kDecay || !InFV((*v.t_xe)[b], (*v.t_ye)[b], (*v.t_ze)[b]))continue;
    const int bid = (*v.t_id)[b];
    const double pB = (*v.t_pEnd)[b];
    int lead = -1;
    for(size_t i = 0; i < v.t_id->size(); i++)
      if((*v.t_parent)[i] == bid && (*v.t_pdg)[i] == 211 && (lead < 0 || (*v.t_p)[i] > (*v.t_p)[lead]))lead = i;
    rr = (lead >= 0 && pB > 0) ? (*v.t_p)[lead]/pB : -1;
    rcls = RClass(rr);
    nInt[rcls]++;
    if(lead < 0)continue;
    if(s_uid->empty() || bh_rr->empty()){ cat[rcls][0]++; continue; }
    const int lid = (*v.t_id)[lead];
    // hits del primario, ordenados a lo largo de la traza
    std::vector<Hit> h;
    double rrMax = -1; for(float x : *bh_rr)rrMax = std::max(rrMax, (double)x);
    for(size_t i = 0; i < bh_rr->size(); i++){
      if((*bh_rr)[i] < 0)continue;
      h.push_back({rrMax - (*bh_rr)[i], (*bh_dEdx)[i], TVector3((*bh_x)[i], (*bh_y)[i], (*bh_z)[i]), (*bh_trueID)[i]});
    }
    std::sort(h.begin(), h.end(), [](const Hit& a, const Hit& c){ return a.s < c.s; });
    nBeamHits = nLeadHits = 0; double sLastBeam = -1;
    for(const Hit& x : h){ if(x.tid == bid){ nBeamHits++; sLastBeam = std::max(sLastBeam, x.s); } if(x.tid == lid)nLeadHits++; }
    // qué hace el primario con la interacción: 1 no es el haz (< kMinHits del haz), 2 acaba en ella (haz sí, superviviente
    // no), 3 sigue con el superviviente. En 2: cuánta traza queda después del último hit del haz, y si el superviviente
    // está reconstruido aparte (otra partícula de la beam slice o fuera)
    if(nBeamHits < kMinHits){ cat[rcls][1]++; continue; }
    if(nLeadHits < kMinHits){
      cat[rcls][2]++;
      vAfter[rcls].push_back(h.back().s - sLastBeam);
      if((*v.t_recoN)[lead] > 0){ if((*v.t_recoInSlice)[lead])leadSlice[rcls]++; else leadOut[rcls]++; }
      continue;
    }
    cat[rcls][3]++;
    nCont[rcls]++;
    sT = sLastBeam; len = h.empty() ? 0 : h.back().s;
    angle = -1;
    { TVector3 db((*v.t_dxe)[b], (*v.t_dye)[b], (*v.t_dze)[b]), dl((*v.t_dx)[lead], (*v.t_dy)[lead], (*v.t_dz)[lead]);
      if(db.Mag() > 0.5 && dl.Mag() > 0.5)angle = db.Angle(dl)*180./M_PI; }
    // 1. codo en la transición, máximo del control (tramo del haz a más de 2L de s_T) y máximo de toda la traza
    kinkT = Kink(h, sT);
    kinkCtrl = -1; kinkMax = -1; sKinkMax = -1;
    for(const Hit& x : h){
      const double k = Kink(h, x.s);
      if(k < 0)continue;
      if(k > kinkMax){ kinkMax = k; sKinkMax = x.s; }
      if(x.s < sT - 2*kL)kinkCtrl = std::max((double)kinkCtrl, k);
    }
    if(kinkT >= 0){ nKinkOK[rcls]++; vKinkT[rcls].push_back(kinkT); for(int t = 0; t < 6; t++)if(kinkT > kThr[t])passT[rcls][t]++; }
    if(kinkCtrl >= 0){ nCtrlOK[rcls]++; vKinkCtrl[rcls].push_back(kinkCtrl); for(int t = 0; t < 6; t++)if(kinkCtrl > kThr[t])passCtrl[rcls][t]++; }
    // criterio reco: el codo máximo de la traza supera el umbral y está a < 5 cm de la transición
    if(kinkMax >= 0)for(int t = 0; t < 6; t++)if(kinkMax > kThr[t] && std::abs(sKinkMax - sT) < 5.)maxNearT[rcls][t]++;
    if(kinkT >= 0 && angle >= 0 && rcls >= 2){
      int ab = 0; while(ab < 4 && angle >= angBins[ab + 1])ab++;
      vKinkByAng[ab].push_back(kinkT);
    }
    // búsqueda combinada: en la transición (a < kFound cm) y en el tramo del haz lejos de ella (control: s en
    // [kL, s_T - 2L]; falsos por traza y longitud del control)
    { bool found[3][4] = {{false}}, fake[3][4] = {{false}};
      const double c0 = kL, c1 = sT - 2*kL;
      const bool hasCtrl = c1 - c0 > 20.;
      for(const Hit& x : h){
        const bool nearT = std::abs(x.s - sT) < kFound, inCtrl = hasCtrl && x.s > c0 && x.s < c1;
        if(!nearT && !inCtrl)continue;
        const double k = Kink(h, x.s), d = MaxdEdx(h, x.s);
        if(k < 0 || d < 0)continue;
        for(int a = 0; a < 3; a++)for(int q = 0; q < 4; q++)if(k > kCK[a] && d > kCD[q]){ if(nearT)found[a][q] = true; if(inCtrl)fake[a][q] = true; }
      }
      for(int a = 0; a < 3; a++)for(int q = 0; q < 4; q++){ if(found[a][q])cFound[rcls][a][q]++; if(hasCtrl && fake[a][q])cFake[rcls][a][q]++; }
      if(hasCtrl){ nCtrlTrk[rcls]++; ctrlLen[rcls] += c1 - c0; } }
    // 2. dE/dx en la transición y 30 cm antes
    dEdxT = MaxdEdx(h, sT); dEdxCtrl = sT > 30 + kVtxWin ? MaxdEdx(h, sT - 30) : -1;
    if(dEdxT > 0)vdT[rcls].push_back(dEdxT); if(dEdxCtrl > 0)vdC[rcls].push_back(dEdxCtrl);
    // 3. hijas de Pandora del primario que empiezan a mitad de traza
    nMid = nMidTrue = 0;
    for(size_t j = 1; j < s_uid->size(); j++){
      if((*s_parent)[j] != (*s_uid)[0])continue;
      const double d = Dist((*s_x)[j], (*s_y)[j], (*s_z)[j], (*s_xe)[0], (*s_ye)[0], (*s_ze)[0]);
      const double d2 = Dist((*s_xe)[j], (*s_ye)[j], (*s_ze)[j], (*s_xe)[0], (*s_ye)[0], (*s_ze)[0]);
      if(std::min(d, d2) < kMidDist)continue;
      nMid++;
      const int ti = r.Index((*s_trueID)[j]);
      if(ti >= 0 && (*v.t_parent)[ti] == bid)nMidTrue++;
    }
    // hijas true visibles de la interacción (salvo el superviviente), para comparar
    nMidTrueAll = 0;
    for(size_t i = 0; i < v.t_id->size(); i++)if((*v.t_parent)[i] == bid && (int)i != lead && (*v.t_nHits2)[i] >= kMinHits && std::abs((*v.t_pdg)[i]) != 2112 && (*v.t_pdg)[i] != 22)nMidTrueAll++;
    if(nMid > 0)midAny[rcls]++; if(nMidTrue > 0)midTrue[rcls]++;
    tsv.Fill();
  }

  printf("\n0. Qué hace el primario de la beam slice con la interacción (hits del plano de colección, por bh_trueID):\n");
  printf("  clase       interacc.  sin primario   primario no es el haz   acaba en la interacción   sigue con el superviviente\n");
  for(int c = 1; c < 5; c++){
    const double t = nInt[c];
    printf("  %-10s %8ld   %5.1f %%        %5.1f %%                %5.1f %%                   %5.1f %%\n", kRName[c], nInt[c],
           100.*cat[c][0]/t, 100.*cat[c][1]/t, 100.*cat[c][2]/t, 100.*cat[c][3]/t);
  }
  printf("  Cuando acaba en la interacción: traza después del último hit del haz (mediana, cm) y superviviente reconstruido aparte (en la beam slice / fuera):\n");
  for(int c = 1; c < 5; c++)printf("  %-10s %6.1f cm   %5.1f %% / %5.1f %%\n", kRName[c], Med(vAfter[c]), cat[c][2] ? 100.*leadSlice[c]/cat[c][2] : 0., cat[c][2] ? 100.*leadOut[c]/cat[c][2] : 0.);
  printf("\nInteracciones inelásticas del pi+ del haz en el FV, por clase de r; el primario sigue con el superviviente (>= %d hits del haz y del superviviente):\n", kMinHits);
  for(int c = 1; c < 5; c++)printf("  %-10s %6ld interacciones, continúa %5ld (%4.1f %%)\n", kRName[c], nInt[c], nCont[c], nInt[c] ? 100.*nCont[c]/nInt[c] : 0.);
  printf("\n1. Codo (ángulo entre las direcciones de %.0f cm de hits a cada lado, grados), en las que continúan:\n", kL);
  printf("  clase       con codo calculable: transición / control    mediana: transición / control (máx. del tramo del haz)\n");
  for(int c = 1; c < 5; c++)printf("  %-10s %10ld / %6ld %30.2f / %.2f\n", kRName[c], nKinkOK[c], nCtrlOK[c], Med(vKinkT[c]), Med(vKinkCtrl[c]));
  printf("\n  Fracción por encima del umbral: en la transición | máximo del control | codo máximo de la traza a < 5 cm de la transición\n  clase     ");
  for(int t = 0; t < 6; t++)printf("  > %4.0f°        ", kThr[t]); printf("\n");
  for(int c = 1; c < 5; c++){
    printf("  %-10s", kRName[c]);
    for(int t = 0; t < 6; t++)printf("  %.2f|%.2f|%.2f", nKinkOK[c] ? 1.*passT[c][t]/nKinkOK[c] : 0., nCtrlOK[c] ? 1.*passCtrl[c][t]/nCtrlOK[c] : 0., nKinkOK[c] ? 1.*maxNearT[c][t]/nKinkOK[c] : 0.);
    printf("\n");
  }
  printf("\n  Codo reco en la transición frente al ángulo true haz-superviviente (r >= 0,5): mediana [16 %%, 84 %%]\n");
  for(int a = 0; a < 5; a++){
    auto& x = vKinkByAng[a]; std::sort(x.begin(), x.end());
    if(x.empty()){ printf("    true %3.0f-%3.0f°: -\n", angBins[a], angBins[a + 1]); continue; }
    printf("    true %3.0f-%3.0f°: %5zu  %.2f [%.2f, %.2f]\n", angBins[a], angBins[a + 1], x.size(), x[x.size()/2], x[size_t(0.16*x.size())], x[size_t(0.84*x.size())]);
  }
  printf("\n2. dE/dx máximo a < %.0f cm (MeV/cm), mediana: transición / 30 cm antes\n", kVtxWin);
  for(int c = 1; c < 5; c++)printf("  %-10s %.2f / %.2f\n", kRName[c], Med(vdT[c]), Med(vdC[c]));
  printf("\n3. Hijas de Pandora del primario que empiezan a mitad de traza (> %.0f cm del final): con alguna / con alguna hija true de la interacción\n", kMidDist);
  for(int c = 1; c < 5; c++)printf("  %-10s %5.1f %% / %5.1f %%\n", kRName[c], nCont[c] ? 100.*midAny[c]/nCont[c] : 0., nCont[c] ? 100.*midTrue[c]/nCont[c] : 0.);
  printf("\n4. Búsqueda combinada (codo > K y dE/dx máximo a < %.0f cm > D en el mismo punto). Encontrada: algún punto a < %.0f cm de la\n", kVtxWin, kFound);
  printf("   transición pasa (sobre las que continúan). Falsa: algún punto del tramo del haz lejos de la transición pasa (sobre las que tienen\n");
  printf("   > 20 cm de control); entre paréntesis, falsas por metro de control.\n");
  for(int c = 1; c < 5; c++){
    printf("  %-10s (control: %ld trazas, %.0f cm de media)\n", kRName[c], nCtrlTrk[c], nCtrlTrk[c] ? ctrlLen[c]/nCtrlTrk[c] : 0.);
    for(int a = 0; a < 3; a++){
      printf("    K > %.0f°:", kCK[a]);
      for(int q = 0; q < 4; q++)printf("   D > %2.0f: %.2f / %.2f (%.2f/m)", kCD[q], nCont[c] ? 1.*cFound[c][a][q]/nCont[c] : 0., nCtrlTrk[c] ? 1.*cFake[c][a][q]/nCtrlTrk[c] : 0.,
                                         ctrlLen[c] > 0 ? 100.*cFake[c][a][q]/ctrlLen[c] : 0.);
      printf("\n");
    }
  }
  fout.cd(); tsv.Write(); fout.Close();
  printf("\nÁrbol sv en %s\n", output);
  return 0;
}
