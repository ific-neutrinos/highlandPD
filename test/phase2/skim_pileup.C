// Pileup del haz en el truth skim: ¿de dónde vienen los signal muons cuya K tiene "otra raíz" (KOrigin 6: la cadena
// de madres no acaba en la partícula del haz que dispara)? Nota técnica, cap. de truth (TODO de pileup).
//
//   root -l -b -q 'test/phase2/skim_pileup.C+("truth_skim_v100s.root", "skim_pileup.root")'
//
// Para cada signal muon (t_sigMu) con KOrigin 6: la raíz (PDG, momento, Generation, Origin), su tiempo respecto al
// haz que dispara (t de la raíz - t del haz, ns), su posición (x, y) al cruzar z = 0 (la cara del TPC) extrapolando
// en recta desde su inicio, la distancia en (x, y) a la del haz que dispara, si la K es hija directa de la raíz, y
// si el mu está reconstruido y en la beam slice. Además, el pileup en general: cuántas partículas primarias del haz
// distintas de la que dispara (primarias: ParentID 0, Origin 4) hay por evento, con su PDG y su tiempo, y cuántas
// llegan al TPC (final en z > 0). Ojo: el módulo no guarda todo el pileup, solo las primarias con algún descendiente
// guardado (SaveTrueAncestors las añade con Generation -999): es el pileup que deja algo, no el total. Resumen por pantalla y árboles "mu" (un signal muon de KOrigin 6 por entrada) y "pu" (una
// primaria de pileup por entrada) en el fichero de salida.
#include "skim_common.h"
#include "TFile.h"
#include "TTree.h"
#include <algorithm>
#include <map>
#include <string>

namespace skpu {
  // (x, y) al cruzar z = z0 en recta desde (x, y, z) con dirección (dx, dy, dz); false si no llega (dz <= 0 y z < z0)
  bool AtZ(double x, double y, double z, double dx, double dy, double dz, double z0, double& xo, double& yo){
    if(z >= z0){ xo = x; yo = y; return true; }
    if(dz <= 1e-6)return false;
    const double s = (z0 - z)/dz;
    xo = x + s*dx; yo = y + s*dy;
    return true;
  }
  const char* Name(int pdg){
    switch(pdg){ case 211: return "pi+"; case -211: return "pi-"; case 2212: return "p"; case 321: return "K+"; case -321: return "K-";
                 case -13: return "mu+"; case 13: return "mu-"; case -11: return "e+"; case 11: return "e-"; case 22: return "gamma";
                 case 2112: return "n"; case 130: return "K0L"; default: return "otro"; }
  }
  template<class M> void PrintMap(const char* title, const M& m, long tot){
    std::vector<std::pair<long, std::string>> v; for(auto& kv : m)v.push_back({kv.second, kv.first});
    std::sort(v.rbegin(), v.rend());
    printf("  %s:", title); for(auto& x : v)printf("  %s %ld (%.0f %%)", x.second.c_str(), x.first, tot ? 100.*x.first/tot : 0.); printf("\n");
  }
  const char* TBin(double dt){
    const double a = std::abs(dt);
    if(a < 1e4)return a < 10 ? "|dt| < 10 ns" : a < 100 ? "10-100 ns" : a < 1000 ? "0,1-1 us" : "1-10 us";
    // > 10 us: con signo, 10-100 us y luego por ms (la ventana de lectura es [-1, 3) ms; 26-nov-2026)
    static std::string s;
    if(a < 1e5)s = dt < 0 ? "-(10-100) us" : "10-100 us";
    else{ const long k = (long)std::floor(dt/1e6); s = "[" + std::to_string(k) + "," + std::to_string(k+1) + ") ms"; }
    return s.c_str();
  }
}

int skim_pileup(const char* input, const char* output = "skim_pileup.root"){
  using namespace skim; using namespace skpu;
  SkimReader r(input);
  TFile fout(output, "RECREATE");
  // árbol mu: signal muons de KOrigin 6
  Int_t rootPDG, rootGen, rootOrigin, kFromRoot, kParentPDG, muReco, muInSlice, rootReachesTPC, run, subrun, event;
  Float_t rootP, dt, xAt0, yAt0, dxyBeam, muP, kStartZ, muStartZ;
  TTree tmu("mu", "signal muons con KOrigin 6 (otra raíz)");
  tmu.Branch("run", &run); tmu.Branch("subrun", &subrun); tmu.Branch("event", &event);
  tmu.Branch("rootPDG", &rootPDG); tmu.Branch("rootGen", &rootGen); tmu.Branch("rootOrigin", &rootOrigin); tmu.Branch("rootP", &rootP);
  tmu.Branch("dt", &dt); tmu.Branch("xAt0", &xAt0); tmu.Branch("yAt0", &yAt0); tmu.Branch("dxyBeam", &dxyBeam); tmu.Branch("rootReachesTPC", &rootReachesTPC);
  tmu.Branch("kFromRoot", &kFromRoot); tmu.Branch("kParentPDG", &kParentPDG); tmu.Branch("kStartZ", &kStartZ);
  tmu.Branch("muP", &muP); tmu.Branch("muStartZ", &muStartZ); tmu.Branch("muReco", &muReco); tmu.Branch("muInSlice", &muInSlice);
  // árbol pu: primarias del haz que no disparan. Ojo: los minitrees hasta v1_00 no guardan el tiempo de las
  // partículas true (el módulo no rellena Position[3]): dt sale 0 o sin sentido; se avisa al final.
  Int_t puPDG, puReaches;
  Float_t puP, puDt, puX0, puY0, puEndZ;
  TTree tpu("pu", "primarias del haz distintas de la que dispara (ParentID 0, Origin 4)");
  tpu.Branch("pdg", &puPDG); tpu.Branch("p", &puP); tpu.Branch("dt", &puDt); tpu.Branch("xAt0", &puX0); tpu.Branch("yAt0", &puY0);
  tpu.Branch("endZ", &puEndZ); tpu.Branch("reachesTPC", &puReaches);

  long nTimeSet = 0;   // partículas con un tiempo distinto de 0 y de -999 (¿guarda el módulo el tiempo?)
  long nEv = 0, nEvBeam = 0, nSig = 0, nSig6 = 0, nPU = 0, nPUReach = 0, nEvPU = 0, nEvPUReach = 0;
  std::map<std::string, long> mRoot, mT, mKPar, mPUpdg, mPUt, mReco;
  std::map<int, long> nByKOrigin;
  for(Long64_t e = 0; r.Next(e); e++){
    const Event& v = r.ev;
    nEv++;
    if(!v.isMC)continue;
    const int b = v.beamIdx;
    if(b >= 0)nEvBeam++;
    double bx0 = -999, by0 = -999;
    if(b >= 0)AtZ((*v.t_x)[b], (*v.t_y)[b], (*v.t_z)[b], (*v.t_dx)[b], (*v.t_dy)[b], (*v.t_dz)[b], 0., bx0, by0);
    const double tBeam = b >= 0 ? (*v.t_t)[b] : 0;
    for(float t : *v.t_t)if(t != 0 && std::abs(t + 999) > 1e-3)nTimeSet++;
    // pileup en general
    int nThis = 0, nThisReach = 0;
    for(size_t i = 0; i < v.t_id->size(); i++){
      if((*v.t_parent)[i] != 0 || (*v.t_origin)[i] != 4 || (int)i == b)continue;   // primarias (ParentID 0) del haz
      puPDG = (*v.t_pdg)[i]; puP = (*v.t_p)[i]; puDt = b >= 0 ? (*v.t_t)[i] - tBeam : -1e9; puEndZ = (*v.t_ze)[i];
      double x0, y0; const bool ok = AtZ((*v.t_x)[i], (*v.t_y)[i], (*v.t_z)[i], (*v.t_dx)[i], (*v.t_dy)[i], (*v.t_dz)[i], 0., x0, y0);
      puX0 = ok ? x0 : -999; puY0 = ok ? y0 : -999;
      puReaches = puEndZ > 0;
      tpu.Fill();
      nPU++; nThis++; if(puReaches){ nPUReach++; nThisReach++; }
      mPUpdg[Name(puPDG)]++; mPUt[b >= 0 ? TBin(puDt) : "sin haz"]++;
    }
    if(nThis)nEvPU++; if(nThisReach)nEvPUReach++;
    // signal muons
    for(size_t i = 0; i < v.t_id->size(); i++){
      if(!(*v.t_sigMu)[i])continue;
      nSig++; nByKOrigin[(*v.t_kOrigin)[i]]++;
      if((*v.t_kOrigin)[i] != 6)continue;
      nSig6++;
      bool cut; const int root = r.Root(i, cut);
      const int k = r.Parent(i);
      run = v.run; subrun = v.subrun; event = v.event;
      rootPDG = (*v.t_pdg)[root]; rootGen = (*v.t_gen)[root]; rootOrigin = (*v.t_origin)[root]; rootP = (*v.t_p)[root];
      dt = b >= 0 ? (*v.t_t)[root] - tBeam : -1e9;
      double x0, y0; const bool ok = AtZ((*v.t_x)[root], (*v.t_y)[root], (*v.t_z)[root], (*v.t_dx)[root], (*v.t_dy)[root], (*v.t_dz)[root], 0., x0, y0);
      xAt0 = ok ? x0 : -999; yAt0 = ok ? y0 : -999;
      dxyBeam = (ok && b >= 0 && bx0 > -900) ? std::hypot(x0 - bx0, y0 - by0) : -1;
      rootReachesTPC = (*v.t_ze)[root] > 0;
      kFromRoot = k >= 0 && (*v.t_parent)[k] == (*v.t_id)[root];
      kParentPDG = k >= 0 ? (r.Parent(k) >= 0 ? (*v.t_pdg)[r.Parent(k)] : 0) : 0;
      kStartZ = k >= 0 ? (*v.t_z)[k] : -999;
      muP = 1000.*(*v.t_p)[i]; muStartZ = (*v.t_z)[i];
      muReco = (*v.t_recoN)[i] > 0; muInSlice = (*v.t_recoInSlice)[i];
      tmu.Fill();
      mRoot[std::string(Name(rootPDG)) + (rootGen == 0 && rootOrigin == 4 ? "" : Form(" (Gen %d, Origin %d)", rootGen, rootOrigin))]++;
      mT[b >= 0 ? TBin(dt) : "sin haz"]++;
      mKPar[kFromRoot ? "K hija directa de la raíz" : Form("K de una %s", Name(kParentPDG))]++;
      mReco[!muReco ? "mu sin reconstruir" : (muInSlice ? "mu en la beam slice" : "mu fuera de la beam slice")]++;
    }
  }

  printf("\nEntradas %ld (con el haz que dispara en la truth %ld)\n", nEv, nEvBeam);
  printf("Signal muons %ld; por KOrigin:", nSig); for(auto& kv : nByKOrigin)printf("  %d: %ld", kv.first, kv.second); printf("\n");
  printf("\n== Signal muons con otra raíz (KOrigin 6): %ld (%.1f %%) ==\n", nSig6, nSig ? 100.*nSig6/nSig : 0.);
  PrintMap("raíz", mRoot, nSig6);
  if(nTimeSet)PrintMap("tiempo de la raíz respecto al haz que dispara", mT, nSig6);
  PrintMap("madre de la K", mKPar, nSig6);
  PrintMap("reconstrucción del mu", mReco, nSig6);
  if(tmu.GetEntries()){
    tmu.Draw("dxyBeam", "dxyBeam>=0", "goff");
    std::vector<double> d(tmu.GetV1(), tmu.GetV1() + tmu.GetSelectedRows()); std::sort(d.begin(), d.end());
    if(!d.empty())printf("  distancia en (x, y) a z = 0 a la del haz que dispara: mediana %.1f cm, p10 %.1f, p90 %.1f (%zu)\n",
                         d[d.size()/2], d[d.size()/10], d[9*d.size()/10], d.size());
    printf("  la raíz llega al TPC (final en z > 0): %lld de %lld\n", tmu.GetEntries("rootReachesTPC"), tmu.GetEntries());
  }
  printf("\n== Pileup en general: primarias del haz distintas de la que dispara (solo las guardadas: con algún descendiente guardado) ==\n");
  printf("  %ld en %ld entradas (%.2f por entrada); eventos con alguna %ld (%.1f %%), con alguna que llega al TPC %ld (%.1f %%)\n",
         nPU, nEv, nEv ? double(nPU)/nEv : 0., nEvPU, nEv ? 100.*nEvPU/nEv : 0., nEvPUReach, nEv ? 100.*nEvPUReach/nEv : 0.);
  PrintMap("PDG", mPUpdg, nPU);
  if(nTimeSet)PrintMap("tiempo respecto al haz que dispara", mPUt, nPU);
  else printf("  AVISO: ninguna partícula true tiene tiempo (0 o -999 en todas): el minitree no lo guarda (módulo hasta v1_00);\n"
              "  no hay estudio de tiempo. Hace falta rellenar Position[3] y PositionEnd[3] en CreateHighlandMiniTree.\n");
  fout.cd(); tmu.Write(); tpu.Write(); fout.Close();
  printf("\nÁrboles mu y pu en %s\n", output);
  return 0;
}
