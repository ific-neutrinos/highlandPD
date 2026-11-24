// Ventana del mu de K+ -> mu+ nu en el truth skim (nota técnica, cap. de truth): momento true de los mu+ hijos de
// desintegración de una K+ (no cósmicos), según el momento final de la K (en reposo = 0 en Geant4, o en vuelo), y la
// fracción dentro de [230, 241] MeV/c (la de IsSignalMuon) y de [210, 260] MeV/c. Una K de 50 MeV/c da un mu de
// 235,5 +- ~26 MeV/c: la ventana estrecha solo deja K de muy poco momento (la aceptancia usaba 50 MeV/c).
//
//   root -l -b -q 'test/phase2/skim_muwindow.C+("truth_skim_v100s.root", "skim_muwindow.root")'
//
// Salida: resumen por pantalla y árbol "mu" (un mu+ de desintegración de K+ por entrada: muP, kPEnd en MeV/c,
// kOrigin, sigMu) para dibujar el espectro (en escala log) en local.
#include "skim_common.h"
#include "TFile.h"
#include "TTree.h"
#include <map>

int skim_muwindow(const char* input, const char* output = "skim_muwindow.root"){
  using namespace skim;
  SkimReader r(input);
  TFile fout(output, "RECREATE");
  Float_t muP, kPEnd; Int_t kOrigin, sigMu, kProcEnd;
  TTree tmu("mu", "mu+ de desintegración de una K+ (no cósmicos)");
  tmu.Branch("muP", &muP); tmu.Branch("kPEnd", &kPEnd); tmu.Branch("kOrigin", &kOrigin); tmu.Branch("sigMu", &sigMu); tmu.Branch("kProcEnd", &kProcEnd);

  // clases del momento final de la K (MeV/c)
  const double kEdges[] = {0, 1e-3, 10, 50, 100, 200, 1e9};
  const char* kName[] = {"en reposo (0)", "0-10", "10-50", "50-100", "100-200", "> 200"};
  const int nK = 6;
  long n[nK] = {0}, inNarrow[nK] = {0}, inWide[nK] = {0}, nSig = 0, nNarrowNotSig = 0;
  for(Long64_t e = 0; r.Next(e); e++){
    const Event& v = r.ev;
    if(!v.isMC)continue;
    for(size_t i = 0; i < v.t_id->size(); i++){
      if((*v.t_pdg)[i] != -13 || (*v.t_procStart)[i] != kDecay || (*v.t_origin)[i] == 2)continue;
      const int k = r.Parent(i);
      if(k < 0 || (*v.t_pdg)[k] != 321 || (*v.t_procEnd)[k] != kDecay)continue;
      muP = 1000.*(*v.t_p)[i]; kPEnd = 1000.*(*v.t_pEnd)[k]; kOrigin = (*v.t_kOrigin)[i]; sigMu = (*v.t_sigMu)[i]; kProcEnd = (*v.t_procEnd)[k];
      tmu.Fill();
      int c = 0; while(c + 1 < nK && !(kPEnd < kEdges[c + 1]))c++;
      if(kPEnd <= 0)c = 0;
      n[c]++;
      const bool narrow = muP > 230 && muP < 241, wide = muP > 210 && muP < 260;
      if(narrow)inNarrow[c]++; if(wide)inWide[c]++;
      if(sigMu)nSig++; else if(narrow)nNarrowNotSig++;
    }
  }
  long tot = 0; for(int c = 0; c < nK; c++)tot += n[c];
  printf("\nmu+ de desintegración de K+ (no cósmicos): %ld; signal muons (IsSignalMuon) %ld\n", tot, nSig);
  printf("  p final de la K (MeV/c)     mu+     en [230, 241]     en [210, 260]\n");
  for(int c = 0; c < nK; c++)
    printf("  %-22s %8ld   %7ld (%5.1f %%)   %7ld (%5.1f %%)\n", kName[c], n[c], inNarrow[c], n[c] ? 100.*inNarrow[c]/n[c] : 0., inWide[c], n[c] ? 100.*inWide[c]/n[c] : 0.);
  printf("  mu+ en [230, 241] que no son signal muon (K en vuelo): %ld\n", nNarrowNotSig);
  fout.cd(); tmu.Write(); fout.Close();
  printf("Árbol mu en %s\n", output);
  return 0;
}
