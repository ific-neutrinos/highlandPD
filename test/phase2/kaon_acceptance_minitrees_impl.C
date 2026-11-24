// Implementación de kaon_acceptance_minitrees.C (se compila con ACLiC).
//
// Aceptancia física (solo truth, sin reconstrucción) de la señal mu+ de K+ en reposo para la medida A: de las K+
// producidas en la primera interacción inelástica del haz que dispara, ¿cuántas acaban dando un mu+ de K+ -> mu+ nu
// en reposo, contenido? Sin cortes en la K: el mapa en (T_K, theta_K) dice qué parte de la producción se ve.
//
// Por evento (MC con truth del haz):
//   - Haz que dispara B (beam->BeamParticle->TrueObject): especie, dónde y cómo acaba. Denominador: B acaba en una
//     interacción inelástica (no Decay) con el vértice en el volumen fiducial kFV.
//   - K+ de la primera interacción: hijas directas de B (medida A estricta) y, para la opción a, hijas de las
//     partículas de la especie de B que descienden de B solo por partículas de esa especie (el haz que sobrevive a
//     interacciones casi elásticas): nivel (número de pasos) y la fracción mínima de momento conservada en cada paso
//     (p inicial de la hija / p final de la madre), para elegir el umbral r después.
//   - Destino de cada K, siguiendo las K+ -> K+ (la hija K+ de mayor momento de una K+ que acaba en inelástico):
//     se para (MomentumEnd < 50 MeV/c) y da mu+ nu (mu+ de 230-241 MeV/c), se para y da otro canal (pi+ pi0 aparte:
//     pi+ de 195-215 MeV/c), se desintegra en vuelo, interacciona sin dejar K+, o sale del volumen activo.
//   - Para mu+ nu: si el mu+ acaba dentro del volumen activo y a > kMuMargin cm de sus caras (Michel incluido).
// Salida: resumen por pantalla (por especie del haz: denominador, K+ por interacción, destinos, mapa de la fracción
// que da mu+ nu contenido en (T_K, theta_K)) y un árbol plano "kprod" (una entrada por K) en el fichero de salida.
#include "kaon_fate_minitrees_impl.C"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

namespace kacc {
  using kfate::InputFiles; using kfate::TreeName;

  const double kActive[6] = {-360., 360., 0., 607.5, 0., 695.};   // como el módulo
  // Volumen fiducial del vértice de la interacción: a > 20 cm de las caras en x e y, y 0 < z < 600 cm. Sin corte
  // en z > 30 cm (donde empieza la traza del haz reconstruida) hasta evaluar el efecto del SCE con datos
  // (decisión del usuario, 24-nov-2026; antes era 30 < z < 600).
  const double kFV[6] = {-340., 340., 20., 587.5, 0., 600.};
  const double kMuMargin = 20.;          // cm, mu+ contenido con su Michel
  const double kKaonMass = 0.493677;     // GeV
  const TVector3 kAxisDir = TVector3(-0.18, -0.20, 0.96).Unit();   // dirección del haz si falta la true

  bool In(const Float_t* x, const double* box, double margin = 0){
    return x[0] > box[0] + margin && x[0] < box[1] - margin && x[1] > box[2] + margin && x[1] < box[3] - margin && x[2] > box[4] + margin && x[2] < box[5] - margin;
  }
  TVector3 Dir(const Float_t* d){ TVector3 v(d[0], d[1], d[2]); return v.Mag() > 0.5 ? v.Unit() : TVector3(0, 0, 0); }

  enum Fate { kRestMuNu = 0, kRestPiPi0, kRestOther, kFlightDecay, kInelastic, kExits, kUnknown, kNFates };
  const char* kFateName[kNFates] = {"en reposo -> mu nu", "en reposo -> pi+ pi0", "en reposo -> otro", "desintegración en vuelo",
                                    "interacciona sin dejar K+", "sale del volumen activo", "sin información"};

  // Bins del mapa
  const std::vector<double> kTBins = {0, 100, 200, 300, 400, 600, 800, 1200, 1e9};   // MeV
  const std::vector<double> kThBins = {0, 15, 30, 45, 60, 90, 180};                  // grados
  int Bin(const std::vector<double>& b, double x){ for(size_t i = 0; i + 1 < b.size(); i++)if(x >= b[i] && x < b[i+1])return i; return -1; }

  struct Row {   // una entrada del árbol kprod
    Int_t entry, beamPDG, beamEndProc, level, fate, muContained, muInAV, kStopInAV, nKK, vtxInFV;
    Float_t beamPEnd, vtxX, vtxY, vtxZ, kT, kP, kTheta, rMin, kEndX, kEndY, kEndZ, muEndX, muEndY, muEndZ;
  };

  struct Acceptance {
    TFile* fout = nullptr; TTree* tk = nullptr; Row r;
    long nEvents = 0, nMC = 0, nNoBeam = 0;
    std::map<int, long> nBeam, nInel, nInelFV, nDecay, nExit;   // por PDG del haz
    std::map<int, std::map<int, long>> nFate;                   // haz -> destino (nivel 0, vértice en FV)
    std::map<int, long> nK0, nK0FV;                             // K+ de nivel 0 (A estricta), con el vértice en FV
    std::map<int, std::vector<std::vector<long>>> prod, good;   // haz -> [T][theta]: producidas y mu nu contenido (nivel 0, FV)
    long noDaughters = 0;                                       // K+ finales sin hijas guardadas (comprobación de la truth)

    void Open(const char* out){
      fout = TFile::Open(out, "RECREATE");
      tk = new TTree("kprod", "K+ de la primera interacción del haz que dispara (truth)");
      tk->Branch("entry", &r.entry); tk->Branch("beamPDG", &r.beamPDG); tk->Branch("beamEndProc", &r.beamEndProc); tk->Branch("beamPEnd", &r.beamPEnd);
      tk->Branch("vtxX", &r.vtxX); tk->Branch("vtxY", &r.vtxY); tk->Branch("vtxZ", &r.vtxZ); tk->Branch("vtxInFV", &r.vtxInFV);
      tk->Branch("level", &r.level); tk->Branch("rMin", &r.rMin);
      tk->Branch("kT", &r.kT); tk->Branch("kP", &r.kP); tk->Branch("kTheta", &r.kTheta);
      tk->Branch("fate", &r.fate); tk->Branch("nKK", &r.nKK); tk->Branch("kStopInAV", &r.kStopInAV);
      tk->Branch("kEndX", &r.kEndX); tk->Branch("kEndY", &r.kEndY); tk->Branch("kEndZ", &r.kEndZ);
      tk->Branch("muInAV", &r.muInAV); tk->Branch("muContained", &r.muContained);
      tk->Branch("muEndX", &r.muEndX); tk->Branch("muEndY", &r.muEndY); tk->Branch("muEndZ", &r.muEndZ);
    }

    // destino de una K+ siguiendo las K+ -> K+
    void FollowK(const AnaTrueParticleB* k, const std::map<int, AnaTrueParticleB*>& byID){
      r.nKK = 0;
      const AnaTrueParticleB* cur = k;
      for(int depth = 0; depth < 20; depth++){
        if(cur->ProcessEnd == AnaTrueParticleB::Decay)break;
        // K+ hija de mayor momento (K+ -> K+)
        const AnaTrueParticleB* next = nullptr;
        for(int d : static_cast<const AnaTrueParticlePD*>(cur)->Daughters){
          auto it = byID.find(d);
          if(it != byID.end() && it->second->PDG == 321 && (!next || it->second->Momentum > next->Momentum))next = it->second;
        }
        if(!next)break;
        cur = next; r.nKK++;
      }
      r.kEndX = cur->PositionEnd[0]; r.kEndY = cur->PositionEnd[1]; r.kEndZ = cur->PositionEnd[2];
      r.kStopInAV = In(cur->PositionEnd, kActive);
      r.muInAV = r.muContained = 0; r.muEndX = r.muEndY = r.muEndZ = -999;
      const auto& dau = static_cast<const AnaTrueParticlePD*>(cur)->Daughters;
      if(dau.empty() && cur->ProcessEnd == AnaTrueParticleB::Decay)noDaughters++;
      if(!r.kStopInAV){ r.fate = kExits; return; }
      if(cur->ProcessEnd == AnaTrueParticleB::Decay){
        if(cur->MomentumEnd > 1e-6){ r.fate = kFlightDecay; return; }   // en reposo = momento final cero (Geant4), como IsSignalMuon
        const AnaTrueParticleB* mu = nullptr; bool pipi = false;
        for(int d : dau){
          auto it = byID.find(d); if(it == byID.end())continue;
          const AnaTrueParticleB* x = it->second;
          if(x->PDG == -13 && x->Momentum > 0.230 && x->Momentum < 0.241)mu = x;
          if(x->PDG == 211 && x->Momentum > 0.195 && x->Momentum < 0.215)pipi = true;
        }
        if(mu){
          r.fate = kRestMuNu;
          r.muEndX = mu->PositionEnd[0]; r.muEndY = mu->PositionEnd[1]; r.muEndZ = mu->PositionEnd[2];
          r.muInAV = In(mu->PositionEnd, kActive); r.muContained = In(mu->PositionEnd, kActive, kMuMargin);
        }
        else r.fate = pipi ? kRestPiPi0 : kRestOther;
        return;
      }
      if(cur->ProcessEnd == AnaTrueParticleB::Unknown){ r.fate = kUnknown; return; }
      r.fate = kInelastic;
    }

    void Process(AnaSpillPD* spill, Long64_t entry){
      nEvents++;
      const AnaEventInfoPD* info = static_cast<AnaEventInfoPD*>(spill->EventInfo);
      if(!info || !info->IsMC)return;
      nMC++;
      const AnaBeamPD* beam = static_cast<AnaBeamPD*>(spill->Beam);
      const AnaTrueParticleB* B = (beam && beam->BeamParticle) ? static_cast<AnaTrueParticleB*>(beam->BeamParticle->TrueObject) : nullptr;
      if(!B){ nNoBeam++; return; }
      std::map<int, AnaTrueParticleB*> byID;
      for(AnaTrueParticleB* t : spill->TrueParticles)byID[t->ID] = t;
      auto itB = byID.find(B->ID); if(itB != byID.end())B = itB->second;   // la de la lista, con sus hijas
      const int pdg = B->PDG;
      nBeam[pdg]++;
      const bool inAV = In(B->PositionEnd, kActive);
      if(!inAV){ nExit[pdg]++; }
      else if(B->ProcessEnd == AnaTrueParticleB::Decay){ nDecay[pdg]++; }
      else {
        nInel[pdg]++;
        if(In(B->PositionEnd, kFV))nInelFV[pdg]++;
      }
      TVector3 bdir = Dir(B->DirectionEnd); if(bdir.Mag() == 0)bdir = kAxisDir;
      // K+ de la primera interacción (nivel 0) y de la opción a (nivel >= 1), recorriendo las partículas de la especie de B
      struct Node { const AnaTrueParticleB* p; int level; double rmin; };
      std::vector<Node> stack{{B, 0, 1e9}};
      while(!stack.empty()){
        const Node nd = stack.back(); stack.pop_back();
        for(int d : static_cast<const AnaTrueParticlePD*>(nd.p)->Daughters){
          auto it = byID.find(d); if(it == byID.end())continue;
          const AnaTrueParticleB* x = it->second;
          if(x->PDG == pdg && nd.level < 5){
            const double ratio = nd.p->MomentumEnd > 0 ? x->Momentum/nd.p->MomentumEnd : -1;
            stack.push_back({x, nd.level + 1, std::min(nd.rmin, ratio)});
            continue;
          }
          if(x->PDG != 321 || x->ProcessStart == AnaTrueParticleB::Decay)continue;
          if(pdg == 321 && nd.level == 0)continue;   // haz K+: su K+ hija es la K del haz tras K+ -> K+, no producción
          r.entry = entry; r.beamPDG = pdg; r.beamEndProc = B->ProcessEnd; r.beamPEnd = B->MomentumEnd;
          r.vtxX = x->Position[0]; r.vtxY = x->Position[1]; r.vtxZ = x->Position[2]; r.vtxInFV = In(x->Position, kFV);
          r.level = nd.level; r.rMin = nd.level ? nd.rmin : 1;
          r.kP = x->Momentum; r.kT = 1000.*(std::sqrt(x->Momentum*x->Momentum + kKaonMass*kKaonMass) - kKaonMass);
          const TVector3 kd = Dir(x->Direction);
          const TVector3 ref = nd.level ? (Dir(nd.p->DirectionEnd).Mag() > 0 ? Dir(nd.p->DirectionEnd) : bdir) : bdir;
          r.kTheta = kd.Mag() > 0 ? kd.Angle(ref)*180./M_PI : -1;
          FollowK(x, byID);
          tk->Fill();
          if(nd.level == 0){
            nK0[pdg]++;
            if(r.vtxInFV){
              nK0FV[pdg]++; nFate[pdg][r.fate]++;
              auto& P = prod[pdg]; auto& G = good[pdg];
              if(P.empty()){ P.assign(kTBins.size() - 1, std::vector<long>(kThBins.size() - 1, 0)); G = P; }
              const int it2 = Bin(kTBins, r.kT), ith = Bin(kThBins, r.kTheta);
              if(it2 >= 0 && ith >= 0){ P[it2][ith]++; if(r.fate == kRestMuNu && r.muContained)G[it2][ith]++; }
            }
          }
        }
      }
    }

    void Print() const {
      std::cout << std::endl << "Entradas " << nEvents << " (MC " << nMC << ", sin truth del haz " << nNoBeam << ")" << std::endl;
      std::cout << "Volumen fiducial del vértice: x [" << kFV[0] << ", " << kFV[1] << "], y [" << kFV[2] << ", " << kFV[3] << "], z [" << kFV[4] << ", " << kFV[5]
                << "] cm; mu contenido: a > " << kMuMargin << " cm de las caras" << std::endl;
      for(auto& kv : nBeam){
        const int pdg = kv.first;
        auto g = [&](const std::map<int, long>& m){ auto it = m.find(pdg); return it == m.end() ? 0L : it->second; };
        const long inelFV = g(nInelFV), k0FV = g(nK0FV);
        std::cout << std::endl << Form("== Haz PDG %d: %ld eventos; acaba fuera del volumen %ld, en desintegración %ld, inelástico %ld (en el FV %ld)",
                                       pdg, kv.second, g(nExit), g(nDecay), g(nInel), inelFV) << std::endl;
        if(!k0FV)continue;
        std::cout << Form("   K+ de la primera interacción con el vértice en el FV: %ld (%.4f por interacción inelástica en el FV)", k0FV, inelFV ? double(k0FV)/inelFV : 0.) << std::endl;
        std::cout << "   destino:" << std::endl;
        auto itf = nFate.find(pdg);
        for(int f = 0; f < kNFates; f++){
          long n = 0; if(itf != nFate.end()){ auto x = itf->second.find(f); if(x != itf->second.end())n = x->second; }
          std::cout << Form("      %-28s %6ld (%5.1f %%)", kFateName[f], n, 100.*n/k0FV) << std::endl;
        }
        auto P = prod.at(pdg), G = good.at(pdg);
        std::cout << "   fracción que da mu+ nu en reposo con el mu contenido, por T_K (filas, MeV) y theta_K (columnas, grados); entre paréntesis, producidas:" << std::endl;
        std::cout << "      T \\ theta   ";
        for(size_t j = 0; j + 1 < kThBins.size(); j++)std::cout << Form("   %3.0f-%-3.0f      ", kThBins[j], kThBins[j+1]);
        std::cout << "    todas" << std::endl;
        std::vector<long> colP(kThBins.size() - 1, 0), colG(kThBins.size() - 1, 0);
        for(size_t i = 0; i + 1 < kTBins.size(); i++){
          long rp = 0, rg = 0;
          std::cout << Form("      %4.0f-%-5s", kTBins[i], i + 2 < kTBins.size() ? Form("%.0f", kTBins[i+1]) : "");
          for(size_t j = 0; j + 1 < kThBins.size(); j++){
            rp += P[i][j]; rg += G[i][j]; colP[j] += P[i][j]; colG[j] += G[i][j];
            std::cout << (P[i][j] ? Form("  %5.3f (%5ld)", double(G[i][j])/P[i][j], P[i][j]) : Form("      - (%5ld)", 0L));
          }
          std::cout << (rp ? Form("   %5.3f (%5ld)", double(rg)/rp, rp) : "") << std::endl;
        }
        long tp = 0, tg = 0;
        std::cout << "      todas      ";
        for(size_t j = 0; j + 1 < kThBins.size(); j++){ tp += colP[j]; tg += colG[j]; std::cout << (colP[j] ? Form("  %5.3f (%5ld)", double(colG[j])/colP[j], colP[j]) : Form("      - (%5ld)", 0L)); }
        std::cout << (tp ? Form("   %5.3f (%5ld)", double(tg)/tp, tp) : "") << std::endl;
      }
      std::cout << std::endl << "Comprobación de la truth: K+ que acaban en desintegración sin hijas guardadas: " << noDaughters << std::endl;
      std::cout << "Opción a (nivel >= 1) y el resto, en el árbol kprod (level, rMin)" << std::endl;
    }

    void Close(){ if(fout){ fout->cd(); tk->Write(); WriteCodeVersion(fout, "kaon_acceptance_minitrees"); fout->Close(); } }
  };
}

int kaon_acceptance_minitrees_impl(const char* input, const char* output, Long64_t maxEntries){
  using namespace kacc;
  const std::vector<std::string> files = InputFiles(input);
  if(files.empty()){ std::cerr << "No hay ficheros en " << input << std::endl; return -1; }
  const std::string treeName = TreeName(TString(files[0]).Contains("*") ? std::string("") : files[0]);
  TChain chain(treeName.empty() ? "highlandana/MiniTree" : treeName.c_str());
  for(const std::string& f : files)chain.Add(f.c_str());
  if(chain.LoadTree(0) < 0){ std::cerr << "No hay entradas (¿árbol " << chain.GetName() << "?)" << std::endl; return -1; }
  AnaSpillPD* spill = nullptr;
  chain.SetBranchAddress("Spill", &spill);
  Acceptance a;
  a.Open(output);
  Long64_t entry = 0;
  LoopMiniTree(chain, spill, maxEntries, files.size(), [&a, &entry](AnaSpillPD* s){ a.Process(s, entry++); });
  std::cout << std::endl << "Ficheros: " << files.size() << std::endl;
  a.Print();
  a.Close();
  std::cout << "Árbol kprod en " << output << std::endl;
  return 0;
}
