// Implementación de beam_axis_minitrees.C (se compila con ACLiC).
#include "pdDataClasses.hxx"
#include "minitree_loop.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TVector3.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace beamaxis {

  // Eje del fcl (CreateHighlandMiniTree.fcl): [0] datos, [1] MC
  const TVector3 kFclStart[2] = {TVector3(-28.6, 424.6, 0.96), TVector3(-30.1, 423.0, 0.0)};
  const TVector3 kFclDir[2]   = {TVector3(-0.18, -0.20, 0.96).Unit(), TVector3(-0.18, -0.20, 0.96).Unit()};

  // Corte del primario parecido al de calidad del haz
  const double kMaxStartZ = 100.;   // cm
  const double kMinCos    = 0.9;    // con la dirección de la instrumentación

  std::vector<std::string> InputFiles(const std::string& input){
    std::vector<std::string> files;
    if(TString(input).EndsWith(".txt")){
      std::ifstream in(input);
      std::string line;
      while(std::getline(in, line)){
        line = TString(line).Strip(TString::kBoth).Data();
        if(!line.empty() && line[0] != '#')files.push_back(line);
      }
    }
    else files.push_back(input);
    return files;
  }

  std::string TreeName(const std::string& file){
    if(file.empty())return "";
    std::unique_ptr<TFile> f(TFile::Open(file.c_str()));
    if(!f || f->IsZombie())return "";
    if(f->Get("highlandana/MiniTree"))return "highlandana/MiniTree";
    if(f->Get("MiniTree"))return "MiniTree";
    return "";
  }

  bool Valid(const Float_t* v){ return v[0] > -900 && v[1] > -900 && v[2] > -900; }

  double Quantile(std::vector<double> v, double f){
    if(v.empty())return -999;
    std::sort(v.begin(), v.end());
    return v[std::min(v.size()-1, size_t(f*v.size()))];
  }

  // Distancia de un punto a la recta del eje
  double DistToAxis(const TVector3& p, const TVector3& start, const TVector3& dir){
    return (p - start).Cross(dir).Mag();
  }

  // Una medida del eje: punto de cruce con z = 0, dirección y, si hay, distancia al eje del fcl
  struct Axis {
    std::vector<double> x0, y0, ux, uy, uz, angle, dist, zref;
    void Add(const TVector3& p, const TVector3& d, int isMC, bool withDist){
      const TVector3 u = d.Unit();
      if(std::fabs(u.Z()) < 1e-3)return;
      const TVector3 c = p + u*((0. - p.Z())/u.Z());   // cruce con z = 0
      x0.push_back(c.X()); y0.push_back(c.Y());
      ux.push_back(u.X()); uy.push_back(u.Y()); uz.push_back(u.Z());
      angle.push_back(u.Angle(kFclDir[isMC])*180./M_PI);
      zref.push_back(p.Z());
      if(withDist)dist.push_back(DistToAxis(p, kFclStart[isMC], kFclDir[isMC]));
    }
    static std::string Q(const std::vector<double>& v, const char* fmt){
      return std::string(Form(fmt, Quantile(v,.16))) + " / " + Form(fmt, Quantile(v,.5)) + " / " + Form(fmt, Quantile(v,.84));
    }
    void Print(const std::string& title) const {
      std::cout << "  " << title << ": " << x0.size() << std::endl;
      if(x0.empty())return;
      std::cout << "    z del punto medido (cm)            " << Q(zref, "%.1f") << std::endl;
      std::cout << "    x en z = 0 (cm)                    " << Q(x0, "%.1f") << std::endl;
      std::cout << "    y en z = 0 (cm)                    " << Q(y0, "%.1f") << std::endl;
      std::cout << "    dirección x                        " << Q(ux, "%.3f") << std::endl;
      std::cout << "    dirección y                        " << Q(uy, "%.3f") << std::endl;
      std::cout << "    dirección z                        " << Q(uz, "%.3f") << std::endl;
      // dirección de las medianas, normalizada: la que iría en el fcl
      const TVector3 m(Quantile(ux,.5), Quantile(uy,.5), Quantile(uz,.5));
      std::cout << Form("    dirección mediana normalizada      [%.3f, %.3f, %.3f]", m.Unit().X(), m.Unit().Y(), m.Unit().Z()) << std::endl;
      std::cout << "    ángulo con la dirección del fcl (°) " << Q(angle, "%.2f") << std::endl;
      if(!dist.empty())
        std::cout << "    inicio: distancia al eje del fcl (cm) " << Q(dist, "%.1f")
                  << Form("   (> 50 cm: %.2f %%)", 100.*std::count_if(dist.begin(), dist.end(), [](double d){ return d > 50; })/dist.size()) << std::endl;
    }
  };

  // Medidas para un momento nominal y datos o MC
  // Distancia al eje del fcl del final de los primarios = haz (MC) o con el corte, por tramos de z
  // del final: si la reco (sin corregir el SCE) se separa del eje a lo largo de la traza
  const double kZBins[] = {0, 50, 100, 200, 300, 450, 700};
  const int kNZ = 6;

  struct Group {
    long nEvents = 0;
    Axis instr, primBeam, primCut, truth;
    std::vector<double> endDist[kNZ];
  };

  struct BeamAxis {
    std::map<std::pair<int,int>, Group> groups;   // (momento nominal, MC)

    void Process(AnaSpillPD* spill){
      const AnaEventInfoPD* info = static_cast<AnaEventInfoPD*>(spill->EventInfo);
      if(!info || spill->Bunches.empty())return;
      const int isMC = info->IsMC ? 1 : 0;
      Group& g = groups[{int(std::lround(info->NominalBeamMom)), isMC}];
      g.nEvents++;

      const AnaBeamPD* beam = static_cast<AnaBeamPD*>(spill->Beam);
      const AnaParticleB* bp = (beam && beam->BeamParticle) ? beam->BeamParticle : nullptr;
      TVector3 instrDir;
      bool hasInstr = false;
      if(bp && Valid(bp->PositionEnd) && Valid(bp->DirectionEnd)){
        instrDir = TVector3(bp->DirectionEnd[0], bp->DirectionEnd[1], bp->DirectionEnd[2]);
        if(instrDir.Mag() > 0){
          hasInstr = true;
          g.instr.Add(TVector3(bp->PositionEnd[0], bp->PositionEnd[1], bp->PositionEnd[2]), instrDir, isMC, false);
        }
      }

      // primario de la beam slice: el módulo guarda primero el primario
      const AnaParticlePD* primary = nullptr;
      for(AnaParticleB* pb : static_cast<AnaBunchPD*>(spill->Bunches[0])->Particles){
        AnaParticlePD* p = static_cast<AnaParticlePD*>(pb);
        if(p->BeamOrigin){ primary = p; break; }
      }
      if(primary && Valid(primary->PositionStart) && Valid(primary->DirectionStart)){
        const TVector3 ps(primary->PositionStart[0], primary->PositionStart[1], primary->PositionStart[2]);
        const TVector3 ds(primary->DirectionStart[0], primary->DirectionStart[1], primary->DirectionStart[2]);
        if(ds.Mag() > 0){
          const AnaTrueParticleB* bt = (isMC && bp) ? static_cast<AnaTrueParticleB*>(bp->TrueObject) : nullptr;
          const AnaTrueParticleB* pt = static_cast<AnaTrueParticleB*>(primary->TrueObject);
          const bool isBeam = bt && pt && pt->ID == bt->ID;
          const bool passCut = hasInstr && ps.Z() < kMaxStartZ && ds.Unit().Dot(instrDir.Unit()) > kMinCos;
          if(isBeam)g.primBeam.Add(ps, ds, isMC, true);
          if(passCut)g.primCut.Add(ps, ds, isMC, true);
          if((isMC ? isBeam : passCut) && Valid(primary->PositionEnd)){
            const TVector3 pe(primary->PositionEnd[0], primary->PositionEnd[1], primary->PositionEnd[2]);
            for(int iz = 0; iz < kNZ; iz++)if(pe.Z() >= kZBins[iz] && pe.Z() < kZBins[iz+1])
              g.endDist[iz].push_back(DistToAxis(pe, kFclStart[isMC], kFclDir[isMC]));
          }
        }
      }

      // trayectoria true del haz al cruzar z = 0 (interpolando entre los dos puntos de alrededor)
      if(isMC && bp && bp->TrueObject){
        const AnaTrueParticlePD* bt = static_cast<const AnaTrueParticlePD*>(bp->TrueObject);
        const std::vector<AnaTrueTrajPointPD>& tr = bt->TrjPoints;
        for(size_t i = 1; i < tr.size(); i++){
          const double z1 = tr[i-1].Position[2], z2 = tr[i].Position[2];
          if(!(z1 < 0 && z2 >= 0))continue;
          const double f = (0. - z1)/(z2 - z1);
          TVector3 p, m;
          for(int k = 0; k < 3; k++){
            p[k] = tr[i-1].Position[k] + f*(tr[i].Position[k] - tr[i-1].Position[k]);
            m[k] = tr[i-1].Momentum[k] + f*(tr[i].Momentum[k] - tr[i-1].Momentum[k]);
          }
          if(m.Mag() > 0)g.truth.Add(p, m, isMC, false);
          break;
        }
      }
    }

    void Print() const {
      for(const auto& kg : groups){
        const int isMC = kg.first.second;
        const Group& g = kg.second;
        std::cout << std::endl << "===== " << (isMC ? "MC" : "Datos") << ", momento nominal " << kg.first.first
                  << " GeV: " << g.nEvents << " eventos con bunch =====" << std::endl;
        std::cout << Form("Eje del fcl: inicio [%.2f, %.2f, %.2f], dirección [%.3f, %.3f, %.3f]",
                          kFclStart[isMC].X(), kFclStart[isMC].Y(), kFclStart[isMC].Z(),
                          kFclDir[isMC].X(), kFclDir[isMC].Y(), kFclDir[isMC].Z()) << std::endl;
        std::cout << "Cuantiles 16 / 50 / 84 %:" << std::endl;
        g.instr.Print("traza de la instrumentación del haz (extremo final)");
        if(isMC)g.primBeam.Print("primario de la beam slice = partícula true del haz (inicio)");
        g.primCut.Print(Form("primario de la beam slice con z inicial < %.0f cm y cos con la instrumentación > %.1f (inicio)", kMaxStartZ, kMinCos));
        if(isMC)g.truth.Print("trayectoria true del haz en z = 0");
        std::cout << "  final del primario (" << (isMC ? "= partícula true del haz" : "con el corte")
                  << "): distancia al eje del fcl (cm), cuantiles 16/50/84/95 %, por z del final" << std::endl;
        for(int iz = 0; iz < kNZ; iz++){
          const std::vector<double>& v = g.endDist[iz];
          std::cout << Form("    z en [%3.0f, %3.0f) cm: %6zu   %5.1f / %5.1f / %5.1f / %5.1f", kZBins[iz], kZBins[iz+1], v.size(),
                            Quantile(v,.16), Quantile(v,.5), Quantile(v,.84), Quantile(v,.95)) << std::endl;
        }
      }
    }
  };
}

int beam_axis_minitrees_impl(const char* input, Long64_t maxEntries){
  using namespace beamaxis;
  const std::vector<std::string> files = InputFiles(input);
  if(files.empty()){ std::cerr << "No hay ficheros en " << input << std::endl; return -1; }
  const std::string treeName = TreeName(TString(files[0]).Contains("*") ? std::string("") : files[0]);
  TChain chain(treeName.empty() ? "highlandana/MiniTree" : treeName.c_str());
  for(const std::string& f : files)chain.Add(f.c_str());
  std::cout << "Ficheros: " << files.size() << std::endl;
  if(chain.LoadTree(0) < 0){ std::cerr << "No hay entradas (¿árbol " << chain.GetName() << "?)" << std::endl; return -1; }

  AnaSpillPD* spill = nullptr;
  chain.SetBranchAddress("Spill", &spill);
  BeamAxis ba;
  const Long64_t n = LoopMiniTree(chain, spill, maxEntries, files.size(), [&ba](AnaSpillPD* s){ ba.Process(s); });
  std::cout << "Entradas leídas: " << n << std::endl;
  ba.Print();
  return 0;
}
