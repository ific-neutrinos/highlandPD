#include "PIDTables.hxx"
#include "pdDataFiles.hxx"
#include "Parameters.hxx"

#include "TFile.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TSpline.h"

#include <cstdlib>
#include <iostream>

namespace{

//data for range-momentum conversion, muons
//http://pdg.lbl.gov/2012/AtomicNuclearProperties/MUON_ELOSS_TABLES/muonloss_289.pdf divided by LAr density for cm

float Range_grampercm_new[29] = {
  9.833E-1/1.396, 1.786E0/1.396, 3.321E0/1.396, 6.598E0/1.396, 1.058E1/1.396, 3.084E1/1.396, 4.250E1/1.396, 6.732E1/1.396,
  1.063E2/1.396,  1.725E2/1.396, 2.385E2/1.396, 4.934E2/1.396, 6.163E2/1.396, 8.552E2/1.396, 1.202E3/1.396, 1.758E3/1.396,
  2.297E3/1.396,  4.359E3/1.396, 5.354E3/1.396, 7.298E3/1.396, 1.013E4/1.396, 1.469E4/1.396, 1.910E4/1.396, 3.558E4/1.396,
  4.326E4/1.396,  5.768E4/1.396, 7.734E4/1.396, 1.060E5/1.396, 1.307E5/1.396};

float KE_MeV_new[29]= {
  10,    14,    20,    30,    40,     80,     100,    140,    200,   300,
  400,   800,   1000,  1400,  2000,   3000,   4000,   8000,   10000, 14000,
  20000, 30000, 40000, 80000, 100000, 140000, 200000, 300000, 400000};

//data for range-momentum conversion, protons

double Range_gpercm_P_Nist[31]={
  1.887E-1/1.396, 3.823E-1/1.396, 6.335E-1/1.396, 1.296/1.396,   2.159/1.396,   7.375/1.396,   1.092E1/1.396, 2.215E1/1.396,
  3.627E1/1.396,  5.282E1/1.396,  7.144E1/1.396,  9.184E1/1.396, 1.138E2/1.396, 1.370E2/1.396, 1.614E2/1.396, 1.869E2/1.396,
  2.132E2/1.396,  2.403E2/1.396,  2.681E2/1.396,  2.965E2/1.396, 3.254E2/1.396, 3.548E2/1.396, 3.846E2/1.396, 4.148E2/1.396,
  4.454E2/1.396,  7.626E2/1.396,  1.090E3/1.396,  1.418E3/1.396, 1.745E3/1.396, 2.391E3/1.396, 3.022E3/1.396};

double KE_MeV_P_Nist[31]={
  10,   15,   20,   30,   40,   80,   100, 150,
  200,  250,  300,  350,  400,  450,  500, 550,
  600,  650,  700,  750,  800,  850,  900, 950,
  1000, 1500, 2000, 2500, 3000, 4000, 5000};

  // Get an object from a file, detached from it, or exit if it is missing
  template<class T> T* GetFromFile(TFile& file, const char* name){
    T* obj = dynamic_cast<T*>(file.Get(name));
    if(!obj){
      std::cerr << "PIDTables. " << name << " not found in " << file.GetName() << std::endl;
      std::exit(1);
    }
    return obj;
  }
}

//********************************************************************
PIDTables::PIDTables(){
//********************************************************************

  // dE/dx templates
  {
    std::string fname = pd::DataFile("pdUtils", ND::params().GetParameterS("pdUtils.PID.dEdxTemplatesFile"));
    TFile file(fname.c_str(), "READ");
    const std::pair<int,const char*> names[] = {{13,"dedx_range_mu"}, {211,"dedx_range_pi"}, {321,"dedx_range_ka"}, {2212,"dedx_range_pro"}};
    for(const auto& n : names){
      TProfile* p = GetFromFile<TProfile>(file, n.second);
      p->SetDirectory(nullptr);
      _dEdxTemplates[n.first] = p;
    }
  }

  // kinetic energy vs range
  {
    std::string fname = pd::DataFile("pdUtils", ND::params().GetParameterS("pdUtils.PID.KEvsRangeFile"));
    TFile file(fname.c_str(), "READ");
    for(const char* particle : {"muon", "pion", "kaon", "proton"})
      _KEvsRange[particle] = GetFromFile<TGraph>(file, particle);
  }

  // NIST tables (same construction as the former globals in pdAnalysisUtils)
  _muonKEvsRGraph  .reset(new TGraph(29, Range_grampercm_new, KE_MeV_new));
  _muonKEvsR       .reset(new TSpline3("KEvsRS", _muonKEvsRGraph.get()));
  _muonRvsKEGraph  .reset(new TGraph(29, KE_MeV_new, Range_grampercm_new));
  _muonRvsKE       .reset(new TSpline3("RvsKES", _muonRvsKEGraph.get()));
  _protonRvsKEGraph.reset(new TGraph(31, KE_MeV_P_Nist, Range_gpercm_P_Nist));
  _protonRvsKE     .reset(new TSpline3("RvsKE_P_S", _protonRvsKEGraph.get()));
}

//********************************************************************
PIDTables::~PIDTables(){
//********************************************************************

  for(auto& t : _dEdxTemplates) delete t.second;
  for(auto& g : _KEvsRange)     delete g.second;
}

//********************************************************************
const TProfile* PIDTables::GetdEdxTemplate(int pdg) const{
//********************************************************************

  auto it = _dEdxTemplates.find(pdg);
  return (it == _dEdxTemplates.end()) ? nullptr : it->second;
}

//********************************************************************
const TGraph& PIDTables::GetKEvsRange(const std::string& particle) const{
//********************************************************************

  auto it = _KEvsRange.find(particle);
  if(it == _KEvsRange.end()){
    std::cerr << "PIDTables::GetKEvsRange(). No kinetic energy vs range for " << particle << std::endl;
    std::exit(1);
  }
  return *it->second;
}
