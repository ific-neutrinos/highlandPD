#ifndef PIDTables_h
#define PIDTables_h

#include <map>
#include <memory>
#include <string>

class TProfile;
class TGraph;
class TSpline3;

/// Tables that depend on the particle hypothesis, shared through
/// pd::services().PID():
///  - dE/dx vs residual range templates per PDG (13, 211, 321, 2212), for
///    the chi2 PID (file pdUtils.PID.dEdxTemplatesFile).
///  - kinetic energy vs range per particle ("muon", "pion", "kaon",
///    "proton"), for the dE/dx likelihoods (file pdUtils.PID.KEvsRangeFile).
///  - NIST range <-> kinetic energy tables of muons and protons in LAr, for
///    the range momentum and the CSDA range (in the code).
/// Everything is loaded once, when the service is created.
class PIDTables{
public:

  PIDTables();
  ~PIDTables();

  /// dE/dx vs residual range template, nullptr if there is none for this PDG
  const TProfile* GetdEdxTemplate(int pdg) const;

  /// Kinetic energy vs range: "muon", "pion", "kaon" or "proton"
  const TGraph& GetKEvsRange(const std::string& particle) const;

  /// NIST tables: muon kinetic energy vs range, and range vs kinetic energy
  /// of muons and protons (MeV, cm)
  const TSpline3& GetMuonKEvsRange()    const {return *_muonKEvsR;}
  const TSpline3& GetMuonRangeVsKE()    const {return *_muonRvsKE;}
  const TSpline3& GetProtonRangeVsKE()  const {return *_protonRvsKE;}

private:

  PIDTables(const PIDTables&) = delete;
  PIDTables& operator=(const PIDTables&) = delete;

  std::map<int, TProfile*> _dEdxTemplates;
  std::map<std::string, TGraph*> _KEvsRange;

  std::unique_ptr<TGraph>   _muonKEvsRGraph, _muonRvsKEGraph, _protonRvsKEGraph;
  std::unique_ptr<TSpline3> _muonKEvsR, _muonRvsKE, _protonRvsKE;
};

#endif
