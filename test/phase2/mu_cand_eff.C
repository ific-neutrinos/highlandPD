// Reparto de la eficiencia de la selección de mu+ de K+ en reposo: de todos los mu+ true (árbol "sig" de
// mu_candidates_minitrees.C, K en reposo) a los que pasan la cadena de cortes de su población
// (InSliceSteps / OutSliceSteps de mu_cand_common.h, sobre el árbol "cand").
//
//   root -l -b -q 'test/phase2/mu_cand_eff.C+("mu_cand_big19.root")'
//
// Cada mu+ se enlaza con su candidato por la entrada del minitree y el p de rango del PFParticle emparejado
// (el de más hits). Pasos: contenido, con PFParticle, traza, completo (con hits), p de rango en la ventana
// del módulo (si no: partido, < 180, o unido a otra cosa, > 300), cortes geométricos del módulo y cada corte
// de la cadena. Columnas: todos, eventos con beam slice (mu en ella o fuera) y eventos sin beam slice.
// Señal como en mu_candidates_minitrees.C (mu+ de una K+ que desciende del haz que dispara; incluye las
// K+ de interacciones secundarias y la K+ del haz tras K+ -> K+, ver kaon_origin_minitrees.C).
#include "mu_cand_common.h"
#include <map>

int mu_cand_eff(const char* in){
  using namespace mucut;
  const std::vector<C> cands = Load(in);
  if(cands.empty())return 1;
  // candidatos de señal por entrada
  std::multimap<int, const C*> sigCand;
  for(const C& c : cands)if(c.group == 0)sigCand.insert({c.entry, &c});

  TFile f(in);
  TTree* ts = (TTree*)f.Get("sig");
  if(!ts){ std::cerr << "No hay árbol sig" << std::endl; return 1; }
  Int_t entry, hasSlice, matched, isTrack, inSlice, complete, contained; Float_t kMomEnd, pRange;
  ts->SetBranchAddress("entry", &entry); ts->SetBranchAddress("hasSlice", &hasSlice); ts->SetBranchAddress("matched", &matched);
  ts->SetBranchAddress("isTrack", &isTrack); ts->SetBranchAddress("inSlice", &inSlice); ts->SetBranchAddress("complete", &complete);
  ts->SetBranchAddress("contained", &contained); ts->SetBranchAddress("kMomEnd", &kMomEnd); ts->SetBranchAddress("pRange", &pRange);
  // desde el 21-nov el árbol sig tiene también los mu+ de K en reposo que no son la señal de antes (isSignalOld)
  Int_t isSignalOld = 1;
  if(ts->GetBranch("isSignalOld"))ts->SetBranchAddress("isSignalOld", &isSignalOld);

  const std::vector<Step> stIn = InSliceSteps(), stOut = OutSliceSteps();
  // filas: nombre; columnas: 0 todos, 1 evento con BS y mu en BS, 2 evento con BS y mu fuera, 3 evento sin BS
  std::vector<std::string> rows;
  std::map<std::string, std::array<double, 4>> n;
  auto Count = [&](const std::string& r, int col){
    if(!n.count(r)){ rows.push_back(r); n[r] = {0, 0, 0, 0}; }
    n[r][0]++; n[r][col]++;
  };
  long notFound = 0;
  for(Long64_t i = 0; i < ts->GetEntries(); i++){
    ts->GetEntry(i);
    if(!isSignalOld || kMomEnd >= 50)continue;
    // columna: hasta tener PFParticle se usa solo hasSlice (col 1 = con BS mientras no se sepa dónde está el mu)
    int col = hasSlice ? 1 : 3;
    Count("mu+ de K+ en reposo", col);
    if(!contained)continue;
    Count("contenido en el volumen activo", col);
    if(!matched)continue;
    if(hasSlice)col = inSlice ? 1 : 2;
    Count("con PFParticle emparejado", col);
    if(!isTrack)continue;
    Count("es una traza", col);
    if(!complete)continue;
    Count("completo (con hits)", col);
    if(pRange < 100){ Count("  (p de rango < 100: muy partido)", col); continue; }
    if(pRange < 180){ Count("  (p de rango 100-180: partido)", col); continue; }
    if(pRange > 300){ Count("  (p de rango > 300: unido a otra cosa)", col); continue; }
    Count("p de rango en [180, 300]", col);
    // candidato en cand (los cortes geométricos del módulo ya están en Load)
    const C* c = nullptr;
    auto range = sigCand.equal_range(entry);
    for(auto it = range.first; it != range.second; ++it)if(std::abs(it->second->pRange - pRange) < 0.01)c = it->second;
    if(!c){ notFound++; continue; }
    Count("cortes geométricos del módulo", col);
    const std::vector<Step>& st = c->inSlice ? stIn : stOut;
    bool ok = true;
    for(size_t k = 0; k < st.size() && ok; k++){
      if(!st[k].pass(*c)){ ok = false; break; }
      Count(std::string(c->inSlice ? "[BS] " : "[fuera] ") + st[k].name, col);
    }
  }
  const auto& tot = n["mu+ de K+ en reposo"];
  std::cout << std::endl << "Eficiencia de la selección (mu+ de K+ en reposo); % sobre el total de cada columna" << std::endl
            << Form("  %-46s %14s %22s %22s %18s", "", "todos", "con BS, mu en BS", "con BS, mu fuera", "sin BS") << std::endl;
  for(const std::string& r : rows){
    std::cout << Form("  %-46s", r.c_str());
    for(int k = 0; k < 4; k++){
      const double den = k == 2 ? tot[1] : tot[k];   // con BS: ambas columnas sobre los eventos con BS
      std::cout << Form(" %7.0f (%5.1f %%)", n[r][k], den > 0 ? 100*n[r][k]/den : 0.);
    }
    std::cout << std::endl;
  }
  std::cout << "(las columnas 'con BS' se dan sobre todos los mu de eventos con beam slice; candidatos no encontrados en cand: " << notFound << ")" << std::endl;
  return 0;
}
