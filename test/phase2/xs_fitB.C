// Medida B: ajuste del pico de p de rango del mu de K+ en reposo por especie del haz y población, con la salida de
// secondaryKaonProductionXS (árbol ana), 23-nov-2026. Mide la señal sin la truth: el fondo lo da el ajuste.
//
//   root -l -b -q 'test/phase2/xs_fitB.C+("xs_v100.root", "<directorio de figuras>")'
//
// Espectro: candidatos con todos los cortes de su cadena salvo la ventana de p de rango. Desde el 24-nov la ventana
// es el último corte y el microtree guarda el nivel de esa muestra (sel_fitlevel). En los microtrees anteriores la
// ventana iba antes del Michel: se toma "llega hasta la ventana" (nivel >= 6 dentro, >= 8 fuera) y el corte de
// Michel por hits se aplica aquí con sus variables (como MuCandMichelHitsCut). Fondo: con la vista normal (microtrees
// con cand_savemode), de todos los eventos sin escalar; antes, de los eventos normales escalado por especie.
// Modelo y pseudoexperimentos: los de mu_cand_fit.C (gaussiana + exp(pol2), máxima verosimilitud de Poisson). El
// "dato" es el MC: mu de K en reposo de todos los eventos más el fondo de los eventos normales, escalado por especie.
// El pico contiene todos los mu de K en reposo, también los del pileup: el resultado se compara con dos verdades, todos
// los mu de K (lo que mide el ajuste) y la señal B (cadena del haz que dispara, origen 0-3). La diferencia es el pileup.
// Especie: la true del haz (beam_truespecies o, en ficheros anteriores, beam_truepdg), agrupada como los Cherenkov.
#include "mu_cand_fit.C"
#include <map>

namespace fitB {
  const int NMAX = 300;
  const char* kSpName[3] = {"e/mu/pi", "K", "p"};
  int SpeciesFromPDG(int pdg){ const int a = std::abs(pdg); return (a == 11 || a == 13 || a == 211) ? 0 : a == 321 ? 1 : a == 2212 ? 2 : -1; }
  const char* Br(TTree* t, const char* nuevo, const char* viejo){ return t->GetBranch(nuevo) ? nuevo : viejo; }
}

int xs_fitB(const char* in, const char* outdir, int nToys = 300){
  using namespace fitB; using namespace mufit;
  TFile f(in);
  TTree* t = (TTree*)f.Get("ana");
  if(!t){ std::cerr << "No hay árbol ana" << std::endl; return 1; }
  Int_t ncand, fulldump, beamTruePDG = 0, beamTrueSp = -1, pop[NMAX], level[NMAX], complete[NMAX], nmich[NMAX];
  Int_t tpdg[NMAX], tproc[NMAX], torig[NMAX], tpar[NMAX], tkor[NMAX];
  Float_t tmom[NMAX], prange[NMAX], cnnmich[NMAX];
  t->SetBranchStatus("*", 0);
  auto B = [&](const char* n, void* a){ t->SetBranchStatus(n, 1); t->SetBranchAddress(n, a); };
  B("ncand", &ncand); B("fulldump", &fulldump); B("cand_pop", pop); B("cand_level", level); B("cand_complete", complete);
  B("cand_prange", prange); B("cand_ar_nmichel10", nmich); B("cand_cnn_michel", cnnmich);
  const bool hasSp = t->GetBranch("beam_truespecies");
  if(hasSp)B("beam_truespecies", &beamTrueSp); else B("beam_truepdg", &beamTruePDG);
  B(Br(t, "cand_truepdg", "cand_tpdg"), tpdg); B(Br(t, "cand_trueproc", "cand_tprocstart"), tproc); B(Br(t, "cand_trueorigin", "cand_torigin"), torig);
  B(Br(t, "cand_trueparentpdg", "cand_tparentpdg"), tpar); B(Br(t, "cand_truekorigin", "cand_tkorigin"), tkor); B(Br(t, "cand_truemom", "cand_tmom"), tmom);

  // Escala del fondo POR ESPECIE: el volcado completo depende de la especie (con haz de K solo ~40 % de los eventos son
  // normales; con pi ~91 %, con p ~98 %), así que el fondo de los eventos normales se escala a todos los de su especie
  const Long64_t n = t->GetEntries();
  double nAll[3] = {0}, nNorm[3] = {0};
  for(Long64_t i = 0; i < n; i++){
    t->GetEntry(i);
    const int sp = hasSp ? beamTrueSp : SpeciesFromPDG(beamTruePDG);
    if(sp < 0 || sp > 2)continue;
    nAll[sp]++; if(!fulldump)nNorm[sp]++;
  }
  const bool normalView = t->GetBranch("cand_savemode");
  double scSp[3]; for(int sp = 0; sp < 3; sp++)scSp[sp] = (normalView || !nNorm[sp]) ? 1. : nAll[sp]/nNorm[sp];
  // muestra del ajuste: sel_fitlevel si existe (ventana la última, Michel ya aplicado); si no, la lógica de antes
  const bool newOrder = t->GetBranch("sel_fitlevel");
  int fitLevel[2] = {6, 8};
  if(newOrder){ Int_t fl[2]; t->SetBranchStatus("sel_fitlevel", 1); t->SetBranchAddress("sel_fitlevel", fl); t->GetEntry(0);
    fitLevel[0] = fl[0]; fitLevel[1] = fl[1]; t->ResetBranchAddress(t->GetBranch("sel_fitlevel")); }

  // rango del ajuste por población: dentro 150-400; fuera 150-350 (fuera, el módulo solo guarda completos los
  // candidatos en su ventana de 150-350 MeV/c)
  const double x0[2] = {150, 150}, x1[2] = {400, 350}, bw = 5;
  const char* popName[2] = {"en la beam slice", "fuera de la beam slice"};
  // [especie][población]: mu de K (todos), señal B, fondo que no es de K
  std::map<int, TH1D*> hK, hS, hB;
  auto key = [](int sp, int p){ return 10*sp + p; };
  for(int sp = 0; sp < 3; sp++)for(int p = 0; p < 2; p++){
    const int nb = int((x1[p] - x0[p])/bw + 0.5);
    hK[key(sp, p)] = new TH1D(Form("hK_%d_%d", sp, p), "", nb, x0[p], x1[p]);
    hS[key(sp, p)] = new TH1D(Form("hS_%d_%d", sp, p), "", nb, x0[p], x1[p]);
    hB[key(sp, p)] = new TH1D(Form("hB_%d_%d", sp, p), "", nb, x0[p], x1[p]);
  }
  for(Long64_t i = 0; i < n; i++){
    t->GetEntry(i);
    const int sp = hasSp ? beamTrueSp : SpeciesFromPDG(beamTruePDG);
    if(sp < 0 || sp > 2)continue;
    for(int c = 0; c < ncand && c < NMAX; c++){
      const int p = pop[c];
      if(!complete[c] || level[c] < fitLevel[p])continue;
      const bool michel = nmich[c] < -900 ? cnnmich[c] > 0.02 : nmich[c] >= 2;
      if(!newOrder && !michel)continue;
      if(prange[c] < x0[p] || prange[c] >= x1[p])continue;
      const bool muK = tpdg[c] == -13 && tproc[c] == 2 && torig[c] != 2 && tmom[c] > 230 && tmom[c] < 241 && (tpar[c] == 321 || tpar[c] == 0);
      if(muK){
        hK[key(sp, p)]->Fill(prange[c]);
        if(tkor[c] >= 0 && tkor[c] <= 3)hS[key(sp, p)]->Fill(prange[c]);
      }
      else if(normalView || !fulldump)hB[key(sp, p)]->Fill(prange[c], scSp[sp]);
    }
  }

  gSystem->mkdir(outdir, true);
  gStyle->SetOptStat(0);
  const std::string pdf = std::string(outdir) + "/xs_fitB.pdf";
  TCanvas cv("cvB", "", 900, 650);
  cv.Print((pdf + "[").c_str());
  TRandom3 rnd(4321);
  std::cout << "Eventos " << n << (normalView ? "; vista normal: fondo de todos los eventos, sin escalar" : "; fondo escalado por especie:")
            << (normalView ? "" : Form(" e/mu/pi x%.3f, K x%.3f, p x%.3f", scSp[0], scSp[1], scSp[2]))
            << (newOrder ? "; muestra del ajuste: sel_fitlevel" : "; muestra del ajuste: nivel de la ventana + Michel (microtree anterior)") << std::endl;
  for(int sp = 0; sp < 3; sp++)for(int p = 0; p < 2; p++){
    TH1D* K = hK[key(sp, p)]; TH1D* S = hS[key(sp, p)]; TH1D* Bk = hB[key(sp, p)];
    const int nb = K->GetNbinsX();
    TH1D hT(*K); hT.SetName(Form("hT_%d_%d", sp, p)); hT.Add(Bk);
    TH1D hD(hT); hD.SetName(Form("hD_%d_%d", sp, p));
    for(int b = 1; b <= nb; b++){ hD.SetBinContent(b, std::round(hT.GetBinContent(b))); hD.SetBinError(b, std::sqrt(hD.GetBinContent(b))); }
    const double nK = K->Integral(), nS = S->Integral();
    std::cout << std::endl << Form("== haz %s, %s (%.0f-%.0f MeV/c): mu de K %.0f (señal B %.0f, pileup y otros %.0f), fondo %.0f ==",
                                   kSpName[sp], popName[p], x0[p], x1[p], nK, nS, nK - nS, Bk->Integral()) << std::endl;
    if(nK < 5){ std::cout << "  demasiado pocos para ajustar" << std::endl; continue; }
    const FitRes r = Fit(&hD, x0[p], x1[p], bw);
    std::cout << Form("  ajuste: N = %.0f ± %.0f (media %.1f, sigma %.1f); frente a todos los mu de K %+.0f %%, frente a la señal B %+.0f %%",
                      r.n, r.en, r.mean, r.sigma, 100*(r.n/nK - 1), nS > 0 ? 100*(r.n/nS - 1) : 0.) << std::endl;
    // figura
    TF1* fm = Model(Form("fmB_%d_%d", sp, p), x0[p], x1[p], bw);
    {
      const double tot = hD.Integral();
      fm->SetParameters(0.1*tot, 228, 9, std::log(std::max(1., hD.GetBinContent(1))), -1, 0);
      fm->SetParLimits(0, 0, 10*tot); fm->SetParLimits(1, 210, 245); fm->SetParLimits(2, 4, 20);
      hD.Fit(fm, "LQR0");
    }
    cv.Clear();
    hD.SetTitle(Form("haz %s, %s;p por rango (hip#acute{o}tesis #mu, MeV/c);candidatos / %g MeV/c", kSpName[sp], popName[p], bw));
    hD.SetMarkerStyle(20); hD.SetMinimum(0); hD.Draw("E");
    Bk->SetLineColor(kAzure + 1); Bk->SetLineWidth(2); Bk->Draw("hist same");
    K->SetLineColor(kOrange + 7); K->SetLineWidth(2); K->Draw("hist same");
    S->SetLineColor(kRed + 1); S->SetLineWidth(2); S->Draw("hist same");
    fm->SetLineColor(kBlack); fm->SetLineWidth(2); fm->Draw("same");
    TLegend leg(0.55, 0.62, 0.89, 0.89); leg.SetBorderSize(0); leg.SetFillStyle(0);
    leg.AddEntry(&hD, "MC (S + B)", "pe"); leg.AddEntry(K, "#mu de K (todos)", "l"); leg.AddEntry(S, "se#tilde{n}al B", "l");
    leg.AddEntry(Bk, "fondo que no es de K", "l"); leg.AddEntry(fm, "ajuste", "l"); leg.Draw();
    cv.Print(pdf.c_str());
    // pseudoexperimentos
    for(double fac : {1., 6.}){
      std::vector<double> nn; int fails = 0;
      for(int it = 0; it < nToys; it++){
        TH1D hToy(hT); hToy.SetName(Form("toyB_%d_%d_%g_%d", sp, p, fac, it));
        for(int b = 1; b <= nb; b++){ const double v = rnd.Poisson(fac*hT.GetBinContent(b)); hToy.SetBinContent(b, v); hToy.SetBinError(b, std::sqrt(v)); }
        const FitRes rt = Fit(&hToy, x0[p], x1[p], bw);
        if(!rt.ok || rt.en <= 0){ fails++; continue; }
        nn.push_back(rt.n/fac);
      }
      std::sort(nn.begin(), nn.end());
      auto q = [&](double x){ return nn.empty() ? -999. : nn[std::min<size_t>(nn.size() - 1, size_t(x*nn.size()))]; };
      const double med = q(0.5), sig = (q(0.84) - q(0.16))/2;
      std::cout << Form("  toys x%-2g: N mediana %.0f (sesgo frente a los mu de K %+.1f %%, frente a la señal B %+.1f %%), error relativo %.1f %%; fallos %d de %d",
                        fac, med, 100*(med/nK - 1), nS > 0 ? 100*(med/nS - 1) : 0., med > 0 ? 100*sig/med : 0., fails, nToys) << std::endl;
    }
  }
  cv.Print((pdf + "]").c_str());
  std::cout << std::endl << "Figuras en " << pdf << std::endl;
  return 0;
}
