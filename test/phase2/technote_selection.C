// Figuras del capítulo de selección de la nota técnica: para cada población (0 en la beam slice, 1 fuera) y cada
// corte de su cadena (secondaryKaonProductionXSSelection), la variable del corte para los candidatos que han pasado
// todos los cortes anteriores (cand_level >= índice del corte - 1), apilada por la categoría true, con la línea del
// corte. Sobre el microtree "ana" del análisis secondaryKaonProductionXS, con las DrawingTools de highland.
//
//   root -l 'test/phase2/technote_selection.C("xs_v100.root", "<technote>/4_selection/figs")'
//
// Con waitplots = true (abajo) se para en cada figura: mirarla y hacer doble clic en el canvas (o pulsar una tecla con
// el canvas activo) para pasar a la siguiente. Para generarlas todas sin parar, waitplots = false o root -b.
//
// categ: candmu por defecto (signal muons en capas excluyentes: A, B que no es A, C que no es B; microtrees desde
// el commit de las capas), o candsigC/B/A (señal de una sola medida). En microtrees anteriores candmu es la
// clasificación antigua por origen. Imprime las purezas al final de cada cadena.

bool waitplots = true;
bool logplots  = true;   // eje y logarítmico en todas las figuras (si no, el de cada corte: SelCut::logy)

struct SelCut {
  int index;              // posición del corte en la cadena de su población (1 = el corte de población)
  std::string name;       // para el nombre del fichero
  std::string var;        // expresión del microtree
  std::string title;      // eje x
  int nbins; double xmin, xmax;
  double lo, hi;          // corte: se queda [lo, hi]; -999 = sin límite
  bool logy;
};

static void DrawCut(DrawingTools& draw, TTree* mc, const std::string& categ, int pop, const SelCut& c, const std::string& outdir){
  const std::string popTag = pop == 0 ? "in" : "out";
  const std::string sel = Form("cand_pop==%d && cand_level>=%d", pop, c.index - 1);
  draw.SetTitleX(c.title);
  const bool logy = logplots || c.logy;
  // el orden importa: SetLogY(true) no hace nada si el mínimo es <= 0, y SetMinY(0) no hace nada en escala log
  if(logy){ draw.SetMinY(0.5); draw.SetLogY(true); }
  else    { draw.SetLogY(false); draw.SetMinY(0); }
  draw.Draw(mc, c.var, c.nbins, c.xmin, c.xmax, categ, sel, "", "UNDER OVER");
  if(c.lo > -998 && c.hi > -998)draw.DrawCutLinesVertical(c.lo, c.hi, true);
  else if(c.lo > -998)draw.DrawCutLineVertical(c.lo, true, "r");
  else if(c.hi > -998)draw.DrawCutLineVertical(c.hi, true, "l");
  TLatex tt; tt.SetNDC(); tt.SetTextSize(0.04);
  tt.DrawLatex(0.15, 0.94, Form("#bf{DUNE:ProtoDUNE-SP} simulation, %s the beam slice", pop == 0 ? "in" : "out of"));
  gPad->RedrawAxis(); gPad->Update();
  if(waitplots && !gROOT->IsBatch()){ gPad->Update(); gPad->WaitPrimitive(); }
  gPad->Print(Form("%s/sel_%s_%02d_%s.pdf", outdir.c_str(), popTag.c_str(), c.index, c.name.c_str()));
}

int technote_selection(const char* in, const char* outdir, const char* categ = "candmu"){
  TFile* f = TFile::Open(in);
  if(!f || f->IsZombie()){ std::cerr << "No se puede abrir " << in << std::endl; return -1; }
  TTree* mc = (TTree*)f->Get("ana");
  gSystem->mkdir(outdir, true);

  DrawingTools draw(in);
  draw.ChangeToProtoDUNEStyle();
  draw.SetStackFillStyle(1001);
  draw.SetLegendPos("r");
  gStyle->SetOptStat(0);

  // la señal abajo del todo en el stack: SetCategorySignal pasa el código al final de la lista de tipos (se dibuja
  // abajo). Se llama de la capa más externa a la más interna, para que quede A abajo, encima B (no A) y encima
  // C (no B); en candsigX, la señal abajo y encima los mu de K de otro origen
  const std::string cat(categ);
  if(cat == "candmu"){ draw.SetCategorySignal(cat, 3); draw.SetCategorySignal(cat, 2); draw.SetCategorySignal(cat, 1); }
  else if(cat.rfind("candsig", 0) == 0){ draw.SetCategorySignal(cat, 2); draw.SetCategorySignal(cat, 1); }

  // distancia mínima de los dos extremos a los huecos entre APAs (z = 230 y 462 cm), como muCandCuts::GapDist
  const std::string gap = "TMath::Min(TMath::Min(abs(cand_stop[][2]-230),abs(cand_stop[][2]-462)),"
                          "TMath::Min(abs(cand_start[][2]-230),abs(cand_start[][2]-462)))";

  // cadenas de secondaryKaonProductionXSSelection::DefineSteps (el corte 1 es la población: no se dibuja)
  const std::vector<SelCut> in_ = {
    {2, "faces",    "TMath::Min(cand_faces_stop,cand_faces_start)", "Distance of the ends to the faces (cm)", 40, 0, 200, 35, -999, false},
    {3, "vertex",   "cand_vtx_min", "Distance of the ends to the beam vertex (cm)", 42, -5, 100, 15, -999, false},
    {4, "chi2mu",   "cand_chi2mu", "#chi^{2}_{#mu}/ndf", 40, 0, 40, -999, 6, false},
    {5, "dedxstop", "cand_dedx_stop", "dE/dx in the last 5 cm (MeV/cm)", 40, 0, 10, 3, -999, false},
    {6, "gaps",     gap, "Distance of the ends to the APA gaps in z (cm)", 40, 0, 100, 5, -999, false},
    {7, "prange",   "cand_prange", "Momentum from range, #mu hypothesis (MeV/c)", 40, 100, 400, 205, 240, false},
    {8, "michel",   "cand_ar_nmichel10", "Michel hits within 10 cm of the stop", 20, 0, 20, 2, -999, false},
  };
  const std::vector<SelCut> out_ = {
    {2, "clearcosmic", "cand_clearcosmic", "Pandora IsClearCosmic", 2, 0, 2, -999, 0.5, false},
    {3, "axis",     "cand_axis_min", "Distance of the closest end to the beam axis (cm)", 30, 0, 150, -999, 60, false},
    {4, "chi2mu",   "cand_chi2mu", "#chi^{2}_{#mu}/ndf", 40, 0, 40, -999, 6, false},
    {5, "dedxstop", "cand_dedx_stop", "dE/dx in the last 5 cm (MeV/cm)", 40, 0, 10, 3, -999, false},
    {6, "chi2p",    "cand_chi2p", "#chi^{2}_{p}/ndf", 40, 0, 400, 80, -999, false},
    {7, "dedxstart","cand_dedx_start", "dE/dx in the first 5 cm (MeV/cm)", 40, 0, 6, 1.8, -999, false},
    {8, "gaps",     gap, "Distance of the ends to the APA gaps in z (cm)", 40, 0, 100, 5, -999, false},
    {9, "prange",   "cand_prange", "Momentum from range, #mu hypothesis (MeV/c)", 40, 100, 400, 205, 240, false},
    {10, "michel",  "cand_ar_nmichel10", "Michel hits within 10 cm of the stop", 20, 0, 20, 2, -999, false},
  };

  for(const auto& c : in_)DrawCut(draw, mc, categ, 0, c, outdir);
  for(const auto& c : out_)DrawCut(draw, mc, categ, 1, c, outdir);

  std::cout << std::endl << "== Purezas al final de la cadena (" << categ << ") ==" << std::endl;
  draw.PrintPurities(mc, categ, Form("cand_pop==0 && cand_level>=%zu", in_.size() + 1));
  draw.PrintPurities(mc, categ, Form("cand_pop==1 && cand_level>=%zu", out_.size() + 1));
  return 0;
}
