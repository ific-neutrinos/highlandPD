// Validación de la variación del mapa de SCE sin copias (SCEScaling).
//
// Compara, en puntos aleatorios del detector, dos formas de obtener los
// desplazamientos y el campo E de un mapa variado:
//  - referencia: el método antiguo, una copia del mapa con cada nodo
//    escalado (guardado en TH3F, float) y los splines reconstruidos;
//  - nuevo: el mapa nominal con un SCEScaling (factor global, o factores
//    por nodo aplicados con la base cardinal del spline).
//
// Uso (con el entorno de highland y highlandPD cargado):
//   root -l -b -q test/sce/run_test_sce_voxel.C
//
// Devuelve el número de comprobaciones que fallan (0 = todo bien).

#include "SpaceCharge.hxx"
#include "TRandom3.h"
#include "TMath.h"
#include <iostream>
#include <vector>

namespace{

// Acceso a los mapas internos para construir la referencia (método antiguo)
class RefSCE: public SpaceCharge{
public:
  // Escala cada nodo del mapa de calibración y reconstruye los splines
  void ScaleNodes(const SCEScaling& scale){
    for(int k = 0; k < 12; k++){
      TH3F* h = CalSCEhistograms[k];
      int dv = (k < 6) ? 1 : 0;
      for(int ix = 1; ix <= h->GetNbinsX(); ix++)
        for(int iy = 1; iy <= h->GetNbinsY(); iy++)
          for(int iz = 1; iz <= h->GetNbinsZ(); iz++){
            int bin = h->GetBin(ix,iy,iz);
            double f = scale.IsVoxel() ? scale.NodeFactor(dv,bin) : scale.Global;
            h->SetBinContent(bin, h->GetBinContent(bin)*f);
          }
    }
    ClearSplines();
    InitializeSplines();
  }
};

struct Stats{
  const char* name;
  double maxabs = 0, maxrel = 0, maxval = 0;
  void Add(double ref, double val){
    double d = TMath::Abs(ref-val);
    double m = TMath::Max(TMath::Abs(ref), TMath::Abs(val));
    if(d > maxabs) maxabs = d;
    if(m > 0 && d/m > maxrel && m > 1e-6*maxval) maxrel = d/m;
    if(m > maxval) maxval = m;
  }
  // tolerancia: relativa al mayor valor visto (los valores cercanos a 0 no cuentan)
  bool Print(double tol) const{
    bool ok = maxabs <= tol*maxval;
    printf("  %-28s |dif| máx %-12.4g valor máx %-10.4g dif/valor máx %-10.3g %s\n",
           name, maxabs, maxval, maxval > 0 ? maxabs/maxval : 0., ok ? "OK" : "FALLA");
    return ok;
  }
};

// Compara la referencia (mapa copiado y escalado, consulta nominal) con el
// nominal consultado con el SCEScaling
int Compare(const char* title, const SpaceCharge& nominal, const RefSCE& ref,
            const SCEScaling& scale, int npoints, double tol){
  TRandom3 rnd(12345);
  Stats pos{"desplazamiento (cm)"}, ef{"campo E relativo"}, efp{"campo E en el bin"};
  for(int i = 0; i < npoints; i++){
    double x = rnd.Uniform(-360, 360), y = rnd.Uniform(5, 600), z = rnd.Uniform(0, 695);
    TVector3 p(x, y, z);
    int tpc = (x < 0) ? 1 : 2;
    TVector3 a = ref.GetCalPosOffsets(p, tpc), b = nominal.GetCalPosOffsets(p, tpc, scale);
    for(int c = 0; c < 3; c++) pos.Add(a[c], b[c]);
    a = ref.GetCalEfieldOffsets(p, tpc); b = nominal.GetCalEfieldOffsets(p, tpc, scale);
    for(int c = 0; c < 3; c++) ef.Add(a[c], b[c]);
    efp.Add(ref.GetEFieldXAtPoint(x,y,z), nominal.GetEFieldXAtPoint(x,y,z,scale));
    efp.Add(ref.GetEFieldYAtPoint(x,y,z), nominal.GetEFieldYAtPoint(x,y,z,scale));
    efp.Add(ref.GetEFieldZAtPoint(x,y,z), nominal.GetEFieldZAtPoint(x,y,z,scale));
  }
  std::cout << title << " (" << npoints << " puntos)" << std::endl;
  int nfail = 0;
  nfail += !pos.Print(tol);
  nfail += !ef.Print(tol);
  nfail += !efp.Print(tol);
  return nfail;
}

} // namespace

int test_sce_voxel(int npoints = 5000){

  SpaceCharge nominal;
  nominal.Initialize();

  // Factores por nodo: un parámetro por nodo, sigma 5 %, variaciones gaussianas
  SCEVoxelSigmas sigmas;
  std::vector<Float_t> vars;
  TRandom3 rnd(4321);
  for(int dv = 0; dv < 2; dv++){
    const TH3F* h = nominal.GetCalHistogram(dv == 0 ? 6 : 0);
    int nbins = (h->GetNbinsX()+2)*(h->GetNbinsY()+2)*(h->GetNbinsZ()+2);
    sigmas.Index[dv].assign(nbins, -1);
    sigmas.Sigma[dv].assign(nbins, 0.);
    for(int ix = 1; ix <= h->GetNbinsX(); ix++)
      for(int iy = 1; iy <= h->GetNbinsY(); iy++)
        for(int iz = 1; iz <= h->GetNbinsZ(); iz++){
          int bin = h->GetBin(ix,iy,iz);
          sigmas.Index[dv][bin] = vars.size();
          sigmas.Sigma[dv][bin] = 0.05;
          vars.push_back(rnd.Gaus(0,1));
        }
  }

  int nfail = 0;

  // 1. Por vóxel: factores distintos en cada nodo
  {
    SCEScaling scale;
    scale.Voxels = &sigmas;
    scale.ToyVariations = &vars[0];
    RefSCE ref; ref.Initialize(); ref.ScaleNodes(scale);
    nfail += Compare("Vóxel, factor distinto en cada nodo", nominal, ref, scale, npoints, 1e-5);
  }

  // 2. Global: un único factor
  {
    SCEScaling scale;
    scale.Global = 1 + 0.08f*1.7f;
    RefSCE ref; ref.Initialize(); ref.ScaleNodes(scale);
    nfail += Compare("Global, factor 1.136", nominal, ref, scale, npoints, 1e-5);
  }

  // 3. Vóxel con el mismo factor en todos los nodos == global con ese factor
  {
    std::vector<Float_t> same(vars.size(), 1.7f);
    SCEVoxelSigmas s8 = sigmas;
    for(int dv = 0; dv < 2; dv++)
      for(auto& sg : s8.Sigma[dv]) sg = 0.08;
    SCEScaling vox; vox.Voxels = &s8; vox.ToyVariations = &same[0];
    SCEScaling glo; glo.Global = 1 + 0.08f*1.7f;
    TRandom3 r(999);
    Stats st{"vóxel uniforme vs global"};
    for(int i = 0; i < npoints; i++){
      TVector3 p(r.Uniform(-360,360), r.Uniform(5,600), r.Uniform(0,695));
      int tpc = (p.X() < 0) ? 1 : 2;
      TVector3 a = nominal.GetCalPosOffsets(p, tpc, glo), b = nominal.GetCalPosOffsets(p, tpc, vox);
      for(int c = 0; c < 3; c++) st.Add(a[c], b[c]);
      a = nominal.GetCalEfieldOffsets(p, tpc, glo); b = nominal.GetCalEfieldOffsets(p, tpc, vox);
      for(int c = 0; c < 3; c++) st.Add(a[c], b[c]);
    }
    std::cout << "Vóxel con todos los factores iguales frente a global (" << npoints << " puntos)" << std::endl;
    nfail += !st.Print(1e-9);
  }

  std::cout << (nfail == 0 ? "RESULTADO: OK" : "RESULTADO: FALLA") << " (" << nfail << ")" << std::endl;
  return nfail;
}
