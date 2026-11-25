// Signal muon candidates (cand_truesignalmuon) por tipo (traza / shower) y nivel de la cadena, y fondo, en el árbol ana de secondaryKaonProductionXS (25-nov-2026)
void xs_showers(const char* f){ TFile* F=TFile::Open(f); TTree* t=(TTree*)F->Get("ana");
 const int N=100; Int_t nc, pop[N], lvl[N], sh[N], sig[N], kor[N], sm[N]; Float_t pr[N], pur[N], dst[N], c2[N];
 t->SetBranchStatus("*",0); auto B=[&](const char* n, void* a){t->SetBranchStatus(n,1); t->SetBranchAddress(n,a);};
 B("ncand",&nc); B("cand_pop",pop); B("cand_level",lvl); B("cand_isshower",sh); B("cand_truesignalmuon",sig); B("cand_truekorigin",kor); B("cand_savemode",sm); B("cand_prange",pr); B("cand_truepur",pur); B("cand_dedx_stop",dst); B("cand_chi2mu",c2);
 // h[pop][isshower][level]
 double h[2][2][12]={0}, hb[2][2][12]={0}; double noDedx[2][2]={0}, nsig[2][2]={0}, sm1[2][2]={0};
 for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
  for(int j=0;j<nc&&j<N;j++){ int p=pop[j]; if(p<0||p>1) continue; int s=sh[j]?1:0; int l=std::min(std::max(lvl[j],0),11);
    if(sig[j]==1){ h[p][s][l]++; nsig[p][s]++; if(dst[j]<0) noDedx[p][s]++; if(sm[j]==1) sm1[p][s]++; }
    else for(int k=0;k<=l;k++) hb[p][s][k]++; } }
 const char* pn[2]={"dentro","fuera"};
 for(int p=0;p<2;p++){ cout<<"== "<<pn[p]<<endl;
  for(int s=0;s<2;s++){ cout<<(s?"  shower":"  traza ")<<Form(": signal muons %5.0f (sin dE/dx en la parada %4.0f, en resumen %4.0f); llegan al nivel >= L:",nsig[p][s],noDedx[p][s],sm1[p][s])<<endl<<"    sig ";
   for(int L=0;L<12;L++){ double c=0; for(int k=L;k<12;k++) c+=h[p][s][k]; cout<<Form(" %5.0f",c);} cout<<endl<<"    fondo";
   for(int L=0;L<12;L++) cout<<Form(" %5.0f",hb[p][s][L]); cout<<endl; } }
}
