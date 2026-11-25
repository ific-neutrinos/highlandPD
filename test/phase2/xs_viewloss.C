// Uso: root -l -b -q 'xs_viewloss.C("xs_v100s_fix.root")' (25-nov-2026)
// Pérdida de signal muons por la vista normal frente a la completa (árbol truth de secondaryKaonProductionXS)
void xs_viewloss(const char* f){ TFile* F=TFile::Open(f); TTree* t=(TTree*)F->Get("truth");
 Int_t al[2], lf[2], kor; t->SetBranchAddress("accum_level",al); t->SetBranchAddress("truemu_level_full",lf); t->SetBranchAddress("truemu_korigin",&kor);
 const int fitL[2]={8,10}, cntL[2]={9,11}; // 1 + kFitLevel, 1 + kCountLevel
 const char* ln[3]={"A","B no A","C no B"}; const char* bn[2]={"dentro","fuera"};
 // n[layer][branch][stage][view]; stages: 0 cand. emparejado (>=1), 1 población (>=2), 2 muestra del ajuste, 3 todos los cortes
 double n[3][2][4][2]={0}, tot[3]={0}, both[3][2][2]={0}; double lostNoCand[3][2]={0}, lostLower[3][2]={0}, gained[3][2]={0};
 int maxLvl=0;
 for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
   int L = (kor==0||kor==1)?0:(kor==2||kor==3)?1:2; tot[L]++;
   bool passN=false, passF=false;
   for(int b=0;b<2;b++){ int th[4]={1,2,fitL[b],cntL[b]}; maxLvl=std::max(maxLvl,std::max(al[b],lf[b]));
     for(int s=0;s<4;s++){ if(al[b]>=th[s]) n[L][b][s][0]++; if(lf[b]>=th[s]) n[L][b][s][1]++; }
     if(lf[b]>=cntL[b] && al[b]<cntL[b]){ if(al[b]==0) lostNoCand[L][b]++; else lostLower[L][b]++; }
     if(al[b]>=cntL[b] && lf[b]<cntL[b]) gained[L][b]++;
     if(al[b]>=cntL[b]) passN=true; if(lf[b]>=cntL[b]) passF=true; }
   if(passN) both[L][0][0]++; if(passF) both[L][0][1]++; }
 cout<<"Entradas "<<t->GetEntries()<<" (signal muons, p true [230, 241]); nivel máximo "<<maxLvl<<endl;
 for(int L=0;L<3;L++) cout<<"  "<<ln[L]<<": "<<tot[L]<<endl;
 const char* sn[4]={"candidato emparejado","en la población","muestra del ajuste B","todos los cortes"};
 for(int b=0;b<2;b++){ cout<<endl<<"== Rama "<<bn[b]<<": normal / completa (pérdida %)"<<endl;
   for(int s=0;s<4;s++){ cout<<Form("  %-22s",sn[s]);
     for(int L=0;L<3;L++){ double a=n[L][b][s][0], c=n[L][b][s][1]; cout<<Form("  %s %5.0f / %5.0f (%5.1f %%)",ln[L],a,c,c?100*(c-a)/c:0.); }
     double a=0,c=0; for(int L=0;L<2;L++){a+=n[L][b][s][0];c+=n[L][b][s][1];} cout<<Form("  | B (A+B no A) %5.0f / %5.0f (%5.1f %%)",a,c,c?100*(c-a)/c:0.)<<endl; }
   cout<<"  todos los cortes, perdidos: sin candidato en la normal / con nivel menor; ganados en la normal:"<<endl;
   for(int L=0;L<3;L++) cout<<Form("    %-7s %4.0f / %4.0f; %4.0f",ln[L],lostNoCand[L][b],lostLower[L][b],gained[L][b])<<endl; }
 cout<<endl<<"== Las dos ramas (todos los cortes en alguna): normal / completa"<<endl;
 for(int L=0;L<3;L++){ double a=both[L][0][0],c=both[L][0][1]; cout<<Form("  %-7s %5.0f / %5.0f (pérdida %5.1f %%; eficiencia %.3f / %.3f)",ln[L],a,c,c?100*(c-a)/c:0.,a/tot[L],c/tot[L])<<endl; }
}
