#include "QPIXUtils.hxx"

#include "QPIXDataClasses.hxx"
#include "CategoryManager.hxx"
#include "pdAnalysisUtils.hxx"

//********************************************************************
void QPIXUtils::AddCustomCategories(){
//********************************************************************

  std::string types[] = {"B8", "hep", "Ar39"    , "Ar42"    , "Bi214"   , "Co60"    , "K40"     , "K42"     , "Kr85"    , "Pb214"   , "Po210"   , "Rn222"   , NAMEOTHER};
  int codes[]         = {8   , 9    , 1000180390, 1000180420, 1000832140, 1000270600, 1000190400, 1000190420, 1000360850, 1000822140, 1000842100, 1000862220, CATOTHER};
  int colors[]        = {2   , 208  , 3         , 4         , 94        , 7         , 6         , 31        ,  92       , 8         , 1         , 9         , COLOTHER};
  const int NPART = sizeof(types)/sizeof(types[0]);
  
  std::reverse(types,  types  + NPART);
  std::reverse(codes,  codes  + NPART);
  std::reverse(colors, colors + NPART);

  anaUtils::_categ->AddCategory("qpix", NPART, types, codes, colors);
}


//********************************************************************
void QPIXUtils::FillCustomCategories(AnaEventPD* event){
//********************************************************************

  AnaEventInfoQPIX* EventInfo = static_cast<AnaEventInfoQPIX*>(event->EventInfo);
  if(!EventInfo->IsBackground){
    if(EventInfo->SolarChain == AnaEventInfoQPIX::SolarChainEnum::B8)
      anaUtils::_categ->SetCode("qpix", 8, CATOTHER);
    if(EventInfo->SolarChain == AnaEventInfoQPIX::SolarChainEnum::hep)
      anaUtils::_categ->SetCode("qpix", 9, CATOTHER);
  }
  else{
    AnaTrueParticleB* TruePart = event->TrueParticles[0];
    anaUtils::_categ->SetCode("qpix", TruePart->PDG, CATOTHER);
  }
}

//********************************************************************
void QPIXUtils::AddQPIXVariables(OutputManager& output){
//********************************************************************
  
  AddVarI(    output, is_background,    "is this a background event?");
  AddVarI(    output, has_potasium ,    "has this event potasium?");
  AddVarI(    output, has_gamma    ,    "has this event a delayed gamma?");
  AddVarF(    output, nu_E         ,    "neutrino energy");
  AddVar4VF(  output, nu_pos       ,    "neutrino position");
  AddVarF(    output, e_E          ,    "electron energy");
  AddVarFixVI(output, nphotons     ,    "detected photons per plane", 7);
  AddVarFixVI(output, wf_total     ,    "PDS total wf"              , 1024);
  AddVarFixVI(output, wf_plane_0   ,    "PDS plane 0 wf"            , 1024);
  AddVarFixVI(output, wf_plane_1   ,    "PDS plane 1 wf"            , 1024);
  AddVarFixVI(output, wf_plane_2   ,    "PDS plane 2 wf"            , 1024);
  AddVarFixVI(output, wf_plane_3   ,    "PDS plane 3 wf"            , 1024);
  AddVarFixVI(output, wf_plane_4   ,    "PDS plane 4 wf"            , 1024);
  AddVarFixVI(output, wf_plane_5   ,    "PDS plane 5 wf"            , 1024);
}

//********************************************************************
void QPIXUtils::FillQPIXVariables(OutputManager& output, AnaEventPD* event){
//********************************************************************

  AnaEventInfoQPIX* EventInfo = static_cast<AnaEventInfoQPIX*>(event->EventInfo);
  
  output.FillVar(               is_background, (int)EventInfo->IsBackground);
  output.FillVar(               has_potasium , (int)EventInfo->HasPotasium);
  output.FillVar(               has_gamma    , (int)EventInfo->HasGamma);
  output.FillVectorVarFromArray(nphotons     , EventInfo->DetectedPhotons, 7);
  for(int i = 0; i < 1024; i++){
    output.FillVectorVar(       wf_total     , EventInfo->PDSwf[0][i]    , i);
    output.FillVectorVar(       wf_plane_0   , EventInfo->PDSwf[1][i]    , i);
    output.FillVectorVar(       wf_plane_1   , EventInfo->PDSwf[2][i]    , i);
    output.FillVectorVar(       wf_plane_2   , EventInfo->PDSwf[3][i]    , i);
    output.FillVectorVar(       wf_plane_3   , EventInfo->PDSwf[4][i]    , i);
    output.FillVectorVar(       wf_plane_4   , EventInfo->PDSwf[5][i]    , i);
    output.FillVectorVar(       wf_plane_5   , EventInfo->PDSwf[6][i]    , i);
  }
  //if(EventInfo->IsBackground)return;
  AnaTrueVertex* vertex = static_cast<AnaTrueVertex*>(event->TrueVertices[0]);
  if(!vertex)return;
  output.FillVar(               nu_E         , vertex->NuEnergy);
  output.FillVectorVarFromArray(nu_pos       , vertex->Position, 4);
  output.FillVar(               e_E          , vertex->LeptonMom);
}

