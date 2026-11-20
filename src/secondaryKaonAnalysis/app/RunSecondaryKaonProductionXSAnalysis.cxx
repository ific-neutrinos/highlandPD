#include "secondaryKaonProductionXSAnalysis.hxx"
#include "AnalysisLoop.hxx"

int main(int argc, char *argv[]){
  secondaryKaonProductionXSAnalysis* ana = new secondaryKaonProductionXSAnalysis();
  AnalysisLoop loop(ana, argc, argv); 
  loop.Execute();
}
