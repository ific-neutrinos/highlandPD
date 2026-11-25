#include "secondaryKaonPairXSAnalysis.hxx"
#include "AnalysisLoop.hxx"

int main(int argc, char *argv[]){
  secondaryKaonPairXSAnalysis* ana = new secondaryKaonPairXSAnalysis();
  AnalysisLoop loop(ana, argc, argv); 
  loop.Execute();
}
