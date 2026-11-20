#include "beamWeightsAnalysis.hxx"
#include "AnalysisLoop.hxx"

int main(int argc, char *argv[]){
  beamWeightsAnalysis* ana = new beamWeightsAnalysis();
  AnalysisLoop loop(ana, argc, argv); 
  loop.Execute();
}
