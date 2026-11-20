#include "brokenTracksAnalysis.hxx"
#include "AnalysisLoop.hxx"

int main(int argc, char *argv[]){
  brokenTracksAnalysis* ana = new brokenTracksAnalysis();
  AnalysisLoop loop(ana, argc, argv); 
  loop.Execute();
}
