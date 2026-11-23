#include "ToyBoxMuCand.hxx"
#include <typeinfo>

//********************************************************************
ToyBoxMuCand::ToyBoxMuCand():ToyBoxPD(){
//********************************************************************

  Candidates.clear();
  Population.clear();
  Level.clear();
  NPreselected = 0;
  HasVertex = false;
  for(int i = 0; i < 3; i++)VertexPos[i] = -999;
  HasTrueVertex = false;
  for(int i = 0; i < 3; i++)TrueVertexPos[i] = -999;
}

//********************************************************************
void ToyBoxMuCand::Reset(){
//********************************************************************

  if (!_ResetCheckDone){
    if( typeid(*this) !=  typeid(ToyBoxMuCand)){
      std::cerr << "ERROR in ToyBoxMuCand::Reset(). Either this mandatory method is not implemented "
                << "by the derived class '" << typeid(*this).name() << "' "
                << "or ToyBoxMuCand::Reset() is called from the Reset method of the derived class. "
                << "Please correct any of these bugs. " << std::endl;
      exit(1);
    }
    _ResetCheckDone=true;
  }
}

//********************************************************************
void ToyBoxMuCand::ResetBase(){
//********************************************************************

  ToyBoxPD::ResetBase();
  Candidates.clear();
  Population.clear();
  Level.clear();
  NPreselected = 0;
  HasVertex = false;
  for(int i = 0; i < 3; i++)VertexPos[i] = -999;
  HasTrueVertex = false;
  for(int i = 0; i < 3; i++)TrueVertexPos[i] = -999;
  SoftReset();  // just reset internal stuff
}
