#include "pdServiceManager.hxx"
#include "SpaceCharge.hxx"
#include "Calorimetry.hxx"

#include "TH1.h"
#include <iostream>

namespace{
  /// Histograms created by the services must not be attached to the
  /// current ROOT directory (it could be the output file).
  class NoAddDirectory{
  public:
    NoAddDirectory():_status(TH1::AddDirectoryStatus()){TH1::AddDirectory(kFALSE);}
    ~NoAddDirectory(){TH1::AddDirectory(_status);}
  private:
    bool _status;
  };
}

//********************************************************************
pdServiceManager& pdServiceManager::Get(){
//********************************************************************

  // Never deleted on purpose: destroying ROOT objects during static
  // destruction, after ROOT's own cleanup, is not safe.
  static pdServiceManager* instance = new pdServiceManager();
  return *instance;
}

//********************************************************************
const SpaceCharge& pdServiceManager::SCE(){
//********************************************************************

  if(!_sce){
    NoAddDirectory guard;
    std::cout << "pdServiceManager::SCE(). Creating nominal space charge map" << std::endl;
    _sce.reset(new SpaceCharge());
    _sce->Initialize();
  }
  return *_sce;
}

//********************************************************************
const Calorimetry& pdServiceManager::Calo(){
//********************************************************************

  if(!_cal){
    const SpaceCharge& sce = SCE();
    NoAddDirectory guard;
    std::cout << "pdServiceManager::Calo(). Creating calorimetry" << std::endl;
    _cal.reset(new Calorimetry(sce));
  }
  return *_cal;
}

//********************************************************************
void pdServiceManager::Print() const{
//********************************************************************

  std::cout << " -------- List of pdServices --------------------------------------" << std::endl;
  std::cout << "  SCE  (nominal space charge map): " << (_sce ? "created" : "not used") << std::endl;
  std::cout << "  Calo (calorimetry):              " << (_cal ? "created" : "not used") << std::endl;
  std::cout << " ------------------------------------------------------------------" << std::endl;
}

//********************************************************************
pdServiceManager& pd::services(){
//********************************************************************

  return pdServiceManager::Get();
}
