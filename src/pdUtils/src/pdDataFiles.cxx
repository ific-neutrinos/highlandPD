#include "pdDataFiles.hxx"

#include "TSystem.h"
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <iostream>

namespace{

  // $<PKG>ROOT/data, or an empty string if the variable is not set
  std::string PackageDataDir(const std::string& package){
    std::string var = package;
    std::transform(var.begin(), var.end(), var.begin(), [](unsigned char c){return std::toupper(c);});
    var += "ROOT";
    const char* root = std::getenv(var.c_str());
    if(!root){
      std::cerr << "pd::DataFile(). Environment variable " << var << " is not set "
                << "(source the highlandPD setup script)" << std::endl;
      std::exit(1);
    }
    return std::string(root) + "/data";
  }

  // <HIGHLANDPD_DATADIR>/<package>, or an empty string if not set
  std::string ExternalDataDir(const std::string& package){
    const char* ext = std::getenv("HIGHLANDPD_DATADIR");
    if(!ext || std::string(ext).empty()) return "";
    return std::string(ext) + "/" + package;
  }

  bool Exists(const std::string& path){
    return !gSystem->AccessPathName(path.c_str()); // AccessPathName returns false if it exists
  }
}

//********************************************************************
std::string pd::DataDir(const std::string& package){
//********************************************************************

  std::string ext = ExternalDataDir(package);
  if(!ext.empty() && Exists(ext)) return ext;

  std::string dir = PackageDataDir(package);
  if(!Exists(dir)){
    std::cerr << "pd::DataDir(). Data directory " << dir << " does not exist" << std::endl;
    std::exit(1);
  }
  return dir;
}

//********************************************************************
std::string pd::DataFile(const std::string& package, const std::string& file){
//********************************************************************

  std::string ext = ExternalDataDir(package);
  if(!ext.empty() && Exists(ext + "/" + file)) return ext + "/" + file;

  std::string path = PackageDataDir(package) + "/" + file;
  if(!Exists(path)){
    std::cerr << "pd::DataFile(). Data file " << file << " of " << package << " not found: " << path;
    if(!ext.empty()) std::cerr << " (nor in " << ext << ")";
    std::cerr << std::endl;
    std::exit(1);
  }
  return path;
}
