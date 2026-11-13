#ifndef pdDataFiles_h
#define pdDataFiles_h

#include <string>

/// Single place where highlandPD resolves the location of its data files
/// (SCE maps, calibrations, PID tables, systematic and correction sources).
///
/// By default the files of package <pkg> are in $<PKG>ROOT/data, e.g.
/// $PDUTILSROOT/data. If the environment variable HIGHLANDPD_DATADIR is
/// set, <HIGHLANDPD_DATADIR>/<pkg> is looked at first, so the data files
/// can live outside the source tree (an external data area, git-lfs...).
///
/// A missing environment variable, directory or file is a fatal error with
/// a clear message (exit code 1).
namespace pd{

  /// Full path of a data file of a package, e.g. DataFile("pdUtils", "ke_vs_range.root")
  std::string DataFile(const std::string& package, const std::string& file);

  /// Data directory of a package, e.g. for BinnedParams: DataDir("pdSystematics")
  std::string DataDir(const std::string& package);
}

#endif
