#ifndef pdServiceManager_h
#define pdServiceManager_h

#include <memory>

class SpaceCharge;
class Calorimetry;

/// Owns the services shared by corrections, systematics and analyses
/// (nominal space charge map, calorimetry). Each service is created the
/// first time it is requested and lives until the end of the program.
/// Services are handed out as const: anything that changes per event or
/// per toy (lifetime, ModBox parameters, varied SCE) is passed as an
/// argument (see CaloContext in Calorimetry.hxx).
///
/// Access it through pd::services(), e.g.
///   const Calorimetry& cal = pd::services().Calo();
class pdServiceManager{
public:

  static pdServiceManager& Get();

  /// Nominal space charge map
  const SpaceCharge& SCE();

  /// Calorimetry, using the nominal space charge map
  const Calorimetry& Calo();

  /// Print the services created so far
  void Print() const;

private:

  pdServiceManager() = default;
  pdServiceManager(const pdServiceManager&) = delete;
  pdServiceManager& operator=(const pdServiceManager&) = delete;

  std::unique_ptr<SpaceCharge> _sce;
  std::unique_ptr<Calorimetry> _cal;
};

namespace pd{
  /// Access to the shared services. See pdServiceManager.
  pdServiceManager& services();
}

#endif
