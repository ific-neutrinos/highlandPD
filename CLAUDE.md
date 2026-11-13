# highlandPD: contexto para trabajar en este repo

Paquete ProtoDUNE de HighLAND. Depende de `../highland` (HighLAND base, ver su `CLAUDE.md`).
Remoto: `github.com/ific-neutrinos/highlandPD`. Rama de integración: `develop`.

Este fichero recoge la auditoría hecha el 13-nov-2026 y el plan acordado. El usuario escribe en español.

## Forma de trabajar (acordada)

- Trabajar en una rama creada desde `develop`, por ejemplo `cleanup/highlandPD`. No tocar `develop` ni hacer merge ni push sin que el usuario lo pida. El usuario no quiere subir nada ni abrir PR hasta tener todos los cambios hechos (fases completas), para no dejar cosas a medias en el remoto.
- Antes de crear la rama: `git status`. Si hay cambios locales sin commitear o sin subir, preguntar al usuario qué hacer con ellos.
- Commits pequeños y temáticos. Compilar después de cada uno.
- **Fase 1 (limpieza): ningún cambio de comportamiento.** Solo borrar o mover ficheros, `.gitignore`, CMake equivalente y quitar prints de depuración. Los cambios que alteran resultados de física van en commits separados y solo después de tener el test de regresión.
- Todo lo que se elimine debe quedar recuperable. Antes de borrar, crear un tag de archivo (p. ej. `archive/pre-cleanup`) o una rama de archivo.

## Compilar

```bash
source ../highland/scripts/setup.sh     # entorno de highland (ROOT, HIGHLANDPATH, ...)
source scripts/setup.sh                 # entorno de highlandPD
source scripts/INSTALL.sh               # primera vez (crea build/, compila e instala en lib/ bin/ include/)
cmake --build build && cmake -P build/cmake_install.cmake   # recompilar
```

`scripts/DOALL.sh` limpia y recompila highland y highlandPD desde cero.
`CMakeLists.txt` compila **todos** los `src/*/src/*.cxx` en una sola `libhighlandPD.so`, y luego 14 ejecutables `*.exe` con la función `highlandpd_add_executable`. Los ejecutables solo dependen de la librería, así que se puede compilar en paralelo (`cmake --build build -j8`).

## Qué hay

- **Núcleo:** `pdEventModel` (`pdDataClasses`), `pdIO` (`PDSPAnalyzerTreeConverter`, `MichelRemovingTreeConverter`, `pdMiniTreeConverter`, `pdBaseConverter`), `pdUtils` (`Calorimetry`, `SpaceCharge`, `pdAnalysisUtils`, `standardPDTree`), `pdCorrections` (13 correcciones), `pdSystematics` (variaciones y pesos), `pdBaseAnalysis`.
- **Análisis:** `secondaryKaonAnalysis` (incluye `coherentfit/`, que no se compila), `dEdxCalibrationAnalysis`, `StoppingProtonAnalysis`, `BrokenTracksAnalysis`, `BeamWeightsAnalysis`, `pdExampleAnalysis`.
- **No se compila:** `secondaryKaonAnalysis/src/coherentfit`. `CreateMiniTreePD.exe` y `CoherentFit.exe` están comentados en el CMake.
- **Solo en local, ignorado por git:** `src/secondaryKaonAnalysis/macros/tesis/` (macros de la tesis del usuario: no borrarlas del disco).

## Recuperar lo eliminado

- Tag `archive/pre-cleanup` (local, sin subir): estado de `develop` (`cc1a931`) antes de la limpieza.
- Copia en disco: `~/Public/highland/backup/highlandPD/`, con la misma estructura de carpetas que el repo. Incluye `LArSoftReader/`, `QPIXAnalysis/`, `unusedCodeToRevisit/`, `macros/tesis/` (con los cambios locales sin commitear que tenía el usuario), el `.pb` del CNN y los ficheros sueltos borrados.

## Fase 1: limpieza (sin cambios de comportamiento)

Hecho en la rama `cleanup/highlandPD` (13-nov-2026), sin push:

- `.gitignore` (productos de compilación, autosaves de emacs, `nohup.out`).
- Quitados del repo: `LArSoftReader/`, `QPIXAnalysis/` (decisión del usuario: sale de highlandPD), `unusedCodeToRevisit/` (también borrado del disco) y `macros/tesis/` (sigue en disco), `nohup.out`, `coherentfit/oldfits`, dos autosaves de emacs, `dEdxCalibrationAnalysis/macros/x.C` (idéntico a `yz.C` y definía `yz()`), los 4 módulos CMake de `/cvmfs` y el modelo CNN `.pb` de MCC11 junto con el parámetro `pdUtils.CNN.PBfile`, que no leía nadie. Los ficheros borrados siguen en el historial de git, así que el `.git` no adelgaza.
- Borrados de disco (no estaban en git): `CMakeLists_backup.txt`, `scripts/setup_backup.sh`, `bin/CoherentFit.exe` (binario viejo).
- Añadido `pdUtils/data/dedx_mpv_profiles.root`, que estaba sin seguir; ningún código lo referencia todavía.
- Quitados los prints de depuración del bucle de eventos: `MIGUEEEE...` en `GetdEdxLikelihood`, el de cada hit en `dEdxLikelihood` y los contadores de `DaughterFinderCorrection::Apply`. El resto de `std::cout` se queda: métodos `Print()`, mensajes de error y mensajes de `Initialize`/`Finalize`. Pasar los errores a `std::cerr` y corregir los textos de `Calorimetry` (los de YZ dicen "X correction") queda para más adelante.
- Rutas `/home/migue/` de `pdUtils/data/likelihood_microtree.C` y `secondaryKaonAnalysis/macros/selfsysttest.C` pasadas a argumentos.
- CMake: función para los ejecutables, sin la cadena de `add_dependencies`, sin `-ansi`, `cmake_minimum_required(3.12)`. `ENABLE_EXPORTS` mantiene el `-rdynamic` que CMake 3.0 añadía siempre (CMP0065). Flags comprobados iguales a los de antes.

Pendiente:

1. **Ficheros de datos de las utilidades** (`pdUtils/data/`, varios MB): decidir cómo gestionarlos (git-lfs, área externa...). El usuario quiere pensarlo más adelante.

## Test de regresión

`test/regression/` (override `regression.parameters.dat`, `run_regression.sh`, `compare.C`). Usar minitrees **sin filtrar**: los de candidatos (`minitree_candidates_*`) no traen las hijas de los candidatos y `secondaryKaon.exe` casca en `EventBoxKaon.cxx:52` (`Daughters[0]` con `Daughters` vacío; el bug ya existía antes de la limpieza).

```bash
test/regression/run_regression.sh <dir_salida> <mc.root|-> <data.root|-> 10000
root -l -b -q 'test/regression/compare.C+("ref/mc.root","<dir_salida>/mc.root")'
```

- MC: `~/Documents/Work/DUNE/analysis/PDSP/secondary_kaon_analysis/files/minitrees/mc/6GeV_prod4a_00_minitree_2023-01-27.root` (5 GB, 47.906 eventos). Con 10.000 eventos: ~3 min, 3.371 eventos guardados. Datos: pendiente de descarga.
- Fase 1 validada (13-nov-2026): con 10.000 eventos de MC, la salida de `archive/pre-cleanup` y la de la rama son idénticas en todos los árboles (`ana`, `all_syst`, `truth`, `header`). En el log solo cambian los prints de depuración quitados.
- El override activa también `RecombinationVariation`, el único sistemático del análisis de kaones que usa `Calorimetry`.
- Referencia actual: `~/Public/highland/backup/regression_ref/recomb/mc.root`, generada con `archive/pre-cleanup`. La de `fase1/` es la anterior, sin recombinación. Cada bug de la Fase 2 debe compararse contra la referencia anterior, y las diferencias deben explicarse por el arreglo.

## Servicios de calorimetría y SCE (en curso, 13-nov-2026)

Decidido con el usuario (primero el bloque de diseño, antes que los bugs sueltos):

- **`pdServiceManager`** en `pdUtils`, con acceso global `pd::services()` al estilo de `ND::params()`. Tiene que ser global porque las correcciones y los sistemáticos no tienen acceso a `ana()`. Usa métodos con tipo (`SCE()`, `Calorimetry()`), creación perezosa (nunca en la inicialización estática), objetos `const` y lista de servicios en el log.
- **Auditoría:** todas las `Calorimetry` y `SpaceCharge` tienen la misma configuración base. Solo cambia estado por evento o por toy: lifetime (`LifetimeVariation`, `SCEVariation`), ModBox A/B (`RecombinationVariation`) y qué SCE usar. Basta con una instancia nominal de cada una, y ese estado viaja como contexto en cada llamada.
- **SCE en los toys sin copias del mapa.** La interpolación (spline natural `"b2e2"` en un eje y bilineal en los otros) es lineal en los valores de los nodos. Variación global: `D_toy = f·D_nominal`. Por vóxel: funciones base del spline (`B_i(x)`, calculadas una vez por eje) por los factores del toy, que ya guarda `ToyExperiment`. Hay que mantener los dos modos (por vóxel era el objetivo; nunca se midieron incertidumbres por vóxel, por eso se usa la global). El sistemático de SCE cuesta hoy ~2 GB y ~10 s (100 mapas de ~20 MB).
- **Validación:** el servicio debe dar salida idéntica bit a bit. El SCE global por factor no saldrá idéntico (hoy el mapa escalado se guarda en `TH3F` float); se acepta con tolerancia, enseñando antes al usuario el tamaño de las diferencias.
- Pasos: (1) servicio con SCE y Calorimetry nominales + contexto, migrando las clases vivas; (2) SCE global por factor; (3) SCE por vóxel; (4) unificar `SCEVariation` y `SCEGeometricVariation`.
- **Paso 1 hecho** (commit `9f8d1d9`): `pdServiceManager` con `pd::services().SCE()` y `pd::services().Calo()`. El método se llama `Calo()` y no `Calorimetry()` porque un método con el nombre de la clase la oculta dentro de `pdServiceManager`. `Calorimetry` es inmutable, recibe un `CaloContext` (`NominalContext()` + cambios) y lee sus constantes de `pdUtils.parameters.dat`. Test idéntico. Las clases que solo usan las calibraciones (`SCECorrection`, `HitPositionSCECorrection`, `HitPitchSCECorrection`, `LifetimeCorrection`, `LifetimeVariation`, `SCEVariation`, `dQdxXCalibration`) solo se han validado compilando: no hay minitree para `dEdxCalibration.exe`.
- Hecho: quitadas las 6 clases muertas (`CalorimetryCalibration`, `RecombinationMCCorrection`, `dQdx{Cal,Norm,XCal,YZCal}Variation`) y sus `.dat`. Copia en el backup.

Encontrado en la auditoría, pendiente:

- `LifetimeCorrection` (solo en las calibraciones) nunca llama a `SetLifetime`: aplica la lifetime de MC (35 ms) también a datos. Los sistemáticos solo corren en MC, donde la lifetime es fija, así que ahí no afecta; la duda es la corrección sobre datos. Preguntar si las calibraciones se corren sobre datos.
- Campo E: `Calorimetry` usa 0,4867 kV/cm y `RecombinationDataCorrection` usa 0,553 (con sus propios ModBox). El usuario no recuerda si es a propósito. Al pasarlos a parámetros, mantener los dos valores con nombres distintos para no cambiar resultados.
- `dQdxXCalibration` lee `$HIGHLANDPDPATH/yz_corrections_mc_lifetime_mypc.root`, un fichero personal que no está en el repo.

## Fase 2 (después del test de regresión): bugs que cambian resultados

Cada uno en su propio commit, comparando micro-trees antes y después.

- `pdSystematics/src/SCEVariation.hxx:33` y `SCEGeometricVariation.hxx:30`: `SpaceCharge* _sce[100]` indexado por `toy_index`. Con `NumberOfToys > 100` se escribe fuera del array. Además crea 100 `SpaceCharge` al inicializar.
- *(Código muerto, ya eliminado.)* `pdSystematics/src/dQdxXCalVariation.cxx:43` y `dQdxYZCalVariation.cxx:44`: `dQdx * (Cal*(1+σ·v)) / Cal`. La calibración se cancela y queda un escalado global, idéntico a `dQdxNormVariation`. Confirmar con el usuario qué se pretendía (¿variación por bins en X o en YZ?).
- *(Código muerto, ya eliminado. `dEdxVariation`, que sí se usa, tiene el mismo patrón: `return` en el bucle de hits y `Chi2ndf` del muón.)* `pdSystematics/src/dQdxCalVariation.cxx:40`: `return` dentro del bucle de hits si `GetBinValues` falla. Las partículas restantes quedan sin variar y la actual a medias. En la línea 53, `Chi2ndf` se toma del muón; en las demás variaciones se toma del protón.
- `pdUtils/src/pdAnalysisUtils.cxx:48`: `TFile` global abierto en la inicialización estática con `getenv("PDUTILSROOT")` sin comprobar. Si la variable no existe, casca al cargar la librería.
- `pdUtils/src/pdAnalysisUtils.cxx:1174,1230,1320,1368`: `TFile::Open("ke_vs_range.root")` en cada llamada a `GetdEdxLikelihood*`, unas 10 veces por candidato desde `secondaryKaonAnalysis/src/kaonTree.cxx`.
- `secondaryKaonAnalysis/src/secondaryKaonAnalysis.cxx:145`: el índice `kdQdxCalibration` se registra con `dEdxVariation`. Resuelto: `dEdxVariation` sustituyó a `dQdxCalVariation` y heredó su índice. Solo queda renombrar el índice.
- **Rendimiento:** hay unas 20 instancias de `Calorimetry` y `SpaceCharge` (una por corrección o sistemático, más una por toy en SCE). Cada `Calorimetry` crea su propio `SpaceCharge` y relee los mismos ficheros (entre ellos `mc_prod4_dQdx_YZCal.txt`, de 4,4 MB). Hay que pasar a una instancia compartida.
- `Calorimetry.cxx`: constantes cableadas (E-field, ModBox A/B, v_drift, APA x) que deberían ir a `pdUtils.parameters.dat`.
- 37 `getenv` sin comprobar.

## Fase 4 (después de reestructurar highland): reestructuración de PD

- Separar el núcleo `pd*` (una librería) de los análisis (una librería y un ejecutable cada uno, opcionales vía CMake).
- Factorizar el copy-paste:
  - Las 4 variaciones `dQdx*` son casi idénticas.
  - El bloque que recalcula Chi2PID está repetido en 9 ficheros.
  - `secondaryKaonAnalysis` vs `secondaryKaonXSAnalysis`: 73 % iguales.
  - `CreateMiniTreeCandidates` vs `CreateMiniTreeEventCandidates`: 76 % iguales.
- `pdBaseAnalysis` debería registrar las correcciones y sistemáticos estándar, activables por parámetros. `pdExampleAnalysis` no usa sistemáticos.
- `EventBoxPD::FillLongTracks` usa `anaUtils::GetAllTracksUsingDet`, que hace `static_cast<AnaTrackB*>` sobre `AnaParticlePD`. Ver el `CLAUDE.md` de highland.
- `pdDataClasses`:
  - Miembros marcados como prueba u obsoletos: `Distance_to_closest_particle` (l.224), `forced_daughter`/`forced_daughter_matched` (l.241-242), `AnaHitPD::dEdx_calib`.
  - `EmptyEvent` y `HasPandora` "should probably go in DataQuality".
  - `CNNwires` está en `AnaBunchPD` y en `AnaEventPD`.
  - Las posiciones SCE de `AnaParticlePD` son `double`; el resto de la jerarquía usa `Float_t`.
- `AnaBeamPD::BeamMomentum` duplica `BeamParticle->Momentum`: `BeamMomSmearingCorrection` iguala uno al otro.
- Coordinación con highland: `EventBoxId` de highland define `kEventBoxPD`, `kEventBoxKaon` y `kEventBoxdEdx`. Cuando highland permita registrar IDs desde fuera, moverlos aquí.

## Mapa del event model

Mapa interactivo de las 56 clases de datos (highland y highlandPD), con herencia y atributos:
https://claude.ai/artifact/XNknDD8vLutMiUZ3534GoQ
