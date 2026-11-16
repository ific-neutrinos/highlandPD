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
test/regression/run_regression.sh <dir_mc>   <mc.root> -          0       # todo el MC, ~12-14 min
test/regression/run_regression.sh <dir_data> -         <data.root> 10000  # ~1,5 min
root -l -b -q 'test/regression/compare.C+("ref/mc.root","<dir_salida>/mc.root")'
```

- MC: `~/Documents/Work/DUNE/analysis/PDSP/secondary_kaon_analysis/files/minitrees/mc/6GeV_prod4a_00_minitree_2023-01-27.root` (5 GB, 47.906 eventos). Con 10.000 eventos: ~3 min, 3.371 eventos guardados.
- Datos: `~/Documents/Work/DUNE/analysis/PDSP/secondary_kaon_analysis/files/minitrees/data/6GeV_prod4a_reco2_minitree_2023-01-27.root` (42 GB, 446.262 eventos). Con 10.000 eventos: ~1,7 min, 2.929 eventos guardados. En datos los toys no varían (los sistemáticos solo actúan en MC) y `truth` está vacío.
- Fase 1 validada (13-nov-2026): con 10.000 eventos de MC, la salida de `archive/pre-cleanup` y la de la rama son idénticas en todos los árboles (`ana`, `all_syst`, `truth`, `header`). En el log solo cambian los prints de depuración quitados.
- El override activa también `RecombinationVariation`, el único sistemático del análisis de kaones que usa `Calorimetry`.
- Referencias en `~/Public/highland/backup/regression_ref/`. **Actual: `masses/`** (todo el MC y 10.000 eventos de datos), tras unificar las masas (PDG 2022). Frente a `full/`: cambian el momento por rango (muón ~1e-4 relativo, protón ~1e-7), las likelihoods de dE/dx (muón hasta 4e-3, protón ~6e-6) y `trueendmom_atAPA` (hasta 0,16 MeV); la selección nominal es idéntica y en `all_syst` el corte 2 (ventana en el momento por rango de la hija) cambia en 5 toys de 4 eventos que estaban en el borde. Anteriores:
  - `full/{mc,data}.root` (todo el MC: 16.113 entradas, 62 con variables de toys; datos: 10.000 eventos): tras los pasos 3 y 4 del SCE; frente a `archive/pre-cleanup` en todo el MC, `ana` idéntico y `all_syst` a 1 ULP de float.
  - `recomb/{mc,data}.root`: generadas con `archive/pre-cleanup`. El servicio de SCE y calorimetría da idéntico.
  - `sce_scaling/{mc,data}.root`: tras el paso 2 (10.000 eventos de MC). Frente a `recomb`, `ana` y datos son idénticos y `all_syst` en MC difiere en 1 ULP de float (`compare.C+(a,b,true,1e-6)` pasa). La de `fase1/` es la anterior, sin recombinación. Cada bug de la Fase 2 debe compararse contra la referencia anterior, y las diferencias deben explicarse por el arreglo.

## Servicios de calorimetría y SCE (en curso, 13-nov-2026)

Decidido con el usuario (primero el bloque de diseño, antes que los bugs sueltos):

- **`pdServiceManager`** en `pdUtils`, con acceso global `pd::services()` al estilo de `ND::params()`. Tiene que ser global porque las correcciones y los sistemáticos no tienen acceso a `ana()`. Usa métodos con tipo (`SCE()`, `Calorimetry()`), creación perezosa (nunca en la inicialización estática), objetos `const` y lista de servicios en el log.
- **Auditoría:** todas las `Calorimetry` y `SpaceCharge` tienen la misma configuración base. Solo cambia estado por evento o por toy: lifetime (`LifetimeVariation`, `SCEVariation`), ModBox A/B (`RecombinationVariation`) y qué SCE usar. Basta con una instancia nominal de cada una, y ese estado viaja como contexto en cada llamada.
- **SCE en los toys sin copias del mapa.** La interpolación (spline natural `"b2e2"` en un eje y bilineal en los otros) es lineal en los valores de los nodos. Variación global: `D_toy = f·D_nominal`. Por vóxel: funciones base del spline (`B_i(x)`, calculadas una vez por eje) por los factores del toy, que ya guarda `ToyExperiment`. Hay que mantener los dos modos (por vóxel era el objetivo; nunca se midieron incertidumbres por vóxel, por eso se usa la global). El sistemático de SCE cuesta hoy ~2 GB y ~10 s (100 mapas de ~20 MB).
- **Validación:** el servicio debe dar salida idéntica bit a bit. El SCE global por factor no saldrá idéntico (hoy el mapa escalado se guarda en `TH3F` float); se acepta con tolerancia, enseñando antes al usuario el tamaño de las diferencias.
- Pasos: (1) servicio con SCE y Calorimetry nominales + contexto, migrando las clases vivas; (2) SCE global por factor; (3) SCE por vóxel; (4) unificar `SCEVariation` y `SCEGeometricVariation`.
- **Paso 2 hecho** (commit `5ae2d90`): `SCEScaling` (por ahora solo `Global`) como argumento opcional de las consultas de `SpaceCharge`, y `CaloContext::SCEScale`. Las variaciones calculan el factor con `GetToyScaling(toy)` y ya no hay mapas por toy ni límite de 100 toys. Memoria 2,44 GB -> 364 MB, inicialización ~9 s -> ~1 s. Se quitó la API que modificaba el mapa y la variación por vóxel antigua (`VarySCEMapLocally`, no usada, que necesitaba una copia del mapa por toy); el paso 3 la reimplementa sin copias.
- **Paso 1 hecho** (commit `9f8d1d9`): `pdServiceManager` con `pd::services().SCE()` y `pd::services().Calo()`. El método se llama `Calo()` y no `Calorimetry()` porque un método con el nombre de la clase la oculta dentro de `pdServiceManager`. `Calorimetry` es inmutable, recibe un `CaloContext` (`NominalContext()` + cambios) y lee sus constantes de `pdUtils.parameters.dat`. Test idéntico. Las clases que solo usan las calibraciones (`SCECorrection`, `HitPositionSCECorrection`, `HitPitchSCECorrection`, `LifetimeCorrection`, `LifetimeVariation`, `SCEVariation`, `dQdxXCalibration`) solo se han validado compilando: no hay minitree para `dEdxCalibration.exe`.
- **Pasos 3 y 4 hechos** (commits `1ae7c28`, `56dbba0`):
  - `SCEVariationBase` reúne lo común (mapa del servicio y `GetToyScaling`). Las dos derivadas mantienen efectos distintos a propósito: `SCEGeometricVariation` (kaones) solo cambia la geometría, porque la calorimetría la cubre `dEdxVariation`; `SCECalorimetricVariation` (antes `SCEVariation`, calibraciones) cambia las posiciones de los hits y el pitch de dQ/dx. No unificarlas.
  - Modo con `pdSystematics.SCE.VariationMode` en `src/pdSystematics/parameters/pdSystematics.parameters.dat`: 0 global (`data/SCE.dat`, sigma 8 %), 1 por vóxel (`data/SCE2.dat`, 21.793 cajas con sigma 5 % de prueba). El modo por vóxel usa un lector propio porque `BinnedParams` tiene `NMAXSYSTSOURCEBINS = 50`.
  - Por vóxel, `SpaceCharge` evalúa el spline con la base cardinal (`_basis`) y los factores de cada nodo, sin copias. Validado con `root -l -b -q test/sce/run_test_sce_voxel.C` frente al método antiguo (copia escalada y splines nuevos): diferencias relativas <= 4e-8.
  - Coste del modo por vóxel en `secondaryKaon.exe` (1000 eventos, 10 toys): 52 s frente a 22 s en global; memoria igual.
  - `SCE2.dat` tiene cajas independientes para los dos volúmenes (42.439); los nodos de x = 0 (cátodo) comparten caja. La implementación antigua repartía los factores por índice de bin entre volúmenes, casi seguro sin querer.
- Las variables de toys solo se rellenan en eventos con super candidato (el corte 3 borra los candidatos si no lo hay): por eso el test usa todo el MC.
- Hecho: quitadas las 6 clases muertas (`CalorimetryCalibration`, `RecombinationMCCorrection`, `dQdx{Cal,Norm,XCal,YZCal}Variation`) y sus `.dat`. Copia en el backup.

Encontrado en la auditoría, pendiente:

- `LifetimeCorrection` (solo en las calibraciones) usa siempre la lifetime de MC. Según el usuario, las calibraciones no se corren sobre datos (solo hay un rastro comentado en `dQdxXCalibration.cxx:90`), así que no se cambia.
- Campo E: `Calorimetry` usa 0,4867 kV/cm y `RecombinationDataCorrection` usa 0,553 (con sus propios ModBox). El usuario no recuerda si es a propósito. Al pasarlos a parámetros, mantener los dos valores con nombres distintos para no cambiar resultados.
- `dQdxXCalibration` lee `$HIGHLANDPDPATH/yz_corrections_mc_lifetime_mypc.root`, un fichero personal que no está en el repo.

## Servicios y ficheros de datos

- `pd::services()` (`pdUtils/src/pdServiceManager`): `SCE()`, `Calo()` y `PID()` (`PIDTables`: plantillas de dE/dx por PDG, energía cinética frente a rango por partícula y tablas NIST de muón y protón; un único servicio porque todo depende de la hipótesis de partícula). Se crean la primera vez que se piden, también dentro del bucle de eventos; el guardia restaura `gDirectory` y evita que los histogramas se cuelguen del directorio actual.
- `pd::DataFile(paquete, fichero)` y `pd::DataDir(paquete)` (`pdUtils/src/pdDataFiles`): único punto para las rutas de datos. `$<PKG>ROOT/data`, o antes `$HIGHLANDPD_DATADIR/<paquete>` si está definida (para tener los datos fuera del repo). Error claro y `exit(1)` si falta algo. No usar `getenv` para rutas de datos.
- Nombres de fichero en `pdUtils.parameters.dat`: `pdUtils.SCE.MapFile`, `pdUtils.Calorimetry.{Norm,XCal,YZCal}File`, `pdUtils.PID.{dEdxTemplatesFile,KEvsRangeFile}`.
- En highland y psyche no hay nada equivalente (solo los singletons `ND::params()` y `ND::versioning()`, y `BinnedParams`). `BinnedParams` hace `exit(0)` si no encuentra su fichero: anotar para la reestructuración de highland. El gestor de servicios y `pd::DataFile` podrían subir a highland como algo genérico.
- Masas: constantes `pdAnaUtils::k{Muon,Pion,Kaon,Proton}Mass` (PDG 2022) y `pdAnaUtils::GetParticleHypothesis(pdg, nombre, masa)`. Las funciones que dependen de la hipótesis avisan por `std::cerr` y devuelven -999 (o -1 en `ComputeCSDARange`) si el PDG no está soportado; no usar valores locales.
- Pendiente: `dQdxXCalibration` lee `$HIGHLANDPDPATH/yz_corrections_mc_lifetime_mypc.root`, que no está en el repo. Las macros de `dEdxCalibrationAnalysis` tienen rutas personales de Fermilab.

## Fase 2 (después del test de regresión): bugs que cambian resultados

Cada uno en su propio commit, comparando micro-trees antes y después. Es lo último que queda en highlandPD antes de pasar a highland.

Resuelto:

- `SpaceCharge* _sce[100]` por toy en las variaciones de SCE: ya no hay mapas por toy (servicio de SCE, paso 2).
- `TFile` global abierto con `getenv` en la inicialización estática de `pdAnalysisUtils` y `TFile::Open("ke_vs_range.root")` en cada `GetdEdxLikelihood*`: ahora lo da `pd::services().PID()`.
- Unas 20 instancias de `Calorimetry` y `SpaceCharge`: una compartida de cada una (`pd::services()`).
- Constantes de `Calorimetry` en `pdUtils.parameters.dat`.
- Masas de las partículas unificadas (commit `b3f8885`).
- `dEdxVariation`: la σ y la variación del toy se leen una vez antes del bucle (había un `return` dentro del bucle de hits que dejaba el evento a medias); si no se puede leer, error y `exit(1)`. Sin cambio de resultados: `dEdx.dat` tiene un solo bin y se consulta siempre en 0,5. Dentro de un toy es un escalado coherente de todos los hits (`1 + 0,03·N(0,1)`, una variación por toy).
- `Chi2ndf` tomado del muón: no es un bug. Todas las variaciones y el conversor lo hacen igual, y las cuatro plantillas de dE/dx tienen el mismo binning (325 bins en 0-26 cm), así que el ndf no depende de la hipótesis.
- `getenv`: solo queda uno vivo sin comprobar, en `dQdxXCalibration.cxx:91` (fichero personal, ver arriba). Los demás son de `pdDataFiles` o están comentados.

Pendiente (sin cambio de resultados):

- Pasar los mensajes de error de `std::cout` a `std::cerr`.
- Textos de `Calorimetry`: los de YZ dicen "X correction".

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
