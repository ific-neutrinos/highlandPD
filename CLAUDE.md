# highlandPD: contexto para trabajar en este repo

Paquete ProtoDUNE de HighLAND. Depende de `../highland` (HighLAND base, ver su `CLAUDE.md`).
Remoto: `github.com/ific-neutrinos/highlandPD`. Rama de integración: `develop`. El usuario escribe en español.

Este fichero es corto a propósito (se carga en cada turno). El detalle está en `notes/`: **leer la nota del tema solo cuando haga falta** (no enlazarlas con `@`, que las cargaría siempre). Al apuntar algo nuevo, el detalle va a su nota y aquí solo el resumen. Versión anterior completa: `~/Public/highland/backup/highlandPD/CLAUDE.md.pre-split-20261123` y el historial de git.

## Notas (leer bajo demanda)

- `notes/analisis_produccion.md`: análisis de producción de K+ (`secondaryKaonProductionXS`, medidas A y B): estado detallado del 23-nov, resultados de selección, ajuste, viabilidad de A, macros `mu_cand_*`, `xs_*`, `kaon_*`, nota técnica, definición de la señal.
- `notes/recetas.md`: lo rutinario paso a paso, con los comandos (estado del workflow, collect, summarize, validar, persistent, análisis en FNAL, traer ficheros, lanzar producción, tags). Para el usuario y para mí.
- `notes/grid_produccion.md`: justIN y scripts de `highlandAnalysis/justin/`, workflows (23322 ... 23480), muestras en scratch y persistent, tokens, prestage, receta de mrb y compilación del módulo en FNAL, recuentos de `count_module`, criterio de cósmicos (i).
- `notes/modulo_minitrees.md`: módulo `CreateHighlandMiniTree` (flags del fcl, eventos marcados, hits), miembros nuevos de `pdDataClasses`, cómo se leen los minitrees (`Clone()`, `RedoLinks`).
- `notes/limpieza.md`: Fase 1 y Fase 2 de la limpieza, auditoría del 16-nov, cómo recuperar lo eliminado, **test de regresión** (comandos, muestras, referencias en `~/Public/highland/backup/regression_ref/`, actual `cleanup/`).
- `notes/servicios.md`: `pd::services()` (SCE, calorimetría, PID), SCE en los toys (global y por vóxel), `pd::DataFile`, masas.
- `notes/fase4.md`: reestructuración pendiente de highlandPD tras la de highland.

## Estado actual (25-nov-2026, noche; detalle en las notas)

- **Rama de trabajo:** `phase2/minitrees` en highlandPD y en `../highlandAnalysis`. Push solo cuando lo pida el usuario (lo hace él). highlandPD: lo posterior a `2757af1` (19 commits a 25-nov noche, hasta `2941b84` y siguientes) está solo en local; FNAL necesita `ca6c30a` (lista xrootd en `check_production_minitree`) y recompilar para `secondaryKaonPairXS.exe` (solo sin macros corriendo: comparten `libhighlandPD.so`). Nunca más de 3-4 procesos de ROOT por dunegpvm.
- **Tags:** `dunesw_v09_56_00d00` (tarball del 23480). **`dunesw_v09_56_00d01`**: highland `50e0989`, highlandPD `5d7f3f9`, highlandAnalysis `5cf5a24` (hits de la traza forzada en 100-600 MeV/c y tres planos, tiempo true); tarball **`0808ff92`** (con el arreglo de `pack_runtime` `4473a18`; el 23524 falló sin `libhighland.so`). Cambios del módulo para una futura `d02`: `notes/modulo_minitrees.md`.
- **Producción `d01` (MC 6 GeV, FermiGrid, 5 ficheros/trabajo):** `v1_00` = **workflow 23527, completo** (7.307/7.307, 0 fallos, 73.070 eventos, 1,763 MB/evento), recogido en `pdsp/mc/20261125v100`, copia a persistent `pdsp/dunesw_v09_56_00d01_claude/minitrees/mc/6GeV/v1_00/{minitrees,stats}` en marcha; será la **muestra de referencia** (sustituye al 23480, que se guarda hasta validarlo). **Lanzados:** `v1_01` 23534, `v1_02` 23535, `v1_03` 23536, `v1_04` 23537, `v1_05` 23538 (recoger con `collect.sh <wf> mc 20261125v1<v>`). `v1_06`-`v1_09` en cinta: **prestage de uno en uno** (decisión del usuario), empezando por `v1_06`; luego `check_locality.sh` y lanzar con el mismo tarball. Estructura de persistent: `proyecto/minitrees/mc/momento/dataset`.
- **Validación de `v1_00` en marcha (FNAL):** `count_module_minitrees.C` y `check_production_minitree.C` sobre `mc_v100d01.txt` (URLs xrootd). Después: truth skim (`run_skim_parallel.sh`) y `skim_compare.C` en local frente al skim del 23480 (`truth_skim_v100s.root`, única copia); luego `run_xs_parallel.sh` (XS y `EXE=secondaryKaonPairXS.exe`) y `kaon_acceptance`. Orden y porqués: último mensaje de plan en `notes/analisis_produccion.md` y recetas 4-7.
- **Con los microtrees del 23527, en local:** `xs_viewloss`, `xs_showers`, `xs_nminus1` (eje a 80 cm por confirmar), `technote_effpur`, `technote_cutflow`, `technote_selection`, figuras de truth, y `compare_pair_xs.C` (parejas de la tesis frente a la nuestra con la misma señal; cortes de la tesis por argumento, por defecto los de `macros/tesis/xs/selection.C`, a confirmar con el usuario; emparejamiento sin probar). Regenerar con eso los números de la nota.
- **Señal:** signal muon = mu+ de K+ -> mu+ nu de una K+ **en reposo**. Capas: C cualquiera, B KOrigin 0-3 (cascada del haz que dispara), A 0-1 (primera interacción). `muCandUtils::IsSignalMuon`, `IsSignalOf`, `CaloHits` (usa los hits de la traza forzada si la traza no tiene). Opción a (pi del haz que sobrevive, ~9-12 % de A irreducible): decisión aplazada al error estadístico de A.
- **Análisis:** `secondaryKaonProductionXS` (dos poblaciones; N-1 hecho: vértice y `IsClearCosmic` útiles) y `secondaryKaonPairXS` (selección de parejas de la tesis recuperada y comparable: cualquier especie, correcciones apagadas, misma señal). El análisis de la tesis (`secondaryKaon*`) está congelado.
- **Nota técnica:** `~/Documents/Work/DUNE/analysis/PDSP/claude/technote` (git propio, en inglés): `1_intro`, `3_truth`, `4_selection` (N-1, `sel_effpur`, cut flow en %). El usuario la está revisando y pasa feedback. Números siempre de los logs. Pendiente en el módulo (d02): primarias del pileup del haz, BI con varias trazas en datos.
- **Recetas de lo rutinario:** `notes/recetas.md` (8: varias versiones con el mismo tarball; tmux con nombre y cerrarlas al acabar).
- **Listas de minitrees en FNAL, siempre con URLs xrootd** (receta 4: `ls ... | sed 's|^/pnfs/|root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/|'`), nunca rutas `/pnfs` ni patrones con `*` sobre `/pnfs` (leerían por NFS). Se nos pasó el 24 y el 25-nov. Comandos para pegar en FNAL: cada orden en una sola línea (la terminal parte las largas).

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

- **Núcleo:** `pdEventModel` (`pdDataClasses`), `pdIO` (`PDSPAnalyzerTreeConverter`, `MichelRemovingTreeConverter`, `pdBaseConverter`; `pdMiniTreeConverter` no lo usa nadie: los minitrees se leen con `HighlandMiniTreeConverter` de highland), `pdUtils` (`Calorimetry`, `SpaceCharge`, `pdAnalysisUtils`, `standardPDTree`), `pdCorrections` (11 correcciones), `pdSystematics` (variaciones y pesos), `pdBaseAnalysis`.
- **Análisis:** `secondaryKaonAnalysis` (análisis de la tesis `secondaryKaon`, `secondaryKaonProductionXS` y `secondaryProton`; incluye `coherentfit/`, que no se compila), `dEdxCalibrationAnalysis`, `stoppingProtonAnalysis`, `brokenTracksAnalysis`, `beamWeightsAnalysis`, `pdExampleAnalysis`.
- **Convención de nombres (la de highland, 20-nov-2026):** paquetes y clases de análisis y selección en minúscula (`beamWeightsAnalysis`, `beamWeightsSelection`); fichero de parámetros `<paquete>.parameters.dat` y prefijo de los parámetros = nombre del paquete (highland busca `<prefijo>.parameters.dat` en `$<PREFIJO>ROOT/parameters`: con otro prefijo no los encuentra); ejecutable = nombre del análisis sin `Analysis` (`secondaryKaon.exe`, `beamWeights.exe`, `pdExample.exe`); herramientas en mayúscula como en highland (`CreateMiniTreeKaon.exe`). `secondaryKaonXS` se llama ahora `secondaryKaonProductionXS` (es el análisis de producción).
- **No se compila:** `secondaryKaonAnalysis/src/coherentfit`. `CreateMiniTreePD.exe` y `CoherentFit.exe` están comentados en el CMake.
- **Solo en local, ignorado por git:** `src/secondaryKaonAnalysis/macros/tesis/` (macros de la tesis del usuario: no borrarlas del disco).

## Mapa del event model

Mapa interactivo de las 56 clases de datos (highland y highlandPD), con herencia y atributos:
https://claude.ai/artifact/XNknDD8vLutMiUZ3534GoQ
