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

## Estado actual (24-nov-2026, tarde; detalle en las notas)

- **Rama de trabajo:** `phase2/minitrees` en highlandPD y en `../highlandAnalysis`. Push solo cuando lo pida el usuario (el 24-nov lo hizo él hasta `c0c0a92`; lo posterior está solo en local). Tag `dunesw_v09_56_00d00` en local (highland `50e0989`, highlandPD `ca2e848`, highlandAnalysis `bc6e41d`): tarball del 23480. **Tag `dunesw_v09_56_00d01` (25-nov, en FNAL): highland `50e0989`, highlandPD `5d7f3f9`, highlandAnalysis `5cf5a24`**: tarball de `v1_00`-`v1_09` (hits de la traza forzada en 100-600 MeV/c y tres planos, tiempo true), subido; `v1_00`: el workflow 23524 falló (tarball sin `libhighland.so`, arreglado en highlandAnalysis `4473a18`); relanzado como **workflow 23527** (25-nov, tarball `0808ff92`, prueba 23526 bien). El push de los tags lo hace el usuario.
- **Producción:** 23480 (`v1_00`, referencia) recogido en `pdsp/mc/20261124v100s` y copiado a persistent; truth skim y `xs_v100s.root` en local (`.../files/minitrees/mc/`). **Plan acordado (24-nov), solo MC:** (1) acabar las macros de dunegpvm05 (`count_module`, `kaon_origin`, `kaon_acceptance`, `pimu`); (2) push (highlandPD hasta `0ae5460`, highlandAnalysis `04cbc5e`: tiempo true en el módulo); (3) relanzar el análisis sobre el 23480 (`xs_v100s_fix.root`) y decidir con él, antes de `v1_01`: hits de la traza forzada de los showers candidatos (`SaveForcedTrackHits`) y cuánta señal pierde el resumen (`truemu_level_full`); (4) módulo, tarball, prueba en la grid, tag `dunesw_v09_56_00d01`; (5) `v1_00`-`v1_09` por bloques: **`v1_00` se relanza** con el mismo tarball (el módulo cambia), y el nuevo será la referencia; el 23480 se guarda hasta validarlo. Nunca más de 3-4 procesos de ROOT por dunegpvm (el 24-nov se colapsó la 05); no recompilar highlandPD con macros corriendo (comparten `libhighlandPD.so`).
- **Señal (decidida el 24-nov):** signal muon = mu+ de K+ -> mu+ nu de una K+ **en reposo** (momento final de la K exactamente 0 en Geant4; el canal se identifica por p_mu, no es un corte). Señal en capas: C cualquier signal muon, B la K viene de la cascada del haz que dispara (KOrigin 0-3), A de la primera interacción (0-1). En el código: `muCandUtils::IsSignalMuon`, `IsSignalOf`; categoría `candmu` en capas excluyentes (1 A, 2 B no A, 3 C no B) y `candsigC/B/A` (`0e2fef3`, `5dade5b`). Un mu con la madre sin guardar no es señal (`5166129`).
- **Análisis y selección:** `secondaryKaonProductionXS`, dos poblaciones. Cambios del 24-nov (detalle en `notes/analisis_produccion.md`): vista normal de los eventos con volcado (sin escalar el fondo), extremo A referido a `RRMax`, variantes del chi2 de mu, Michel en hits propios, ventana de p de rango al final de la cadena (`sel_fitlevel`, `sel_countlevel`). Revisión corte a corte dentro de la beam slice hecha: umbrales en el óptimo (barrido conjunto con media muestra: no mejora). Abiertos: rama de showers, escala del rango (pico en 224 MeV/c), separación pi/mu (fondo dominante), población de fuera. Figura de mérito: eff x pur; para B, el error del ajuste (`xs_fitB.C`). Macros del skim en local: `test/phase2/skim_*.C`. El análisis de la tesis (`secondaryKaon*`) está congelado.
- **Nota técnica:** `~/Documents/Work/DUNE/analysis/PDSP/claude/technote` (git propio, en inglés). `1_intro`, `3_truth` (ventana del mu y pileup del 23480 añadidos el 24-nov, `d3767be`) y borrador de `4_selection` (con la vista normal). Números siempre de los logs. Pendientes en `notes/analisis_produccion.md`: rerun de la aceptancia y del reparto por especie sobre el 23480, estudio del pi del haz que sobrevive, tiempo del pileup (con `v1_01`), versión de Geant4 de dunesw v09_56.
- **Recetas de lo rutinario:** `notes/recetas.md`.

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
