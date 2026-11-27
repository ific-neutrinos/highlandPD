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

## Estado actual (27-nov-2026, mañana; detalle en las notas)

- **Rama de trabajo:** `phase2/minitrees` en highlandPD y `../highlandAnalysis`. Push solo cuando lo pida el usuario (lo hace él). **Sin subir:** highlandPD desde `dc819c9` (macros `measB_*`, notas; lo de antes ya está en el remoto) y highlandAnalysis `ed27a81`, `a07b356`, `7b42292` (`to_persistent.sh`: solo token, sin proxy, token y vault propios, avance al saltar). Los `xs_v10Xd01b` se hicieron con highlandPD hasta `dc819c9`. Nunca más de 3-4 procesos de ROOT por dunegpvm; no recompilar con macros corriendo.
- **Tags:** `dunesw_v09_56_00d01` (highland `50e0989`, highlandPD `5d7f3f9`, highlandAnalysis `5cf5a24`), tarball `0808ff92`. Cambios del módulo para `d02`: `notes/modulo_minitrees.md`.
- **Producción `d01` (MC 6 GeV):** `v1_00`-`v1_05` (23527, 23534-23538) recogidos y validados (366.460 eventos); copia a persistent en marcha (máquina dedicada, receta 5; `v1_01`, `v1_02` ya estaban). `v1_06` 23551 y `v1_07` 23552 recogidos y validados (`pdsp/mc/20261126v10<v>`, 116.120 y 69.840 eventos); `v1_08` = 23562 recogido y validado (`pdsp/mc/20261127v108`, 5.835 ficheros, 58.350 eventos, 0 fallos); `v1_09` = 23563 recogiéndose.
- **`xs_v10Xd01b.root`** (`v1_00`-`v1_05`, con `beam_quality_sce`) en local, juntados en `xs_v100to105d01b.root` (mismo orden de eventos que los `xs_*d01`). Capítulo 5 de la nota rehecho con la calidad SCE (27-nov). Pendiente: el mismo análisis sobre `v1_06`-`v1_09` (pseudo-datos).
- **En local** (`~/Documents/Work/DUNE/analysis/PDSP/claude/files/minitrees/mc`): `xs_v10Xd01.root` y `xs_v100to105d01.root`, `pair_v10Xd01.root` (los seis), `truth_skim_v100d01.root`, `opt_flat3_v100to105d01.root` (candidatos), BDT final dentro `opt_bdt_final/`.
- **Selección: BDT** (decisión del usuario). Dentro, 10 variables (`as_mean5,k_dedx,k_n,cnn_michel,michel_own,nmich10,cosbeam,dedx_stop,ar_charge10,intpeak`, umbral 0,19): B 4,2 % / A 7,7-8,0 % (cortes 6,1 / 12,2). Fuera, punto dulce 13 variables (+ `axis,clear,cosmiccrit,slice_min,ar_charge10`): 4,9 / 10,5 % (`cosmiccrit` aparcado). Sin variables ligadas al momento (no esculpe el fondo). Falta describirlo en el capítulo de selección de la nota y validar sus entradas en datos.
- **Medida B** (nota, `5_measB`): por partícula del haz que pasa la calidad del haz (con SCE: `beam_quality_sce`, desde el 27-nov en todo el capítulo), solo en la beam slice de momento; numerador = mu de K de la cascada del haz en argón; resultado por pi+ con la fracción de pi del grupo e/mu/pi (83 % con la calidad). Ajuste: forma del mu de K + forma del pi de **K+ -> pi+ pi0 en reposo** (pico en 160-180 MeV/c, era el "hueco" del fondo) del MC + exp(pol2). Con el BDT: sesgo +3,5 % (0,4 sigma, mismo MC), error ~6 % en datos (e/mu/pi), ~12 % (K).
- **Denominador de B por pi+** (opción b del usuario): fracción de pi medida en el TPC (`measB_beampid.C`; e por CNN EM > 0,8, mu por `beam_reach`, variable nueva del 27-nov). Correr el análisis otra vez en FNAL (v1_00-v1_09, `xs_*d01c`) para tenerla.
- **Medida C** (parejas K -> mu para la elástica de K que paran): mu del BDT + corte en p + madre (de Pandora o forced mother, `cand_fmc_*`/`cand_fmk_*`, `310615b`) frente al buscador de parejas de la tesis; `measC_compare.C`, `measC_fm.C`. Necesita la pasada nueva de FNAL.
- **Pasada de FNAL acordada (27-nov), cuando el usuario diga que 23562 (`v1_08`) y 23563 (`v1_09`) han acabado:** (1) recoger, resumir y validar `v1_08`/`v1_09` y hacer las listas xrootd de `v1_06`-`v1_09`; (2) push (usuario) de highlandPD hasta `d6a3539` o posterior y de highlandAnalysis; en FNAL pull y recompilar; (3) `secondaryKaonProductionXS` sobre los 10 bloques -> `xs_v1XXd01c.root` (beam_reach, forced mother, dispersiones, interacción del haz), repartido en máquinas, máx. 3-4 ROOT por dunegpvm; (4) `secondaryKaonPairXS` sobre `v1_04` (vacío) y `v1_06`-`v1_09`.
- **Siguiente:** (1) hecho: calidad SCE en B; (2) medida simulada de B con `v1_06`-`v1_09` como pseudo-datos (cierre, pulls); (3) medida C (¿etiqueta del mu + K al inicio o parejas?, con la truth nueva); (4) explorar recuperar pileup y eventos sin beam slice (actividad en la región de entrada del haz).
- **Decisiones del usuario (26-nov):** B solo en la beam slice por ahora; el p se queda; selección = BDT; señal de A con el haz en el FV; opción a sí (marcada); `cosmiccrit` aparcado.
- **Señal:** signal muon = mu+ de K+ -> mu+ nu de una K+ en reposo. Capas: C cualquiera, B KOrigin 0-3, A 0-1. Pileup del haz (KOrigin 6) fuera de tiempo; con el BI con varias trazas (d02) podría usarse (A x1,7).
- **Nota técnica:** `~/Documents/Work/DUNE/analysis/PDSP/claude/technote` (git propio, en inglés). `4_selection` con `v1_00`-`v1_05` (cortes); `5_measB` empezado (muestra, origen, definición, calidad del haz, composición, modelo de ajuste, pico Kpipi, versión BDT). Números siempre de los logs. Figuras apiladas con la señal abajo.
- **Recetas:** `notes/recetas.md`. Listas de minitrees en FNAL con URLs xrootd y sin glob. Copia a persistent: receta 5 (token y vault propios, sin proxy). Comandos para FNAL en una sola línea.

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
