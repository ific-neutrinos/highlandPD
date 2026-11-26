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

## Estado actual (26-nov-2026, tarde; detalle en las notas)

- **Rama de trabajo:** `phase2/minitrees` en highlandPD y en `../highlandAnalysis`. Push solo cuando lo pida el usuario (lo hace él). Solo en local: highlandPD desde `fadb0e1` (26-nov: eje a 80 cm `df3cec0`, truth de la medida C `675d80c`, macros `opt_*`, notas) y highlandAnalysis `ed27a81` (`to_persistent.sh` solo con token, sin el proxy de justIN). Nunca más de 3-4 procesos de ROOT por dunegpvm; no recompilar en FNAL con macros corriendo (comparten `libhighlandPD.so`).
- **Tags:** `dunesw_v09_56_00d01`: highland `50e0989`, highlandPD `5d7f3f9`, highlandAnalysis `5cf5a24`; tarball **`0808ff92`**. Cambios del módulo para `d02`: `notes/modulo_minitrees.md`.
- **Producción `d01` (MC 6 GeV):** `v1_00` 23527 (referencia, validada: igual que el 23480; truth idéntica), `v1_01`-`v1_05` 23534-23538 (todo bien, 366.460 eventos en total, `pdsp/mc/20261125v10<v>`). `v1_00` y parte de `v1_01` en persistent; copia de `v1_01`-`v1_05` en marcha (tmux `persist`, en serie; con `BEARER_TOKEN_FILE=/tmp/bt_interactive_u<uid>` y `X509_USER_PROXY=/tmp/no_existe`, receta 5). **`v1_06` = 23551, `v1_07` = 23552 corriendo** (recoger en `pdsp/mc/20261126v10<v>`); prestage de `v1_08` y `v1_09` en marcha.
- **En local** (`~/Documents/Work/DUNE/analysis/PDSP/claude/files/minitrees/mc`): `xs_v10<v>d01.root` (`v1_00`-`v1_05`) y `xs_v100to105d01.root` (juntado); `pair_v100d01`, `pair_v101d01`, `pair_v102d01`; `truth_skim_v100d01.root`; `opt_flat2_v100to105d01.root`. En FNAL: `pair_v103`, `v105` corriendo, `pair_v104d01.root` hecho (traer y comprobar entradas por evento frente a `v1_01`/`v1_02`).
- **Pendiente en FNAL (cuando acaben los `pair`):** `git pull` + recompilar (eje a 80 cm, truth de C) y volver a correr `secondaryKaonProductionXS` sobre `v1_00`-`v1_05` (`run_xs_parallel.sh`, una máquina por bloque).
- **Selección:** cortes actuales en su óptimo (reoptimizarlos no gana). Un **BDT sin variables ligadas al momento** (sin `trmean`, `chi2*`) da dentro B 4,0 % / A 7,5 % (cortes 6,1 / 12,2) y fuera B 4,3 / A 9,8 % (11,6 / 31), sin esculpir el fondo; duplica la señal antes de la ventana. Lo que aporta: inicio/K (carga media por hit al inicio del mu, idea del usuario; `k_dedx`) y Michel; fuera, las de cósmicos. **BDT reducido** de 8 variables: dentro 4,6 / 8,5 %. **Buscando el punto dulce** (`opt_fwd.py`, selección hacia delante; `opt_fwd_in.log`, `opt_fwd_out.log`): el usuario prefiere un BDT pequeño con variables físicas, validable una a una en datos.
- **Plan acordado:** punto dulce del BDT → **medida simulada, empezando por B**: `v1_00`-`v1_05` como MC (entrenamiento, eficiencias, plantillas), `v1_06`-`v1_09` como pseudo-datos; cierre, pulls, estabilidad, dentro frente a fuera → datos de 6 GeV → otras energías (3, 7 GeV). Luego A (opción a pendiente: con error estadístico de A ~12 % con cortes, el 9-12 % irreducible pesa) y C (¿etiqueta del mu + K al inicio, o parejas? estudiar con la truth nueva). Medida B por especie: BI, ajuste del pico por población, eficiencia, restar el pileup (C no B), pureza del PID.
- **Señal:** signal muon = mu+ de K+ -> mu+ nu de una K+ **en reposo**. Capas: C cualquiera, B KOrigin 0-3, A 0-1. Pileup del haz (KOrigin 6, ~40 % de los signal muons): fuera de tiempo, en [-3, 2) ms; con el BI podrían usarse (A x1,7, B x1,66; `skim_pileup_layers.C`), a ver en datos.
- **Nota técnica:** `~/Documents/Work/DUNE/analysis/PDSP/claude/technote` (git propio, en inglés). `4_selection` regenerado con `v1_00`-`v1_05` (`99bc899`); falta reflejar el eje a 80 cm (figuras y cut flow tras volver a correr) y el BDT. Números siempre de los logs.
- **Recetas:** `notes/recetas.md`. Listas de minitrees en FNAL con URLs xrootd y sin glob (`ls $D | grep ... | sed`; con 11.626 ficheros el glob falla). Comandos para pegar en FNAL: cada orden en una sola línea.

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
