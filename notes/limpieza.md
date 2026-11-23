# Limpieza de highlandPD: Fases 1 y 2, auditoría y test de regresión

Movido de CLAUDE.md el 23-nov-2026 (sin cambios de texto).

## Recuperar lo eliminado

- Tag `archive/pre-cleanup` (local, sin subir): estado de `develop` (`cc1a931`) antes de la limpieza.
- Copia en disco: `~/Public/highland/backup/highlandPD/`, con la misma estructura de carpetas que el repo. Incluye `LArSoftReader/`, `QPIXAnalysis/`, `unusedCodeToRevisit/`, `macros/tesis/` (con los cambios locales sin commitear que tenía el usuario), el `.pb` del CNN y los ficheros sueltos borrados.

## Fase 1: limpieza (sin cambios de comportamiento)

Hecho en la rama `cleanup/highlandPD` (13-nov-2026), sin push:

- `.gitignore` (productos de compilación, autosaves de emacs, `nohup.out`).
- Quitados del repo: `LArSoftReader/`, `QPIXAnalysis/` (decisión del usuario: sale de highlandPD), `unusedCodeToRevisit/` (también borrado del disco) y `macros/tesis/` (sigue en disco), `nohup.out`, `coherentfit/oldfits`, dos autosaves de emacs, `dEdxCalibrationAnalysis/macros/x.C` (idéntico a `yz.C` y definía `yz()`), los 4 módulos CMake de `/cvmfs` y el modelo CNN `.pb` de MCC11 junto con el parámetro `pdUtils.CNN.PBfile`, que no leía nadie. Los ficheros borrados siguen en el historial de git, así que el `.git` no adelgaza.
- Borrados de disco (no estaban en git): `CMakeLists_backup.txt`, `scripts/setup_backup.sh`, `bin/CoherentFit.exe` (binario viejo).
- Añadido `pdUtils/data/dedx_mpv_profiles.root`, que estaba sin seguir; ningún código lo referencia todavía.
- Quitados los prints de depuración del bucle de eventos: `MIGUEEEE...` en `GetdEdxLikelihood`, el de cada hit en `dEdxLikelihood` y los contadores de `DaughterFinderCorrection::Apply`. El resto de `std::cout` se queda: métodos `Print()`, mensajes de error y mensajes de `Initialize`/`Finalize`. (Los errores ya van por `std::cerr`.)
- Rutas `/home/migue/` de `pdUtils/data/likelihood_microtree.C` y `secondaryKaonAnalysis/macros/selfsysttest.C` pasadas a argumentos.
- CMake: función para los ejecutables, sin la cadena de `add_dependencies`, sin `-ansi`, `cmake_minimum_required(3.12)`. `ENABLE_EXPORTS` mantiene el `-rdynamic` que CMake 3.0 añadía siempre (CMP0065). Flags comprobados iguales a los de antes.

Pendiente:

1. **Ficheros de datos de las utilidades** (`pdUtils/data/`, varios MB): decidir cómo gestionarlos (git-lfs, área externa...). El usuario quiere pensarlo más adelante.

## Test de regresión

`test/regression/` (override `regression.parameters.dat`, `run_regression.sh`, `compare.C`). Usar minitrees **sin filtrar**: los de candidatos (`minitree_candidates_*`) no traen las hijas de los candidatos. Ya no hacen caer `secondaryKaon.exe` (el `Daughters[0]` de `EventBoxKaon` se arregló en `389e350`; comprobado con el de MC, 1.503 eventos), pero la selección no sirve para comparar. `check.sh <ref> <salida> <mc> <data>` corre MC y datos en paralelo y compara.

```bash
test/regression/run_regression.sh <dir_mc>   <mc.root> -          0       # todo el MC, ~8-9 min
test/regression/run_regression.sh <dir_data> -         <data.root> 10000  # ~1,5 min
root -l -b -q 'test/regression/compare.C+("ref/mc.root","<dir_salida>/mc.root")'
```

- MC: `~/Documents/Work/DUNE/analysis/PDSP/secondary_kaon_analysis/files/minitrees/mc/6GeV_prod4a_00_minitree_2023-01-27.root` (5 GB, 47.906 eventos). Con 10.000 eventos: ~3 min, 3.371 eventos guardados.
- Datos: `~/Documents/Work/DUNE/analysis/PDSP/secondary_kaon_analysis/files/minitrees/data/6GeV_prod4a_reco2_minitree_2023-01-27.root` (42 GB, 446.262 eventos). Con 10.000 eventos: ~1,7 min, 2.929 eventos guardados. En datos los toys no varían (los sistemáticos solo actúan en MC) y `truth` está vacío.
- Fase 1 validada (13-nov-2026): con 10.000 eventos de MC, la salida de `archive/pre-cleanup` y la de la rama son idénticas en todos los árboles (`ana`, `all_syst`, `truth`, `header`). En el log solo cambian los prints de depuración quitados.
- El override activa también `RecombinationVariation`, el único sistemático del análisis de kaones que usa `Calorimetry`.
- Referencias en `~/Public/highland/backup/regression_ref/`. **Actual: `cleanup/`** (Fase 2 de `phase2.md`): igual que `fixA/` salvo que ya no están las ramas `*_vtx_michelscore` y `*_vtx_nhits` de `kaonTree`. Anterior: `fixA/` (16-nov), tras los arreglos A1–A4 y A6 de la auditoría y quitar `*_hit_dedx_cal`. Frente a `lklwin/`: MC con la misma selección, cambian `beam_quality` (677 entradas, cortes de 6 GeV) y `bestcandidateparticlereduced` (30, ahora 3210); en datos cambian también los χ² del mejor candidato y su hija en 11 eventos (hasta un 18 %), sin cambio de selección. Anteriores:
  - `lklwin/` (todo el MC y 10.000 eventos de datos), tras aplicar el filtro 1 < dE/dx < 100 en todas las likelihoods. Frente a `masses/` solo cambia `bestcandidate_lkl_muon` en 5 eventos de MC (hits fuera de la ventana con densidad casi nula para el muón; para protón y kaón ya valía 0 y se saltaban); datos idénticos. Anteriores:
  - `masses/`: tras unificar las masas (PDG 2022). Frente a `full/`: cambian el momento por rango (muón ~1e-4 relativo, protón ~1e-7), las likelihoods de dE/dx (muón hasta 4e-3, protón ~6e-6) y `trueendmom_atAPA` (hasta 0,16 MeV); la selección nominal es idéntica y en `all_syst` el corte 2 (ventana en el momento por rango de la hija) cambia en 5 toys de 4 eventos que estaban en el borde. Anteriores:
  - `full/{mc,data}.root` (todo el MC: 16.113 entradas, 62 con variables de toys; datos: 10.000 eventos): tras los pasos 3 y 4 del SCE; frente a `archive/pre-cleanup` en todo el MC, `ana` idéntico y `all_syst` a 1 ULP de float.
  - `recomb/{mc,data}.root`: generadas con `archive/pre-cleanup`. El servicio de SCE y calorimetría da idéntico.
  - `sce_scaling/{mc,data}.root`: tras el paso 2 (10.000 eventos de MC). Frente a `recomb`, `ana` y datos son idénticos y `all_syst` en MC difiere en 1 ULP de float (`compare.C+(a,b,true,1e-6)` pasa). La de `fase1/` es la anterior, sin recombinación. Cada bug de la Fase 2 debe compararse contra la referencia anterior, y las diferencias deben explicarse por el arreglo.

## Fase 2 de la limpieza (después del test de regresión): bugs que cambian resultados

No confundir con `phase2.md` (ampliar los minitrees para los análisis de kaones), que es otro plan, paralelo a esta limpieza.

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

- Todas las `GetdEdxLikelihood*` filtran los hits con 1 < dE/dx < 100 (antes solo `GetdEdxLikelihoodFreeRange`). `ComputeAveragedEdxOverResRange` devuelve -999 si ningún hit pasa los cortes (antes NaN). `IsStoppingInFV` se llama ahora `HasEndInFV`. `ComputeTruncatedMean` no ordena a propósito (media de la parte central de la traza). `GetBeamTrueParticle(spill)` y `GetTrueBeamParticle(event)` se dejan separadas por decisión del usuario.
- `pdAnalysisUtils` limpiado sin cambio de resultados: fuera las funciones sin uso (`ComputeKineticEnergy`, `ComputeDistanceToVertex`, `isBeamLike`, `ComputeTrackLengthFromHitPosition`) y el código comentado; `Chi2PID_UpToRR` es ahora `Chi2PID(part, pdg, maxRR)`; las cuatro `GetdEdxLikelihood*` comparten un helper; `GetTrueParticle(event, ID)` usa `anaUtils::GetTrueParticleByID`; arreglados accesos fuera de rango en `ComputeResidualRange` y `EstimateHitsDirection` con menos de 2 hits.
- Mensajes de error por `std::cerr` (los últimos estaban en `SCEVariationBase` y en las utilidades madre-hija de `pdAnalysisUtils`) y textos de `Calorimetry` corregidos. Los `std::cout` que quedan son `Print()` y mensajes de inicialización o finalización. `coherentfit/` (no se compila) no se ha tocado.

## Auditoría del 16-nov-2026 (hecha)

Compila desde cero sin errores ni warnings. Ejecutables probados con 2.000 eventos de MC y datos. `test/sce` estaba roto desde `33b9a51` y se arregló en `bba9295`.
Toys con σ = 0 (copia de `pdSystematics/data` con las σ a cero vía `$HIGHLANDPD_DATADIR`): en todo el MC la selección de los 161.130 toys es idéntica a la nominal y las variables difieren ≤ 1,5e-5.

Hecho:
- Bugs que cambian resultados (un commit cada uno, referencia `fixA/`): χ² recalculados en `RecombinationDataCorrection` y `dEdxMCCorrection` (helper `pdAnaUtils::RecomputeChi2PID`/`CopyChi2PID`, también en las variaciones); cortes de calidad del haz de 6 y 7 GeV; categoría de señal del mejor candidato (`SetIsSignal`); `secondaryKaonXSAnalysis` (hoy `secondaryKaonProductionXSAnalysis`) lee su `xsAnalysis.ApplyRecombinationCorrection`; `pdExampleAnalysis` rellena los hits de las hijas; `UndoSystematic` de `LifetimeVariation` y `SCECalorimetricVariation`.
- Sin cambio de resultados (test idéntico a `lklwin/`): bugs latentes (desbordamiento de `size()-1` en los cortes de kaones, `dau_mom`, `return` silencioso de `RecombinationVariation`, etiqueta del paso XS); código muerto (clases sin uso, cortes desactivados, `SpaceCharge` forward, funciones de nietas de `standardPDTree`, `PDCounters`); ~170 líneas de código comentado; parámetros que no se leían; ficheros de datos sin uso (copia en el backup); warnings; macros y scripts sin rutas personales (`$PDUTILSROOT`, argumentos, `$MICROTREEDIR` en `load.C`).
- Decidido con el usuario: los filtros de los minitrees de candidatos no importan (el análisis usa los minitrees enteros); el signo de la corrección de posición del SCE se deja.

Aplazado por decisión del usuario (16-nov-2026), no tocar hasta que lo pida:
- Calibraciones: el efecto de `SCECalorimetricVariation` no llega a `dQdx_elife` (lo que rellenan las calibraciones X/YZ) porque la llamada a la lifetime está comentada; y el dE/dx que recalcula `LifetimeVariation` parte de `hit.dQdx`, que no recoge la variación. Validar cuando haya un minitree de calibración.
- Resuelto (23-nov): `secondaryKaonProductionXSAnalysis` define la señal de la medida A (`IsProductionSignal`: la primera K+ de la cadena, subiendo las K+ -> K+, es hija directa del haz que dispara; igual que `kaon_origin_minitrees`), decisión del usuario. Es solo producción: el canal de desintegración entra en la eficiencia. Con cadenas cortadas (minitrees antiguos) la K no es señal.

Encontrado al renombrar (20-nov-2026), sin arreglar (ninguno lo causa el renombre; `ana` idéntico antes y después):
- Resuelto (23-nov, hotfix `7cda084`, `8732278`): `secondaryProton` y `secondaryKaonProductionXS` daban un `accum_level` del árbol `truth` distinto en cada ejecución. `FinalizeConfiguration` (copiado de `secondaryKaon`, que tiene una rama por candidato) llamaba a `GetAccumCutLevel(i)` con `i` = índice del candidato, con una sola rama. Ahora usa la rama 0; `ana` idéntico. Valgrind no lo detectaba.
- Resuelto (23-nov, `368313f`): `pdExample.exe` no arrancaba porque volvía a añadir `seltrk_csdarange_prot`. Ahora pasa la partícula del haz al relleno estándar. Ojo: `pdExample.exe` no sobrescribe un `-o` que ya existe.
- `beamWeights` y `brokenTracks` solo leen el formato antiguo de pduneana (`pduneana/beamana`): no corren sobre los minitrees de highland.
- Calibraciones (aparcadas por el usuario): con el prefijo de parámetros arreglado ya los leen, pero `dEdxCalibration` y `dQdxYZCalibration` acaban con corrupción de memoria en el bucle de eventos y `dQdxXCalibration` en `Initialize` (fichero personal `yz_corrections_mc_lifetime_mypc.root`).

Rama `cleanup/highlandPD` subida y PR a `develop` abierta (16-nov-2026), pendiente de revisión y merge. Revisado `phase2.md` (16-nov): la rama no quita nada del event model que ese plan necesite (solo `PDCounters`, sin uso, y comentarios). Conviene hacer el merge de esta rama antes de crear la de `phase2.md` desde `develop`.

