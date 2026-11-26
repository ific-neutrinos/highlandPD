# Recetas: lo rutinario, paso a paso

Comandos para copiar y pegar. Los detalles y el porqué, en `grid_produccion.md` y en la cabecera de cada script
(`highlandAnalysis/justin/*.sh`, `highlandPD/test/phase2/*.C`). Escrito el 24-nov-2026.

Convenciones: `<wf>` número del workflow, `<fecha>` AAAAMMDD con un sufijo que diga qué es (p. ej. `20261124v100s`),
`USER` = miagarc. Todo lo largo, **dentro de tmux** (si se cierra Apptainer, `nohup` no basta).

## 0. Entrar en FNAL

```bash
ssh dunegpvmXX.fnal.gov            # una dunegpvm para cada token (ver 0b)
tmux new -s <nombre>               # o: tmux attach -t <nombre>
# Apptainer SL7 como siempre, y dentro:
cd /exp/dune/app/users/miagarc/larsoft/dunesw_v09_56_00d00
source setup_larsoft.sh            # token de justIN: samweb, justin, compilar
```

0b. Tokens: `source setup_larsoft.sh` (por defecto, `PD_TOKEN=justin`) para justIN y samweb;
`PD_TOKEN=interactive source setup_larsoft.sh` solo para escribir en persistent (paso 5). samweb siempre con
`--no-cert-auth`.

## 1. ¿Ha acabado el workflow?

```bash
justin show-files --workflow-id <wf> \
  | grep -oE '\b(finding|unallocated|allocated|outputting|processed|notfound|failed)\b' | sort | uniq -c
justin show-workflows --workflow-id <wf>
```

Acabado: todo `processed` (algún `failed`/`notfound` se mira aparte) y el workflow en `finished`. Los trabajos
(`show-jobs`) no sirven de criterio: justIN lanza de más.

## 2. Recoger la salida

```bash
cd <highlandAnalysis>/justin
./collect.sh <wf> mc <fecha>
```

Mueve (sin copiar) de `/pnfs/dune/scratch/users/miagarc/fnal/<wf>/` a `/pnfs/dune/scratch/users/miagarc/pdsp/mc/<fecha>/`.
Usar siempre una carpeta nueva: no sobrescribe, y avisa si un fichero ya existe.

## 3. ¿Ha ido bien? (resumen)

```bash
./summarize.sh /pnfs/dune/scratch/users/miagarc/pdsp/mc/<fecha> lists/filelist_mc_6GeV_v1_00.txt
```

Mirar: ficheros bien / mal / sin stats, s/evento, `VmHWM` (memoria máxima; si roza `RSS_MIB`, subirlo la próxima
vez) y MB/evento. Si falla alguno:

```bash
RETRY_LIST=lists/retry_<wf>.txt ./summarize.sh /pnfs/dune/scratch/users/miagarc/pdsp/mc/<fecha> lists/filelist_mc_6GeV_v1_00.txt
# rehacer el tarball con esa lista (paso 8) y relanzar:
LIST=retry_<wf>.txt FILES_PER_JOB=1 SITES=US_FNAL-FermiGrid ./submit.sh mc <N> $INPUT_TAR_DIR_LOCAL
```

## 4. Lista de las salidas y validación

Las macros se corren sobre la lista de **salidas** (minitrees), no sobre `filelist_*.txt` (que son las entradas de
SAM). **Siempre con URLs xrootd**, no con rutas `/pnfs`: las macros leen lo que diga la lista, y con rutas `/pnfs`
leerían por el montaje NFS de dCache. Por xrootd, ~7 eventos/s por macro (una pasada entera, 1,5-2 h):

```bash
D=/pnfs/dune/scratch/users/miagarc/pdsp/mc/<fecha>; ls $D | grep '^minitree_.*\.root$' | sed "s|^|$D/|; s|^/pnfs/|root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/|" > mc_<nombre>.txt
# sin glob: con 11.626 ficheros (v1_03) "ls .../minitree_*.root" da "Argument list too long"
head -1 mc_<nombre>.txt; wc -l mc_<nombre>.txt   # root://fndcadoor.fnal.gov:1094/pnfs/fnal.gov/usr/dune/...; nº de ficheros
```

Se puede correr sobre scratch sin esperar a persistent. Cuando ya esté en persistent, rehacer la lista con
`persistent` en lugar de `scratch` (el scratch se purga). El aviso `Plugin version SecClnt ... secztn` de xrootd es
inofensivo.

Macros que se corren sobre cada muestra nueva (dentro de tmux; las de una sola pasada, con `&` y su log):

Orden (25-nov; como mucho 3-4 procesos de ROOT a la vez, cada `run_*_parallel.sh` son 3):

1. Validación del contenido: `count_module` (eventos marcados, beam slice, primario = haz; igual que la muestra
   anterior si la reco no cambió) y, con un módulo nuevo, `check_production_minitree` (hits de la traza forzada solo
   dentro de la ventana, tiempo true). Si fallan, se paran los workflows de las demás versiones.
2. Truth skim y `skim_compare.C` en local frente al skim de la producción anterior de los mismos eventos (truth
   idéntica: los resultados de truth anteriores siguen valiendo).
3. Análisis: `secondaryKaonProductionXS` (microtree de referencia) y `secondaryKaonPairXS` (comparación con la tesis).
4. Figuras de la nota: `kaon_acceptance` (y `kaon_origin`, `pimu_endcharge` si hacen falta).

```bash
root -l -b -q 'test/phase2/count_module_minitrees.C("mc_<nombre>.txt")' > count_<nombre>.log 2>&1 &
root -l -b -q 'test/phase2/check_production_minitree.C("mc_<nombre>.txt")' > check_<nombre>.log 2>&1 &
test/phase2/run_skim_parallel.sh mc_<nombre>.txt truth_skim_<nombre>.root      # truth skim, ~40 kB/evento (3 trozos)
test/phase2/run_xs_parallel.sh mc_<nombre>.txt xs_<nombre>.root                # análisis (3 trozos)
EXE=secondaryKaonPairXS.exe test/phase2/run_xs_parallel.sh mc_<nombre>.txt pair_<nombre>.root
root -l -b -q 'test/phase2/kaon_acceptance_minitrees.C("mc_<nombre>.txt","kaon_acceptance_<nombre>.root")' > kaon_acceptance_<nombre>.log 2>&1 &
root -l -b -q 'test/phase2/kaon_origin_minitrees.C("mc_<nombre>.txt")' > kaon_origin_<nombre>.log 2>&1 &
root -l -b -q 'test/phase2/pimu_endcharge.C("mc_<nombre>.txt")' > pimu_endcharge_<nombre>.log 2>&1 &
```

En local, con los microtrees: `skim_compare.C`, `xs_viewloss`, `xs_showers`, `xs_nminus1`, `technote_effpur`,
`technote_cutflow`, `technote_selection`, figuras de truth (`technote_truth_figs.C`) y `compare_pair_xs.C`.

Traer a local los `.root` (salvo los trozos) y los `.log`.

Validación (en highlandPD, con su entorno):

```bash
cd <highlandPD>
source ../highland/scripts/setup.sh && source scripts/setup.sh
root -l -b -q 'test/phase2/count_module_minitrees.C("mc_<nombre>.txt")' > count_<nombre>.log
```

Mirar que no haya incoherencias entre flags y contenido, y que los números sean como los de la muestra anterior.

## 5. Copiar a persistent (el scratch se purga)

**Organización (decidida con el usuario el 25-nov-2026, desde el tag `d01`):** un proyecto por tag del código, y
dentro, separados minitrees y microtrees, MC y datos, momento del haz y dataset:

```
/pnfs/dune/persistent/users/miagarc/pdsp/<proyecto>/
    README.txt                                  # commits del tag, workflows, definiciones de SAM
    minitrees/<mc|data>/<momento>/<dataset>/minitrees/   # minitree_*.root
    minitrees/<mc|data>/<momento>/<dataset>/stats/       # stats_*.txt
    microtrees/<análisis>/<mc|data>/<momento>/<dataset>/ # salidas de highland (xs_*.root, ...)
```

`<proyecto>` = `dunesw_v09_56_00d01_claude` (tag + sufijo); `<momento>` = `6GeV`; `<dataset>` = el bloque de SAM
abreviado (`v1_00` ... `v1_09` para `PDSPProd4a_MC_6GeV_gen_datadriven_reco1_sce_datadriven_v1_XX`; en datos, el
runset). Las carpetas antiguas (`pdsp/mc/20261119big`, `20261120v100`, `20261124v100s`) se quedan como están hasta
borrarlas. El scratch sigue por fecha (`collect.sh`).

En una dunegpvm con el token interactivo:

```bash
PD_TOKEN=interactive source setup_larsoft.sh
cd <highlandAnalysis>/justin
P=/pnfs/dune/persistent/users/miagarc/pdsp/dunesw_v09_56_00d01_claude/minitrees/mc/6GeV/v1_00
./to_persistent.sh /pnfs/dune/scratch/users/miagarc/pdsp/mc/<fecha> $P/minitrees
./to_persistent.sh /pnfs/dune/scratch/users/miagarc/pdsp/mc/<fecha> $P/stats 'stats_*.txt'
```

Se puede relanzar tal cual si se corta (salta lo ya copiado). Los fallos van a `to_persistent_failed.txt`.
Ritmo: ~110 ficheros/min (una muestra de 7.300 ficheros, ~1 h por tipo de fichero). **Nunca** `cp`/`rsync` por NFS a
`/pnfs`: solo `to_persistent.sh` (ifdh).

## 6. Correr el análisis de highland sobre la muestra

```bash
cd <highlandPD>
git pull                                                   # si hay commits nuevos (push hecho desde aquí)
source ../highland/scripts/setup.sh && source scripts/setup.sh
cmake --build build -j8 && cmake -P build/cmake_install.cmake   # recompilar tras cada pull
test/phase2/run_xs_parallel.sh mc_<nombre>.txt xs_<nombre>.root 10
```

Avance de los trozos:

```bash
for f in xs_<nombre>_trozos/xs_[0-9][0-9].log; do echo "$f: $(grep -a -o 'entry: [0-9]* of [0-9]* ([0-9]*%)' $f | tail -1)"; done
```

Al acabar escribe `Hecho: xs_<nombre>.root`. Comprobar que tiene las ramas esperadas:

```bash
root -l -b -q -e 'TFile f("xs_<nombre>.root"); ((TTree*)f.Get("ana"))->Print("candmu*")'
```

## 7. Traer ficheros al portátil

Desde el portátil:

```bash
scp miagarc@dunegpvmXX.fnal.gov:<ruta>/xs_<nombre>.root ~/Documents/Work/DUNE/analysis/PDSP/claude/files/minitrees/mc/
```

## 8. Lanzar una producción nueva (bloque de MC o datos)

Varias versiones a la vez con el mismo tarball (25-nov, `v1_01`-`v1_09`): las listas ya están dentro del tarball, así que
no hace falta rehacerlo; lanzar solo las que `check_locality.sh` da en disco (si están en cinta los trabajos esperan y
acaban por wall time). Bucle con `DRYRUN=1` primero:

```bash
source tardir.sh
for v in 01 02 03; do
  L=filelist_mc_6GeV_v1_$v.txt; N=$(( ($(wc -l < lists/$L) + 4) / 5 ))
  DRYRUN=1 LIST=$L FILES_PER_JOB=5 SITES=US_FNAL-FermiGrid ./submit.sh mc $N $INPUT_TAR_DIR_LOCAL
done
```

`workflows.txt` no guarda la versión: apuntarla a mano. Recoger cada una en su carpeta (`collect.sh <wf> mc <fecha>v1<v>`)
y a persistent en `.../6GeV/v1_<v>/{minitrees,stats}`. En dunegpvm, con el entorno de dunesw, `sed` avisa de
`libselinux.so.1: no version information available` (por `dune_oslibs` en `LD_LIBRARY_PATH`): inofensivo
(`env -u LD_LIBRARY_PATH` lo quita). Las tmux no se cierran solas: una por tarea con nombre y `tmux kill-session -t <n>`
al acabar (mirar antes con `tmux capture-pane -pt <n> -S -30` que no sigue corriendo, sobre todo los prestage).

1. Lista de entrada: `./make_filelist.sh <definición SAM> <N> lists/filelist_<algo>.txt` (en tmux: una llamada a
   samweb por fichero; `SAMWEB_OPTS=--no-cert-auth`).
2. Sacar de cinta: `samweb --no-cert-auth prestage-dataset --defname=<definición>` (en tmux; los errores de DNS del
   final no importan). Estado: `./check_locality.sh lists/filelist_<algo>.txt 10` (todo `ONLINE_AND_NEARLINE`).
3. Si cambió el código, el fcl o las listas: recompilar el módulo (`mrb z && mrbsetenv && fixcmake && checkenv && mrb i -j8`)
   y rehacer el tarball: `./pack_runtime.sh lists pdsp_runtime.tar.gz`, subirlo con `justin-cvmfs-upload` y guardar la
   ruta (`INPUT_TAR_DIR_LOCAL`, en `tardir.sh`).
4. Lanzar: `LIST=filelist_<algo>.txt FILES_PER_JOB=5 SITES=US_FNAL-FermiGrid ./submit.sh mc <N trabajos> $INPUT_TAR_DIR_LOCAL`
   (N = ficheros / 5, redondeando hacia arriba; `DRYRUN=1` delante para ver el comando sin lanzarlo). Apunta el
   workflow en `workflows.txt`.
5. Después: pasos 1 a 6.

## 9. Tag del código de una producción

Sobre los commits con los que se hizo el tarball (no los posteriores), en highland, highlandPD y highlandAnalysis:

```bash
git tag -a dunesw_v09_56_00dXX <commit> -m "<qué producción>"
git push origin dunesw_v09_56_00dXX
```
