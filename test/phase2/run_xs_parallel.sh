#!/bin/bash
# Corre secondaryKaonProductionXS.exe (u otro análisis: EXE=secondaryKaonPairXS.exe) sobre una lista de minitrees partida en N trozos en paralelo y junta las salidas.
#
#   test/phase2/run_xs_parallel.sh <lista.txt> <salida.root> [N trozos, por defecto 10]
#
# Con el entorno de highland y highlandPD cargado. Lanzarlo dentro de tmux: tarda lo que tarde un trozo (~0,3 s por
# evento: v1_00 entero en 10 trozos, ~40 min). Deja los trozos y sus logs en <salida>_trozos/ y los junta con hadd en
# <salida.root> (los árboles ana y truth se concatenan; config y header quedan repetidos, uno por trozo).
set -u
LIST=${1:?uso: run_xs_parallel.sh <lista.txt> <salida.root> [N]}
OUT=${2:?uso: run_xs_parallel.sh <lista.txt> <salida.root> [N]}
N=${3:-3}   # no más de 3-4 procesos de ROOT por dunegpvm (el 24-nov se colapsó la 05)
EXE=${EXE:-secondaryKaonProductionXS.exe}
[[ -n "${HIGHLANDPDPATH:-}" ]] || { echo "Carga antes el entorno de highland y highlandPD"; exit 1; }
DIR=${OUT%.root}_trozos
mkdir -p "$DIR"
grep -v '^\s*#' "$LIST" | grep -v '^\s*$' > "$DIR/lista_entera.txt"
split -n l/$N -d -a 2 "$DIR/lista_entera.txt" "$DIR/lista_"
for l in "$DIR"/lista_[0-9][0-9]; do
  k=${l##*_}
  rm -f "$DIR/xs_$k.root"
  "$HIGHLANDPDPATH/bin/$EXE" -o "$DIR/xs_$k.root" "$l" > "$DIR/xs_$k.log" 2>&1 &
done
echo "$N trozos lanzados (logs en $DIR); esperando..."
wait
bad=0
for l in "$DIR"/lista_[0-9][0-9]; do
  k=${l##*_}
  [[ -s "$DIR/xs_$k.root" ]] || { echo "trozo $k sin salida: mira $DIR/xs_$k.log"; bad=1; }
done
(( bad )) && exit 1
rm -f "$OUT"
hadd -f "$OUT" "$DIR"/xs_[0-9][0-9].root > "$DIR/hadd.log" 2>&1 && echo "Hecho: $OUT" || echo "Falló hadd: mira $DIR/hadd.log"
