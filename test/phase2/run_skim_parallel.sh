#!/bin/bash
# Corre truth_skim_minitrees.C sobre una lista de minitrees partida en N trozos en paralelo y junta las salidas.
#
#   test/phase2/run_skim_parallel.sh <lista.txt> <salida.root> [N trozos, por defecto 10]
#
# Con el entorno de highland y highlandPD cargado. Lanzarlo dentro de tmux. Compila la macro una vez antes (con
# una entrada) para que los trozos no compilen a la vez. Deja los trozos y sus logs en <salida>_trozos/ y los junta
# con hadd en <salida.root> (el TNamed "version" queda repetido, uno por trozo). ~40 kB por evento: v1_00 entero,
# ~3 GB.
set -u
LIST=${1:?uso: run_skim_parallel.sh <lista.txt> <salida.root> [N]}
OUT=${2:?uso: run_skim_parallel.sh <lista.txt> <salida.root> [N]}
N=${3:-3}   # no más de 3-4 procesos de ROOT por dunegpvm
[[ -n "${HIGHLANDPDPATH:-}" ]] || { echo "Carga antes el entorno de highland y highlandPD"; exit 1; }
MACRO=$(cd "$(dirname "$0")" && pwd)/truth_skim_minitrees.C
DIR=${OUT%.root}_trozos
mkdir -p "$DIR"
grep -v '^\s*#' "$LIST" | grep -v '^\s*$' > "$DIR/lista_entera.txt"
split -n l/$N -d -a 2 --additional-suffix=.txt "$DIR/lista_entera.txt" "$DIR/lista_"   # .txt: la macro lo lee como lista
echo "Compilando la macro..."
root -l -b -q "$MACRO(\"$DIR/lista_00.txt\", \"$DIR/compila.root\", 1)" > "$DIR/compila.log" 2>&1
grep -q "Árbol skim en" "$DIR/compila.log" || { echo "Falló la compilación o la prueba con una entrada: mira $DIR/compila.log"; exit 1; }
rm -f "$DIR/compila.root"
for l in "$DIR"/lista_[0-9][0-9].txt; do
  k=${l##*_}; k=${k%.txt}
  rm -f "$DIR/skim_$k.root"
  root -l -b -q "$MACRO(\"$l\", \"$DIR/skim_$k.root\")" > "$DIR/skim_$k.log" 2>&1 &
done
echo "$N trozos lanzados (logs en $DIR); esperando..."
wait
bad=0
for l in "$DIR"/lista_[0-9][0-9].txt; do
  k=${l##*_}; k=${k%.txt}
  grep -q "Árbol skim en" "$DIR/skim_$k.log" && [[ -s "$DIR/skim_$k.root" ]] || { echo "trozo $k sin salida: mira $DIR/skim_$k.log"; bad=1; }
done
(( bad )) && exit 1
rm -f "$OUT"
hadd -f "$OUT" "$DIR"/skim_[0-9][0-9].root > "$DIR/hadd.log" 2>&1 && echo "Hecho: $OUT" || echo "Falló hadd: mira $DIR/hadd.log"
