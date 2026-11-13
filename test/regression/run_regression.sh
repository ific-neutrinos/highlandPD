#!/bin/bash
# Test de regresión de secondaryKaon.exe sobre un minitree de MC y otro de datos.
#
# Uso (con el entorno de highland y highlandPD ya cargado):
#   test/regression/run_regression.sh <dir_salida> <minitree_mc> <minitree_data> [n_eventos]
#
# Escribe <dir_salida>/{mc,data}.root y sus logs. Para comparar dos ejecuciones:
#   root -l -b -q 'test/regression/compare.C+("ref/mc.root","new/mc.root")'
#
# El ejecutable se toma de $HIGHLANDPDPATH/bin, así que el resultado depende
# de qué setup.sh se haya cargado (el repo o un worktree de otra versión).

if [[ -z "${HIGHLANDPDPATH}" ]]; then
    echo "Carga primero el entorno: source ../highland/scripts/setup.sh && source scripts/setup.sh"
    exit 1
fi
if [[ $# -lt 3 ]]; then
    echo "Uso: $0 <dir_salida> <minitree_mc> <minitree_data> [n_eventos]"
    exit 1
fi

OUT="$(mkdir -p "$1" && cd "$1" && pwd)"
MC="$2"
DATA="$3"
NEV="${4:-0}"   # 0 = todos los eventos

# El override vive junto a este script
PARAMS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/regression.parameters.dat"
EXE="${HIGHLANDPDPATH}/bin/secondaryKaon.exe"

NOPT=""
[[ "${NEV}" -gt 0 ]] && NOPT="-n ${NEV}"

status=0
for sample in mc data; do
    if [[ ${sample} == mc ]]; then input="${MC}"; else input="${DATA}"; fi
    echo "== ${sample}: ${input}"
    start=$(date +%s)
    # Se ejecuta desde el directorio de salida para no dejar ficheros sueltos
    (cd "${OUT}" && "${EXE}" -p "${PARAMS}" ${NOPT} -o "${OUT}/${sample}.root" "${input}") \
        > "${OUT}/${sample}.log" 2>&1
    rc=$?
    echo "   exit ${rc}, $(( $(date +%s) - start )) s, log en ${OUT}/${sample}.log"
    [[ ${rc} -ne 0 ]] && status=1
done
exit ${status}
