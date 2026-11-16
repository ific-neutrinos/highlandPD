#!/bin/bash
# Test de regresión completo: corre todo el MC y 10.000 eventos de datos en
# paralelo y compara cada salida con la referencia.
#
# Uso (con el entorno de highland y highlandPD ya cargado):
#   test/regression/check.sh <dir_referencia> <dir_salida> <minitree_mc> <minitree_data> [tolerancia]
#
# <dir_referencia> debe contener mc.root y data.root. La tolerancia relativa
# es opcional (0 = igualdad exacta). Devuelve 0 si todo coincide.

if [[ -z "${HIGHLANDPDPATH}" ]]; then
    echo "Carga primero el entorno: source ../highland/scripts/setup.sh && source scripts/setup.sh"
    exit 1
fi
if [[ $# -lt 4 ]]; then
    echo "Uso: $0 <dir_referencia> <dir_salida> <minitree_mc> <minitree_data> [tolerancia]"
    exit 1
fi

REF="$(cd "$1" && pwd)"
mkdir -p "$2"
OUT="$(cd "$2" && pwd)"
MC="$3"
DATA="$4"
TOL="${5:-0}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

"${HERE}/run_regression.sh" "${OUT}" "${MC}" - 0 &
"${HERE}/run_regression.sh" "${OUT}" - "${DATA}" 10000 &
wait

status=0
for sample in mc data; do
    echo "== ${sample}"
    # ACLiC compila compare.C en el directorio de salida
    res=$(cd "${OUT}" && root -l -b -q -e 'gSystem->SetBuildDir("aclic",true)' \
        "${HERE}/compare.C+(\"${REF}/${sample}.root\",\"${OUT}/${sample}.root\",true,${TOL})" 2>&1 \
        | grep -E "^\s+\[|RESULTADO")
    echo "${res}"
    echo "${res}" | grep -q "RESULTADO: idénticos" || status=1
done
exit ${status}
