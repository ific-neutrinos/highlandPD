# Limpia y recompila desde cero highland y highlandPD.
# Uso (desde cualquier sitio): source scripts/DOALL.sh

# directorio que contiene highland/ y highlandPD/
TOPDIR="$(cd -P -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd -P)"
echo $TOPDIR
cd $TOPDIR/highland

source scripts/cleanup.sh
source scripts/INSTALL.sh

cd $TOPDIR/highlandPD

source scripts/cleanup.sh
source scripts/INSTALL.sh
