#
set echo
cd $cwd
#
# --- Usage:  ./Make.com >& Make.log
#
# --- make hycom with TYPE from this directory's name (src_*_$TYPE).
# --- set ARCH to the correct value for this machine.
# --- assumes dimensions.h is correct for $TYPE.
#
setenv ARCH sp3GPFS
#
setenv TYPE `echo $cwd | awk -F"_" '{print $NF}'`
#
if (! -e ../config/${ARCH}_${TYPE}) then
  echo "ARCH = " $ARCH "  TYPE = " $TYPE "  is not supported"
  exit 1
else
  make hycom ARCH=$ARCH TYPE=$TYPE
endif
