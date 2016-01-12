#!/bin/csh
#
set echo
#
# --- sample directory softlinks for topography and source.
#
setenv R ATLb2.00
setenv N 01
#
touch   regional.grid.a
/bin/rm *.[ab]
#
/bin/ln -s ../topo/regional.grid.a   regional.grid.a
/bin/ln -s ../topo/regional.grid.b   regional.grid.b
#
/bin/ln -s ../topo/depth_${R}_${N}.a regional.depth.a
/bin/ln -s ../topo/depth_${R}_${N}.b regional.depth.b
#
if (! -e src) then
  /bin/ln -s ../../ALL/sample/src_2.1.00 src
endif
