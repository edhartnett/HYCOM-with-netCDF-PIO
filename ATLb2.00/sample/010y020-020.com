#
set echo
#
# --- sample transport sections from a sequence of archive files
#
setenv D  ../expt_01.0/data
setenv S  ../expt_01.0/data
setenv E  010
setenv YF 0020
setenv YL 0020
#
touch regional.depth.a regional.depth.b
if (-z regional.depth.a) then
  /bin/rm regional.depth.a
  /bin/ln -s ../topo/depth_ATLb2.00_01.a regional.depth.a
endif
if (-z regional.depth.b) then
  /bin/rm regional.depth.b
  /bin/ln -s ../topo/depth_ATLb2.00_01.b regional.depth.b
endif
#
touch regional.grid.a regional.grid.b
if (-z regional.grid.a) then
  /bin/rm regional.grid.a
  /bin/ln -s ../topo/regional.grid.a .
endif
if (-z regional.grid.b) then
  /bin/rm regional.grid.b
  /bin/ln -s ../topo/regional.grid.b .
endif
#
touch   ${D}/trans.${YF}-${YL}.d
/bin/rm ${D}/trans.${YF}-${YL}.d
#
./src/transport <<E-o-D
${D}/trans.${YF}-${YL}.d
${D}/archv.${YF}_016_00.b
${D}/archv.${YL}_346_00.b
 000	  'iexpt ' = experiment number x10
   0	  'yrflag' = days in year flag (0=360,  1=366,  2=366J1, 3=actual)
  30.0	  'tranfq' = number of days between archive input
  57	  'idm   ' = longitudinal array size
  52	  'jdm   ' = latitudinal  array size
  22	  'kdm   ' = number of layers
   3	  'ntrans' = number of transport sections
Cross-Equator
  25	  'if    ' = first index of the transport section, longitude
  55	  'il    ' = last  index of the transport section, longitude
  11	  'jf    ' = first index of the transport section, latitude
  11	  'jl    ' = last  index of the transport section, latitude
Along-Equator, 39W
  30	  'if    ' = first index of the transport section, longitude
  30	  'il    ' = last  index of the transport section, longitude
  09	  'jf    ' = first index of the transport section, latitude
  13	  'jl    ' = last  index of the transport section, latitude
Along-Equator, 19W
  40	  'if    ' = first index of the transport section, longitude
  40	  'il    ' = last  index of the transport section, longitude
  09	  'jf    ' = first index of the transport section, latitude
  13	  'jl    ' = last  index of the transport section, latitude
E-o-D
