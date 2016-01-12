#!/bin/csh
#
set echo
#
# --- merge all layers into one layer from a single HYCOM mean archive file.
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
# --- D,y,d select the archive files.
#
setenv D ../expt_01.0/data
#
foreach y ( 0020 )
    touch      ${D}/archMN.${y}_01L.b
    /bin/rm -f ${D}/archMN.${y}_01L.[ab]
    ../../ALL/archive/src/mrgl_archv <<E-o-D
${D}/archMN.${y}.a
${D}/archMN.${y}_01L.a
000	'iexpt ' = experiment number x10 (000=from archive file)
  0	'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3-actual)
 57	'idm   ' = longitudinal array size
 52	'jdm   ' = latitudinal  array size
 22	'kdmold' = original number of layers
  1	'kdmnew' = target   number of layers
 25.0	'thbase' = reference density (sigma units)
 22	'laybot' = last layer in next combination of layers
 27.89	'sigma ' = layer  1  density (sigma units)
E-o-D
end
