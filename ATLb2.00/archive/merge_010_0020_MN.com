#!/bin/csh
#
set echo
#
# --- merge layers from a single HYCOM mean archive file.
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
    touch      ${D}/archMN.${y}_12L.b
    /bin/rm -f ${D}/archMN.${y}_12L.[ab]
    ../../ALL/archive/src/mrgl_archv <<E-o-D
${D}/archMN.${y}.a
${D}/archMN.${y}_12L.a
000	'iexpt ' = experiment number x10 (000=from archive file)
  0	'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3-actual)
 57	'idm   ' = longitudinal array size
 52	'jdm   ' = latitudinal  array size
 22	'kdmold' = original number of layers
 12	'kdmnew' = target   number of layers
 25.0	'thbase' = reference density (sigma units)
  1	'laybot' = last layer in next combination of layers
  3	'laybot' = last layer in next combination of layers
 21.0	'sigma ' = layer  2  density (sigma units)
  5	'laybot' = last layer in next combination of layers
 22.0	'sigma ' = layer  3  density (sigma units)
  7	'laybot' = last layer in next combination of layers
 23.6	'sigma ' = layer  4  density (sigma units)
  9	'laybot' = last layer in next combination of layers
 25.0	'sigma ' = layer  5  density (sigma units)
 11	'laybot' = last layer in next combination of layers
 25.95	'sigma ' = layer  6  density (sigma units)
 13	'laybot' = last layer in next combination of layers
 26.66	'sigma ' = layer  7  density (sigma units)
 15	'laybot' = last layer in next combination of layers
 27.12	'sigma ' = layer  8  density (sigma units)
 17	'laybot' = last layer in next combination of layers
 27.45	'sigma ' = layer  9  density (sigma units)
 19	'laybot' = last layer in next combination of layers
 27.70	'sigma ' = layer 10  density (sigma units)
 20	'laybot' = last layer in next combination of layers
 22	'laybot' = last layer in next combination of layers
 27.89	'sigma ' = layer 12  density (sigma units)
E-o-D
end
