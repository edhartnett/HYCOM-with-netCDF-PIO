#!/bin/csh
#
set echo
#
# --- remap single archive file, to interface depths from archv2data2t.
#
setenv R ATLb2.00
setenv T 01
setenv E  010
setenv X 01.0
setenv Y 0020
setenv A ""
#
cd /net/ajax/data/$user/hycom/${R}/expt_${X}/data/meanstd
#cd /scr/$user/hycom/${R}/expt_${X}/data/meanstd
#
touch regional.depth.a regional.depth.b
if (-z regional.depth.a) then
  /bin/rm regional.depth.a
  /bin/ln -s ~/hycom/${R}/topo/depth_${R}_${T}.a regional.depth.a
endif
if (-z regional.depth.b) then
  /bin/rm regional.depth.b
  /bin/ln -s ~/hycom/${R}/topo/depth_${R}_${T}.b regional.depth.b
endif
#
touch regional.grid.a regional.grid.b
if (-z regional.grid.a) then
  /bin/rm regional.grid.a
  /bin/ln -s ~/hycom/${R}/topo/regional.grid.a .
endif
if (-z regional.grid.b) then
  /bin/rm regional.grid.b
  /bin/ln -s ~/hycom/${R}/topo/regional.grid.b .
endif
#
/bin/rm -f archMN.${Y}${A}_NEW.a archMN.${Y}${A}_NEW.b
~wallcraf/hycom/ALL/archive/src/remapi_archv <<E-o-D
archMN.${Y}${A}.a
archMN.${Y}${A}_NEW.a
${E}_archMN.${Y}${A}_3td.a
${E}	  'iexpt ' = experiment number x10 (000=from archive file)
   1	  'yrflag' = days in year flag (0=360,  1=366,  2=366J1, 3=actual)
  57	  'idm   ' = longitudinal array size
  52	  'jdm   ' = latitudinal  array size
  22	  'kdmold' = original number of layers
  14	  'kdmnew' = target   number of layers
  25.0	  'thbase' = reference density (sigma units)
   1.00   'sigma ' = layer  1 isopycnal target density (sigma units)
   2.00   'sigma ' = layer  2 isopycnal target density (sigma units)
   3.00   'sigma ' = layer  3 isopycnal target density (sigma units)
   4.00   'sigma ' = layer  4 isopycnal target density (sigma units)
   5.00   'sigma ' = layer  5 isopycnal target density (sigma units)
   6.00   'sigma ' = layer  6 isopycnal target density (sigma units)
   7.00   'sigma ' = layer  7 isopycnal target density (sigma units)
   8.00   'sigma ' = layer  8 isopycnal target density (sigma units)
   9.00   'sigma ' = layer  9 isopycnal target density (sigma units)
  10.00   'sigma ' = layer 10 isopycnal target density (sigma units)
  11.00   'sigma ' = layer 11 isopycnal target density (sigma units)
  12.00   'sigma ' = layer 12 isopycnal target density (sigma units)
  13.00   'sigma ' = layer 13 isopycnal target density (sigma units)
  14.00   'sigma ' = layer 14 isopycnal target density (sigma units)
E-o-D
