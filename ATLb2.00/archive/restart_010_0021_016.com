#!/bin/csh
#
set echo
#
# --- Form a HYCOM restart file from a HYCOM archive file.
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
# --- data directory
#
setenv D ../expt_01.0/data
#
# ---  input archive file
# ---  input restart file
# --- output restart file
#
../../ALL/archive/src/archv2restart <<E-o-D
${D}/archv.0021_016_00.a
${D}/restart_021.a
${D}/restart_021_016.a
000	'iexpt ' = experiment number x10 (000=from archive file)
  0	'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3-actual)
 57	'idm   ' = longitudinal array size
 52	'jdm   ' = latitudinal  array size
 22	'kdm   ' = number of layers
 25.0	'thbase' = reference density (sigma units)
5760.0	'baclin' = baroclinic time step (seconds), int. divisor of 86400
E-o-D
