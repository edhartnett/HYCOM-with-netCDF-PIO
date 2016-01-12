#!/bin/csh
#
set echo
#
# --- interpolate to 3-d z-levels from a single HYCOM archive file.
# --- z-levels, via linear interpolation, at Levitus depths.
#
# --- output is netCDF.
# --- this is an example, customize it for your datafile needs.
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
# --- optional title and institution.
#
setenv CDF_TITLE	"HYCOM ATLb2.00"
#setenv CDF_INST 	"RSMAS"
setenv CDF_INST 	"Naval Research Laboratory"
#
# --- D,y,d select the archive files.
#
setenv D ../expt_01.0/data
#
foreach y ( 0021 )
  foreach d ( 016 )
    setenv CDF051  ${D}/archv.${y}_${d}_00_bot.nc
    setenv CDF021  ${D}/archv.${y}_${d}_00_mlt.nc
    setenv CDF022  ${D}/archv.${y}_${d}_00_3di.nc
    setenv CDF033  ${D}/archv.${y}_${d}_00_3zt.nc
    setenv CDF034  ${D}/archv.${y}_${d}_00_3zs.nc
    setenv CDF035  ${D}/archv.${y}_${d}_00_3zr.nc
    setenv CDF037  ${D}/archv.${y}_${d}_00_3zu.nc
    setenv CDF038  ${D}/archv.${y}_${d}_00_3zv.nc
    setenv CDF039  ${D}/archv.${y}_${d}_00_3zc.nc
    setenv CDF040  ${D}/archv.${y}_${d}_00_3zw.nc
    /bin/rm $CDF051  $CDF021  $CDF022
    /bin/rm $CDF033  $CDF034  $CDF035  $CDF037  $CDF038  $CDF039
    /bin/rm $CDF040
    ../../ALL/archive/src/archv2ncdf3z <<E-o-D
${D}/archv.${y}_${d}_00.a
netCDF
 000	'iexpt ' = experiment number x10 (000=from archive file)
   0	'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3-actual)
  57	'idm   ' = longitudinal array size
  52	'jdm   ' = latitudinal  array size
  22	'kdm   ' = number of layers
  25.0	'thbase' = reference density (sigma units)
   0	'smooth' = smooth the layered fields (0=F,1=T)
   1	'iorign' = i-origin of plotted subregion
   1	'jorign' = j-origin of plotted subregion
   0	'idmp  ' = i-extent of plotted subregion (<=idm; 0 implies idm)
   0	'jdmp  ' = j-extent of plotted subregion (<=jdm; 0 implies jdm)
   1	'itype ' = interpolation type (0=sample,1=linear)
  33	'kz    ' = number of depths to sample
   0.0	'z     ' = sample depth  1
  10.0	'z     ' = sample depth  2
  20.0	'z     ' = sample depth  3
  30.0	'z     ' = sample depth  4
  50.0	'z     ' = sample depth  5
  75.0	'z     ' = sample depth  6
 100.0	'z     ' = sample depth  7
 125.0	'z     ' = sample depth  8
 150.0	'z     ' = sample depth  9
 200.0	'z     ' = sample depth 10
 250.0	'z     ' = sample depth 11
 300.0	'z     ' = sample depth 12
 400.0	'z     ' = sample depth 13
 500.0	'z     ' = sample depth 14
 600.0	'z     ' = sample depth 15
 700.0	'z     ' = sample depth 16
 800.0	'z     ' = sample depth 17
 900.0	'z     ' = sample depth 18
1000.0	'z     ' = sample depth 19
1100.0	'z     ' = sample depth 20
1200.0	'z     ' = sample depth 21
1300.0	'z     ' = sample depth 22
1400.0	'z     ' = sample depth 23
1500.0	'z     ' = sample depth 24
1750.0	'z     ' = sample depth 25
2000.0	'z     ' = sample depth 26
2500.0	'z     ' = sample depth 27
3000.0	'z     ' = sample depth 28
3500.0	'z     ' = sample depth 29
4000.0	'z     ' = sample depth 30
4500.0	'z     ' = sample depth 31
5000.0	'z     ' = sample depth 32
5500.0	'z     ' = sample depth 33
  51	'botio ' = bathymetry  I/O unit (0 no I/O)
   0	'mltio ' = mix.l.thk.  I/O unit (0 no I/O)
   0.0	'tempml' = temperature jump across mixed-layer (degC,  0 no I/O)
   0.0	'densml' =     density jump across mixed-layer (kg/m3, 0 no I/O)
 -22	'infio ' = intf. depth I/O unit (0 no I/O, <0 label with layer #)
  40	'wvlio ' = w-velocity  I/O unit (0 no I/O)
  37	'uvlio ' = u-velocity  I/O unit (0 no I/O)
  38	'vvlio ' = v-velocity  I/O unit (0 no I/O)
  39	'splio ' = speed       I/O unit (0 no I/O)
  33	'temio ' = temperature I/O unit (0 no I/O)
  34	'salio ' = salinity    I/O unit (0 no I/O)
  35	'tthio ' = density     I/O unit (0 no I/O)
E-o-D
  end
end
