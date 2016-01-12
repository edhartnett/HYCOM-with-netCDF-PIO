#!/bin/csh
#
set echo
#
# --- form meridional stream function from a HYCOM mean archive file.
#
# --- output can be formatted or unformatted (BINARY)
# --- or use archv2ncdfsf for netCDF output.
#
# --- output is HYCOM .a files, converted to "raw" .A files.
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
# --- D,y,d select the archive files.
#
setenv D ../expt_01.0/data
#
foreach y ( 0020 )
    setenv FOR021  ${D}/archMN.${y}_sfi.d
    setenv FOR022  ${D}/archMN.${y}_sfl.d
    setenv FOR023  ${D}/archMN.${y}_sfz.d
    /bin/rm $FOR021  $FOR022  $FOR023
    ../../ALL/archive/src/archv2datasf <<E-o-D
${D}/archMN.${y}.a
(1p5g14.5)
 000	'iexpt ' = experiment number x10 (000=from archive file)
   0	'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3-actual)
  57	'idm   ' = longitudinal array size
  52	'jdm   ' = latitudinal  array size
  22	'kdm   ' = number of layers
   1	'iorign' = i-origin of plotted subregion
   1	'jorign' = j-origin of plotted subregion
   0	'idmp  ' = i-extent of plotted subregion (<=idm; 0 implies idm)
   0	'jdmp  ' = j-extent of plotted subregion (<=jdm; 0 implies jdm)
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
 -21	'infio ' = intf. depth I/O unit (0 no I/O, <0 label with layer #)
 -22	'sfnlio' = layer strfn I/O unit (0 no I/O, <0 label with layer #))
  23	'sfn3io' =     z strfn I/O unit (0 no I/O)
E-o-D
end
