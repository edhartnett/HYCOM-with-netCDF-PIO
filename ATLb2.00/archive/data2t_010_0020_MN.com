#!/bin/csh
#
set echo
#
# --- extract 2-d fields along temperature surfaces
# --- from a single HYCOM mean archive file.
# --- configured for 22 layers.
#
# --- output can be formatted, unformatted (BINARY), .[ab] (HYCOM).
# --- or use archv2ncdf2d for netCDF output.
#
# --- output is HYCOM .a files, optionally converted to "raw" .A files.
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
    setenv FOR051A ${D}/archMN.${y}_bot.a
    setenv FOR051  ${D}/archMN.${y}_bot.b
    setenv FOR031A ${D}/archMN.${y}_3tl.a
    setenv FOR031  ${D}/archMN.${y}_3tl.b
    setenv FOR032A ${D}/archMN.${y}_3td.a
    setenv FOR032  ${D}/archMN.${y}_3td.b
    setenv FOR033A ${D}/archMN.${y}_3tt.a
    setenv FOR033  ${D}/archMN.${y}_3tt.b
    setenv FOR034A ${D}/archMN.${y}_3ts.a
    setenv FOR034  ${D}/archMN.${y}_3ts.b
    setenv FOR035A ${D}/archMN.${y}_3tr.a
    setenv FOR035  ${D}/archMN.${y}_3tr.b
    /bin/rm $FOR051
    /bin/rm $FOR051A
    /bin/rm $FOR031  $FOR032  $FOR033  $FOR034  $FOR035
    /bin/rm $FOR031A $FOR032A $FOR033A $FOR034A $FOR035A
    ../../ALL/archive/src/archv2data2t <<E-o-D
${D}/archMN.${y}.a
HYCOM
000	'iexpt ' = experiment number x10 (000=from archive file)
  0	'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3-actual)
 57	'idm   ' = longitudinal array size
 52	'jdm   ' = latitudinal  array size
 40	'itest ' = longitudinal test point (optional, default 0)
 40	'jtest ' = latitudinal  test point (optional, default 0)
 22	'kdm   ' = number of layers
 25.0	'thbase' = reference density (sigma units)
  0	'smooth' = smooth fields before plotting (0=F,1=T)
  1	'iorign' = i-origin of plotted subregion
  1	'jorign' = j-origin of plotted subregion
  0	'idmp  ' = i-extent of plotted subregion (<=idm; 0 implies idm)
  0	'jdmp  ' = j-extent of plotted subregion (<=jdm; 0 implies jdm)
  4	'ktemp ' = number of temperature surfaces to sample
  4.0	'tsur  ' = sample temperaure
  3.6	'tsur  ' = sample temperaure
  2.8	'tsur  ' = sample temperaure
  2.4	'tsur  ' = sample temperaure
  0	'botio ' = bathymetry       I/O unit (0 no I/O)
 31	'layio ' = layer k   i.dep. I/O unit (0 no I/O)
 32	'depio ' = layer k   thick. I/O unit (0 no I/O)
 33	'temio ' = layer k   temp   I/O unit (0 no I/O)
 34	'salio ' = layer k   saln.  I/O unit (0 no I/O)
 35	'tthio ' = layer k   dens,  I/O unit (0 no I/O)
E-o-D
####
#### --- convert HYCOM .a files to RAW files (no padding, spval=1.e10).
#### --- comment this out if you don't need RAW files.
####
###    foreach t ( bot 3tl 3td 3tt 3ts 3tr )
###      if (-e ${D}/archMN.${y}_${t}.a) then
###        /bin/rm -f ${D}/archMN.${y}_${t}.A
###        ../../ALL/bin/hycom2raw ${D}/archMN.${y}_${t}.a 57 52 1.e10 ${D}/archMN.${y}_${t}.A
####       /bin/rm -f ${D}/archMN.${y}_${t}.a
###      endif
###    end
end
