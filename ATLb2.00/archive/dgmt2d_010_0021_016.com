#!/bin/csh
#
set echo
#
# --- extract 2-d fields from a single HYCOM archive file.
# --- configured to output only layers 9-12
#
# --- output can be formatted, unformatted (BINARY), .[ab] (HYCOM).
# --- or use archv2ncdf2d for netCDF output.
#
# --- output is formatted (longitude, latitude, value), 
# --- suitable for input as irregularly spaced data into GMT.
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
foreach y ( 0021 )
  foreach d ( 016 )
    foreach k ( 09 10 11 12 )
      setenv FOR051  ${D}/archv.${y}_${d}_00_bot.d
      setenv FOR021  ${D}/archv.${y}_${d}_00_flx.d
      setenv FOR022  ${D}/archv.${y}_${d}_00_ice.d
      setenv FOR023  ${D}/archv.${y}_${d}_00_fsd.d
      setenv FOR024  ${D}/archv.${y}_${d}_00_mld.d
      setenv FOR025  ${D}/archv.${y}_${d}_00_bld.d
      setenv FOR031  ${D}/archv.${y}_${d}_00_2di_${k}.d
      setenv FOR032  ${D}/archv.${y}_${d}_00_2dh_${k}.d
      setenv FOR033  ${D}/archv.${y}_${d}_00_2dt_${k}.d
      setenv FOR034  ${D}/archv.${y}_${d}_00_2ds_${k}.d
      setenv FOR035  ${D}/archv.${y}_${d}_00_2dr_${k}.d
      setenv FOR036  ${D}/archv.${y}_${d}_00_2dn_${k}.d
      setenv FOR037  ${D}/archv.${y}_${d}_00_2du_${k}.d
      setenv FOR038  ${D}/archv.${y}_${d}_00_2dv_${k}.d
      setenv FOR039  ${D}/archv.${y}_${d}_00_2dc_${k}.d
      setenv FOR040  ${D}/archv.${y}_${d}_00_2dw_${k}.d
      /bin/rm $FOR051
      /bin/rm $FOR021  $FOR022  $FOR023  $FOR024  $FOR025
      /bin/rm $FOR031  $FOR032  $FOR033  $FOR034  $FOR035
      /bin/rm $FOR036  $FOR037  $FOR038  $FOR039  $FOR040
      ../../ALL/archive/src/archv2data2d <<E-o-D
${D}/archv.${y}_${d}_00.a
(2f10.4,f18.8)
000	'iexpt ' = experiment number x10 (000=from archive file)
  0	'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3-actual)
 57	'idm   ' = longitudinal array size
 52	'jdm   ' = latitudinal  array size
 22	'kdm   ' = number of layers
 25.0	'thbase' = reference density (sigma units)
  0	'smooth' = smooth fields before plotting (0=F,1=T)
  0	'mthin ' = mask thin layers from plots   (0=F,1=T)
  1	'iorign' = i-origin of plotted subregion
  1	'jorign' = j-origin of plotted subregion
  0	'idmp  ' = i-extent of plotted subregion (<=idm; 0 implies idm)
  0	'jdmp  ' = j-extent of plotted subregion (<=jdm; 0 implies jdm)
 51	'botio ' = bathymetry       I/O unit (0 no I/O)
 21	'flxio ' = surf. heat flux  I/O unit (0 no I/O)
  0	'empio ' = surf. evap-pcip  I/O unit (0 no I/O)
  0	'ttrio ' = surf. temp trend I/O unit (0 no I/O)
  0	'strio ' = surf. saln trend I/O unit (0 no I/O)
  0	'icvio ' = ice coverage     I/O unit (0 no I/O)
  0	'ithio ' = ice thickness    I/O unit (0 no I/O)
  0	'ictio ' = ice temperature  I/O unit (0 no I/O)
 23	'sshio ' = sea surf. height I/O unit (0 no I/O)
  0	'bsfio ' = baro. strmfn.    I/O unit (0 no I/O)
  0	'uvmio ' = mix. lay. u-vel. I/O unit (0 no I/O)
  0	'vvmio ' = mix. lay. v-vel. I/O unit (0 no I/O)
  0	'spmio ' = mix. lay. speed  I/O unit (0 no I/O)
 25	'bltio ' = bnd. lay. thick. I/O unit (0 no I/O)
 24	'mltio ' = mix. lay. thick. I/O unit (0 no I/O)
  0	'sstio ' = mix. lay. temp.  I/O unit (0 no I/O)
  0	'sssio ' = mix. lay. saln.  I/O unit (0 no I/O)
  0	'ssdio ' = mix. lay. dens.  I/O unit (0 no I/O)
-$k	'kf    ' = first output layer (=0 end output; <0 label with layer #)
 $k	'kl    ' = last  output layer
 37	'uvlio ' = layer k   u-vel. I/O unit (0 no I/O)
 38	'vvlio ' = layer k   v-vel. I/O unit (0 no I/O)
  0	'splio ' = layer k   speed. I/O unit (0 no I/O)
 40	'iwvio ' = layer k   i.vel. I/O unit (0 no I/O)
 31	'infio ' = layer k   i.dep. I/O unit (0 no I/O)
 32	'thkio ' = layer k   thick. I/O unit (0 no I/O)
 33	'temio ' = layer k   temp   I/O unit (0 no I/O)
 34	'salio ' = layer k   saln.  I/O unit (0 no I/O)
 35	'tthio ' = layer k   dens,  I/O unit (0 no I/O)
  0	'sfnio ' = layer k  strmfn. I/O unit (0 no I/O)
  0	'kf    ' = first output layer (=0 end output; <0 label with layer #)
E-o-D
    end
  end
end
