#!/bin/csh
#
set echo
#
# --- extract 2-d fields from a single HYCOM archive file.
# --- configured for 22 layers.
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
    setenv CDF021  ${D}/archv.${y}_${d}_00_flx.nc
    setenv CDF022  ${D}/archv.${y}_${d}_00_ice.nc
    setenv CDF023  ${D}/archv.${y}_${d}_00_fsd.nc
    setenv CDF024  ${D}/archv.${y}_${d}_00_mix.nc
    setenv CDF031  ${D}/archv.${y}_${d}_00_3di.nc
    setenv CDF032  ${D}/archv.${y}_${d}_00_3dh.nc
    setenv CDF033  ${D}/archv.${y}_${d}_00_3dt.nc
    setenv CDF034  ${D}/archv.${y}_${d}_00_3ds.nc
    setenv CDF035  ${D}/archv.${y}_${d}_00_3dr.nc
    setenv CDF036  ${D}/archv.${y}_${d}_00_3dn.nc
    setenv CDF037  ${D}/archv.${y}_${d}_00_3du.nc
    setenv CDF038  ${D}/archv.${y}_${d}_00_3dv.nc
    setenv CDF039  ${D}/archv.${y}_${d}_00_3dc.nc
    setenv CDF040  ${D}/archv.${y}_${d}_00_3dw.nc
    /bin/rm $CDF051
    /bin/rm $CDF021  $CDF022  $CDF023  $CDF024
    /bin/rm $CDF031  $CDF032  $CDF033  $CDF034  $CDF035
    /bin/rm $CDF036  $CDF037  $CDF038  $CDF039  $CDF040
    ../../ALL/archive/src/archv2ncdf2d <<E-o-D
${D}/archv.${y}_${d}_00.a
netCDF
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
 24	'uvmio ' = mix. lay. u-vel. I/O unit (0 no I/O)
 24	'vvmio ' = mix. lay. v-vel. I/O unit (0 no I/O)
  0	'spmio ' = mix. lay. speed  I/O unit (0 no I/O)
 24	'bltio ' = bnd. lay. thick. I/O unit (0 no I/O)
 24	'mltio ' = mix. lay. thick. I/O unit (0 no I/O)
 24	'sstio ' = mix. lay. temp.  I/O unit (0 no I/O)
 24	'sssio ' = mix. lay. saln.  I/O unit (0 no I/O)
 24	'ssdio ' = mix. lay. dens.  I/O unit (0 no I/O)
  1	'kf    ' = first output layer (=0 end output; <0 label with layer #)
 22	'kl    ' = last  output layer
 37	'uvlio ' = layer k   u-vel. I/O unit (0 no I/O)
 38	'vvlio ' = layer k   v-vel. I/O unit (0 no I/O)
 39	'splio ' = layer k   speed. I/O unit (0 no I/O)
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
