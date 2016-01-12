#!/bin/csh
#
set echo
#
# --- extract 2-d fields from one year's archives.
# --- configured for monthly archives and 16 layers.
#
# --- output can be formatted, unformatted (BINARY), .[ab] (HYCOM), or netCDF.
#
# --- output is HYCOM .a files, converted to "raw" .A files.
# --- all fields are merged into a single annual file.
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
setenv D ../expt_02.4/data
#
foreach y ( 0031 )
  foreach d ( 016 046 076 106 136 166 196 226 256 286 316 346 )
    setenv FOR051A ${D}/archv.${y}_${d}_00_bot.a
    setenv FOR051  ${D}/archv.${y}_${d}_00_bot.b
    setenv FOR021A ${D}/archv.${y}_${d}_00_flx.a
    setenv FOR021  ${D}/archv.${y}_${d}_00_flx.b
    setenv FOR022A ${D}/archv.${y}_${d}_00_ice.a
    setenv FOR022  ${D}/archv.${y}_${d}_00_ice.b
    setenv FOR023A ${D}/archv.${y}_${d}_00_fsd.a
    setenv FOR023  ${D}/archv.${y}_${d}_00_fsd.b
    setenv FOR024A ${D}/archv.${y}_${d}_00_mix.a
    setenv FOR024  ${D}/archv.${y}_${d}_00_mix.b
    setenv FOR031A ${D}/archv.${y}_${d}_00_3di.a
    setenv FOR031  ${D}/archv.${y}_${d}_00_3di.b
    setenv FOR032A ${D}/archv.${y}_${d}_00_3dh.a
    setenv FOR032  ${D}/archv.${y}_${d}_00_3dh.b
    setenv FOR033A ${D}/archv.${y}_${d}_00_3dt.a
    setenv FOR033  ${D}/archv.${y}_${d}_00_3dt.b
    setenv FOR034A ${D}/archv.${y}_${d}_00_3ds.a
    setenv FOR034  ${D}/archv.${y}_${d}_00_3ds.b
    setenv FOR035A ${D}/archv.${y}_${d}_00_3dr.a
    setenv FOR035  ${D}/archv.${y}_${d}_00_3dr.b
    setenv FOR036A ${D}/archv.${y}_${d}_00_3dn.a
    setenv FOR036  ${D}/archv.${y}_${d}_00_3dn.b
    setenv FOR037A ${D}/archv.${y}_${d}_00_3du.a
    setenv FOR037  ${D}/archv.${y}_${d}_00_3du.b
    setenv FOR038A ${D}/archv.${y}_${d}_00_3dv.a
    setenv FOR038  ${D}/archv.${y}_${d}_00_3dv.b
    setenv FOR039A ${D}/archv.${y}_${d}_00_3dm.a
    setenv FOR039  ${D}/archv.${y}_${d}_00_3dm.b
    /bin/rm $FOR051
    /bin/rm $FOR051A
    /bin/rm $FOR021  $FOR022  $FOR023  $FOR024
    /bin/rm $FOR021A $FOR022A $FOR023A $FOR024A
    /bin/rm $FOR031  $FOR032  $FOR033  $FOR034  $FOR035
    /bin/rm $FOR031A $FOR032A $FOR033A $FOR034A $FOR035A
    /bin/rm $FOR036  $FOR037  $FOR038  $FOR039
    /bin/rm $FOR036A $FOR037A $FOR038A $FOR039A
    ../../ALL/archive/src/archv2data2d <<E-o-D
${D}/archv.${y}_${d}_00.a
HYCOM
000	'iexpt ' = experiment number x10 (000=from archive file)
  0	'yrflag' = days in year flag (0=360J16,1=366J16,2=366J01,3-actual)
 57	'idm   ' = longitudinal array size
 52	'jdm   ' = latitudinal  array size
 16	'kdm   ' = number of layers
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
 -1	'kf    ' = first output layer (=0 end output; <0 label with layer #)
 16	'kl    ' = last  output layer
 37	'uvlio ' = layer k   u-vel. I/O unit (0 no I/O)
 38	'vvlio ' = layer k   v-vel. I/O unit (0 no I/O)
  0	'splio ' = layer k   speed. I/O unit (0 no I/O)
 31	'infio ' = layer k   i.dep. I/O unit (0 no I/O)
 32	'thkio ' = layer k   thick. I/O unit (0 no I/O)
 33	'temio ' = layer k   temp   I/O unit (0 no I/O)
 34	'salio ' = layer k   saln.  I/O unit (0 no I/O)
 35	'tthio ' = layer k   dens,  I/O unit (0 no I/O)
  0	'sfnio ' = layer k  strmfn. I/O unit (0 no I/O)
  0	'kf    ' = first output layer (=0 end output; <0 label with layer #)
E-o-D
#
# --- convert HYCOM .a files to RAW files (no padding, spval=1.e10).
# --- comment this out if you don't need RAW files.
#
    foreach t ( bot flx ice fsd mix 3di 3dh 3dt 3ds 3dr 3dn 3du 3dv 3dm )
      if (-e ${D}/archv.${y}_${d}_00_${t}.a) then
        /bin/rm -f ${D}/archv.${y}_${d}_00_${t}.A
        ../../ALL/bin/hycom2raw ${D}/archv.${y}_${d}_00_${t}.a 57 52 1.e10 ${D}/archv.${y}_${d}_00_${t}.A
        /bin/rm -f ${D}/archv.${y}_${d}_00_${t}.a
      endif
    end
  end
end
#
# --- merge individual files into yearly files.
# --- comment this out if you don't need yearly files.
#
cd ${D}
foreach t ( bot flx ice fsd mix 3di 3dh 3dt 3ds 3dr 3dn 3du 3dv 3dm )
  if (-e archv.${y}_016_00_${t}.A) then
    cat archv.${y}_???_00_${t}.A >! archv.${y}_${t}.A
  endif
  if (-e archv.${y}_016_00_${t}.a) then
    cat archv.${y}_???_00_${t}.a >! archv.${y}_${t}.a
  endif
  if (-e archv.${y}_016_00_${t}.b) then
    cat archv.${y}_???_00_${t}.b >! archv.${y}_${t}.b
  endif
# /bin/rm -f archv.${y}_???_00_${t}.A
# /bin/rm -f archv.${y}_???_00_${t}.a
# /bin/rm -f archv.${y}_???_00_${t}.b
end
