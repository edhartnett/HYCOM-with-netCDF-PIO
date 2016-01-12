#!/bin/csh
#
set echo
#
# --- extract 2-d fields from a single HYCOM archive file.
# --- configured for 22 layers.
#
# --- output can be formatted, unformatted (BINARY), .[ab] (HYCOM).
# --- or use archv2ncdf2d for netCDF output.
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
foreach y ( 0021 )
  foreach d ( 016 )
    setenv FOR051A ${D}/archv.${y}_${d}_00_bot.a
    setenv FOR051  ${D}/archv.${y}_${d}_00_bot.b
    setenv FOR041A ${D}/archv.${y}_${d}_00_flx.a
    setenv FOR041  ${D}/archv.${y}_${d}_00_flx.b
    setenv FOR042A ${D}/archv.${y}_${d}_00_emp.a
    setenv FOR042  ${D}/archv.${y}_${d}_00_emp.b
    setenv FOR043A ${D}/archv.${y}_${d}_00_bft.a
    setenv FOR043  ${D}/archv.${y}_${d}_00_bft.b
    setenv FOR044A ${D}/archv.${y}_${d}_00_bfs.a
    setenv FOR044  ${D}/archv.${y}_${d}_00_bfs.b
    setenv FOR045A ${D}/archv.${y}_${d}_00_bfa.a
    setenv FOR045  ${D}/archv.${y}_${d}_00_bfa.b
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
    setenv FOR039A ${D}/archv.${y}_${d}_00_3dc.a
    setenv FOR039  ${D}/archv.${y}_${d}_00_3dc.b
    setenv FOR040A ${D}/archv.${y}_${d}_00_3dw.a
    setenv FOR040  ${D}/archv.${y}_${d}_00_3dw.b
    /bin/rm $FOR051
    /bin/rm $FOR051A
    /bin/rm          $FOR022  $FOR023  $FOR024
    /bin/rm          $FOR022A $FOR023A $FOR024A
    /bin/rm $FOR031  $FOR032  $FOR033  $FOR034  $FOR035
    /bin/rm $FOR031A $FOR032A $FOR033A $FOR034A $FOR035A
    /bin/rm $FOR036  $FOR037  $FOR038  $FOR039  $FOR040
    /bin/rm $FOR036A $FOR037A $FOR038A $FOR039A $FOR040A
    /bin/rm $FOR041  $FOR042  $FOR043  $FOR044  $FOR045
    /bin/rm $FOR041A $FOR042A $FOR043A $FOR044A $FOR045A
    ../../ALL/archive/src/archv2data2d <<E-o-D
${D}/archv.${y}_${d}_00.a
HYCOM
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
 41	'flxio ' = surf. heat flux  I/O unit (0 no I/O)
 42	'empio ' = surf. evap-pcip  I/O unit (0 no I/O)
 43	'tbfio ' = temp. bou. flux  I/O unit (0 no I/O)
 44	'sbfio ' = saln. bou. flux  I/O unit (0 no I/O)
 45	'abfio ' = total bou. flux  I/O unit (0 no I/O)
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
 22	'kl    ' = last  output layer
 37	'uvlio ' = layer k   u-vel. I/O unit (0 no I/O)
 38	'vvlio ' = layer k   v-vel. I/O unit (0 no I/O)
  0	'splio ' = layer k   speed. I/O unit (0 no I/O)
  0	'iwvio ' = layer k   i.vel. I/O unit (0 no I/O)
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
    foreach t ( bot flx emp bft bfs bfa ice fsd mix 3di 3dh 3dt 3ds 3dr 3dn 3du 3dv 3dm )
      if (-e ${D}/archv.${y}_${d}_00_${t}.a) then
        /bin/rm -f ${D}/archv.${y}_${d}_00_${t}.A
        ../../ALL/bin/hycom2raw ${D}/archv.${y}_${d}_00_${t}.a 57 52 1.e10 ${D}/archv.${y}_${d}_00_${t}.A
#       /bin/rm -f ${D}/archv.${y}_${d}_00_${t}.a
      endif
    end
  end
end
